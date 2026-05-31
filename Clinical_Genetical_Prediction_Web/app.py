from datetime import datetime
import re
import os
import pickle
import pandas as pd
import numpy as np
from flask import Flask, render_template, request, jsonify, session, redirect, url_for
from functools import wraps
import firebase_admin
from firebase_admin import credentials, auth, firestore, db
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'fallback-key')  # Use environment variable for secret key in production

# Initialize Firebase Admin SDK
# Place your firebase_credentials.json in root directory
cred = credentials.Certificate('Secrets/Config/Admin_SDK.json')
firebase_admin.initialize_app(cred, {
    'databaseURL': 'https://dna-scan-default-rtdb.europe-west1.firebasedatabase.app'  # Update with your Realtime DB URL
})

# Initialize Firestore
firestore_client = firestore.client()



# Authentication decorator
def login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if 'user_id' not in session:
            return jsonify({'error': 'Authentication required'}), 401
        return f(*args, **kwargs)
    return decorated_function

def _build_row(col_list, features):
    """Zero-init from exact pkl cols → fill via aliases + direct match + variant_* fallback."""
    FIELD_MAP = {
        'age':'Diagnosis Age','fraction_genome_altered':'Fraction Genome Altered',
        'aneuploidy_score':'Aneuploidy Score','buffa_hypoxia_score':'Buffa Hypoxia Score',
        'ragnum_hypoxia_score':'Ragnum Hypoxia Score','winter_hypoxia_score':'Winter Hypoxia Score',
        'msi_mantis_score':'MSI MANTIS Score','msisensor_score':'MSIsensor Score',
        'mutation_count':'Mutation Count','tmb_nonsynonymous':'TMB (nonsynonymous)',
        'sex_encoded':'Sex_encoded','race_encoded':'Race_encoded',
        'total_mutations':'total_mutations','unique_genes':'unique_genes',
        'missense_mutations':'missense_mutations','silent_mutations':'silent_mutations',
        'high_impact_mutations':'high_impact_mutations','tmb':'tmb',
        'unique_chromosomes':'unique_chromosomes',
        'n_moderate_impact':'n_MODERATE_impact','n_low_impact':'n_LOW_impact',
        'n_modifier_impact':'n_MODIFIER_impact','n_high_impact':'n_HIGH_impact',
        'mean_impact_ord':'mean_impact_ord','max_impact_ord':'max_impact_ord',
        'mean_vaf':'mean_VAF','max_vaf':'max_VAF','std_vaf':'std_VAF',
        'mean_t_depth':'mean_t_depth','mean_n_depth':'mean_n_depth',
        'mean_t_alt_count':'mean_t_alt_count',
        'mean_pp_ord':'mean_pp_ord','mean_pp_score':'mean_pp_score',
        'n_probably_damaging':'n_probably_damaging','n_possibly_damaging':'n_possibly_damaging',
        'n_pp_benign':'n_pp_benign','mean_sift_ord':'mean_sift_ord',
        'mean_sift_score':'mean_sift_score','n_deleterious_sift':'n_deleterious_sift',
        'n_tolerated_sift':'n_tolerated_sift','has_stop_gained':'has_stop_gained',
        'has_frameshift':'has_frameshift','has_splice':'has_splice',
        'has_utr_variant':'has_utr_variant','has_intronic':'has_intronic',
        'has_downstream':'has_downstream','n_synonymous_csq':'n_synonymous_csq',
        'n_nonsynonymous_csq':'n_nonsynonymous_csq','ratio_nonsyn_syn':'ratio_nonsyn_syn',
        'n_snp':'n_snp','n_indel':'n_indel','snp_fraction':'snp_fraction',
        'n_cosmic_hits':'n_cosmic_hits','n_rare_population':'n_rare_population',
        'n_with_protein_domain':'n_with_protein_domain','domain_coverage_frac':'domain_coverage_frac',
        'mean_ncallers':'mean_ncallers','max_ncallers':'max_ncallers',
    }
    row = {col: 0.0 for col in col_list}
    for col in col_list:
        if col in features:
            row[col] = float(features[col])
    for k, v in FIELD_MAP.items():
        if k in features and v in row:
            row[v] = float(features[k])
    total = float(features.get('total_mutations', 1)) or 1
    for vcol, val in {
        "variant_Missense_Mutation": features.get('missense_mutations', 0),
        "variant_Silent":            features.get('silent_mutations', 0),
        "variant_Nonsense_Mutation": features.get('high_impact_mutations', 0),
        "variant_3'UTR":0,"variant_5'UTR":0,"variant_Intron":0,
        "variant_Frame_Shift_Del":0,"variant_Splice_Site":0,
        "variant_In_Frame_Del":0,"variant_In_Frame_Ins":0,
    }.items():
        if vcol in row:
            row[vcol] = float(val)
    if 'domain_coverage_frac' in row and row['domain_coverage_frac'] == 0.0:
        row['domain_coverage_frac'] = float(features.get('n_with_protein_domain', 0)) / total
    return row


def predict_risk(features: dict) -> dict:
    """Two-stage: genomic → Model 1 → pred_probs → Model 2 → Cox."""
    # Stage 1
    m1_row = _build_row(m1_feature_cols, features)
    X_m1   = pd.DataFrame([m1_row])[m1_feature_cols]
    X_m1_s = scaler_m1.transform(X_m1)
    pred_prob_dict = {}
    for i, model in enumerate(variant_models):
        prob = (float(variant_mean_probs[i]) if model is None
                else float(model.predict_proba(X_m1_s)[0, 1]))
        pred_prob_dict[f'pred_prob_{top_variants[i]}'] = prob

    # Stage 2
    combined = {**features, **pred_prob_dict}
    m2_row   = _build_row(all_feature_cols_m2, combined)
    for k, v in pred_prob_dict.items():
        if k in m2_row:
            m2_row[k] = v
    X_m2   = pd.DataFrame([m2_row])[all_feature_cols_m2]
    X_m2_s = scaler_m2.transform(X_m2)
    X_pca  = pca.transform(X_m2_s)

    stage_enc   = gb_stage.predict(X_pca[:, :25])[0]
    stage_label = stage_encoder.inverse_transform([stage_enc])[0]
    stage_conf  = float(np.max(gb_stage.predict_proba(X_pca[:, :25])[0]))
    variant_probs = {top_variants[i]: round(v, 4) for i, v in enumerate(pred_prob_dict.values())}

    cox_df     = pd.DataFrame(X_pca[:, :15], columns=[f'PC_{i}' for i in range(15)])
    risk_score = float(cph.predict_partial_hazard(cox_df).values[0])

    median_survival = None
    try:
        surv = cph.predict_survival_function(cox_df)
        times, probs_arr = surv.index.values, surv.iloc[:, 0].values
        idx = (probs_arr <= 0.5).argmax() if (probs_arr <= 0.5).any() else -1
        if idx != -1 and not np.isnan(times[idx]):
            median_survival = float(times[idx])
    except Exception:
        pass

    risk_category = ("LOW RISK"      if risk_score <= risk_percentiles[0] else
                     "MODERATE RISK" if risk_score <= risk_percentiles[1] else "HIGH RISK")
    return {
        'risk_score':             round(risk_score, 4),
        'risk_category':          risk_category,
        'predicted_stage':        stage_label,
        'stage_confidence':       round(stage_conf * 100, 2),
        'variant_probabilities':  variant_probs,
        'median_survival_months': round(median_survival, 1) if median_survival else None,
        'median_survival_years':  round(median_survival / 12, 1) if median_survival else None,
    }


@app.route('/')
def home():
    """Landing page presenting the Cox model and its capabilities."""
    return render_template('home.html')

@app.route('/profile_setup')
@login_required
def profile_setup():
    return render_template('profile_setup.html')


@app.route('/login')
def login_page():
    """Login page with Firebase auth integration"""
    return render_template('login.html')

@app.route('/dashboard')
def dashboard():
    if 'user_id' not in session:
        return redirect(url_for('login_page'))
    
    user_id = session['user_id']
    doc = firestore_client.collection('users').document(user_id).get()
    
    if doc.exists:
        profile = doc.to_dict()
        # Only require name and email (not full medical profile) for dashboard access
        if not profile.get('name'):
            return redirect(url_for('profile_setup'))
    else:
        # Create minimal user document if missing
        firestore_client.collection('users').document(user_id).set({
            'email': session.get('user_email'),
            'created_at': firestore.SERVER_TIMESTAMP
        })
        return redirect(url_for('profile_setup'))
    
    return render_template('dashboard.html', user_email=session.get('user_email'))

@app.route('/api/check_profile', methods=['GET'])
@login_required
def check_profile():
    """Check if user has completed profile setup"""
    user_id = session['user_id']
    doc = firestore_client.collection('users').document(user_id).get()
    if doc.exists:
        profile = doc.to_dict()
        # Check required fields
        required = ['name', 'date_of_birth', 'height', 'weight', 'sex', 'has_cancer_history']
        if all(k in profile and profile[k] not in (None, '') for k in required):
            return jsonify({'profile_complete': True})
    return jsonify({'profile_complete': False})

@app.route('/api/save_profile', methods=['POST'])
@login_required
def save_profile():
    """Save or update user profile information"""
    try:
        data = request.get_json()
        user_id = session['user_id']
        
        # Validate inputs
        name = data.get('name', '').strip()
        if not name or len(name) < 2:
            return jsonify({'error': 'Valid name is required'}), 400
        
        dob = data.get('date_of_birth')  # Expect YYYY-MM-DD
        if not dob or not re.match(r'^\d{4}-\d{2}-\d{2}$', dob):
            return jsonify({'error': 'Valid date of birth (YYYY-MM-DD) required'}), 400
        
        # Optional: validate age > 0 and < 120
        try:
            birth_date = datetime.strptime(dob, '%Y-%m-%d')
            today = datetime.now()
            age = today.year - birth_date.year - ((today.month, today.day) < (birth_date.month, birth_date.day))
            if age < 0 or age > 120:
                return jsonify({'error': 'Invalid age'}), 400
        except:
            return jsonify({'error': 'Invalid date format'}), 400
        
        height = float(data.get('height', 0))
        if height <= 0 or height > 300:
            return jsonify({'error': 'Height must be between 1 and 300 cm'}), 400
        
        weight = float(data.get('weight', 0))
        if weight <= 0 or weight > 500:
            return jsonify({'error': 'Weight must be between 1 and 500 kg'}), 400
        
        sex = data.get('sex')
        if sex not in ['Male', 'Female', 'Other']:
            return jsonify({'error': 'Sex must be Male, Female, or Other'}), 400
        
        has_cancer_history = bool(data.get('has_cancer_history', False))
        
        # Update Firestore
        user_ref = firestore_client.collection('users').document(user_id)
        user_ref.set({
            'name': name,
            'date_of_birth': dob,
            'height': height,
            'weight': weight,
            'sex': sex,
            'has_cancer_history': has_cancer_history,
            'email': session.get('user_email'),
            'updated_at': firestore.SERVER_TIMESTAMP
        }, merge=True)
        
        return jsonify({'success': True, 'message': 'Profile saved'})
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/verify_token', methods=['POST'])
def verify_token():
    """Verify Firebase ID token and create session"""
    try:
        data = request.get_json()
        id_token = data.get('id_token')
        
        # Verify the ID token
        decoded_token = auth.verify_id_token(id_token)
        user_id = decoded_token['uid']
        user_email = decoded_token.get('email', '')
        
        # Store in session
        session['user_id'] = user_id
        session['user_email'] = user_email
        
        # Create/update user document in Firestore
        user_ref = firestore_client.collection('users').document(user_id)
        # Prepare user data
        user_data = {
            'email': user_email,
            'last_login': firestore.SERVER_TIMESTAMP,
        }
        # Only set created_at if this is a new user (document doesn't exist yet)
        if not user_ref.get().exists:
            user_data['created_at'] = firestore.SERVER_TIMESTAMP

        user_ref.set(user_data, merge=True)
        return jsonify({'success': True, 'user_id': user_id})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 401


@app.route('/api/get_profile', methods=['GET'])
@login_required
def get_profile():
    """Fetch user profile data (name, etc.)"""
    user_id = session['user_id']
    doc = firestore_client.collection('users').document(user_id).get()
    if doc.exists:
        data = doc.to_dict()
        return jsonify({'success': True, 'profile': data})
    return jsonify({'success': False, 'error': 'Profile not found'}), 404





@app.route('/api/logout', methods=['POST'])
def logout():
    """Clear session"""
    session.clear()
    return jsonify({'success': True})

@app.route('/api/predict', methods=['POST'])
@login_required
def predict():
    """Run prediction and store results in Firestore & Realtime DB"""
    try:
        data = request.get_json()
        
        # Extract & validate core fields
        age               = float(data['age'])
        stage_data        = int(data['stage_data'])
        total_mutations   = int(data['total_mutations'])
        missense_mutations = int(data['missense_mutations'])
        high_impact_mutations = int(data['high_impact_mutations'])

        if not (0 <= age <= 120):
            return jsonify({'error': 'Invalid age'}), 400
        if stage_data not in [0, 1, 2, 3]:
            return jsonify({'error': 'Stage must be 0-3'}), 400
        if total_mutations < 0 or missense_mutations < 0 or high_impact_mutations < 0:
            return jsonify({'error': 'Mutations cannot be negative'}), 400
        if missense_mutations > total_mutations or high_impact_mutations > total_mutations:
            return jsonify({'error': 'Mutation counts exceed total'}), 400

        # Build full feature dict (all notebook fields)
        features = {
            'age':                     age,
            'stage_data':              stage_data,
            'total_mutations':         total_mutations,
            'missense_mutations':      missense_mutations,
            'high_impact_mutations':   high_impact_mutations,
            'unique_genes':            float(data.get('unique_genes', 80)),
            'silent_mutations':        float(data.get('silent_mutations', 30)),
            'tmb':                     float(data.get('tmb', 5.2)),
            'tmb_nonsynonymous':       float(data.get('tmb_nonsynonymous', 4.5)),
            'mutation_count':          float(data.get('mutation_count', total_mutations)),
            'fraction_genome_altered': float(data.get('fraction_genome_altered', 0.15)),
            'aneuploidy_score':        float(data.get('aneuploidy_score', 0.5)),
            'buffa_hypoxia_score':     float(data.get('buffa_hypoxia_score', 0.0)),
            'ragnum_hypoxia_score':    float(data.get('ragnum_hypoxia_score', 0.0)),
            'winter_hypoxia_score':    float(data.get('winter_hypoxia_score', 0.0)),
            'msi_mantis_score':        float(data.get('msi_mantis_score', 0.0)),
            'msisensor_score':         float(data.get('msisensor_score', 0.0)),
            'sex_encoded':             int(data.get('sex_encoded', 0)),
            'race_encoded':            int(data.get('race_encoded', 0)),
        }

        # Make prediction
        result = predict_risk(features)

        # Store under /predictions/{uid}/ — matches RTDB rules
        user_id = session['user_id']
        rt_ref = db.reference(f'/predictions/{user_id}')
        rt_ref.push({
            'userEmail': session.get('user_email'),
            'timestamp': datetime.utcnow().isoformat(),
            'input':     {k: v for k, v in features.items()},
            'result':    result,
        })
        
        
        
        return jsonify({
            'success': True,
            'result': result,
            'message': 'Prediction stored successfully'
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/delete_prediction/<prediction_id>', methods=['DELETE'])
@login_required
def delete_prediction(prediction_id):
    try:
        user_id = session['user_id']
        db.reference(f'/predictions/{user_id}/{prediction_id}').delete()
        return jsonify({'success': True})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/history', methods=['GET'])
@login_required
def get_history():
    try:
        user_id = session['user_id']
        # Read only this user's predictions — aligned with RTDB rules
        predictions_ref = db.reference(f'/predictions/{user_id}')
        all_data = predictions_ref.get() or {}

        history = [
            {
                'id':        key,
                'timestamp': value.get('timestamp'),
                'input':     value.get('input'),
                'result':    value.get('result')
            }
            for key, value in all_data.items()
        ]
        history.sort(key=lambda x: x['timestamp'], reverse=True)
        history = history[:20]
        return jsonify({'success': True, 'history': history})
    except Exception as e:
        print(f"Error: {e}")
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    import joblib
    print("Loading prediction models...")
    model_dir = './saved_models'
    scaler_m1           = joblib.load(f'{model_dir}/scaler_m1.pkl')
    scaler_m2           = joblib.load(f'{model_dir}/scaler_m2.pkl')
    pca                 = joblib.load(f'{model_dir}/pca.pkl')
    stage_encoder       = joblib.load(f'{model_dir}/stage_encoder.pkl')
    gb_stage            = joblib.load(f'{model_dir}/gb_stage.pkl')
    variant_models      = joblib.load(f'{model_dir}/variant_models.pkl')
    variant_mean_probs  = joblib.load(f'{model_dir}/variant_mean_probs.pkl')
    top_variants        = joblib.load(f'{model_dir}/top_variants.pkl')
    cph                 = joblib.load(f'{model_dir}/cox_model.pkl')
    risk_percentiles    = joblib.load(f'{model_dir}/risk_percentiles.pkl')
    m1_feature_cols     = joblib.load(f'{model_dir}/m1_feature_cols.pkl')
    all_feature_cols    = joblib.load(f'{model_dir}/all_feature_cols.pkl')
    all_feature_cols_m2 = joblib.load(f'{model_dir}/all_feature_cols_m2.pkl')
    race_encoder        = joblib.load(f'{model_dir}/race_encoder.pkl')
    print(f"Models loaded — M1:{len(m1_feature_cols)} cols  M2:{len(all_feature_cols_m2)} cols")
    app.run(debug=True, host='0.0.0.0', port=5000)