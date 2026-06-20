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
    """Zero-init from exact pkl cols, fill via alias map + direct name match."""
    FIELD_MAP = {
        'age': 'Diagnosis Age', 'fraction_genome_altered': 'Fraction Genome Altered',
        'aneuploidy_score': 'Aneuploidy Score', 'buffa_hypoxia_score': 'Buffa Hypoxia Score',
        'ragnum_hypoxia_score': 'Ragnum Hypoxia Score', 'winter_hypoxia_score': 'Winter Hypoxia Score',
        'msi_mantis_score': 'MSI MANTIS Score', 'msisensor_score': 'MSIsensor Score',
        'mutation_count': 'Mutation Count', 'tmb_nonsynonymous': 'TMB (nonsynonymous)',
        'sex_encoded': 'Sex_encoded', 'race_encoded': 'Race_encoded',
        'mean_vaf': 'mean_VAF', 'max_vaf': 'max_VAF', 'std_vaf': 'std_VAF',
        'n_high_impact': 'n_HIGH_impact', 'n_moderate_impact': 'n_MODERATE_impact',
        'n_low_impact': 'n_LOW_impact', 'n_modifier_impact': 'n_MODIFIER_impact',
    }
    row = {col: 0.0 for col in col_list}
    for col in col_list:
        if col in features:
            row[col] = float(features[col])
    for k, v in FIELD_MAP.items():
        if k in features and v in row:
            row[v] = float(features[k])

    total = float(features.get('total_mutations', 1)) or 1
    if 'domain_coverage_frac' in row and row['domain_coverage_frac'] == 0.0:
        row['domain_coverage_frac'] = float(features.get('n_with_protein_domain', 0)) / total
    if 'snp_fraction' in row and row['snp_fraction'] == 0.0:
        n_snp, n_indel = float(features.get('n_snp', 0)), float(features.get('n_indel', 0))
        row['snp_fraction'] = n_snp / (n_snp + n_indel) if (n_snp + n_indel) else 0.0
    if 'ratio_nonsyn_syn' in row and row['ratio_nonsyn_syn'] == 0.0:
        n_syn = float(features.get('n_synonymous_csq', 0))
        row['ratio_nonsyn_syn'] = float(features.get('n_nonsynonymous_csq', 0)) / n_syn if n_syn else 0.0
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

        age = float(data['age'])
        total_mutations = float(data['total_mutations'])

        if not (0 <= age <= 120):
            return jsonify({'error': 'Invalid age'}), 400
        if total_mutations < 0:
            return jsonify({'error': 'Total mutations cannot be negative'}), 400

        def g(key, default):
            return float(data.get(key, default))

        features = {
            # Clinical
            'age':                     age,
            'sex_encoded':             int(data.get('sex_encoded', 0)),
            'race_encoded':            int(data.get('race_encoded', 0)),
            'fraction_genome_altered': g('fraction_genome_altered', 0.15),
            'aneuploidy_score':        g('aneuploidy_score', 0.5),
            'mutation_count':          g('mutation_count', total_mutations),
            'tmb_nonsynonymous':       g('tmb_nonsynonymous', 4.5),
            'buffa_hypoxia_score':     g('buffa_hypoxia_score', 0.0),
            'ragnum_hypoxia_score':    g('ragnum_hypoxia_score', 0.0),
            'winter_hypoxia_score':    g('winter_hypoxia_score', 0.0),
            'msi_mantis_score':        g('msi_mantis_score', 0.0),
            'msisensor_score':         g('msisensor_score', 0.0),
            # Genomic — burden, calling, VAF
            'total_mutations':         total_mutations,
            'unique_genes':            g('unique_genes', 80.0),
            'unique_chromosomes':      g('unique_chromosomes', 12.0),
            'tmb':                     g('tmb', 150.0),
            'mean_vaf':                g('mean_vaf', 0.28),
            'max_vaf':                 g('max_vaf', 0.55),
            'std_vaf':                 g('std_vaf', 0.12),
            'mean_t_depth':            g('mean_t_depth', 110.0),
            'mean_n_depth':            g('mean_n_depth', 105.0),
            'mean_t_alt_count':        g('mean_t_alt_count', 28.0),
            # Functional impact
            'mean_impact_ord':         g('mean_impact_ord', 1.8),
            'n_high_impact':           g('n_high_impact', 5.0),
            'n_moderate_impact':       g('n_moderate_impact', 90.0),
            'n_low_impact':            g('n_low_impact', 30.0),
            'n_modifier_impact':       g('n_modifier_impact', 25.0),
            # Protein damage prediction
            'mean_pp_ord':             g('mean_pp_ord', 1.1),
            'mean_pp_score':           g('mean_pp_score', 0.72),
            'n_probably_damaging':     g('n_probably_damaging', 35.0),
            'n_possibly_damaging':     g('n_possibly_damaging', 20.0),
            'n_pp_benign':             g('n_pp_benign', 30.0),
            'mean_sift_ord':           g('mean_sift_ord', 1.4),
            'mean_sift_score':         g('mean_sift_score', 0.04),
            'n_deleterious_sift':      g('n_deleterious_sift', 55.0),
            'n_tolerated_sift':        g('n_tolerated_sift', 25.0),
            # Variant consequence flags
            'has_stop_gained':         int(data.get('has_stop_gained', 0)),
            'has_frameshift':          int(data.get('has_frameshift', 0)),
            'has_splice':              int(data.get('has_splice', 0)),
            'has_utr_variant':         int(data.get('has_utr_variant', 0)),
            'has_intronic':            int(data.get('has_intronic', 0)),
            'has_downstream':          int(data.get('has_downstream', 0)),
            # Sequence composition
            'n_synonymous_csq':        g('n_synonymous_csq', 28.0),
            'n_nonsynonymous_csq':     g('n_nonsynonymous_csq', 90.0),
            'n_snp':                   g('n_snp', 140.0),
            'n_indel':                 g('n_indel', 10.0),
            # Database & annotation
            'n_cosmic_hits':           g('n_cosmic_hits', 3.0),
            'n_rare_population':       g('n_rare_population', 135.0),
            'n_with_protein_domain':   g('n_with_protein_domain', 80.0),
            'mean_ncallers':           g('mean_ncallers', 4.0),
            'max_ncallers':            g('max_ncallers', 5.0),
        }

        result = predict_risk(features)

        user_id = session['user_id']
        rt_ref = db.reference(f'/predictions/{user_id}')
        rt_ref.push({
            'userEmail': session.get('user_email'),
            'timestamp': datetime.utcnow().isoformat(),
            'input':     features,
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