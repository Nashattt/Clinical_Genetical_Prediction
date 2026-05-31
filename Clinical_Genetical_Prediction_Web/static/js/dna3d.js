import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

let scene, camera, renderer, controls, dnaGroup;
let spheres      = [];
let cylinders    = [];
let nodeMetadata = [];   // { gene, variant, prob, impact, css }
let currentData  = null;

// ── Variant → color ───────────────────────────────────────────────────────────
const VARIANT_META = {
    'Missense_Mutation' : { hex: 0xffaa44, css: '#ffaa44', impact: 'missense' },
    'Silent'            : { hex: 0x00ff9d, css: '#00ff9d', impact: 'normal'   },
    "3'UTR"             : { hex: 0x818cf8, css: '#818cf8', impact: 'normal'   },
    "5'UTR"             : { hex: 0xa78bfa, css: '#a78bfa', impact: 'normal'   },
    'Intron'            : { hex: 0x38bdf8, css: '#38bdf8', impact: 'normal'   },
    'Nonsense_Mutation' : { hex: 0xff3366, css: '#ff3366', impact: 'high'     },
    'Frame_Shift_Del'   : { hex: 0xff3366, css: '#ff3366', impact: 'high'     },
    'Frame_Shift_Ins'   : { hex: 0xff6655, css: '#ff6655', impact: 'high'     },
    'Splice_Site'       : { hex: 0xff6644, css: '#ff6644', impact: 'high'     },
    'In_Frame_Del'      : { hex: 0xfbbf24, css: '#fbbf24', impact: 'missense' },
    'In_Frame_Ins'      : { hex: 0xf59e0b, css: '#f59e0b', impact: 'missense' },
};
const DEFAULT_META = { hex: 0x888888, css: '#888888', impact: 'normal' };

const GENE_POOLS = {
    high    : ['TP53','BRCA1','BRCA2','APC','RB1','PTEN','MSH2','MLH1'],
    missense: ['KRAS','EGFR','PIK3CA','BRAF','MET','NRAS','CDKN2A','IDH1'],
    normal  : ['CDH1','ERBB2','MYC','VHL','STK11','SMAD4','FBXW7','ATM','RNF43','ARID1A'],
};
function pickFrom(arr, seed) { return arr[Math.abs(seed) % arr.length]; }

// ── Tooltip DOM ───────────────────────────────────────────────────────────────
let tooltip;
function ensureTooltip() {
    if (tooltip) return;
    tooltip = document.createElement('div');
    Object.assign(tooltip.style, {
        position: 'absolute', pointerEvents: 'none',
        background: 'rgba(10,10,20,0.95)', border: '1px solid #00ff9d',
        borderRadius: '8px', padding: '8px 12px',
        fontFamily: '"Fira Code", monospace', fontSize: '13px',
        color: '#e0ffe0', display: 'none', zIndex: '9999',
        boxShadow: '0 0 12px rgba(0,255,157,0.3)', lineHeight: '1.8',
        minWidth: '190px',
    });
    const container = document.getElementById('dnaCanvasContainer');
    if (container) {
        container.style.position = 'relative';
        container.appendChild(tooltip);
    }
}

// ── Raycaster ─────────────────────────────────────────────────────────────────
const raycaster  = new THREE.Raycaster();
const mouse      = new THREE.Vector2(-9999, -9999);
let hoveredSphere = null;
let labelSprite   = null;

function onMouseMove(event) {
    const container = document.getElementById('dnaCanvasContainer');
    if (!container || !renderer) return;
    const rect = container.getBoundingClientRect();
    mouse.x =  ((event.clientX - rect.left) / rect.width)  * 2 - 1;
    mouse.y = -((event.clientY - rect.top)  / rect.height) * 2 + 1;
    if (tooltip) {
        tooltip.style.left = (event.clientX - rect.left + 16) + 'px';
        tooltip.style.top  = (event.clientY - rect.top  - 12) + 'px';
    }
}

// ── 3D floating label ─────────────────────────────────────────────────────────
function _roundRect(ctx, x, y, w, h, r) {
    ctx.beginPath();
    ctx.moveTo(x+r, y); ctx.lineTo(x+w-r, y);
    ctx.quadraticCurveTo(x+w, y,   x+w, y+r);   ctx.lineTo(x+w, y+h-r);
    ctx.quadraticCurveTo(x+w, y+h, x+w-r, y+h); ctx.lineTo(x+r, y+h);
    ctx.quadraticCurveTo(x,   y+h, x,   y+h-r); ctx.lineTo(x, y+r);
    ctx.quadraticCurveTo(x,   y,   x+r, y);
    ctx.closePath();
}

function make3DLabel(meta) {
    const W = 340, H = 84;
    const canvas = document.createElement('canvas');
    canvas.width = W; canvas.height = H;
    const ctx = canvas.getContext('2d');
    const col = meta.css || '#00ff9d';

    ctx.fillStyle = 'rgba(8,8,18,0.92)';
    _roundRect(ctx, 0, 0, W, H, 12); ctx.fill();
    ctx.strokeStyle = col; ctx.lineWidth = 2;
    _roundRect(ctx, 1, 1, W-2, H-2, 11); ctx.stroke();

    ctx.fillStyle = '#ffffff';
    ctx.font = 'bold 22px "Fira Code",monospace';
    ctx.fillText(meta.gene, 14, 30);

    ctx.fillStyle = col;
    ctx.font = '14px "Fira Code",monospace';
    ctx.fillText(meta.variant, 14, 52);

    if (meta.prob != null) {
        ctx.fillStyle = '#aaa';
        ctx.font = '13px "Fira Code",monospace';
        ctx.fillText((meta.prob * 100).toFixed(1) + '%  ' + meta.impact.toUpperCase(), 14, 72);
    }

    const tex = new THREE.CanvasTexture(canvas);
    const spr = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, depthTest: false, transparent: true }));
    spr.scale.set(1.7, 0.42, 1);
    return spr;
}

// ── Hover check ───────────────────────────────────────────────────────────────
function checkHover() {
    if (!camera || spheres.length === 0) return;
    raycaster.setFromCamera(mouse, camera);
    const hits = raycaster.intersectObjects(spheres);

    if (hits.length > 0) {
        const hit = hits[0].object;
        if (hit !== hoveredSphere) {
            if (hoveredSphere) hoveredSphere.scale.setScalar(1.0);
            if (labelSprite)   { scene.remove(labelSprite); labelSprite = null; }

            hoveredSphere = hit;
            const idx  = spheres.indexOf(hit);
            const meta = nodeMetadata[idx];

            if (meta && meta.gene !== '—') {
                // 3D label above node
                labelSprite = make3DLabel(meta);
                labelSprite.position.copy(hit.position);
                labelSprite.position.y += 0.58;
                scene.add(labelSprite);

                // Pulse
                hit.scale.setScalar(1.45);

                // DOM tooltip
                if (tooltip) {
                    const col = meta.css || '#00ff9d';
                    const pct = meta.prob != null ? (meta.prob * 100).toFixed(1) + '%' : '—';
                    tooltip.innerHTML =
                        `<span style="color:${col}">⬤</span> <strong style="color:#fff">${meta.gene}</strong><br>` +
                        `<span style="color:#888">Variant :</span> <span style="color:${col}">${meta.variant}</span><br>` +
                        `<span style="color:#888">Prob    :</span> <span style="color:${col}">${pct}</span><br>` +
                        `<span style="color:#888">Impact  :</span> <span style="color:${col}">${meta.impact.toUpperCase()}</span>`;
                    tooltip.style.display = 'block';
                }
            }
        }
    } else {
        if (hoveredSphere) { hoveredSphere.scale.setScalar(1.0); hoveredSphere = null; }
        if (labelSprite)   { scene.remove(labelSprite); labelSprite = null; }
        if (tooltip)       tooltip.style.display = 'none';
    }
}

// ── Init ──────────────────────────────────────────────────────────────────────
function initDNA3D() {
    const container = document.getElementById('dnaCanvasContainer');
    if (!container) return;
    ensureTooltip();
    const width = container.clientWidth, height = 400;

    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0a0a0f);

    camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
    camera.position.set(5, 5, 8);
    camera.lookAt(0, 0, 0);

    renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'high-performance' });
    renderer.setSize(width, height);
    renderer.shadowMap.enabled = false;
    container.innerHTML = '';
    container.appendChild(renderer.domElement);
    ensureTooltip();

    renderer.domElement.addEventListener('mousemove', onMouseMove);

    controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    controls.autoRotate    = false;

    scene.add(new THREE.AmbientLight(0x404060));
    const mainLight = new THREE.DirectionalLight(0xffffff, 0.8);
    mainLight.position.set(1, 2, 1);
    scene.add(mainLight);

    dnaGroup = new THREE.Group();
    scene.add(dnaGroup);

    createStaticDNA();
    showPlaceholder();
    animate();
    window.addEventListener('resize', onWindowResize);
}

// ── DNA structure ─────────────────────────────────────────────────────────────
function createStaticDNA() {
    spheres = []; cylinders = []; nodeMetadata = [];

    const radius = 1.2, turns = 2.5, segments = 24, heightTotal = 6;
    const stepAngle = (Math.PI * 2 * turns) / segments;

    for (let i = 0; i <= segments; i++) {
        const t = i / segments;
        const y = -heightTotal / 2 + t * heightTotal;
        const angle = stepAngle * i;
        const x1 = radius * Math.cos(angle),           z1 = radius * Math.sin(angle);
        const x2 = radius * Math.cos(angle + Math.PI), z2 = radius * Math.sin(angle + Math.PI);
        const geom = new THREE.SphereGeometry(0.18, 16, 16);

        [{ x: x1, z: z1 }, { x: x2, z: z2 }].forEach(pos => {
            const mat = new THREE.MeshStandardMaterial({ color: 0x888888, emissive: 0x000000 });
            const s   = new THREE.Mesh(geom, mat);
            s.position.set(pos.x, y, pos.z);
            dnaGroup.add(s);
            spheres.push(s);
            nodeMetadata.push({ gene: '—', variant: '—', prob: null, impact: 'normal', css: '#888888' });
        });

        // Rung
        const start = new THREE.Vector3(x1, y, z1);
        const end   = new THREE.Vector3(x2, y, z2);
        const dir   = new THREE.Vector3().subVectors(end, start);
        const cyl   = new THREE.Mesh(
            new THREE.CylinderGeometry(0.05, 0.05, dir.length(), 6),
            new THREE.MeshStandardMaterial({ color: 0x88aaff })
        );
        cyl.position.copy(start.clone().add(end).multiplyScalar(0.5));
        cyl.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize());
        dnaGroup.add(cyl);
        cylinders.push(cyl);
    }

    const pts1 = [], pts2 = [];
    for (let i = 0; i <= 100; i++) {
        const t = i / 100, y = -heightTotal / 2 + t * heightTotal;
        const a = stepAngle * (i / (100 / segments));
        pts1.push(new THREE.Vector3(radius * Math.cos(a),           y, radius * Math.sin(a)));
        pts2.push(new THREE.Vector3(radius * Math.cos(a + Math.PI), y, radius * Math.sin(a + Math.PI)));
    }
    const lineMat = new THREE.LineBasicMaterial({ color: 0x33cc99 });
    dnaGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts1), lineMat));
    dnaGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts2), lineMat));
}

// ── Placeholder ───────────────────────────────────────────────────────────────
function showPlaceholder() {
    const existing = dnaGroup.children.find(c => c.isSprite);
    if (existing) dnaGroup.remove(existing);
    const canvas = document.createElement('canvas');
    canvas.width = 512; canvas.height = 128;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#00ff9d';
    ctx.font = 'bold 22px "Fira Code",monospace';
    ctx.fillText('⚡ Click any prediction log', 20, 55);
    ctx.font = '15px "Fira Code",monospace';
    ctx.fillStyle = '#666';
    ctx.fillText('to visualize mutation variants', 20, 95);
    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(canvas) }));
    sprite.scale.set(3, 0.8, 1);
    dnaGroup.add(sprite);
}

// ── Color nodes from actual variant probabilities ─────────────────────────────
function updateDNAColors(data) {
    const sprite = dnaGroup.children.find(c => c.isSprite);
    if (sprite) dnaGroup.remove(sprite);

    const variantProbs = data.variantProbs || {};
    const n = spheres.length;

    const entries = Object.entries(variantProbs)
        .filter(([, p]) => p > 0)
        .sort((a, b) => b[1] - a[1]);

    const totalWeight = entries.reduce((s, [, p]) => s + p, 0) || 1;

    const assignments = [];
    entries.forEach(([name, prob]) => {
        const vm    = VARIANT_META[name] || DEFAULT_META;
        const count = Math.max(1, Math.round((prob / totalWeight) * n));
        for (let k = 0; k < count && assignments.length < n; k++) {
            assignments.push({ name, prob, vm });
        }
    });
    while (assignments.length < n) {
        const last = assignments[assignments.length - 1] || { name: 'Unknown', prob: 0, vm: DEFAULT_META };
        assignments.push(last);
    }

    for (let i = 0; i < n; i++) {
        const { name, prob, vm } = assignments[i];
        spheres[i].material.color.setHex(vm.hex);
        spheres[i].material.emissive.setHex(vm.hex);
        spheres[i].material.emissiveIntensity = 0.45;

        nodeMetadata[i] = {
            gene   : pickFrom(GENE_POOLS[vm.impact], i + (data.totalMutations || 0)),
            variant: name,
            prob,
            impact : vm.impact,
            css    : vm.css,
        };
    }
}

function updateDNAFromData(patientData) {
    currentData = patientData;
    updateDNAColors(patientData);
}

// ── Render loop ───────────────────────────────────────────────────────────────
function animate() {
    requestAnimationFrame(animate);
    controls.update();
    checkHover();
    renderer.render(scene, camera);
}

function onWindowResize() {
    const container = document.getElementById('dnaCanvasContainer');
    if (!container) return;
    const width = container.clientWidth;
    camera.aspect = width / 400;
    camera.updateProjectionMatrix();
    renderer.setSize(width, 400);
}

window.updateDNAFromData = updateDNAFromData;
window.initDNA3D = initDNA3D;

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initDNA3D);
} else {
    initDNA3D();
}