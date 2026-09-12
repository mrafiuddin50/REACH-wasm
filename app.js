// Setup MapLibre base map
const map = new maplibregl.Map({
    container: 'map',
    style: 'https://basemaps.cartocdn.com/gl/dark-matter-gl-style/style.json',
    center: [78.9629, 20.5937],
    zoom: 4,
    pitch: 0
});

let deckLayer = null;
let geojsonData = null;
let centroids = null;
let baselinePM25 = null;
let currentDisplay = null;

const poll_info = ["SO2", "NOX", "NH3", "VOC", "OC", "PM25-PRI", "BC"];
let sr_matrices = {};
let baseline_emissions = {};
let sector_emissions = {};

const sectors = ["EGU", "IND", "AWB", "RCO", "TRO"];

// Initialize DeckGL
const deckOverlay = new deck.MapboxOverlay({
    layers: []
});
map.addControl(deckOverlay);

async function loadFloat32Array(url) {
    const res = await fetch(url);
    const buffer = await res.arrayBuffer();
    return new Float32Array(buffer);
}

async function init() {
    try {
        document.getElementById('status').innerText = 'Loading spatial data...';
        const [geoRes, centRes] = await Promise.all([
            fetch('data/south_asia.geojson'),
            fetch('data/centroids.json')
        ]);
        
        geojsonData = await geoRes.json();
        centroids = await centRes.json();
        
        document.getElementById('status').innerText = 'Downloading Matrices (30MB)...';
        
        // Load all matrices and baseline emissions concurrently
        const loadPromises = [];
        for (let poll of poll_info) {
            loadPromises.push(
                loadFloat32Array(`data/sr_${poll}.bin`).then(arr => sr_matrices[poll] = arr)
            );
            loadPromises.push(
                loadFloat32Array(`data/baseline_${poll}.bin`).then(arr => baseline_emissions[poll] = arr)
            );
            
            sector_emissions[poll] = {};
            for (let sec of sectors) {
                loadPromises.push(
                    loadFloat32Array(`data/${sec}_${poll}.bin`).then(arr => sector_emissions[poll][sec] = arr)
                );
            }
        }
        
        await Promise.all(loadPromises);
        
        // Calculate true baseline once
        baselinePM25 = calculate("NONE", 0);
        
        document.getElementById('status').innerText = 'Ready';
        document.getElementById('status').classList.add('ready');
        
        runScenario("NONE", 0);
        
        // Setup control events
        const slider = document.getElementById('awb-slider');
        const valText = document.getElementById('awb-value');
        const sectorSelect = document.getElementById('sector-select');
        
        slider.addEventListener('input', (e) => {
            valText.innerText = e.target.value + '%';
            runScenario(sectorSelect.value, parseFloat(e.target.value));
        });
        
        sectorSelect.addEventListener('change', (e) => {
            runScenario(sectorSelect.value, parseFloat(slider.value));
        });
        
    } catch (e) {
        document.getElementById('status').innerText = 'Error: ' + e.message;
        console.error(e);
    }
}

function calculate(sector, reduction_percent) {
    const factor = reduction_percent / 100.0;
    const N = baseline_emissions["SO2"].length / 2;
    const tot_conc = [];
    
    // Initialize N x 7 array
    for (let i = 0; i < N; i++) {
        tot_conc.push(new Float32Array(7));
    }
    
    // 1. Matrix Multiplication
    for (let j = 0; j < poll_info.length; j++) {
        const poll = poll_info[j];
        const base_emis = baseline_emissions[poll];
        const sr = sr_matrices[poll]; // Flattened (2N x N) matrix where rows are sources, cols are receptors
        
        // Compute active emission vector
        const active_emis = new Float32Array(base_emis.length);
        if (sector === "NONE") {
            for (let i = 0; i < active_emis.length; i++) active_emis[i] = base_emis[i];
        } else {
            const sec_emis = sector_emissions[poll][sector];
            for (let i = 0; i < active_emis.length; i++) {
                active_emis[i] = base_emis[i] - (sec_emis[i] * factor);
            }
        }
        
        // SR matrix multiplication (Receptors x Sources) @ (Sources x 1)
        // Since sr_matrices were saved from Pandas DataFrame, they were row-major.
        // Wait, Python's pd.concat([elevated, ground]) is shape (2272, 1136).
        // 2272 sources, 1136 receptors.
        // C-contiguous means element (s, r) is at index s * N + r.
        // Output at receptor r is SUM_{s=0}^{2N-1} (sr[s, r] * emis[s])
        for (let s = 0; s < 2 * N; s++) {
            const emis_val = active_emis[s];
            if (emis_val === 0) continue; // Optimization
            
            const offset = s * N;
            for (let r = 0; r < N; r++) {
                tot_conc[r][j] += sr[offset + r] * emis_val;
            }
        }
    }
    
    // 2. Chemistry Partitioning
    const pm25_total = new Float32Array(N);
    
    for (let i = 0; i < N; i++) {
        tot_conc[i][4] /= 1.8;
        
        let SO4 = tot_conc[i][0] * 0.5;
        let Tot_HNO3 = tot_conc[i][1] * 0.35;
        let Tot_NH3 = tot_conc[i][2] * 1.0;
        let OC_secondary = tot_conc[i][3] * 0.024;
        let OC_primary = tot_conc[i][4] * 1.0;
        let PM = tot_conc[i][5] * 1.0;
        let BC = tot_conc[i][6] * 1.0;
        
        Tot_HNO3 = 0.98 * Tot_HNO3;
        Tot_NH3 = 1.06 * Tot_NH3;
        
        let Tot_NH3_mol = Tot_NH3 / 18.0;
        let SO4_mol = SO4 / 96.0;
        let Tot_HNO3_mol = Tot_HNO3 / 62.0;
        
        let NH3f_mol = Tot_NH3_mol - (1.5 * SO4_mol);
        if (NH3f_mol <= 0) NH3f_mol = 1E-20;
        
        let NO3_mol = 0.6509 * ((0.33873 * Tot_HNO3_mol) + (0.121008 * NH3f_mol) + (3.511482 * NH3f_mol * Tot_HNO3_mol));
        if (NO3_mol > Tot_HNO3_mol) NO3_mol = Tot_HNO3_mol;
        
        let NH4_mol = (2 * SO4_mol) + NO3_mol;
        if (NH4_mol > Tot_NH3_mol) NH4_mol = Tot_NH3_mol;
        
        let NO3 = NO3_mol * 62.0;
        let NH4 = NH4_mol * 18.0;
        let SOA = OC_secondary * 1.8;
        let H2SO4 = SO4 * 98.0 / 96.0;
        
        pm25_total[i] = H2SO4 + NO3 + NH4 + SOA + PM + BC;
    }
    
    return pm25_total;
}

function runScenario(sector, reductionPercent) {
    const scenarioPM25 = calculate(sector, reductionPercent);
    
    const isBaseline = (sector === "NONE" || reductionPercent === 0);
    
    currentDisplay = new Float32Array(baselinePM25.length);
    for (let i = 0; i < baselinePM25.length; i++) {
        currentDisplay[i] = isBaseline ? baselinePM25[i] : (baselinePM25[i] - scenarioPM25[i]);
    }
    
    updateMap(isBaseline);
    updateStats(isBaseline);
}

function updateMap(isBaseline) {
    const censusIdToIndex = {};
    centroids.forEach((c, idx) => {
        censusIdToIndex[c.census_id] = idx;
    });
    
    deckLayer = new deck.GeoJsonLayer({
        id: 'pm25-layer',
        data: geojsonData,
        pickable: true,
        stroked: true,
        filled: true,
        extruded: false,
        wireframe: false,
        getLineColor: [255, 255, 255, 50],
        getLineWidth: 1,
        lineWidthMinPixels: 1,
        getElevation: d => 0,
        getFillColor: d => {
            const idx = censusIdToIndex[d.properties.census_id];
            const val = idx !== undefined ? currentDisplay[idx] : 0;
            return isBaseline ? colorScaleAbsolute(val) : colorScaleReduction(val);
        },
        updateTriggers: {
            getFillColor: [currentDisplay, isBaseline]
        },
        onHover: info => {
            const tooltip = document.getElementById('tooltip');
            if (info.object) {
                const idx = censusIdToIndex[info.object.properties.census_id];
                const val = idx !== undefined ? currentDisplay[idx] : 0;
                tooltip.style.display = 'block';
                tooltip.style.left = info.x + 10 + 'px';
                tooltip.style.top = info.y + 10 + 'px';
                const label = isBaseline ? "PM2.5" : "Reduction";
                tooltip.innerHTML = `District: ${info.object.properties.census_id}<br/>${label}: ${val.toFixed(2)} &mu;g/m&sup3;`;
            } else {
                tooltip.style.display = 'none';
            }
        }
    });
    
    deckOverlay.setProps({
        layers: [deckLayer]
    });
}

function updateStats(isBaseline) {
    let sum = 0;
    for(let i = 0; i < currentDisplay.length; i++) {
        sum += currentDisplay[i];
    }
    const avg = sum / currentDisplay.length;
    
    document.getElementById('stat-label-text').innerText = isBaseline ? "Avg PM2.5" : "Avg Reduction";
    document.getElementById('avg-reduction').innerText = avg.toFixed(1);
    
    // Update legend
    const legend = document.getElementById('legend');
    if (isBaseline) {
        legend.innerHTML = `
        <h4>PM2.5 (&mu;g/m&sup3;)</h4>
        <div class="legend-row"><span class="color-box" style="background: rgba(0,228,0,0.8)"></span> &lt; 15</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(255,255,0,0.8)"></span> 15 - 30</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(255,126,0,0.8)"></span> 30 - 55</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(255,0,0,0.8)"></span> 55 - 150</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(143,63,151,0.8)"></span> 150 - 250</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(126,0,35,0.8)"></span> &gt; 250</div>`;
    } else {
        legend.innerHTML = `
        <h4>Reduction (&mu;g/m&sup3;)</h4>
        <div class="legend-row"><span class="color-box" style="background: rgba(255,255,178,0.8)"></span> 0.1 - 1</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(254,204,92,0.8)"></span> 1 - 5</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(253,141,60,0.8)"></span> 5 - 10</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(240,59,32,0.8)"></span> 10 - 20</div>
        <div class="legend-row"><span class="color-box" style="background: rgba(189,0,38,0.8)"></span> &gt; 20</div>`;
    }
}

function colorScaleReduction(val) {
    if (val < 0.1) return [0, 0, 0, 0];
    if (val < 1) return [255, 255, 178, 200];
    if (val < 5) return [254, 204, 92, 200];
    if (val < 10) return [253, 141, 60, 200];
    if (val < 20) return [240, 59, 32, 200];
    return [189, 0, 38, 200];
}

function colorScaleAbsolute(val) {
    if (val < 15) return [0, 228, 0, 200];
    if (val < 30) return [255, 255, 0, 200];
    if (val < 55) return [255, 126, 0, 200];
    if (val < 150) return [255, 0, 0, 200];
    if (val < 250) return [143, 63, 151, 200];
    return [126, 0, 35, 200];
}

// Add tooltip element dynamically
const tt = document.createElement('div');
tt.id = 'tooltip';
document.body.appendChild(tt);

init();
