/* docs/script.js */

// Tab Navigation Logic
document.querySelectorAll('.nav-item').forEach(item => {
  item.addEventListener('click', function() {
    // Remove active class from all tabs
    document.querySelectorAll('.nav-item').forEach(nav => nav.classList.remove('active'));
    document.querySelectorAll('.tab-pane').forEach(pane => pane.classList.remove('active'));
    
    // Add active class to clicked tab
    this.classList.add('active');
    const targetId = this.getAttribute('data-target');
    document.getElementById(targetId).classList.add('active');
  });
});

// Helper function to render a DataTable from a CSV
function loadTable(csvUrl, tableId) {
  Papa.parse(csvUrl, {
    download: true,
    header: true,
    skipEmptyLines: true,
    complete: function(results) {
      if (results.data.length === 0) return;
      
      const columns = Object.keys(results.data[0]).map(key => ({
        title: key,
        data: key,
        defaultContent: ""
      }));
      
      $(`#${tableId}`).DataTable({
        data: results.data,
        columns: columns,
        pageLength: 10,
        scrollX: true,
        destroy: true
      });
    }
  });
}

// Load Tables
document.addEventListener("DOMContentLoaded", function() {
  loadTable('data/STable2_17_FDR_hits_complete.csv', 'table-fdr');
  loadTable('data/STable1_full_MR_screen.csv', 'table-full');
  loadTable('data/STable4_mediation_integrated_evidence.csv', 'table-mediation');
  
  // Render Plotly Graphs
  renderVolcanoPlot();
  renderForestPlot();
});

// Render Volcano Plot from STable1_full_MR_screen.csv
function renderVolcanoPlot() {
  Papa.parse('data/STable1_full_MR_screen.csv', {
    download: true,
    header: true,
    dynamicTyping: true,
    skipEmptyLines: true,
    complete: function(results) {
      const data = results.data;
      
      // Separate into significant (FDR < 0.05) and non-significant
      const sig = data.filter(d => d.FDR < 0.05);
      const nonsig = data.filter(d => d.FDR >= 0.05 || d.FDR === null);
      
      const traceSig = {
        x: sig.map(d => d.OR ? Math.log(d.OR) : null),
        y: sig.map(d => d.pvalue ? -Math.log10(d.pvalue) : null),
        mode: 'markers',
        type: 'scatter',
        name: 'FDR < 0.05',
        marker: { color: '#D55E00', size: 8, opacity: 0.8, line: {width:0.5, color:'white'} },
        text: sig.map(d => `${d.protein} - ${d.cancer_outcome}<br>FDR: ${d.FDR}`),
        hoverinfo: 'text'
      };
      
      const traceNonSig = {
        x: nonsig.map(d => d.OR ? Math.log(d.OR) : null),
        y: nonsig.map(d => d.pvalue ? -Math.log10(d.pvalue) : null),
        mode: 'markers',
        type: 'scatter',
        name: 'Non-significant',
        marker: { color: '#bbbbbb', size: 8, opacity: 0.6 },
        text: nonsig.map(d => `${d.protein} - ${d.cancer_outcome}`),
        hoverinfo: 'text'
      };
      
      const layout = {
        xaxis: { title: 'ln(OR) = β', zeroline: true },
        yaxis: { title: '−log₁₀(p)', zeroline: false },
        margin: { l: 50, r: 20, t: 30, b: 50 },
        paper_bgcolor: 'white',
        plot_bgcolor: 'white'
      };
      
      Plotly.newPlot('plot-volcano', [traceNonSig, traceSig], layout, {responsive: true});
    }
  });
}

// Render Forest Plot from STable_master_evidence.csv
function renderForestPlot() {
  Papa.parse('data/STable_master_evidence.csv', {
    download: true,
    header: true,
    dynamicTyping: true,
    skipEmptyLines: true,
    complete: function(results) {
      let data = results.data.filter(d => d.mr_or !== null && d.mr_or !== undefined);
      
      // Sort by tier_short, cancer_mr, mr_or
      data.sort((a, b) => {
        if (a.tier_short !== b.tier_short) return String(a.tier_short).localeCompare(String(b.tier_short));
        if (a.cancer_mr !== b.cancer_mr) return String(a.cancer_mr).localeCompare(String(b.cancer_mr));
        return a.mr_or - b.mr_or;
      });
      
      const labels = data.map(d => `${d.protein} (${d.tier_short} — ${d.cancer_mr})`);
      const ors = data.map(d => d.mr_or);
      const or_lo = data.map(d => d.mr_or_lo);
      const or_hi = data.map(d => d.mr_or_hi);
      
      const trace = {
        x: ors,
        y: labels,
        mode: 'markers',
        type: 'scatter',
        marker: { color: '#D55E00', size: 10 },
        error_x: {
          type: 'data',
          symmetric: false,
          array: or_hi.map((hi, i) => hi - ors[i]),
          arrayminus: ors.map((or, i) => or - or_lo[i]),
          color: '#333'
        },
        text: data.map(d => `<b>${d.protein}</b><br>Cancer: ${d.cancer_mr}<br>Tier: ${d.tier_short}<br>OR: ${d.mr_or.toFixed(3)}`),
        hoverinfo: 'text'
      };
      
      const layout = {
        xaxis: { title: 'Odds ratio (95% CI)' },
        yaxis: { title: '', automargin: true },
        margin: { l: 250, r: 20, t: 30, b: 50 },
        height: Math.max(500, labels.length * 25),
        paper_bgcolor: 'white',
        plot_bgcolor: 'white'
      };
      
      Plotly.newPlot('plot-forest', [trace], layout, {responsive: true});
    }
  });
}
