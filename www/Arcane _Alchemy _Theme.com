/* === Arcane Alchemy Theme === */
body {
  background-color: #0d0d1a; /* Deep dark blue-violet background */
  color: #d9e3ff; /* Electric blue text */
  font-family: 'Consolas', 'Courier New', monospace;
}

h1, h2, h3, h4 {
  color: #a4c8ff; /* Bright blue for headers */
  text-shadow: 0 0 6px #7aaeff;
}

a, .btn, .shiny-download-link {
  background-color: #1e1e3f;
  color: #ffffff;
  border: 1px solid #4da6ff;
  text-shadow: 0 0 5px #4da6ff;
}

a:hover, .btn:hover {
  background-color: #003366;
  color: #ffffff;
}

.shiny-input-container {
  color: #d9e3ff;
}

input, select, .form-control {
  background-color: #1b1b2e;
  color: #e0eaff;
  border: 1px solid #4da6ff;
}

.nav-tabs > li > a {
  background-color: #1b1b2e;
  color: #99ccff;
  border: 1px solid #4da6ff;
}

.nav-tabs > li.active > a,
.nav-tabs > li.active > a:focus,
.nav-tabs > li.active > a:hover {
  background-color: #0f0f22;
  color: #ffffff;
  border: 1px solid #7aaeff;
  box-shadow: 0 0 6px #4da6ff;
}

.table {
  background-color: #1a1a2e;
  color: #d0eaff;
  border-color: #4da6ff;
}

.table thead {
  background-color: #0f0f22;
  color: #7aaeff;
}

#heatmap_plot, #volcano_plot, #pcaplot, #umapplot, #rf_roc_plot {
  background-color: #0f0f22;
  border: 1px solid #4da6ff;
  box-shadow: 0 0 10px #4da6ff;
}

/* Notification popup glow */
.shiny-notification {
  background-color: #000019;
  color: #d9e3ff;
  border-left: 5px solid #4da6ff;
  box-shadow: 0 0 12px #4da6ff;
}
