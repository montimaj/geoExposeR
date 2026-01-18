# ============================================================================
# Generate Automated Package Dependency/Workflow Graph for gwArsenicR
# ============================================================================
# This script creates a workflow diagram showing how different packages
# integrate within the gwArsenicR analysis pipeline.
#
# Required packages: DiagrammeR, DiagrammeRsvg, rsvg
# ============================================================================

# Install required packages if not available
required_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ============================================================================
# JOSS Publication Figure - Narrow Vertical Layout for Single-Column PDF
# ============================================================================
# Optimized for JOSS single-column format (~6.5 inch width, narrow aspect ratio)

joss_workflow <- grViz("
digraph gwArsenicR_joss {
  
  # Graph attributes - narrow vertical layout
  graph [
    rankdir = TB,
    fontname = 'Arial',
    fontsize = 10,
    nodesep = 0.2,
    ranksep = 0.4,
    bgcolor = 'white',
    compound = true,
    newrank = true,
    splines = ortho
  ]
  
  # Node defaults - compact
  node [
    fontname = 'Arial',
    fontsize = 8,
    style = 'filled,rounded',
    shape = 'box',
    margin = '0.08,0.04',
    penwidth = 1.0,
    width = 1.0,
    height = 0.4
  ]
  
  # Edge defaults
  edge [
    fontname = 'Arial',
    fontsize = 7,
    color = '#546E7A',
    penwidth = 1.0,
    arrowsize = 0.5
  ]
  
  # ==== INPUT DATA (2x2 grid) ====
  subgraph cluster_input {
    label = 'Input Data'
    labelloc = 't'
    fontname = 'Arial Bold'
    fontsize = 9
    fontcolor = '#1565C0'
    style = 'rounded,filled'
    color = '#1976D2'
    fillcolor = '#E3F2FD'
    margin = 8
    
    usgs [label = 'USGS Arsenic\\nRasters', fillcolor = 'white', color = '#1976D2']
    epa [label = 'EPA Water\\nSystem Data', fillcolor = 'white', color = '#1976D2']
    wells [label = 'Private Well\\nData', fillcolor = 'white', color = '#1976D2']
    health [label = 'Health\\nOutcomes', fillcolor = 'white', color = '#1976D2']
    
    {rank = same; usgs; epa}
    {rank = same; wells; health}
  }
  
  # ==== DATA INGESTION ====
  subgraph cluster_ingest {
    label = 'Data Ingestion  [terra, sf, exactextractr]'
    labelloc = 't'
    fontname = 'Arial Bold'
    fontsize = 9
    fontcolor = '#2E7D32'
    style = 'rounded,filled'
    color = '#388E3C'
    fillcolor = '#E8F5E9'
    margin = 8
    
    ingest [label = 'prepare_usgs_data()\\nprepare_epa_data()\\nadd_well_percentages()\\nformat_health_data()', fillcolor = 'white', color = '#388E3C', width = 2.2, height = 0.7]
  }
  
  # ==== PROBABILITY INTEGRATION ====
  subgraph cluster_integrate {
    label = 'Probability Integration  [data.table]'
    labelloc = 't'
    fontname = 'Arial Bold'
    fontsize = 9
    fontcolor = '#E65100'
    style = 'rounded,filled'
    color = '#F57C00'
    fillcolor = '#FFF3E0'
    margin = 8
    
    convert [label = 'convert_epa_to_multinomial()\\ncreate_weighted_prob_matrix()', fillcolor = 'white', color = '#F57C00', width = 2.2, height = 0.5]
  }
  
  # ==== MULTIPLE IMPUTATION ====
  subgraph cluster_impute {
    label = 'Multiple Imputation  [Amelia, mice]'
    labelloc = 't'
    fontname = 'Arial Bold'
    fontsize = 9
    fontcolor = '#6A1B9A'
    style = 'rounded,filled'
    color = '#7B1FA2'
    fillcolor = '#F3E5F5'
    margin = 8
    
    impute [label = 'impute_arsenic_exposure()\\nimpute_additional_variables()', fillcolor = 'white', color = '#7B1FA2', width = 2.2, height = 0.5]
  }
  
  # ==== STATISTICAL MODELING ====
  subgraph cluster_model {
    label = 'Statistical Modeling  [lme4, broom.mixed]'
    labelloc = 't'
    fontname = 'Arial Bold'
    fontsize = 9
    fontcolor = '#B71C1C'
    style = 'rounded,filled'
    color = '#C62828'
    fillcolor = '#FFEBEE'
    margin = 8
    
    model [label = 'regression_analysis()\\npool_estimates_by_term()', fillcolor = 'white', color = '#C62828', width = 2.2, height = 0.5]
  }
  
  # ==== OUTPUT ====
  subgraph cluster_output {
    label = 'Outputs'
    labelloc = 't'
    labeljust = 'l'
    fontname = 'Arial Bold'
    fontsize = 9
    fontcolor = '#37474F'
    style = 'rounded,filled'
    color = '#455A64'
    fillcolor = '#ECEFF1'
    margin = 8
    
    output [label = 'Imputed Datasets\\nPooled Estimates (beta, SE, CI, p)\\nSensitivity Analysis Results', fillcolor = 'white', color = '#455A64', width = 2.4, height = 0.6]
  }
  
  # ==== EDGES - Simple vertical flow ====
  usgs -> ingest [tailport = s, headport = nw]
  epa -> ingest [tailport = s, headport = n]
  wells -> ingest [tailport = s, headport = n]
  health -> ingest [tailport = s, headport = ne]
  
  ingest -> convert
  convert -> impute
  impute -> model
  model -> output
}
")

# ============================================================================
# Save Figures
# ============================================================================

# Export SVG from DiagrammeR
svg_joss <- export_svg(joss_workflow)
svg_content <- as.character(svg_joss)

# Save as PNG - 600 DPI for JOSS (3.5 inch width * 600 = 2100 px, 5.5 inch height * 600 = 3300 px)
rsvg_png(charToRaw(svg_content), "workflow.png", width = 2100, height = 3300)
cat("Saved: workflow.png (600 DPI)\n")

cat("\n✓ All figures generated successfully!\n")