# ============================================================================
# Generate the geoExposeR Workflow Figure (v1.1.0) for the SoftwareX Paper
# ============================================================================
# This script draws the geoExposeR v1.1.0 analysis pipeline and the R
# packages used at each stage. Every geoExposeR function named in the figure
# is defined in R/ (checked below), and every package listed is called in
# that stage.
#
# Outputs:
#   man/figures/workflow.png  README/pkgdown figure, 1800 px wide, which draws
#                             the model_type x model_family engine grid itself.
#   paper/workflow.pdf        Paper figure, vector drawing (SoftwareX: EPS/PDF
#                             with embedded fonts or text as graphics).
#   paper/workflow.png        Paper figure, bitmapped line drawing, 3600 px wide
#                             at 1000 dpi (SoftwareX minimum: 1000 dpi, 3543 px
#                             for a single column). The paper lists the
#                             regression engines in its Table 1.
# The paper/ manuscript folder is not part of the repository, so the paper
# figures are written only when a local paper/ folder exists.
#
# Usage (from the package root): Rscript scripts/generate_workflow_figure.R
# Required packages: DiagrammeR, DiagrammeRsvg, rsvg, png
# ============================================================================

# Install required packages if not available
required_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg", "png")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Package root (the directory with DESCRIPTION), whether the script is run
# from the root or from scripts/
root <- if (file.exists("DESCRIPTION")) "." else ".."
stopifnot(file.exists(file.path(root, "DESCRIPTION")))
paper_dir <- file.path(root, "paper")

# ============================================================================
# SoftwareX Publication Figure - Portrait Layout for a Single Column
# ============================================================================
# The natural (PDF) width is about 283 pt. Printed at 90 mm the 8 pt body
# text is about 7.2 pt and the smallest text (7 pt) about 6.3 pt; at 140 mm
# they are about 11 pt and 9.8 pt. Colors come from the
# Okabe-Ito palette (distinguishable with color-vision deficiency), and every
# stage is labeled in text, so color is never the only cue. The optional
# data-ingestion stage uses a dashed outline and dashed arrows.
#
# Note: DiagrammeR::grViz() replaces every single quote with a double quote,
# so labels use HTML entities (e.g. &#8217; for an apostrophe) instead.

text_font <- "Arial"
# Courier New sets the layout metrics in viz.js (0.6 em monospace advance);
# the SVG is re-mapped to Menlo, which has the same advance but heavier strokes.
mono_font <- "Courier New"
mono_render <- "Menlo, DejaVu Sans Mono, Courier New, monospace"

# viz.js (Graphviz 2.40) truncates POINT-SIZE to whole points
fs_title <- 9 # stage titles and the main-function label
fs_body <- 8 # descriptions
fs_code <- 7 # function/argument names (Menlo looks larger than Arial)
fs_small <- 7 # package lists and secondary lines

# Okabe-Ito palette: border color and a light tint for the stage background
pal <- list(
  ingest = c(border = "#009E73", fill = "#E0F3EC"), # bluish green
  files = c(border = "#56B4E9", fill = "#E3F2FB"), # sky blue
  integrate = c(border = "#E69F00", fill = "#FDF0D5"), # orange
  impute = c(border = "#CC79A7", fill = "#F7E6EF"), # reddish purple
  model = c(border = "#0072B2", fill = "#DCEAF4"), # blue
  output = c(border = "#6E6E6E", fill = "#ECECEC"), # grey
  main = c(border = "#3C3C3C", fill = "#FAFAFA") # neutral
)

# ---- HTML-label helpers ----------------------------------------------------
nbsp <- "&#160;"
code <- function(x, size = fs_code, bold = FALSE) {
  x <- if (bold) paste0("<B>", x, "</B>") else x
  sprintf('<FONT FACE="%s" POINT-SIZE="%s">%s</FONT>', mono_font, size, x)
}
txt <- function(x, size = fs_body, bold = FALSE) {
  x <- if (bold) paste0("<B>", x, "</B>") else x
  sprintf('<FONT POINT-SIZE="%s">%s</FONT>', size, x)
}
# One row: geoExposeR function (monospace) + short description. Small
# non-breaking spaces pad the text away from the panel edge and between the
# columns, and leave room for glyphs (arrows, Greek letters) that the viz.js
# width estimate under-measures.
pad <- txt(nbsp, fs_small)
fun_row <- function(fun, desc) {
  sprintf(
    paste0('<TR><TD ALIGN="LEFT" VALIGN="BOTTOM">%s%s%s</TD>',
           '<TD ALIGN="LEFT" VALIGN="BOTTOM">%s%s</TD></TR>'),
    pad, code(fun), pad, txt(desc), pad
  )
}
# One full-width row spanning both columns
wide_row <- function(content, align = "CENTER") {
  sprintf('<TR><TD COLSPAN="2" ALIGN="%s">%s</TD></TR>', align, content)
}
# White panel holding the rows of a stage (the tinted stage box is rounded)
panel <- function(rows, border) {
  sprintf(
    paste0(
      '<TABLE BORDER="0" COLOR="%s" BGCOLOR="white" ',
      'CELLBORDER="0" CELLSPACING="2" CELLPADDING="0">%s</TABLE>'
    ),
    border, paste(rows, collapse = "")
  )
}
# A stage: bold title with its package list on the same line, then a panel
stage_label <- function(title, pkgs, rows, border) {
  heading <- txt(title, fs_title, bold = TRUE)
  if (!is.null(pkgs)) {
    pkg_txt <- txt(paste0("[", pkgs, "]"), fs_small)
    # A long package list goes on its own line so it does not set the width
    heading <- if (nchar(pkgs) > 45) {
      paste0(heading, "<BR/>", pkg_txt)
    } else {
      paste0(heading, txt(strrep(nbsp, 2), fs_title), pkg_txt)
    }
  }
  sprintf(
    paste0(
      '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="1">',
      "<TR><TD>%s</TD></TR><TR><TD>%s</TD></TR></TABLE>>"
    ),
    heading, panel(rows, border)
  )
}
# An input file: argument name, description and key columns
file_label <- function(arg, desc, cols) {
  sprintf(
    paste0(
      '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="0">',
      "<TR><TD>%s</TD></TR><TR><TD>%s</TD></TR>",
      "<TR><TD>%s</TD></TR></TABLE>>"
    ),
    code(arg, fs_body, bold = TRUE), txt(desc), cols
  )
}
node_attrs <- function(p, dashed = FALSE) {
  sprintf(
    "style = '%s', fillcolor = '%s', color = '%s'",
    if (dashed) "rounded,filled,dashed" else "rounded,filled",
    p[["fill"]], p[["border"]]
  )
}

# ---- Stage contents (v1.1.0) ------------------------------------------------
# Exported helpers in R/data-ingestion.R (terra, sf, exactextractr, Hmisc)
ingest_label <- stage_label(
  "Data ingestion (optional)",
  "terra, sf, exactextractr, Hmisc",
  c(
    fun_row("prepare_prob_data()", "GeoTIFFs &#8594; 3-category probabilities"),
    fun_row("prepare_conc_data()", "concentrations &#8594; lognormal"),
    fun_row("add_well_percentages()", "join private-well % column"),
    fun_row("format_health_data()", "FIPS-style IDs, range checks"),
    fun_row("subset_by_geography()", "filter by ID prefix, list, regex"),
    fun_row("validate_prepared_inputs()", "check columns, ranges, IDs")
  ),
  pal$ingest[["border"]]
)

# The three files read by perform_sensitivity_analysis(); population weights
# are columns inside the exposure files (pop_well_col, pwell_col)
prob_label <- file_label(
  "model_prob_csv", "Modeled probabilities",
  paste0(code("id_col, prob_cols,", fs_small), "<BR/>",
         code("pop_well_col", fs_small))
)
conc_label <- file_label(
  "conc_data_csv", "Measured concentrations",
  paste0(code("id_col, pwell_col,", fs_small), "<BR/>",
         txt("lognormal or raw values", fs_small))
)
health_label <- file_label(
  "birth_data_txt", "Health outcomes",
  paste0(code("record_id_col,", fs_small), "<BR/>",
         txt("targets, covariates", fs_small))
)

# R/data-loading.R (data.table I/O; stats::plnorm for the lognormal CDF)
integrate_label <- stage_label(
  "Probability integration",
  "data.table, stats",
  c(
    fun_row("validate_data_alignment()", "align rows by id_col"),
    fun_row("convert_conc_to_multinomial()", "lognormal CDF at cutoffs"),
    fun_row("create_weighted_prob_matrix()",
            "w&#183;P_model + (1 &#8722; w)&#183;P_meas")
  ),
  pal$integrate[["border"]]
)

# R/imputation.R (mice::mice / mice::complete for optional covariates)
impute_label <- stage_label(
  "Multiple imputation",
  "mice",
  c(
    fun_row("impute_exposure()", "ndraws exposure draws per unit"),
    fun_row("impute_additional_variables()", "optional MICE (impute_vars)")
  ),
  pal$impute[["border"]]
)

# R/regression.R: engine grid = model_type (rows) x model_family (columns);
# the engines are stats::lm/glm, lme4::lmer/glmer, nnet::multinom and
# mclogit::mblogit, called in fit_single_model() / fit_multinomial_mixed();
# reformulas splits the random-effects terms for mclogit.
# The paper figure omits the grid because the paper lists it in Table 1.
eng <- function(x) code(x, fs_code)
hdr <- function(x) paste0(pad, txt(x, fs_small, bold = TRUE), pad)
head_cell <- function(x, align = "CENTER") {
  sprintf('<TD ALIGN="%s" SIDES="B" BORDER="1">%s</TD>', align, x)
}
body_cell <- function(x, align = "CENTER") {
  sprintf('<TD ALIGN="%s" VALIGN="BOTTOM">%s</TD>', align, x)
}
engine_row <- function(type, engines) {
  paste0(
    "<TR>", body_cell(paste0(pad, txt(type, fs_small)), "LEFT"),
    paste(vapply(engines, function(e) body_cell(eng(e)), ""), collapse = ""),
    "</TR>"
  )
}
engine_grid <- paste0(
  '<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="1">',
  "<TR><TD></TD>",
  '<TD COLSPAN="3">', code("model_family", fs_code, bold = TRUE), "</TD>",
  "</TR><TR>",
  head_cell(paste0(pad, code("model_type", fs_code, bold = TRUE)), "LEFT"),
  head_cell(hdr("Gaussian")), head_cell(hdr("GLM family")),
  head_cell(hdr("multinomial")),
  "</TR>",
  engine_row("fixed", c("lm", "glm", "multinom")),
  engine_row("random | mixed", c("lmer", "glmer", "mblogit")),
  "</TABLE>"
)

make_model_label <- function(with_grid) {
  stage_label(
    "Regression and pooling",
    "stats, dplyr, lme4, nnet, mclogit, reformulas, broom.mixed",
    c(
      if (with_grid) wide_row(engine_grid),
      fun_row("regression_analysis()", "per target and imputed dataset"),
      fun_row("fit_single_model()",
              if (with_grid) "fit with the engine above"
              else "fit with the engine in Table 1"),
      fun_row("fit_multinomial_mixed()", "mixed multinomial (mblogit)"),
      fun_row("pool_estimates_by_term()", "per exposure term (and y.level)"),
      fun_row("pool_single_estimate()", "Rubin&#8217;s rules")
    ),
    pal$model[["border"]]
  )
}
model_label <- make_model_label(with_grid = FALSE)

output_label <- stage_label(
  "Outputs (per target)",
  NULL,
  c(
    wide_row(paste0(txt(paste0("Pooled &#946;, SE, statistic, 95% CI and ",
                               "p-value for each")), pad)),
    wide_row(paste0(txt(paste0("exposure category vs. reference ",
                               "(by y.level if multinomial)")), pad)),
    wide_row(txt(paste0("data frame + ",
                        code("&lt;target&gt;_pooled_results.csv", fs_small)),
                 fs_small))
  ),
  pal$output[["border"]]
)

main_label <- sprintf(
  "<%s>",
  code("perform_sensitivity_analysis()", fs_title, bold = TRUE)
)

# ---- Check that every geoExposeR function shown is defined in R/ ------------
r_dir <- file.path(root, "R")
if (dir.exists(r_dir)) {
  labels <- c(ingest_label, integrate_label, impute_label, model_label,
              main_label)
  shown <- unique(unlist(regmatches(
    labels, gregexpr("[A-Za-z_.][A-Za-z0-9_.]*(?=\\(\\))", labels, perl = TRUE)
  )))
  src <- unlist(lapply(list.files(r_dir, "\\.R$", full.names = TRUE),
                       readLines, warn = FALSE))
  defined <- sub("^([A-Za-z0-9_.]+)\\s*<-\\s*function.*$", "\\1",
                 grep("^[A-Za-z0-9_.]+\\s*<-\\s*function", src, value = TRUE))
  missing <- setdiff(shown, defined)
  if (length(missing) > 0) {
    stop("Functions shown but not defined in R/: ",
         paste(missing, collapse = ", "))
  }
  cat("All", length(shown), "functions shown are defined in R/\n")
}

# ---- Graph -------------------------------------------------------------------
make_dot <- function(model_label) paste0("
digraph geoExposeR_softwarex {

  graph [
    rankdir = TB,
    fontname = '", text_font, "',
    nodesep = 0.1,
    ranksep = 0.13,
    bgcolor = 'white',
    compound = true,
    newrank = true,
    splines = ortho,
    pad = 0.05
  ]

  node [
    fontname = '", text_font, "',
    shape = 'box',
    margin = '0.07,0.05',
    penwidth = 1.1
  ]

  edge [
    color = '#1A1A1A',
    penwidth = 1.0,
    arrowsize = 0.55
  ]

  # ==== DATA INGESTION (optional) ====
  ingest [label = ", ingest_label, ",
          ", node_attrs(pal$ingest, dashed = TRUE), "]

  # ==== INPUT FILES read by perform_sensitivity_analysis() ====
  prob   [label = ", prob_label, ", ", node_attrs(pal$files), "]
  conc   [label = ", conc_label, ", ", node_attrs(pal$files), "]
  health [label = ", health_label, ", ", node_attrs(pal$files), "]
  {rank = same; prob; conc; health}
  prob, conc, health [height = 0.66]

  # ==== MAIN PIPELINE ====
  subgraph cluster_main {
    label = ", main_label, "
    labelloc = 't'
    fontname = '", text_font, "'
    style = 'rounded,filled'
    color = '", pal$main[["border"]], "'
    fillcolor = '", pal$main[["fill"]], "'
    penwidth = 1.1
    margin = 4

    integrate [label = ", integrate_label, ", ", node_attrs(pal$integrate), "]
    impute    [label = ", impute_label, ", ", node_attrs(pal$impute), "]
    model     [label = ", model_label, ", ", node_attrs(pal$model), "]
  }

  # ==== OUTPUTS ====
  output [label = ", output_label, ", ", node_attrs(pal$output), "]

  # ==== EDGES ====
  ingest -> prob   [style = dashed]
  ingest -> conc   [style = dashed]
  ingest -> health [style = dashed]
  prob   -> integrate [lhead = cluster_main]
  conc   -> integrate [lhead = cluster_main]
  health -> integrate [lhead = cluster_main]
  integrate -> impute
  impute -> model
  model -> output [ltail = cluster_main]
}
")

# ============================================================================
# Save Figures
# ============================================================================


# rsvg returns premultiplied RGBA; composite onto white to get plain RGB
on_white <- function(bitmap) {
  pmin(bitmap[, , 1:3, drop = FALSE] + (1 - as.vector(bitmap[, , 4])), 1)
}

# Graphviz leaves uneven blank space at the edges; crop the SVG to the drawn
# content plus a uniform margin (in pt) so both outputs are balanced.
crop_svg <- function(svg, margin = 2) {
  vb <- as.numeric(strsplit(
    sub('.*viewBox="([^"]+)".*', "\\1", svg), " "
  )[[1]])
  probe <- rsvg(charToRaw(svg), width = 2000)
  ink <- apply(on_white(probe), c(1, 2), min) < 0.98
  rows <- range(which(rowSums(ink) > 0))
  cols <- range(which(colSums(ink) > 0))
  s <- vb[3] / dim(probe)[2] # pt per probe pixel
  x0 <- vb[1] + (cols[1] - 1) * s - margin
  y0 <- vb[2] + (rows[1] - 1) * s - margin
  w <- diff(cols) * s + s + 2 * margin
  h <- diff(rows) * s + s + 2 * margin
  svg <- sub('width="[^"]+pt" height="[^"]+pt"',
             sprintf('width="%.2fpt" height="%.2fpt"', w, h), svg)
  sub('viewBox="[^"]+"',
      sprintf('viewBox="%.2f %.2f %.2f %.2f"', x0, y0, w, h), svg)
}
render_svg <- function(dot) {
  svg <- as.character(export_svg(grViz(dot)))
  # Render monospace text with Menlo (same advance width as the layout font)
  svg <- gsub(
    sprintf('font-family="%s"', mono_font),
    sprintf('font-family="%s"', mono_render),
    svg,
    fixed = TRUE
  )
  crop_svg(svg)
}

# Flatten onto white (no alpha channel) and store the DPI in the pHYs chunk
save_png <- function(svg, png_file, width_px, dpi) {
  bitmap <- rsvg(charToRaw(svg), width = width_px)
  png::writePNG(on_white(bitmap), png_file, dpi = dpi)
  cat(sprintf(
    "Saved: %s (%d x %d px, %d dpi, %.1f x %.1f mm)\n",
    png_file, dim(bitmap)[2], dim(bitmap)[1], dpi,
    dim(bitmap)[2] / dpi * 25.4, dim(bitmap)[1] / dpi * 25.4
  ))
}

# README/pkgdown figure with the engine grid
readme_svg <- render_svg(make_dot(make_model_label(with_grid = TRUE)))
save_png(readme_svg, file.path(root, "man", "figures", "workflow.png"), 1800, 300)

# Paper figure: vector PDF at natural size (1 pt = 1/72 in) and a bitmapped
# line drawing 3.6 in x 1000 dpi = 3600 px wide, height by aspect
if (dir.exists(paper_dir)) {
  paper_svg <- render_svg(make_dot(model_label))
  pdf_file <- file.path(paper_dir, "workflow.pdf")
  rsvg_pdf(charToRaw(paper_svg), pdf_file)
  cat("Saved:", pdf_file, "\n")
  save_png(paper_svg, file.path(paper_dir, "workflow.png"), 3600, 1000)
} else {
  cat("No paper/ folder: skipped the paper figure\n")
}

cat("\nAll figures generated successfully.\n")
