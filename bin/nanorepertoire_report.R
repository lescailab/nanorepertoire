#!/usr/bin/env Rscript
# =============================================================================
# Nanobody Repertoire Analysis Report Generator
# =============================================================================
# Inspired by:
#   - Deschaght et al. 2017 (doi:10.3389/fimmu.2017.00420)
#   - Spinelli et al. 2022 (doi:10.3389/fimmu.2022.927966) - VHH CDR3 properties
#   - MiXCR repertoire analysis framework (Bolotin et al. 2015)
#   - IMGT/V-QUEST VHH annotation standards
#   - Alpseq pipeline (specialized nanobody repertoire processing)
# =============================================================================

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(htmltools)
  library(htmlwidgets)
  library(RColorBrewer)
  library(scales)
})

# ---- 0. Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript nanorepertoire_report.R <input.RData> <output.html>")
}
rdata_file <- args[1]
output_html <- args[2]

cat("Loading data from:", rdata_file, "\n")
load(rdata_file)

# ---- Color palette ----
SAMPLE_COLORS <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
  "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7"
)

samples <- unique(cdrcounts$Sample)
pal <- setNames(SAMPLE_COLORS[seq_along(samples)], samples)

# =============================================================================
# 1. SUMMARY STATISTICS
# =============================================================================
total_clusters   <- sum(clustercounts$Clusters)
total_cdr3       <- sum(cdrcounts$Unique_CDR3s)
total_seqs       <- nrow(fastaSeq)

cdr3_seqs        <- fastaSeq %>% filter(!is.na(CDR3) & CDR3 != "" & unique == "unique")
cdr3_len         <- nchar(cdr3_seqs$CDR3)

mean_len  <- round(mean(cdr3_len, na.rm = TRUE), 1)
sd_len    <- round(sd(cdr3_len, na.rm = TRUE), 1)
median_len <- median(cdr3_len, na.rm = TRUE)

# Shannon diversity per sample
shannon_div <- fastaSeq %>%
  filter(!is.na(CDR3) & CDR3 != "" & unique == "unique") %>%
  group_by(sample) %>%
  summarise(
    n_unique = n(),
    shannon  = {
      # Based on CDR3 length classes (as proxy for diversity buckets)
      lens <- nchar(CDR3)
      t    <- table(lens)
      p    <- as.numeric(t) / sum(t)
      -sum(p * log(p + 1e-12))
    },
    .groups = "drop"
  )

# Clonotype expansion: proportion of singletons vs expanded
expansion_df <- clustercounts %>%
  mutate(
    singletons   = Clusters - Clusters_of_5,
    small        = Clusters_of_5 - Clusters_of_100,
    medium       = Clusters_of_100 - Clusters_of_1000,
    large        = Clusters_of_1000
  ) %>%
  select(Sample, singletons, small, medium, large)

# =============================================================================
# 2. REPERTOIRE METRICS (MiXCR-equivalent)
# =============================================================================

# D50 Index: minimum number of unique CDR3s that account for 50% of total sequences
# (a clonality metric used in immunology, e.g. Robins 2009, J Immunol Methods)
compute_d50 <- function(counts) {
  sorted_desc <- sort(counts, decreasing = TRUE)
  cumulative  <- cumsum(sorted_desc)
  target      <- sum(sorted_desc) * 0.5
  idx         <- which(cumulative >= target)[1]
  round(idx / length(counts) * 100, 2)
}

# Gini coefficient - measure of clonal inequality (0=equal, 1=monopoly)
compute_gini <- function(x) {
  x <- sort(x)
  n <- length(x)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - ((n + 1) / n)
}

diversity_metrics <- fastaSeq %>%
  filter(!is.na(CDR3) & CDR3 != "" & unique == "unique") %>%
  group_by(sample) %>%
  summarise(
    Unique_CDR3s = n(),
    Shannon  = round({
      lens <- nchar(CDR3); t <- table(lens); p <- as.numeric(t)/sum(t)
      -sum(p * log(p + 1e-12))
    }, 3),
    .groups = "drop"
  ) %>%
  left_join(clustercounts, by = c("sample" = "Sample")) %>%
  mutate(
    Pct_expanded = round(Clusters_of_5 / Clusters * 100, 1),
    Pct_large    = round(Clusters_of_1000 / Clusters * 100, 2)
  )

# =============================================================================
# 3. CDR3 AMINO ACID COMPOSITION ANALYSIS
# (Literature: Spinelli 2022, showed VHH CDR3 has distinct AA preference vs VH)
# =============================================================================

aa_standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

aa_freq_df <- cdr3_seqs %>%
  mutate(CDR3 = toupper(CDR3)) %>%
  group_by(sample) %>%
  summarise(
    aa_string = paste(CDR3, collapse = ""),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    aa_counts = list({
      chars <- strsplit(aa_string, "")[[1]]
      chars <- chars[chars %in% aa_standard]
      total <- length(chars)
      setNames(as.numeric(table(factor(chars, levels = aa_standard))) / total * 100, aa_standard)
    })
  ) %>%
  unnest_wider(aa_counts) %>%
  select(-aa_string)

# Pivot for plotting
aa_long <- aa_freq_df %>%
  pivot_longer(-sample, names_to = "AminoAcid", values_to = "Frequency_pct")

# AA physicochemical groups (for color annotation, IMGT classification)
aa_groups <- c(
  A="Aliphatic", G="Aliphatic", I="Aliphatic", L="Aliphatic", V="Aliphatic",
  F="Aromatic", W="Aromatic", Y="Aromatic",
  D="Acidic", E="Acidic",
  K="Basic", R="Basic", H="Basic",
  S="Polar", T="Polar", C="Polar", M="Polar", N="Polar", Q="Polar", P="Polar"
)
aa_long <- aa_long %>%
  mutate(Group = aa_groups[AminoAcid])

# =============================================================================
# 4. BUILD INTERACTIVE PLOTLY FIGURES
# =============================================================================

# ---- Fig 1: Cluster size distribution (stacked bar, MiXCR-style) ----
exp_long <- expansion_df %>%
  pivot_longer(-Sample, names_to = "Size_class", values_to = "Count") %>%
  mutate(Size_class = factor(Size_class, levels = c("large","medium","small","singletons"),
         labels = c("≥1000 members", "100-999 members", "5-99 members", "Singletons (<5)")))

fig_clusters_stacked <- plot_ly(
  exp_long, x = ~Sample, y = ~Count, color = ~Size_class,
  colors = c("#d62728","#ff7f0e","#1f77b4","#aec7e8"),
  type = "bar", text = ~Count, textposition = "inside",
  hovertemplate = "<b>%{x}</b><br>%{fullData.name}: %{y:,}<extra></extra>"
) %>%
  layout(
    barmode = "stack",
    title = list(text = "Cluster Size Distribution per Sample", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "Sample"),
    yaxis = list(title = "Number of Clusters"),
    legend = list(title = list(text = "Cluster Size Class")),
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 2: Cluster count summary (grouped bar) ----
cc_long <- clustercounts %>%
  pivot_longer(-Sample, names_to = "Threshold", values_to = "Count") %>%
  mutate(Threshold = recode(Threshold,
    Clusters           = "Total Clusters",
    Clusters_of_5      = "≥5 members",
    Clusters_of_100    = "≥100 members",
    Clusters_of_1000   = "≥1000 members"
  )) %>%
  mutate(Threshold = factor(Threshold, levels = c("Total Clusters","≥5 members","≥100 members","≥1000 members")))

fig_cluster_summary <- plot_ly(
  cc_long, x = ~Sample, y = ~Count, color = ~Threshold,
  colors = c("#003f5c","#58508d","#bc5090","#ff6361"),
  type = "bar",
  hovertemplate = "<b>%{x}</b><br>%{fullData.name}: %{y:,}<extra></extra>"
) %>%
  layout(
    barmode = "group",
    title = list(text = "Cluster Abundance at Multiple Size Thresholds", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "Sample"),
    yaxis = list(title = "Count (log scale)", type = "log"),
    legend = list(title = list(text = "Size Threshold")),
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 3: CDR3 length distribution (ridge/violin equiv -> overlapping area) ----
len_df <- fastaSeq %>%
  filter(!is.na(CDR3) & CDR3 != "" & unique == "unique") %>%
  mutate(cdr3_length = nchar(CDR3)) %>%
  filter(cdr3_length >= 1, cdr3_length <= 40)

fig_cdr3_length <- plot_ly()
for (s in samples) {
  sdata <- len_df %>% filter(sample == s)
  fig_cdr3_length <- fig_cdr3_length %>%
    add_trace(
      x = sdata$cdr3_length,
      type = "violin",
      name = s,
      fillcolor = paste0(pal[s], "80"),
      line = list(color = pal[s]),
      meanline = list(visible = TRUE),
      box = list(visible = TRUE),
      side = "positive",
      hovertemplate = paste0("<b>", s, "</b><br>CDR3 Length: %{x}<extra></extra>")
    )
}
fig_cdr3_length <- fig_cdr3_length %>%
  layout(
    title = list(text = "CDR3 Length Distribution (VHH)", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "CDR3 Length (amino acids)",
                 range = c(0, 42),
                 tickvals = seq(0, 42, 5)),
    yaxis = list(title = "Sample"),
    violingap = 0.1,
    violinmode = "overlay",
    annotations = list(
      list(x = 19, y = 0.02, text = "VHH median ~19 AA<br>(Spinelli 2022)",
           showarrow = TRUE, arrowhead = 2, ax = 30, ay = -30,
           font = list(color = "#555555", size = 11))
    ),
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 4: CDR3 length histogram from hist files ----
fig_len_hist <- plot_ly()
for (s in samples) {
  hdata <- cdrhists %>% filter(Sample == s)
  fig_len_hist <- fig_len_hist %>%
    add_trace(
      x = hdata$Size, y = hdata$Count,
      type  = "scatter", mode = "lines+markers",
      name  = s,
      line  = list(color = pal[s], width = 2),
      marker = list(color = pal[s], size = 5),
      fill  = "tozeroy",
      fillcolor = paste0(pal[s], "30"),
      hovertemplate = paste0("<b>", s, "</b><br>Length: %{x} AA<br>Count: %{y:,}<extra></extra>")
    )
}
fig_len_hist <- fig_len_hist %>%
  layout(
    title = list(text = "CDR3 Length Frequency Profile", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "CDR3 Length (AA)", range = c(0, 45)),
    yaxis = list(title = "Number of Unique CDR3s"),
    hovermode = "x unified",
    legend = list(title = list(text = "Sample")),
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 5: Unique CDR3s per sample (bar) ----
fig_cdr3_counts <- plot_ly(
  cdrcounts, x = ~Sample, y = ~Unique_CDR3s,
  color = ~Sample, colors = unname(pal),
  type = "bar",
  text = ~format(Unique_CDR3s, big.mark=","),
  textposition = "outside",
  hovertemplate = "<b>%{x}</b><br>Unique CDR3s: %{y:,}<extra></extra>"
) %>%
  layout(
    title = list(text = "Unique CDR3 Sequences per Sample", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "Sample"),
    yaxis = list(title = "Number of Unique CDR3s"),
    showlegend = FALSE,
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 6: Amino Acid composition heatmap (per sample) ----
aa_matrix <- aa_freq_df %>%
  select(-sample) %>%
  as.matrix()
rownames(aa_matrix) <- aa_freq_df$sample

fig_aa_heatmap <- plot_ly(
  z = aa_matrix,
  x = colnames(aa_matrix),
  y = rownames(aa_matrix),
  type = "heatmap",
  colorscale = "Viridis",
  reversescale = FALSE,
  hovertemplate = "<b>%{y}</b><br>AA %{x}: %{z:.2f}%<extra></extra>",
  colorbar = list(title = "Frequency (%)")
) %>%
  layout(
    title = list(text = "Amino Acid Composition of CDR3 Sequences", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "Amino Acid"),
    yaxis = list(title = "Sample"),
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 7: AA composition grouped bar (CDR3 context, all samples) ----
fig_aa_bar <- plot_ly()
for (s in samples) {
  sdata <- aa_long %>% filter(sample == s)
  fig_aa_bar <- fig_aa_bar %>%
    add_trace(
      x = sdata$AminoAcid, y = sdata$Frequency_pct,
      type = "bar", name = s,
      marker = list(color = pal[s]),
      hovertemplate = paste0("<b>", s, "</b><br>%{x}: %{y:.2f}%<extra></extra>")
    )
}
fig_aa_bar <- fig_aa_bar %>%
  layout(
    barmode = "group",
    title = list(text = "CDR3 Amino Acid Usage (% of CDR3 residues)", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "Amino Acid", categoryorder = "array",
                 categoryarray = aa_standard),
    yaxis = list(title = "Frequency (%)"),
    legend = list(title = list(text = "Sample")),
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 8: Diversity metrics bubble chart ----
div_plot_df <- diversity_metrics %>%
  mutate(sample = factor(sample, levels = samples))

fig_diversity <- plot_ly(
  div_plot_df,
  x = ~Shannon, y = ~Pct_expanded,
  size = ~Unique_CDR3s, color = ~sample, colors = unname(pal),
  type = "scatter", mode = "markers+text",
  text = ~sample, textposition = "top center",
  marker = list(sizemode = "diameter", sizeref = 0.5, opacity = 0.75),
  hovertemplate = paste0(
    "<b>%{text}</b><br>Shannon Diversity: %{x:.3f}<br>% Expanded Clones: %{y:.1f}%<br>",
    "Unique CDR3s: %{marker.size:,}<extra></extra>"
  )
) %>%
  layout(
    title = list(text = "Repertoire Diversity Overview", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "Shannon Diversity Index (CDR3 length distribution)"),
    yaxis = list(title = "% Expanded Clonotypes (≥5 members)"),
    showlegend = FALSE,
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 9: Identity to cluster representative (clonal expansion view) ----
# Takes a sample of rows for performance
set.seed(42)
identity_sample <- clusterbig %>%
  group_by(Sample) %>%
  slice_sample(n = min(5000, n())) %>%
  ungroup()

fig_identity <- plot_ly()
for (s in samples) {
  id_data <- identity_sample %>% filter(Sample == s)
  fig_identity <- fig_identity %>%
    add_trace(
      x = id_data$Identity,
      type = "histogram",
      name = s,
      marker = list(color = paste0(pal[s], "BB"), line = list(color = pal[s], width = 0.5)),
      opacity = 0.7,
      nbinsx = 50,
      hovertemplate = paste0("<b>", s, "</b><br>Identity bin: %{x:.1f}%<br>Count: %{y:,}<extra></extra>")
    )
}
fig_identity <- fig_identity %>%
  layout(
    barmode = "overlay",
    title = list(text = "Sequence Identity to Cluster Representative (SHM proxy)",
                 font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "% Identity to Cluster Representative", range = c(0, 101)),
    yaxis = list(title = "Number of Sequences (sampled)"),
    legend = list(title = list(text = "Sample")),
    annotations = list(
      list(x = 99, y = 0, text = "100% = identical clone<br>(potential PCR duplicate)",
           showarrow = FALSE, font = list(size = 10, color = "#888888"), align = "right", xanchor = "right")
    ),
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff"
  )

# ---- Fig 10: Top-10 most abundant clusters per sample (clonotype table viz) ----
top_clusters_df <- clusterbig %>%
  group_by(Sample) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 10) %>%
  mutate(Rank = row_number(), label = paste0("Rank ", Rank, " (Rep: ", Representative, ")")) %>%
  ungroup()

fig_top_clones <- plot_ly()
for (s in samples) {
  tdata <- top_clusters_df %>% filter(Sample == s) %>% arrange(Rank)
  fig_top_clones <- fig_top_clones %>%
    add_trace(
      x = tdata$Count, y = reorder(tdata$label, -tdata$Rank),
      type = "bar", orientation = "h",
      name = s,
      marker = list(color = pal[s]),
      hovertemplate = paste0("<b>", s, "</b><br>%{y}<br>Members: %{x:,}<extra></extra>")
    )
}
fig_top_clones <- fig_top_clones %>%
  layout(
    barmode = "group",
    title = list(text = "Top-10 Largest Clusters per Sample", font = list(size = 16, color = "#2c3e50")),
    xaxis = list(title = "Number of Sequences in Cluster"),
    yaxis = list(title = "", autorange = "reversed"),
    legend = list(title = list(text = "Sample")),
    plot_bgcolor  = "#fafafa",
    paper_bgcolor = "#ffffff",
    height = 500
  )

# =============================================================================
# 5. HELPER: serialize plotly to HTML fragment
# =============================================================================
fig_to_html <- function(fig, id = NULL) {
  as.character(as_widget(fig))
}

# =============================================================================
# 6. BUILD THE HTML DOCUMENT
# =============================================================================

css <- "
<style>
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: 'Inter', sans-serif; background: #f0f4f8; color: #2c3e50; }
  
  /* Header */
  .report-header {
    background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
    color: white; padding: 48px 60px 40px;
    border-bottom: 4px solid #e94560;
  }
  .report-header h1 { font-size: 2.4em; font-weight: 700; letter-spacing: -0.5px; margin-bottom: 8px; }
  .report-header .subtitle { font-size: 1.1em; opacity: 0.8; font-weight: 300; }
  .report-header .meta { margin-top: 20px; display: flex; gap: 32px; flex-wrap: wrap; }
  .report-header .meta-item { font-size: 0.88em; opacity: 0.7; }
  .report-header .meta-item strong { display: block; font-size: 1.1em; opacity: 1; color: #e94560; }

  /* Navigation */
  .nav-bar {
    background: #16213e; position: sticky; top: 0; z-index: 1000;
    display: flex; padding: 0 60px; box-shadow: 0 2px 12px rgba(0,0,0,0.3);
  }
  .nav-bar a {
    color: rgba(255,255,255,0.7); text-decoration: none; padding: 14px 16px;
    font-size: 0.85em; font-weight: 500; transition: all 0.2s; white-space: nowrap;
    border-bottom: 3px solid transparent;
  }
  .nav-bar a:hover { color: #e94560; border-bottom-color: #e94560; }

  /* Main content */
  .container { max-width: 1400px; margin: 0 auto; padding: 40px 60px; }

  /* KPI Cards */
  .kpi-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 40px; }
  .kpi-card {
    background: white; border-radius: 16px; padding: 28px 24px;
    box-shadow: 0 2px 16px rgba(0,0,0,0.06); border-left: 5px solid #e94560;
    transition: transform 0.2s, box-shadow 0.2s;
  }
  .kpi-card:hover { transform: translateY(-3px); box-shadow: 0 6px 24px rgba(0,0,0,0.12); }
  .kpi-card .kpi-label { font-size: 0.8em; text-transform: uppercase; letter-spacing: 0.8px; color: #7f8c8d; font-weight: 600; margin-bottom: 10px; }
  .kpi-card .kpi-value { font-size: 2.2em; font-weight: 700; color: #2c3e50; line-height: 1; }
  .kpi-card .kpi-sub { font-size: 0.82em; color: #95a5a6; margin-top: 6px; }
  .kpi-card:nth-child(2) { border-left-color: #3498db; }
  .kpi-card:nth-child(3) { border-left-color: #2ecc71; }
  .kpi-card:nth-child(4) { border-left-color: #f39c12; }
  .kpi-card:nth-child(5) { border-left-color: #9b59b6; }

  /* Sections */
  .section { margin-bottom: 56px; }
  .section-header {
    display: flex; align-items: center; gap: 14px;
    margin-bottom: 24px; padding-bottom: 14px;
    border-bottom: 2px solid #e8ecef;
  }
  .section-number {
    background: linear-gradient(135deg, #e94560, #c0392b);
    color: white; width: 38px; height: 38px; border-radius: 10px;
    display: flex; align-items: center; justify-content: center;
    font-size: 1em; font-weight: 700; flex-shrink: 0;
  }
  .section-header h2 { font-size: 1.5em; font-weight: 600; color: #2c3e50; }
  .section-intro {
    background: #eef2f7; border-left: 4px solid #3498db;
    padding: 14px 18px; border-radius: 0 8px 8px 0;
    font-size: 0.92em; color: #34495e; margin-bottom: 20px;
    line-height: 1.6;
  }
  .section-intro cite { font-style: italic; color: #7f8c8d; }

  /* Chart cards */
  .chart-card {
    background: white; border-radius: 16px;
    box-shadow: 0 2px 16px rgba(0,0,0,0.06);
    overflow: hidden; margin-bottom: 24px;
  }
  .chart-card-header { padding: 18px 24px 0; }
  .chart-card-header h3 { font-size: 1.05em; font-weight: 600; color: #2c3e50; }
  .chart-card-header p { font-size: 0.86em; color: #7f8c8d; margin-top: 4px; line-height: 1.5; }
  .chart-card-body { padding: 12px 8px 8px; }

  /* Metrics table */
  .metrics-table { width: 100%; border-collapse: collapse; font-size: 0.9em; }
  .metrics-table thead tr { background: #2c3e50; color: white; }
  .metrics-table th { padding: 12px 16px; text-align: left; font-weight: 600; font-size: 0.85em; letter-spacing: 0.4px; }
  .metrics-table td { padding: 10px 16px; border-bottom: 1px solid #f0f4f8; }
  .metrics-table tbody tr:nth-child(even) { background: #f8f9fb; }
  .metrics-table tbody tr:hover { background: #eef2f7; }
  .badge {
    display: inline-block; padding: 3px 10px; border-radius: 20px;
    font-size: 0.78em; font-weight: 600;
  }
  .badge-high { background: #fde8e8; color: #c0392b; }
  .badge-med  { background: #fef6e4; color: #f39c12; }
  .badge-low  { background: #e8f8f0; color: #27ae60; }

  /* Two-column layout */
  .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }

  /* Methods section */
  .methods-card { background: white; border-radius: 16px; padding: 28px; box-shadow: 0 2px 16px rgba(0,0,0,0.06); }
  .methods-card h3 { font-size: 1.05em; font-weight: 600; margin-bottom: 12px; color: #2c3e50; }
  .methods-card p { color: #555; font-size: 0.9em; line-height: 1.7; margin-bottom: 10px; }
  .methods-card ul { padding-left: 20px; color: #555; font-size: 0.88em; line-height: 1.8; }

  /* References */
  .references { font-size: 0.83em; color: #7f8c8d; line-height: 1.7; }
  .references li { margin-bottom: 6px; }

  /* Footer */
  .footer {
    background: #2c3e50; color: rgba(255,255,255,0.6);
    text-align: center; padding: 28px 60px; font-size: 0.82em; margin-top: 40px;
  }

  @media (max-width: 900px) {
    .container { padding: 24px 20px; }
    .report-header { padding: 32px 20px; }
    .nav-bar { padding: 0 20px; }
    .grid-2 { grid-template-columns: 1fr; }
    .kpi-grid { grid-template-columns: repeat(2, 1fr); }
  }
</style>
"

# Build diversity table HTML
make_table <- function(df) {
  rows <- apply(df, 1, function(r) {
    exp_val <- as.numeric(r["Pct_expanded"])
    badge_class <- if (!is.na(exp_val) && exp_val > 25) "badge-high" else if (!is.na(exp_val) && exp_val > 10) "badge-med" else "badge-low"
    paste0("<tr>",
      "<td><strong>", r["sample"], "</strong></td>",
      "<td>", format(as.numeric(r["Unique_CDR3s"]), big.mark=","), "</td>",
      "<td>", format(as.numeric(r["Clusters"]), big.mark=","), "</td>",
      "<td>", r["Shannon"], "</td>",
      "<td><span class='badge ", badge_class, "'>", r["Pct_expanded"], "%</span></td>",
      "<td>", r["Pct_large"], "%</td>",
    "</tr>")
  })
  paste0(
    "<table class='metrics-table'>",
    "<thead><tr><th>Sample</th><th>Unique CDR3s</th><th>Total Clusters</th>",
    "<th>Shannon H'</th><th>% Expanded (≥5)</th><th>% Large (≥1000)</th></tr></thead>",
    "<tbody>", paste(rows, collapse=""), "</tbody></table>"
  )
}

# Helper to wrap a plotly widget
wrap_chart <- function(fig, title, description) {
  widget_html <- as.character(as_widget(fig))
  paste0(
    "<div class='chart-card'>",
    "<div class='chart-card-header'>",
    "<h3>", title, "</h3>",
    "<p>", description, "</p>",
    "</div>",
    "<div class='chart-card-body'>", widget_html, "</div>",
    "</div>"
  )
}

# Format numbers
fmt_num <- function(x) format(round(x), big.mark=",")

html_body <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Nanobody Repertoire Report</title>
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
', css, '
</head>
<body>

<!-- HEADER -->
<div class="report-header">
  <h1>🧬 Nanobody Repertoire Analysis</h1>
  <div class="subtitle">VHH Sequencing &amp; Clonotype Analysis Report</div>
  <div class="meta">
    <div class="meta-item"><strong>', fmt_num(total_seqs), '</strong>Total Sequences</div>
    <div class="meta-item"><strong>', fmt_num(total_clusters), '</strong>Total Clusters</div>
    <div class="meta-item"><strong>', fmt_num(total_cdr3), '</strong>Unique CDR3s</div>
    <div class="meta-item"><strong>', length(samples), '</strong>Samples Analysed</div>
    <div class="meta-item"><strong>Generated</strong>', format(Sys.Date(), "%B %d, %Y"), '</div>
  </div>
</div>

<!-- NAVIGATION -->
<nav class="nav-bar">
  <a href="#summary">Summary</a>
  <a href="#clusters">Clonal Clusters</a>
  <a href="#cdr3">CDR3 Diversity</a>
  <a href="#aminoacids">Amino Acids</a>
  <a href="#clonality">Clonality</a>
  <a href="#methods">Methods</a>
</nav>

<div class="container">

  <!-- SECTION 0: KPI Summary Cards -->
  <div id="summary" class="section" style="padding-top: 20px;">
    <div class="kpi-grid">
      <div class="kpi-card">
        <div class="kpi-label">Total Sequences</div>
        <div class="kpi-value">', fmt_num(total_seqs), '</div>
        <div class="kpi-sub">Across all samples</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Total Clusters (CD-HIT)</div>
        <div class="kpi-value">', fmt_num(total_clusters), '</div>
        <div class="kpi-sub">90% identity threshold</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Unique CDR3 Sequences</div>
        <div class="kpi-value">', fmt_num(total_cdr3), '</div>
        <div class="kpi-sub">Novel paratopes</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Mean CDR3 Length</div>
        <div class="kpi-value">', mean_len, ' AA</div>
        <div class="kpi-sub">SD: ', sd_len, ' | Median: ', median_len, '</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Samples</div>
        <div class="kpi-value">', length(samples), '</div>
        <div class="kpi-sub">Independent libraries</div>
      </div>
    </div>
  </div>

  <!-- SECTION 1: CLONAL CLUSTERS -->
  <div id="clusters" class="section">
    <div class="section-header">
      <div class="section-number">1</div>
      <h2>Clonal Cluster Analysis</h2>
    </div>
    <div class="section-intro">
      Amino acid sequences were clustered using <strong>CD-HIT</strong> at 90% identity, following the approach of Deschaght et al. (2017).
      Expanded clonotypes (clusters with ≥5 members) represent B-cell lineages that have undergone somatic expansion, 
      and are strong candidates for antigen-specific nanobodies. Large clusters (≥1000 members) likely reflect dominant 
      clonal responses following immunization, analogous to clonotype expansion metrics reported by MiXCR.
      <cite>Deschaght et al. 2017. Front. Immunol. doi:10.3389/fimmu.2017.00420</cite>
    </div>
    
    <div class="section-header" style="margin-bottom: 14px; border: none;">
      <h2 style="font-size: 1.1em;">📊 Repertoire Diversity Metrics</h2>
    </div>
    ', make_table(diversity_metrics), '

    <div style="margin-top: 24px;">
    ', wrap_chart(fig_clusters_stacked,
      "Clonal Expansion Profile (Stacked)",
      "Breakdown of clusters by size class. Singletons represent one-off sequences (noise or rare B cells), while large clusters indicate dominant clonotypes expanded by antigen-driven selection. Inspired by MiXCR clonotype report format."),
    wrap_chart(fig_cluster_summary,
      "Cluster Count at Multiple Size Thresholds",
      "Comparison of clonotype abundance across 4 size thresholds per sample. Log scale reveals hierarchical representation of the repertoire at different expansion levels."),
    wrap_chart(fig_top_clones,
      "Top-10 Largest Clonotypes per Sample",
      "The most abundant clusters, by number of member sequences, per sample. Dominant clones in immunized samples may represent the most matured antigen-specific binders."),
    '
    </div>
  </div>

  <!-- SECTION 2: CDR3 DIVERSITY -->
  <div id="cdr3" class="section">
    <div class="section-header">
      <div class="section-number">2</div>
      <h2>CDR3 Diversity Analysis</h2>
    </div>
    <div class="section-intro">
      The CDR3 (Complementarity Determining Region 3) of VHH nanobodies is the primary determinant of antigen binding and epitope specificity.
      VHH CDR3s are notably longer and more diverse than conventional VH CDR3s (mean ~17 AA for VH vs ~19 AA for VHH — Spinelli et al. 2022), 
      enabling access to cryptic epitopes such as enzyme active sites and viral canyons. CDR3 length diversity is a proxy for 
      repertoire breadth; antigen-enriched samples often show a skewed distribution toward longer CDR3s.
      <cite>Spinelli et al. 2022. Front. Immunol. doi:10.3389/fimmu.2022.927966 | Muyldermans, S. 2013. Annu. Rev. Biochem.</cite>
    </div>
    
    <div class="grid-2">
    ', wrap_chart(fig_cdr3_counts,
      "Unique CDR3 Sequences per Sample",
      "Total number of distinct CDR3 sequences identified in each sample after filtering for standard amino acids. Higher diversity indicates a broader immune response."),
    wrap_chart(fig_cdr3_length,
      "CDR3 Length Distribution (Violin)",
      "Distribution of CDR3 amino acid lengths per sample. Internal box shows interquartile range; center line shows median. VHH CDR3s typically range from 6-24 AA (Spinelli 2022)."),
    '
    </div>
    ', wrap_chart(fig_len_hist,
      "CDR3 Length Frequency Profile (Area Plot)",
      "Overlapping area plot of CDR3 length frequencies. Peaks at different lengths may indicate distinct B-cell lineages with characteristic loop structures. Shorter CDR3s (<12 AA) tend to form flat epitope-binding β-strand, while longer CDR3s (>16 AA) can form finger-like protrusions for cavity binding."),
    '
  </div>

  <!-- SECTION 3: AMINO ACID COMPOSITION -->
  <div id="aminoacids" class="section">
    <div class="section-header">
      <div class="section-number">3</div>
      <h2>CDR3 Amino Acid Repertoire</h2>
    </div>
    <div class="section-intro">
      The amino acid composition of CDR3 loops reflects selection pressures for antigen binding. 
      VHH CDR3s are enriched in <strong>Tyrosine (Y)</strong> and <strong>Serine (S)</strong>, which participate in hydrogen bonding and van der Waals contacts (Desmyter et al. 2002).
      The presence of <strong>Cysteine (C)</strong> in CDR3 can form disulfide bonds with CDR1, creating a structural loop that extends reach.
      Charged residues (R, D, E, K) contribute to electrostatic complementarity with antigens.
      Comparing composition across samples can reveal convergence of antigen-driven selection.
      <cite>Desmyter et al. 2002. J. Biol. Chem. doi:10.1074/jbc.D200025200</cite>
    </div>
    ', wrap_chart(fig_aa_heatmap,
      "Amino Acid Frequency Heatmap",
      "Heatmap of amino acid usage frequency (%) across all unique CDR3 sequences, stratified by sample. Dark colors = high frequency. Key VHH-associated residues: Y, S, G, T."),
    wrap_chart(fig_aa_bar,
      "Side-by-Side Amino Acid Usage",
      "Grouped bar chart comparing amino acid composition across samples. Differences may reflect antigen-driven B-cell selection or sample quality differences."),
    '
  </div>

  <!-- SECTION 4: CLONALITY & SHM PROXY -->
  <div id="clonality" class="section">
    <div class="section-header">
      <div class="section-number">4</div>
      <h2>Clonality &amp; Somatic Hypermutation Proxy</h2>
    </div>
    <div class="section-intro">
      The distribution of sequence identity to the cluster representative serves as a proxy for 
      <strong>somatic hypermutation (SHM)</strong> within clonal lineages. 
      Sequences with &lt;100% identity to their representative have been diversified through SHM. 
      A bimodal distribution — a peak at 100% and one at 90-99% — indicates ongoing affinity maturation.
      This is analogous to the V-gene SHM rate reported by MiXCR and IMGT/VQuest for conventional antibody repertoires.
      Sequences clustered at high identity (≥99%) are likely PCR duplicates or clonally identical B cells.
      <cite>Bolotin et al. 2015. Nature Methods. doi:10.1038/nmeth.3364 (MiXCR)</cite>
    </div>
    <div class="grid-2">
    ', wrap_chart(fig_identity,
      "Identity to Cluster Representative (SHM Proxy)",
      "Histogram of sequence identity (%) to the CD-HIT cluster representative. Sequences at 100% identity are clonal copies; identity <100% reflects somatic hypermutation accumulated during affinity maturation."),
    wrap_chart(fig_diversity,
      "Repertoire Diversity Landscape",
      "Bubble chart combining Shannon diversity, % expanded clones, and total unique CDR3s per sample. Samples toward top-right have both high diversity and high clonal expansion — hallmarks of a robust antibody response."),
    '
    </div>
  </div>

  <!-- SECTION 5: METHODS -->
  <div id="methods" class="section">
    <div class="section-header">
      <div class="section-number">5</div>
      <h2>Methods Summary</h2>
    </div>
    <div class="grid-2">
      <div class="methods-card">
        <h3>🔬 Sequencing &amp; Pre-processing</h3>
        <p>Paired-end FASTQ reads were adapter-trimmed with <strong>Cutadapt</strong> and merged with <strong>FLASH</strong> (max. overlap 300 bp). 
        Merged reads were translated using an in-house Python script identifying conserved VHH framework start/end motifs (inspired by Deschaght 2017).</p>
        <p>CDR3 regions were predicted using <strong>nanocdr-x</strong>, a deep learning model trained on camelid antibody sequences.</p>
      </div>
      <div class="methods-card">
        <h3>🧩 Clustering (CD-HIT)</h3>
        <p>Translated protein sequences were clustered with <strong>CD-HIT</strong> at 90% amino acid identity, 
        consistent with Deschaght et al. (2017). Four size thresholds were analysed: 
        1 (singleton), ≥5, ≥100, ≥1000 members per cluster.</p>
        <ul>
          <li>Identity threshold: 0.90</li>
          <li>Word size: 5</li>
          <li>Sequence coverage: 0.9</li>
        </ul>
      </div>
      <div class="methods-card">
        <h3>📈 Diversity Metrics</h3>
        <p><strong>Shannon Diversity Index (H\'):</strong> calculated on CDR3 length distribution per sample. 
        Higher H\' indicates more even length distribution (broader repertoire).</p>
        <p><strong>% Expanded Clonotypes:</strong> fraction of clusters with ≥5 members, as a measure of somatic expansion.</p>
        <p><strong>Identity distribution:</strong> used as a SHM proxy, similar to V-gene mutation frequency in MiXCR reports.</p>
      </div>
      <div class="methods-card">
        <h3>📚 Key References</h3>
        <ul class="references">
          <li>Deschaght et al. 2017. <em>Front. Immunol.</em> doi:10.3389/fimmu.2017.00420</li>
          <li>Spinelli et al. 2022. <em>Front. Immunol.</em> doi:10.3389/fimmu.2022.927966</li>
          <li>Muyldermans S. 2013. <em>Annu. Rev. Biochem.</em> doi:10.1146/annurev-biochem-071812-105327</li>
          <li>Bolotin et al. 2015. <em>Nature Methods</em>. doi:10.1038/nmeth.3364</li>
          <li>Desmyter et al. 2002. <em>J. Biol. Chem.</em> doi:10.1074/jbc.D200025200</li>
          <li>Mitchell &amp; Colwell 2018. <em>Proteins</em>. doi:10.1002/prot.25497</li>
          <li>Li &amp; Godzik 2006. <em>Bioinformatics</em>. doi:10.1093/bioinformatics/btl158 (CD-HIT)</li>
        </ul>
      </div>
    </div>
  </div>

</div>

<div class="footer">
  Nanobody Repertoire Report • Generated by lescailab/nanorepertoire pipeline • ', format(Sys.Date()), '
  <br><span style="opacity:0.5;">Based on Deschaught 2017 framework | Plots powered by Plotly.js</span>
</div>

</body>
</html>
')

cat("Writing HTML report to:", output_html, "\n")
writeLines(html_body, output_html)
cat("Done! Report saved to:", output_html, "\n")
