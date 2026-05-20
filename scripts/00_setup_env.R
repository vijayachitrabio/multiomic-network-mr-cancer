# setup_env.R
# Core MR packages
install.packages(c(
  "TwoSampleMR",       # two-sample MR framework
  "MendelianRandomization",
  "mr.raps",           # overlap-robust MR estimator
  "coloc",             # colocalization
  "susieR"             # fine-mapping (SuSiE)
), repos="https://cloud.r-project.org")

# Visualisation
install.packages(c(
  "ggplot2",
  "ComplexHeatmap",    # pan-cancer heatmap
  "ggraph",            # network diagrams
  "igraph",
  "patchwork",
  "forestplot",
  "circlize"           # colour scales for ComplexHeatmap
), repos="https://cloud.r-project.org")

# Survival analysis (Phase 5 / RAP only)
install.packages(c("survival", "survminer"), repos="https://cloud.r-project.org")

# Data handling and parallelisation
install.packages(c(
  "data.table",
  "tidyverse",
  "furrr",             # parallel map over protein-cancer pairs
  "future"
), repos="https://cloud.r-project.org")

print("Environment setup complete.")
