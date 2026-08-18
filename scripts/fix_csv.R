library(data.table)
library(dplyr)

# Read the first 6 lines of the existing CSV
existing_csv <- fread("results/decode/decode_batch_mr_results.csv", nrows=5) # 5 rows of data + header

# Define the recovered rows
recovered_rows <- data.frame(
  id.exposure = c("eZ7Zru", "FrS02v", "rRXEbn", "pJy4hR", "Z1Bvv4", "zJmT9V", "vYuvC4"),
  id.outcome = c("tm1PKW", "aZkPnA", "NSYMp3", "pwldjT", "BR0UL8", "o3LrZu", "MnFOSV"),
  outcome = c("Breast", "Breast", "Breast", "Breast", "Breast", "Breast", "Breast"),
  exposure = c("FGF5 (deCODE)", "UMOD (deCODE)", "APOE (deCODE)", "KLB (deCODE)", "FGFR4 (deCODE)", "INHBB_Activin_B", "CGREF1 (deCODE)"),
  method = c("Wald ratio", "Wald ratio", "Wald ratio", "Inverse variance weighted", "Wald ratio", "Wald ratio", "Harmonization Failed"),
  nsnp = c(1, 1, 1, 2, 1, 1, NA),
  b = c(-13.50233, -53.77875, 0.09107944, 0.02065292, -125.809, 2.632068, NA),
  se = c(1.804367, 10.71327, 0.02642684, 0.005093133, 35.2898, 0.7940364, NA),
  pval = c(7.256599e-14, 5.171897e-07, 0.0005679512, 5.012299e-05, 0.0003638223, 0.0009171018, NA),
  lo_ci = c(-17.03889, -74.77677, 0.03928284, 0.01067038, -194.977, 1.075757, NA),
  up_ci = c(-9.965775, -32.78073, 0.1428761, 0.03063546, -56.64099, 4.188379, NA),
  or = c(1.367764e-06, 4.407432e-24, 1.095356, 1.020868, 2.300623e-55, 13.90249, NA),
  or_lci95 = c(3.982018e-08, 3.348589e-33, 1.040065, 1.010728, 2.101655e-85, 2.932211, NA),
  or_uci95 = c(4.698065e-05, 5.801087e-15, 1.153587, 1.03111, 2.518429e-25, 65.91588, NA)
)

final_df <- bind_rows(existing_csv, recovered_rows)

# Write a clean combined version
write.csv(final_df, "results/decode/decode_batch_mr_results_combined.csv", row.names=FALSE)
cat("Successfully wrote all cleanly formatted rows to results/decode/decode_batch_mr_results_combined.csv\n")
