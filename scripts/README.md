# Scripts Note

Please note that `83_coloc_abf_palindromic_sensitivity.R` and `86_susie_faithful_ld_stability.R` are the authoritative scripts for the sensitivity sweep. 

Scripts `82_...R`, `84_...R`, and `85_...R` contain a buggy LD builder (omitting the biallelic-SNV restriction, the reference-panel MAF filter, and eigenvalue regularisation) which understates EFNA1 and ATRAID. They are retained strictly for the audit trail. Their numbers must never be cited.
