# R/globals.R
utils::globalVariables(c(
  "seqname","start","end",
  "gene_state_contrasts",
  "avg_ln_neg_dNdS","avg_ln_pos_dNdS",
  "or","contrast","or_lo","or_hi",
  "odds_ratio","or_ci_lower","or_ci_upper",
  "label","significant",
  "enrichment","enrichment_plot","y_lab","p_adj","pos_count","is_inf_enrichment",
  "dNdS_A","dNdS_B","delta"
))
