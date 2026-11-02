# Contrasts two comparisons of the bundled example. Expected values are calculated here from
# the per-gene tables, independently of dnds_contrast().
skip_if_not_installed("orthologr")
skip_if_not_installed("txdbmaker")
skip_if(Sys.which("diamond") == "", "DIAMOND is not on PATH")

use_example_copy()

test_that("contrasts have the right direction and pairing", {
  suppressMessages(suppressWarnings({
    extract_cds(comparison_file = "two_comparisons.txt", output_dir = ".")
    calculate_dnds(comparison_file = "two_comparisons.txt", output_dir = ".", comp_cores = test_cores())
    append_annotations(comparison_file = "two_comparisons.txt", output_dir = ".")
    dnds_contrast(comparison_file = "two_comparisons.txt", contrast_file = "contrast_file.txt", output_dir = ".", make_plots = FALSE)
  }))
  got <- utils::read.delim("contrast/AvB_vs_AvC_global_q_vs_q_dnds_contrast.tsv")
  a <- utils::read.delim("example_AvB/example_AvB_dnds_annot.tsv")
  b <- utils::read.delim("example_AvC/example_AvC_dnds_annot.tsv")
  shared <- merge(a, b, by = "query_id")
  shared <- shared[is.finite(shared$dNdS.x) & is.finite(shared$dNdS.y) & shared$dNdS.x < 10 & shared$dNdS.y < 10, ]
  delta <- shared$dNdS.x - shared$dNdS.y
  expect_equal(got$n, nrow(shared))
  expect_equal(got$mean_delta, mean(delta), tolerance = 1e-8)
  expect_gt(got$mean_delta, 0) # genes with excess nonsynonymous changes against B only
  expect_equal(got$wilcox_p_value, stats::wilcox.test(delta)$p.value, tolerance = 1e-8)
  expect_equal(got$n_pos_A, sum(shared$dNdS.x > 1))
  expect_equal(got$n_pos_B, sum(shared$dNdS.y > 1))
})
