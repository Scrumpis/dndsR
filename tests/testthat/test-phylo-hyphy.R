# The HyPhy engine layer, without HyPhy: real BUSTED JSON from HyPhy 2.5.103 (paths
# removed, see fixtures/hyphy/README.md), and fake executables for each failure.
# Expected values were read from the JSON with Python, not with dndsR.
fixture <- function(name) test_path("fixtures", "hyphy", name)

test_that("BUSTED results are read from real HyPhy JSON with selection evidence", {
  m <- .hyphy_methods$busted
  j <- .hyphy_read_json(fixture("busted_2.5.103_OG0001.json"), m$required)
  v <- m$parse(j)
  expect_identical(v$lrt, 7.307789067929662)
  expect_identical(v$p_value, 0.01294505114975508)
  expect_true(v$constrained_fitted)
  expect_identical(v$lnl_unconstrained, -2548.822260013832)
  expect_identical(v$lnl_constrained, -2552.476154547797)
  expect_identical(c(v$fg_omega_1, v$fg_omega_2, v$fg_omega_3), c(0.9926780342396962, 0, 4.508709733549398))
  expect_identical(c(v$fg_prop_1, v$fg_prop_2, v$fg_prop_3), c(0, 0.5146089747894831, 0.4853910252105169))
  expect_equal(v$tree_length, 0.29896195644701284, tolerance = 1e-12)
  expect_identical(.hyphy_tested_branches(j), c("Node2", "spA", "spB"))
})

test_that("a BUSTED run without evidence has LRT 0, p 0.5, and no constrained fit", {
  m <- .hyphy_methods$busted
  v <- m$parse(.hyphy_read_json(fixture("busted_2.5.103_OG0002.json"), m$required))
  expect_identical(v$lrt, 0)
  expect_identical(v$p_value, 0.5)
  expect_false(v$constrained_fitted)
  expect_true(is.na(v$lnl_constrained))
  # HyPhy's p-value: 0.5 * (1 - pchisq(LRT, 2)), a 50:50 mixture of chi-square 0 and 2.
  expect_equal(0.5 * stats::pchisq(7.307789067929662, 2, lower.tail = FALSE), 0.01294505114975508, tolerance = 1e-10)
})

test_that("missing or broken JSON is an error, not a result", {
  dir <- withr::local_tempdir()
  writeLines('{"test results": {"LRT": 1', file.path(dir, "bad.json"))
  expect_error(.hyphy_read_json(file.path(dir, "bad.json"), "fits"), "could not be parsed")
  writeLines('{"fits": {}}', file.path(dir, "partial.json"))
  expect_error(.hyphy_read_json(file.path(dir, "partial.json"), c("fits", "tested")), "lacks: tested")
  expect_error(.hyphy_read_json(file.path(dir, "none.json"), "fits"), "no JSON")
})

test_that("settings default to HyPhy's, and invalid ones stop", {
  s <- .hyphy_settings("busted", list())
  expect_identical(s[c("srv", "multiple_hits", "rates", "syn_rates", "error_sink", "grid_size", "starting_points")],
                   list(srv = "Yes", multiple_hits = "None", rates = 3L, syn_rates = 3L, error_sink = "No",
                        grid_size = 250L, starting_points = 1L))
  expect_identical(.hyphy_settings("busted", list(multiple_hits = "double+triple"))$multiple_hits, "Double+Triple")
  expect_error(.hyphy_settings("busted", list(srv = "maybe")), "Invalid srv")
  expect_error(.hyphy_settings("busted", list(rates = 0)), "Invalid rates")
  expect_error(.hyphy_settings("busted", list(genetic_code = "Vertebrate mtDNA")), "Invalid genetic_code")
})

test_that("every BUSTED setting is passed to HyPhy explicitly", {
  a <- .hyphy_methods$busted$args(.hyphy_settings("busted", list()), "a.fas", "t.nwk", "o.json")
  expect_identical(a, c("--alignment", "'a.fas'", "--tree", "'t.nwk'", "--branches", "FG", "--code", "Universal",
                        "--srv", "Yes", "--rates", "3", "--syn-rates", "3", "--multiple-hits", "'None'",
                        "--error-sink", "No", "--grid-size", "250", "--starting-points", "1", "--kill-zero-lengths", "Yes",
                        "--output", "'o.json'"))
  no_srv <- .hyphy_methods$busted$args(.hyphy_settings("busted", list(srv = "No")), "a.fas", "t.nwk", "o.json")
  expect_false("--syn-rates" %in% no_srv)
})

# A stand-in for HyPhy. FAKE_MODE: error, truncated, mismatch, or copy (writes a fixture).
fake_hyphy <- function(dir, mode, version = "2.5.103") {
  path <- file.path(dir, paste0("hyphy-", mode))
  writeLines(c(
    "#!/bin/sh",
    sprintf('if [ "$1" = --version ]; then echo "HYPHY %s(MP) for Linux"; exit 0; fi', version),
    'out=""; prev=""; for a in "$@"; do [ "$prev" = --output ] && out=$a; prev=$a; done',
    'echo "$@" > args.txt',
    switch(mode,
      error = c('printf "Error:## ASSERTION FAILED\\nThe input alignment must not contain stop codons\\n" > errors.log', "exit 1"),
      truncated = c('printf "{\\"test results\\": {" > "$out"', "exit 0"),
      mismatch = c(sprintf('sed "s/\\"spB\\":\\"test\\"/\\"spB\\":\\"background\\"/" "%s" > "$out"',
                           normalizePath(fixture("busted_2.5.103_OG0001.json"))), "exit 0"),
      copy = c(sprintf('cp "%s" "$out"', normalizePath(fixture("busted_2.5.103_OG0001.json"))), "exit 0"),
      # HyPhy deleted the clade's stem (dndsR's Node3) as zero-length, and tested the rest.
      zerolen = c("echo '### Deleted 1 zero-length internal branches: `Node3`'",
                  sprintf('sed "s/\\"Node2\\":\\"test\\"/\\"Node2\\":\\"background\\"/" "%s" > "$out"', normalizePath(fixture("busted_2.5.103_OG0001.json"))), "exit 0"),
      # HyPhy deleted a background branch; the tested branches are the requested ones.
      bgzero = c("echo '### Deleted 1 zero-length internal branches: `Node4`'",
                 sprintf('sed "s/\\"Node2\\":\\"test\\"/\\"Node3\\":\\"test\\"/" "%s" > "$out"', normalizePath(fixture("busted_2.5.103_OG0001.json"))), "exit 0")
    )
  ), path)
  Sys.chmod(path, "0755")
  path
}

# One orthogroup whose labelled branches are spA, spB, and Node2, as in the fixture.
fake_inputs <- function(dir, scope = "terminal") {
  ex <- system.file("extdata", "phylo_example", package = "dndsR", mustWork = TRUE)
  dir.create(file.path(dir, "aln"))
  file.copy(file.path(ex, "alignments", "OG0001.fasta"), file.path(dir, "aln"))
  # The fixture's tree named the internal node Node2; this one-row hypothesis keeps the
  # stem of spA and spB as dndsR's Node3, so "copy" mode reports a label mismatch too.
  writeLines("((((spA,spB),(spC,spD)),(spE,spF)),(spG,spH));", file.path(dir, "tree.nwk"))
  writeLines(c("scope\ttaxa", paste0(scope, "\tspA,spB")), file.path(dir, "fg.tsv"))
}

run_fake <- function(mode, version = "2.5.103", ..., scope = "terminal") {
  skip_on_os("windows")
  skip_if_not_installed("ape")
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  fake_inputs(dir, scope)
  res <- phylo_busted(alignment_dir = file.path(dir, "aln"), species_tree_file = file.path(dir, "tree.nwk"),
                      foreground_file = file.path(dir, "fg.tsv"), output_dir = dir, hyphy = fake_hyphy(dir, mode, version), ...)
  list(res = res, dir = dir)
}

test_that("a HyPhy error marks the orthogroup failed, keeps the message, and leaves the family", {
  expect_warning(r <- run_fake("error"), "failed and are not in the multiple-testing family")
  expect_identical(r$res$status, "failed")
  expect_identical(r$res$reason, "hyphy_error")
  expect_identical(r$res$exclusion_basis, "software")
  expect_identical(r$res$stage, "analysis")
  expect_match(r$res$detail, "must not contain stop codons")
  expect_identical(r$res$family_size, 0L)
  expect_true(is.na(r$res$p_adj))
})

test_that("truncated HyPhy output is a failure, not a result", {
  expect_warning(r <- run_fake("truncated"), "failed")
  expect_identical(r$res$reason, "hyphy_output_invalid")
})

test_that("a difference between labelled and tested branches is a failure", {
  expect_warning(r <- run_fake("mismatch"), "failed")
  expect_identical(r$res$reason, "label_mismatch")
  expect_match(r$res$detail, "dndsR labelled spA,spB")
})

test_that("HyPhy is called with the seed, thread count, and settings, and they are recorded", {
  expect_warning(r <- run_fake("copy"), "failed") # the fixture tested Node2 too
  d <- file.path(r$dir, "busted")
  args <- readLines(file.path(d, "hyphy", "OG0001", "args.txt"))
  expect_match(args, "^CPU=1 ENV=RANDOM_SEED=1 busted --alignment alignment.fas --tree tree.nwk --branches FG", fixed = FALSE)
  prov <- jsonlite::read_json(file.path(d, "busted_provenance.json"))
  expect_identical(prov$random_seed, 1L)
  expect_identical(prov$hyphy_cpu, 1L)
  expect_identical(prov$hyphy_version, "2.5.103")
  expect_identical(prov$settings$srv, "Yes")
  expect_identical(prov$species_tree$md5, unname(tools::md5sum(file.path(r$dir, "tree.nwk"))))
  expect_identical(readLines(file.path(d, "hyphy", "OG0001", "tree.nwk")),
                   "(((spA{FG},spB{FG})Node3,(spC,spD)Node4)Node2,(spE,spF)Node5,(spG,spH)Node6);")
})

test_that("an extra start that fails makes the orthogroup fail, naming the seed", {
  expect_warning(r <- run_fake("error", extra_seeds = 2:3), "failed")
  expect_identical(r$res$reason, "hyphy_error")
  expect_error(suppressWarnings(run_fake("error", extra_seeds = c(1L, 2L))), "different from seed")
})

test_that("diagnostics across starts use each model's best likelihood, and never replace the primary", {
  v <- function(lrt, unc, con) list(values = list(lrt = lrt, lnl_unconstrained = unc, lnl_constrained = con))
  # Primary: LRT 2 (unc -100, con -101). Start 2: better alternative (-99.5), no constrained fit.
  # Start 3: better null (-100.8). Best alternative -99.5; best null max(-101, -99.5, -100.8) = -99.5.
  d <- .phylo_start_diagnostics(list(v(2, -100, -101), v(0, -99.5, NA), v(1.6, -100, -100.8)), c(1L, 2L, 3L))
  expect_identical(d$n_starts, 3L)
  expect_identical(c(d$lrt_min_starts, d$lrt_max_starts), c(0, 2))
  expect_identical(d$lnl_unconstrained_best, -99.5)
  expect_identical(d$lnl_null_best, -99.5)
  expect_identical(d$lrt_best_fits, 0)
  expect_identical(d$p_best_fits, 0.5)
  expect_false(d$primary_unconstrained_is_best)
  d2 <- .phylo_start_diagnostics(list(v(2, -100, -101), v(1, -100.2, -100.7)), c(1L, 2L))
  expect_equal(d2$lrt_best_fits, 2 * (-100 - -100.7))
  expect_true(d2$primary_unconstrained_is_best)
  expect_false(d2$primary_null_is_best)
})

test_that("a requested branch HyPhy deletes as zero-length excludes the orthogroup, by policy", {
  r <- run_fake("zerolen", scope = "clade")
  expect_identical(r$res$status, "excluded")
  expect_identical(r$res$reason, "foreground_branch_zero_length")
  expect_identical(r$res$exclusion_basis, "policy")
  expect_match(r$res$detail, "deleted requested branch\\(es\\) Node3")
  expect_identical(r$res$family_size, 0L)
})

test_that("a background branch HyPhy deletes as zero-length is flagged, and the gene tested", {
  r <- run_fake("bgzero", scope = "clade")
  expect_identical(r$res$status, "tested")
  expect_match(r$res$qc_flags, "hyphy_deleted_zero_length_branches:1")
})

test_that("an unvalidated HyPhy version is announced", {
  expect_message(suppressWarnings(run_fake("error", version = "2.5.90")), "HyPhy 2.5.90 has not been checked")
})

test_that("a missing HyPhy stops with installation advice", {
  withr::local_options(dndsR.hyphy = NULL)
  withr::local_envvar(DNDSR_HYPHY = "", PATH = "/nonexistent")
  expect_error(.hyphy_path(), "HyPhy was not found.*bioconda")
})

test_that("the CLI finds the phylogenetic commands and keeps paths and choices as typed", {
  expect_true(all(c("phylo_prepare", "phylo_busted", "phylo_example", "verify_phylo_example") %in% cli_list_commands()))
  p <- parse_dnds_opts(args = c("--species-tree-file", "NA", "--multiple-hits", "Double+Triple", "--srv", "No", "--seed", "7",
                                "--extra-seeds", "2,3", "--prepared-dir", "prepared"))
  expect_identical(p$extra_seeds, c("2", "3"))
  expect_identical(p$prepared_dir, "prepared")
  expect_identical(p$species_tree_file, "NA")
  expect_identical(p$multiple_hits, "Double+Triple")
  expect_identical(p$srv, "No")
  expect_identical(p$seed, 7L)
  expect_error(cli_main(c("phylo_busted", "--foregrund-file", "x")), "--foregrund-file")
})
