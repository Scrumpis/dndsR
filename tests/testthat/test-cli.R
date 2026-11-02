parse <- function(...) dndsR:::parse_dnds_opts(args = c(...))

test_that("values keep their meaning through parsing", {
  expect_identical(parse("--comparison-name", "007")$comparison_name, "007")
  expect_identical(parse("--comparison-name", "NA")$comparison_name, "NA")
  expect_identical(parse("--query-fasta", "a,b.fa")$query_fasta, "a,b.fa")
  expect_identical(parse("--filter-expr", "x %in% c(1,2)")$filter_expr, "x %in% c(1,2)")
  expect_identical(parse("--fdr-method", "none")$fdr_method, "none")
  expect_identical(parse("--sides", "query,subject")$sides, c("query", "subject"))
})

test_that("numbers, including negative ones, are read as numbers", {
  expect_identical(parse("--x-axis-min", "-2")$x_axis_min, -2L)
  expect_identical(parse("--x-axis-min=-2")$x_axis_min, -2L)
  expect_identical(parse("--max-dnds", "Inf")$max_dnds, Inf)
  expect_identical(parse("--max-dnds", "1.5")$max_dnds, 1.5)
  expect_identical(parse("--overwrite", "true")$overwrite, TRUE)
})

test_that("short and long forms are equivalent", {
  expect_identical(parse("-c", "f.txt")$comparison_file, "f.txt")
  expect_identical(parse("--comparison_file", "f.txt")$comparison_file, "f.txt")
  expect_identical(parse("-t", "4")$threads, 4L)
  expect_identical(parse("-o", "out")$output_dir, "out")
  expect_true(parse("-v", "-c", "f.txt")$verbose)
  expect_true(parse("--overwrite", "-t", "2")$overwrite)
})

test_that("an equals sign inside a value does not split it", {
  expect_identical(parse("--filter-expr", "x==1")$filter_expr, "x==1")
  expect_identical(parse("--filter-expr=x==1")$filter_expr, "x==1")
})

test_that("unknown commands and options are errors", {
  expect_error(cli_main(c("nonsense")), "Unknown command")
  expect_error(cli_main(c("nonsense", "--help")), "Unknown command")
  expect_error(cli_main(c("dnds_contrasts", "-c", "x")), "Did you mean: dnds_contrast")
  expect_error(cli_main(c("extract_cds", "--overwite", "-c", "x")), "--overwite")
})

test_that("failed workers stop the run with their messages", {
  skip_on_os("windows")
  expect_error(
    .parallel_lapply(1:3, function(i) if (i == 2L) stop("bad input") else i, threads = 2L, what = "test"),
    "1 of 3 jobs failed.*bad input"
  )
  expect_equal(.parallel_lapply(1:3, function(i) i * 2L, threads = 2L), list(2L, 4L, 6L))
  expect_equal(.parallel_lapply(1:3, function(i) i, threads = 1L), list(1L, 2L, 3L))
})

test_that("a worker that dies without a result stops the run", {
  skip_on_os("windows")
  expect_error(
    .parallel_lapply(1:2, function(i) if (i == 1L) tools::pskill(Sys.getpid(), tools::SIGKILL) else i, threads = 2L, what = "test"),
    "died without returning"
  )
})

test_that("logs record errors", {
  log <- withr::local_tempfile(fileext = ".log")
  expect_error(.dndsr_with_log(log, tag = "t", expr = stop("the reason")), "the reason")
  expect_match(paste(readLines(log), collapse = "\n"), "ERROR: the reason", fixed = TRUE)
})

test_that("the command line exits nonzero on failure", {
  run <- function(...) {
    suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c("-e", shQuote("dndsR::cli_main()"), "--", ...), stdout = FALSE, stderr = FALSE))
  }
  expect_gt(run("nonsense"), 0L)
  expect_gt(run("calculate_dnds", "-c", "does-not-exist.txt"), 0L)
  expect_equal(run("__commands"), 0L)
})
