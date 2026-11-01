# Runs every R code block in docs/tutorials/using-r.md, in order, and stops on any error.
#   Rscript tools/dev/test-r-tutorial.R <path to using-r.md>
lines <- readLines(commandArgs(TRUE)[1])
blocks <- list(); inside <- FALSE; current <- character()
for (line in lines) {
  if (line == "```r") { inside <- TRUE; current <- character() }
  else if (line == "```" && inside) { inside <- FALSE; blocks[[length(blocks) + 1L]] <- current }
  else if (inside) current <- c(current, line)
}
stopifnot(length(blocks) >= 3L)
for (i in seq_along(blocks)) {
  invisible(utils::capture.output(suppressMessages(eval(parse(text = blocks[[i]]), envir = globalenv()))))
  cat("R tutorial block", i, "ok\n")
}
cat("R tutorial OK:", length(blocks), "blocks ran\n")
