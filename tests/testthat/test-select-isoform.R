test_that("the longest isoform is chosen, whatever the input order", {
  got <- .select_longest_isoform(
    tx = c("g1.b", "g1.a", "g2.a"), gene = c("g1", "g1", "g2"), length_nt = c(900L, 1200L, 300L)
  )
  expect_equal(got$tx, c("g1.a", "g2.a"))
  expect_equal(got$n_isoforms, c(2L, 1L))
  shuffled <- .select_longest_isoform(
    tx = c("g2.a", "g1.a", "g1.b"), gene = c("g2", "g1", "g1"), length_nt = c(300L, 1200L, 900L)
  )
  expect_equal(shuffled, got)
})

test_that("a length divisible by 3 beats a longer one that is not", {
  got <- .select_longest_isoform(c("a", "b"), c("g", "g"), c(301L, 300L))
  expect_equal(got$tx, "b")
})

test_that("ties are broken by transcript name", {
  expect_equal(.select_longest_isoform(c("b", "a"), c("g", "g"), c(300L, 300L))$tx, "a")
})

test_that("an invalid longest isoform loses to a shorter valid one", {
  got <- .select_longest_isoform(c("long", "short"), c("g", "g"), c(951L, 906L), valid = c(FALSE, TRUE))
  expect_equal(got$tx, "short")
})

test_that("a gene with no valid isoform still gets a row, so it can be reported", {
  got <- .select_longest_isoform(c("a", "b"), c("g", "g"), c(951L, 906L), valid = c(FALSE, FALSE))
  expect_equal(got$tx, "a")
})
