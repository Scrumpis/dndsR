# Pilot gates of studies L and B (docs/dev/lineage-summary-study.md, lineage-robustness-study.md),
# as fixed before the pilots ran.
#   Rscript --vanilla phylo-study-pilot-gates.R L <root> <pilot units file> <full units file> <budget core-hours>
#   Rscript --vanilla phylo-study-pilot-gates.R B <root> <pilot units file> <full units file> <budget core-hours>
# Cost: per unit and method, core-seconds per gene measured in the pilot (wall seconds x cpus / genes),
# projected to the full units' genes; BUSTED's cost on units where the pilot did not run it is
# taken from the pilot's error-free unit of the same tree.
a <- commandArgs(TRUE); study <- a[1]; root <- a[2]; pilot <- utils::read.delim(a[3], stringsAsFactors = FALSE)
full <- utils::read.delim(a[4], stringsAsFactors = FALSE); budget <- as.numeric(a[5])
ud <- function(u) file.path(root, "units", u)
gate <- function(name, ok, detail) cat(sprintf("[%s] %s: %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
done <- file.exists(file.path(ud(pilot$unit), "DONE"))
gate("1 units complete", all(done), paste(sum(done), "of", nrow(pilot)))
cost <- list()
for (i in which(done)) {
  u <- pilot$unit[i]; ms <- utils::read.delim(file.path(ud(u), "method_seconds.tsv"), stringsAsFactors = FALSE)
  cpus <- utils::read.delim(file.path(ud(u), "timing.tsv"))$cpus[1]
  ms$per_gene <- ms$seconds * cpus / pilot$n[i]
  cost[[u]] <- data.frame(unit = sub("-pilot$", "", u), tree = pilot$tree[i], method = ms$method, per_gene = ms$per_gene)
}
cost <- do.call(rbind, cost)
proj <- 0
for (i in seq_len(nrow(full))) {
  u <- full$unit[i]; ms <- strsplit(full$methods[i], ",")[[1]]
  runs <- if (full$mask[i] == "yes") c(ms, paste0(ms, "_mask")) else ms
  for (m in runs) {
    pg <- cost$per_gene[cost$unit == u & cost$method == m]
    if (!length(pg)) pg <- cost$per_gene[grepl("^B-E0-", cost$unit) & cost$tree == sub("m$", "", full$tree[i]) & cost$method == m]
    if (!length(pg)) pg <- max(cost$per_gene[cost$method == m]) # no measurement on this tree: the largest measured
    proj <- proj + mean(pg) * full$n[i] / 3600
  }
}
print(stats::aggregate(per_gene ~ method + tree, cost, mean))
gate("cost", proj <= budget, sprintf("projected %.0f core-hours for the full units, budget %.0f", proj, budget))
if (study == "L") {
  ok_tested <- vapply(pilot$unit[done], function(u) { r <- utils::read.delim(file.path(ud(u), "branch_omega", "branch_omega_results.tsv")); mean(r$status == "tested") }, numeric(1))
  gate("1 tested share", all(ok_tested >= 0.95), paste(round(ok_tested, 3), collapse = ", "))
  m <- unlist(lapply(pilot$unit[done & pilot$cell %in% c("LG1", "LA1")], function(u) utils::read.delim(file.path(ud(u), "sim", "genes.tsv"))$gene_effect))
  gate("2 gene effect mean", abs(mean(m) - 1) <= 3 * stats::sd(m) / sqrt(length(m)), sprintf("mean %.3f, se %.3f, n %d", mean(m), stats::sd(m) / sqrt(length(m)), length(m)))
  # Only the simulator's gamma_classes definition (sourcing the file would run phylo-sim.R's 5-argument entry point)
  ex <- parse(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "phylo-lineage-sim2.R"))
  eval(ex[[which(vapply(ex, function(e) startsWith(paste(deparse(e), collapse = ""), "gamma_classes <-"), logical(1)))]])
  gate("2 gamma classes mean", abs(mean(gamma_classes(0.5, 0.2, 8L)) - 0.2) < 1e-9, sprintf("%.12f", mean(gamma_classes(0.5, 0.2, 8L))))
  for (u in pilot$unit[done & pilot$cell == "LA1"]) {
    g <- utils::read.delim(file.path(ud(u), "sim", "genes.tsv"), stringsAsFactors = FALSE); cut <- stats::qlnorm(2 / 3, -0.7^2 / 2, 0.7)
    lost <- grepl("t1", g$dropped)
    gate("2 LA1 loss only above cutoff", all(g$gene_effect[lost] > cut), sprintf("%d of %d genes lost t1 (%.1f%%; planned 16.7%%)", sum(lost), nrow(g), 100 * mean(lost)))
  }
  for (u in pilot$unit[done & pilot$cell == "LS0"]) {
    b <- utils::read.delim(file.path(ud(u), "branch_omega", "branch_omega_branches.tsv"), stringsAsFactors = FALSE)
    b <- b[b$in_species_tree & !grepl("+", b$species_tree_branches, fixed = TRUE) & !grepl("high_ds", b$flags), ]
    po <- tapply(seq_len(nrow(b)), b$species_tree_branches, function(i) (sum(b$NdN_free[i]) / sum(b$N_free[i])) / (sum(b$SdS_free[i]) / sum(b$S_free[i])))
    gate("3 LS0 pooled omega within 10% of 0.2", all(abs(po / 0.2 - 1) <= 0.10), paste(names(po), round(po, 3), collapse = " "))
  }
}
if (study == "B") {
  for (u in pilot$unit[done]) {
    cl <- file.path(ud(u), "sim", "corrupted.tar.gz")
    if (grepl("^B-E[12]-", u)) {
      log <- utils::read.delim(pipe(sprintf("tar -xzOf %s corrupted/corrupt_log.tsv", cl)))
      g <- utils::read.delim(file.path(ud(u), "sim", "genes.tsv"))
      gate(paste("2", u, "one segment per alignment"), nrow(log) == nrow(g) && !anyDuplicated(log$alignment), paste(nrow(log), "segments,", nrow(g), "genes"))
    }
    if (grepl("^B-MIX-", u)) {
      log <- utils::read.delim(pipe(sprintf("tar -xzOf %s corrupted/corrupt_log.tsv", cl)))
      g <- utils::read.delim(file.path(ud(u), "sim", "genes.tsv"))
      cat("  mixture corruption per taxon (planned t1 .032 t2 .020 t3 .056 t4 .083 t5 .154 t6 .056 t7 .083 t8 .032):",
          paste(names(table(log$taxon)), round(as.vector(table(log$taxon)) / nrow(g), 3), collapse = " "), "; selected", sum(g$selected), "of", nrow(g), "\n")
    }
    if (grepl("^B-E3-", u)) {
      r <- utils::read.delim(file.path(ud(u), "absrel", "absrel_results.tsv"))
      gate(paste("2", u, "realigned and tested"), mean(r$status == "tested") >= 0.95, sprintf("%.1f%% tested", 100 * mean(r$status == "tested")))
    }
    if (grepl("^B-E4-", u)) {
      b <- utils::read.delim(file.path(ud(u), "branch_omega", "branch_omega_branches.tsv"), stringsAsFactors = FALSE)
      b <- b[b$terminal %in% TRUE, ]
      ok <- vapply(split(b, b$orthogroup), function(x) x$t_m0[x$branch == "t1"] >= 3 * stats::median(x$t_m0[x$branch != "t1"]), logical(1))
      gate(paste("2", u, "paralog t1 branch >= 3x median"), mean(ok) >= 0.9, sprintf("%.0f%% of %d genes", 100 * mean(ok), length(ok)))
    }
  }
}
