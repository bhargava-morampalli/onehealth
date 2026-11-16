# Reference values for onehealth's statistics, computed with R and the vegan package.
#
#     Rscript tests/data/stats_reference/reference.R
#
# Reads cases.tsv (same directory) and writes reference.tsv. Nothing here uses
# onehealth: each value comes from a function of base R or vegan that statisticians
# and ecologists already rely on, or from the definition of the test itself.
#
#   Simpson's index of diversity     vegan::simpson.unb  (Hunter & Gaston's unbiased form)
#   its sampling variance            simulation: the index of 2,000,000 multinomial samples
#   rarefaction                      vegan::rarefy
#   Wilson interval                  stats::prop.test(correct = FALSE)
#   Benjamini-Hochberg q-values      stats::p.adjust(method = "BH")
#   chi-squared statistic            stats::chisq.test(correct = FALSE)
#   its permutation p-value          stats::chisq.test(simulate.p.value = TRUE, B = 1,000,000)
#   lineage enrichment p-value       exact, from the hypergeometric distribution (stats::dhyper)
#
# The simulations use a fixed seed, so a rerun with the same R gives the same file.
# vegan warns that counts without singletons may not be counts of individuals; they are
# here, so that warning is silenced.

suppressPackageStartupMessages(library(vegan))
set.seed(20261116)

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
cases <- read.delim(file.path(here, "cases.tsv"), comment.char = "#", quote = "", stringsAsFactors = FALSE)

out <- character(0)
emit <- function(kind, case, quantity, value) {
  out <<- c(out, paste(kind, case, quantity, ifelse(is.na(value), "NA", trimws(formatC(value, digits = 17, format = "g"))), sep = "\t"))
}
numbers <- function(text) suppressWarnings(as.numeric(strsplit(text, ",")[[1]]))
named <- function(text) {
  parts <- strsplit(text, ";")[[1]]
  cut <- regexpr("=[^=]*$", parts)
  setNames(lapply(substring(parts, cut + 1), numbers), substring(parts, 1, cut - 1))
}

for (i in seq_len(nrow(cases))) {
  kind <- cases$kind[i]; case <- cases$case[i]; values <- cases$values[i]

  if (kind == "counts") {
    x <- numbers(values); n <- sum(x)
    emit(kind, case, "simpson", if (n < 2) NA else suppressWarnings(simpson.unb(x)))
    if (n >= 2) {
      # Sampling variance of the index: draw samples of n isolates from a population
      # with these lineage frequencies and take the variance of the index over them.
      draws <- 2e6; chunk <- 2e5; d <- numeric(0)
      for (k in seq_len(draws / chunk)) {
        m <- rmultinom(chunk, n, x / n)
        d <- c(d, 1 - colSums(m * (m - 1)) / (n * (n - 1)))
      }
      v <- var(d)
      emit(kind, case, "simpson_variance", v)
      emit(kind, case, "simpson_variance_se", sd((d - mean(d))^2) / sqrt(draws))
      for (size in sort(unique(pmin(n, pmax(1, c(1, 2, round(n / 4), round(n / 2), n - 1, n)))))) {
        emit(kind, case, paste0("rarefy_", size), as.numeric(suppressWarnings(rarefy(x, size))))
      }
    }
  }

  if (kind == "proportion") {
    x <- numbers(values)
    ci <- suppressWarnings(prop.test(x[1], x[2], correct = FALSE)$conf.int)
    emit(kind, case, "wilson_low", ci[1]); emit(kind, case, "wilson_high", ci[2])
  }

  if (kind == "pvalues") {
    q <- p.adjust(numbers(values), method = "BH")
    for (k in seq_along(q)) emit(kind, case, paste0("q_", k), q[k])
  }

  if (kind == "features") {
    parts <- named(values); sizes <- parts[["sizes"]]
    for (feature in setdiff(names(parts), "sizes")) {
      present <- parts[[feature]]
      table <- rbind(present, sizes - present)
      emit(kind, paste(case, feature, sep = "|"), "chi2", unname(suppressWarnings(chisq.test(table, correct = FALSE)$statistic)))
      emit(kind, paste(case, feature, sep = "|"), "p_simulated", chisq.test(table, simulate.p.value = TRUE, B = 1e6)$p.value)
    }
  }

  if (kind == "lineages") {
    parts <- named(values); table <- do.call(rbind, parts)
    total <- sum(table); source_size <- colSums(table)
    for (lineage in names(parts)) {
      size <- sum(parts[[lineage]])
      if (size < 5) next
      for (s in seq_along(source_size)) {
        if (source_size[s] == 0) next
        expected <- size * source_size[s] / total
        x <- 0:min(size, source_size[s])
        probability <- dhyper(x, size, total - size, source_size[s])
        as_extreme <- abs(x - expected) >= abs(parts[[lineage]][s] - expected) - 1e-9
        key <- paste(case, lineage, s, sep = "|")
        emit(kind, key, "expected", expected)
        emit(kind, key, "p_exact", sum(probability[as_extreme]))
      }
    }
  }
}

header <- c(
  "# Reference values computed by reference.R; do not edit by hand. See README.md.",
  paste0("# ", R.version.string, ", vegan ", as.character(packageVersion("vegan"))),
  "kind\tcase\tquantity\tvalue"
)
writeLines(c(header, out), file.path(here, "reference.tsv"))
cat("wrote", length(out), "values\n")
