## Shared test scaffolding.
##
## Convention borrowed from the parent project's check_*.R scripts: report the
## measured numbers in a table AND assert a tolerance where one is genuinely
## defensible. A test that only prints is a diagnostic; a test that only
## asserts hides the magnitude of what it checked. These do both.

.TEST_RESULTS <- new.env(parent = emptyenv())
.TEST_RESULTS$rows <- list()

test_init <- function(name) {
  cat("\n", strrep("=", 74), "\n", name, "\n", strrep("=", 74), "\n", sep = "")
  .TEST_RESULTS$current <- name
}

## check(): record and report one assertion.
##   value    -- the measured quantity
##   op       -- "<", "<=", "==", "TRUE"
##   bound    -- threshold (ignored for "TRUE")
check <- function(label, value, op = "<", bound = NA, fmt = "%.3e", note = "") {
  ok <- switch(op,
               "<"    = isTRUE(all(value <  bound)),
               "<="   = isTRUE(all(value <= bound)),
               ">"    = isTRUE(all(value >  bound)),
               ">="   = isTRUE(all(value >= bound)),
               "=="   = isTRUE(all(value == bound)),
               "TRUE" = isTRUE(all(value)),
               stop("unknown op"))
  ## For lower-bound checks the interesting extreme is the minimum.
  rep_val <- if (op %in% c(">", ">=")) min(value) else max(value)
  vs <- if (op == "TRUE") as.character(all(value)) else sprintf(fmt, rep_val)
  bs <- if (op == "TRUE") "" else sprintf(paste0(" ", op, " ", fmt), bound)
  cat(sprintf("  [%s] %-52s %s%s%s\n", if (ok) "PASS" else "FAIL", label, vs, bs,
              if (nzchar(note)) paste0("   (", note, ")") else ""))
  .TEST_RESULTS$rows[[length(.TEST_RESULTS$rows) + 1]] <-
    data.frame(test = .TEST_RESULTS$current, check = label,
               value = if (op == "TRUE") NA_real_ else rep_val,
               bound = if (op == "TRUE") NA_real_ else bound,
               pass = ok, stringsAsFactors = FALSE)
  invisible(ok)
}

test_summary <- function() {
  if (!length(.TEST_RESULTS$rows)) { cat("\nNo checks recorded.\n"); return(invisible(TRUE)) }
  df <- do.call(rbind, .TEST_RESULTS$rows)
  n_fail <- sum(!df$pass)
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat(sprintf("TEST SUMMARY: %d checks, %d passed, %d FAILED\n",
              nrow(df), sum(df$pass), n_fail))
  if (n_fail > 0) {
    cat("\nFailures:\n")
    print(df[!df$pass, c("test", "check", "value", "bound")], row.names = FALSE)
  }
  cat(strrep("=", 74), "\n")
  invisible(n_fail == 0)
}

test_reset <- function() { .TEST_RESULTS$rows <- list(); invisible(NULL) }
