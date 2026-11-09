# X2 Fine–Gray subdistribution hazard models (registered plan §X2). Base `survival` package only.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
suppressMessages(library(survival))
d <- read.csv(file.path(out, "_design_for_R.csv"))
d$ev <- factor(d$status, levels = 0:2, labels = c("censor", "recurrence", "death"))
clin <- c("age_rcs1", "age_rcs2", "sex_Male", "tnm_stage_3_0", "tumor_location_proximal", "chemotherapy_adjuvant_Y")
scores <- c("EMT", "E2F", "Wnt", "p53")
fit_fg <- function(vars, cause, label) {
  fg <- finegray(as.formula(paste("Surv(time, ev) ~", paste(vars, collapse = "+"))), data = d, etype = cause)
  m <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", paste(vars, collapse = "+"))),
             data = fg, weight = fgwt)
  s <- summary(m)$conf.int
  data.frame(model = label, term = rownames(s), sHR = s[, 1], lo = s[, 3], hi = s[, 4],
             p = summary(m)$coefficients[, ncol(summary(m)$coefficients)], row.names = NULL)
}
res <- rbind(
  fit_fg(clin, "recurrence", "FG recurrence: clinical"),
  fit_fg(c(clin, scores), "recurrence", "FG recurrence: clinical + 4 scores"),
  fit_fg(clin, "death", "FG death w/o recurrence: clinical"),
  fit_fg(c(clin, scores), "death", "FG death w/o recurrence: clinical + 4 scores"))
write.csv(res, file.path(out, "fine_gray.csv"), row.names = FALSE)
cat(R.version.string, "\n")
cat("survival", as.character(packageVersion("survival")), "\n")
cat("rows", nrow(res), "\n")
