# X2 deterministic Aalen–Johansen cumulative incidence (tie-aware) via survival::survfit multistate. Base `survival` only.
args <- commandArgs(trailingOnly = TRUE); out <- args[1]
suppressMessages(library(survival))
d <- read.csv(file.path(out, "_aj_input.csv"))
d$ev <- factor(d$status, levels = 0:2, labels = c("censor", "recurrence", "death"))
res <- list(); curves <- list()
for (g in c("all", "stage II", "stage III")) {
  s <- if (g == "all") d else d[d$group == g, ]
  f <- survfit(Surv(time, ev) ~ 1, data = s)
  for (cause in c("recurrence", "death")) {
    j <- match(cause, f$states)
    for (h in c(36, 60)) {
      sm <- summary(f, times = h)
      res[[length(res) + 1]] <- data.frame(group = g, cause = cause, horizon = h, cif = sm$pstate[, j],
                                           lo = sm$lower[, j], hi = sm$upper[, j])
    }
    curves[[length(curves) + 1]] <- data.frame(group = g, cause = cause, time = f$time, cif = f$pstate[, j],
                                               lo = f$lower[, j], hi = f$upper[, j])
  }
}
write.csv(do.call(rbind, res), file.path(out, "_aj_table.csv"), row.names = FALSE)
write.csv(do.call(rbind, curves), file.path(out, "_aj_curves.csv"), row.names = FALSE)
cat("ok\n")
