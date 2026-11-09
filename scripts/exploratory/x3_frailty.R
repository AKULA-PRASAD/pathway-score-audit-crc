# X3 batch frailty sensitivity models (registered plan §X3). Base `survival` package only.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
suppressMessages(library(survival))
d <- read.csv(file.path(out, "_design_for_R.csv"))
d$batch <- factor(d$batch)
vars <- c("age_rcs1", "age_rcs2", "sex_Male", "tnm_stage_3_0", "tumor_location_proximal",
          "chemotherapy_adjuvant_Y", "EMT", "E2F", "Wnt", "p53")
rhs <- paste(vars, collapse = " + ")
coefs <- function(m) {
  b <- coef(m)[vars]; se <- sqrt(diag(vcov(m)))[match(vars, names(coef(m)))]
  lapply(setNames(seq_along(vars), vars), function(i) list(beta = unname(b[i]), se = unname(se[i]),
                                                          HR = unname(exp(b[i]))))
}
m0 <- coxph(as.formula(paste("Surv(time, event) ~", rhs)), data = d)
# gamma frailty, theta estimated (default method selects theta by marginal likelihood)
mg <- coxph(as.formula(paste("Surv(time, event) ~", rhs, "+ frailty(batch, distribution = 'gamma')")), data = d)
theta_hat <- mg$history[[1]]$theta
# profile marginal (corrected) log-likelihood over a theta grid
grid <- sort(unique(c(1e-6, 0.005, 0.01, seq(0.015, 0.05, by = 0.0025), 0.075, 0.1, 0.13, 0.15, 0.2, 0.3, seq(0.35, 0.8, by = 0.025), 1.0)))  # refined after independent review (deviation 6)
prof <- sapply(grid, function(th) {
  m <- coxph(as.formula(paste("Surv(time, event) ~", rhs, "+ frailty(batch, distribution = 'gamma', theta =", th, ")")), data = d)
  m$history[[1]]$c.loglik
})
ll_max <- max(prof, mg$history[[1]]$c.loglik)
inside <- grid[2 * (ll_max - prof) <= qchisq(0.95, 1)]
lrt <- 2 * (mg$history[[1]]$c.loglik - m0$loglik[2])
# gaussian frailty check
mn <- coxph(as.formula(paste("Surv(time, event) ~", rhs, "+ frailty(batch, distribution = 'gaussian')")), data = d)
res <- list(
  no_frailty = list(loglik = m0$loglik[2], coefficients = coefs(m0)),
  gamma = list(theta_hat = theta_hat, marginal_loglik = mg$history[[1]]$c.loglik,
               profile_grid = grid, profile_loglik = prof,
               theta_95_profile_interval_on_grid = c(min(inside), max(inside)),
               lrt_vs_no_frailty = lrt,
               lrt_p_boundary_half_chisq = 0.5 * pchisq(max(lrt, 0), 1, lower.tail = FALSE),
               coefficients = coefs(mg),
               batch_frailties = setNames(as.list(exp(mg$frail)), levels(d$batch))),
  gaussian = list(variance_hat = mn$history[[1]]$theta, coefficients = coefs(mn)),
  software = list(R = R.version.string, survival = as.character(packageVersion("survival"))))
# minimal JSON writer (no extra packages)
tojson <- function(x) {
  if (is.list(x)) {
    if (!is.null(names(x)) && all(names(x) != "")) {
      paste0("{", paste(sprintf('"%s": %s', names(x), sapply(x, tojson)), collapse = ", "), "}")
    } else paste0("[", paste(sapply(x, tojson), collapse = ", "), "]")
  } else if (is.character(x)) {
    if (length(x) == 1) sprintf('"%s"', x) else paste0("[", paste(sprintf('"%s"', x), collapse = ", "), "]")
  } else {
    v <- ifelse(is.finite(x), formatC(x, digits = 10, format = "g"), "null")
    if (length(x) == 1) v else paste0("[", paste(v, collapse = ", "), "]")
  }
}
writeLines(tojson(res), file.path(out, "frailty.json"))
cat("ok\n")
