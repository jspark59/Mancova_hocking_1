## rohwer_example.R
## Reproduces Section 5 and Table 1 of
##   "Variable selection theorem for the multivariate analysis of covariance model"
##   by S.-H. Yoon, P. Busababodhin, and J.-S. Park.
## Requires: heplots (data), car (MANOVA tests), MASS (generalized inverse).
## Run: source("R/rohwer_example.R")

library(heplots); library(car); library(MASS)
data(Rohwer, package = "heplots")

Y  <- as.matrix(Rohwer[, c("SAT", "PPVT", "Raven")])
Lo <- as.numeric(Rohwer$SES == "Lo")
Rohwer$SESna <- Lo * Rohwer$na                       # x_Lo o na

## ---------------- 1. Model checks and the test of H0: B_r = 0 -------------
fit_int <- lm(Y ~ SES * (n + s + ns + na + ss), data = Rohwer)
fit_cs  <- lm(Y ~ SES + n + s + ns + na + ss, data = Rohwer)
fit_full <- lm(Y ~ SES + n + s + ns + na + ss + SESna, data = Rohwer)
fit_red  <- lm(Y ~ SES + s + ns + na + ss + SESna, data = Rohwer)   # omit n

cat("=== Five SES-by-covariate interactions (Pillai) ===\n")
print(anova(fit_cs, fit_int))
cat("\n=== Interactions tested separately (Type II) ===\n")
print(Anova(fit_int))
cat("\n=== Covariates in the common-slope model (Type II) ===\n")
print(Anova(fit_cs))
cat("\n=== Remaining four interactions after adding SES:na ===\n")
print(anova(fit_full, fit_int))
cat("\n=== Test of H0: B_r = 0 (omitting n) ===\n")
print(anova(fit_red, fit_full))

## ---------------- 2. Matrices of the paper --------------------------------
X  <- cbind(1, Hi = 1 - Lo, Lo = Lo)                  # n x 3, rank 2
Zp <- with(Rohwer, cbind(s, ns, na, SESna, ss))       # retained covariates
Zr <- with(Rohwer, cbind(n))                          # omitted covariate
Z  <- cbind(Zp, Zr)
nobs <- nrow(Y); d <- ncol(Y)
p <- ncol(Zp); r <- ncol(Zr); s <- qr(X)$rank; q <- p + r
f <- nobs - s - q                                     # 61

XtXg <- ginv(crossprod(X))
Q   <- diag(nobs) - X %*% XtXg %*% t(X)
A   <- solve(t(Z) %*% Q %*% Z)                        # A = (Z'QZ)^{-1}
A11 <- A[1:p, 1:p]
A22 <- A[(p + 1):q, (p + 1):q, drop = FALSE]
Mp  <- t(Zp) %*% Q %*% Zp
D   <- solve(Mp, t(Zp) %*% Q %*% Zr)
Bhat   <- A %*% t(Z) %*% Q %*% Y
Bhat_r <- Bhat[(p + 1):q, , drop = FALSE]
R   <- Q - Q %*% Z %*% A %*% t(Z) %*% Q               # residual projection, full
Sig_hat <- t(Y) %*% R %*% Y / f                       # unbiased estimator of Sigma
trS <- sum(diag(Sig_hat))
Hm  <- t(Bhat_r) %*% solve(A22, Bhat_r)               # hypothesis matrix

## ---------------- 3. Estimates of tr(Delta) -------------------------------
tSH  <- sum(diag(solve(Sig_hat, Hm)))                 # tr(Sigma_hat^{-1} H)
tr_u <- (f - d - 1) / f * tSH - r * d                 # unbiased moment estimator
tr_pos <- max(0, tr_u)                                # positive-part estimator
## MLE for r = 1, d = 3, with Sigma replaced by Sigma_hat:
## log L(w) = -w/2 + log{ sinh(sqrt(w t)) / sqrt(w t) },  t = tr(Sigma_hat^{-1} H)
loglik <- function(w) {
  x <- sqrt(w * tSH)
  -w / 2 + x + log1p(-exp(-2 * x)) - log(2 * x)
}
tr_mle <- optimize(loglik, c(1e-10, 50), maximum = TRUE)$maximum
cat("\n=== Estimates of tr(Delta) ===\n")
print(round(c(f = f, trSH = tSH, unbiased = tr_u,
              positive_part = tr_pos, MLE = tr_mle), 4))

## ---------------- 4. Table 1: traces of the MSE matrices ------------------
## psi: Lo - Hi difference at na = mean(na); prediction: new Lo child at means
vvec <- function(x0, z0) {
  a <- drop(X %*% XtXg %*% x0); u <- drop(z0 - t(Z) %*% a)
  list(a = a, u = u, v = drop(t(D) %*% u[1:p] - u[(p + 1):q]))
}
abar <- mean(Rohwer$na)
zbar <- colMeans(Z); zbar["SESna"] <- abar
covtr <- function(x0, z0, add = 0) {                  # traces of covariances
  w <- vvec(x0, z0); up <- w$u[1:p]
  list(F = (sum(w$a^2) + drop(t(w$u) %*% A %*% w$u)) * trS + add,
       R = (sum(w$a^2) + drop(t(up) %*% solve(Mp, up))) * trS + add,
       v = w$v)
}
cPsi  <- covtr(c(0, -1, 1), c(0, 0, 0, abar, 0, 0))
cPred <- covtr(c(1, 0, 1), zbar, add = trS)           # new observation adds Sigma
bpF <- trS * sum(diag(A11)); bpR <- trS * sum(diag(solve(Mp)))

tab1 <- t(sapply(c(0, 0.5, 1, 1.5, 2, 3, 4), function(tD) {
  Br <- sqrt(tD / tSH) * Bhat_r                       # B_r = c * Bhat_r
  c(trDelta = tD,
    Bp_F = bpF,       Bp_R = bpR + sum((D %*% Br)^2),
    psi_F = cPsi$F,   psi_R = cPsi$R + sum((t(Br) %*% cPsi$v)^2),
    pred_F = cPred$F, pred_R = cPred$R + sum((t(Br) %*% cPred$v)^2))
}))
cat("\n=== Table 1 ===\n")
print(round(tab1, 3))

equal_traces <- tSH * c(
  Bp   = (bpF - bpR) / sum((D %*% Bhat_r)^2),
  psi  = (cPsi$F - cPsi$R) / sum((t(Bhat_r) %*% cPsi$v)^2),
  pred = (cPred$F - cPred$R) / sum((t(Bhat_r) %*% cPred$v)^2))
cat("\n=== tr(Delta) at which the two traces are equal ===\n")
print(round(equal_traces, 3))

## ---------------- 5. Estimators of Sigma ----------------------------------
Sig_til <- t(Y) %*% (Q - Q %*% Zp %*% solve(Mp, t(Zp)) %*% Q) %*% Y / (nobs - s - p)
cat("\n=== Diagonal elements of Sigma_hat and Sigma_tilde ===\n")
print(round(rbind(full = diag(Sig_hat), reduced = diag(Sig_til)), 3))
