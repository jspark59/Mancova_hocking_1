## formula_check.R
## Additional check (not reported in the paper): exact traces of the MSE
## matrices from the formulas of Section 4 for a small simulated one-way
## MANCOVA design, compared with a Monte Carlo estimate at tr(Delta) = 1.
## Requires: MASS. Run: source("R/formula_check.R")

ginv <- MASS::ginv

## ---------------- Design (deterministic, no random numbers) -------------
ns <- c(6, 8, 10)
g  <- rep(0:2, ns)                     # treatment index 0, 1, 2
n  <- length(g); i <- 1:n
X  <- cbind(1, outer(g, 0:2, "==") * 1) # n x 4, rank 3
s  <- qr(X)$rank
z1 <- cos(i) + 0.5 * g
z2 <- (i %% 7) / 7
z3 <- 0.6 * z1 + 0.4 * sin(2 * i) - 0.3 * g
Zp <- cbind(z1, z2); Zr <- cbind(z3); Z <- cbind(Zp, Zr)
p <- ncol(Zp); r <- ncol(Zr); d <- 2
Sig <- matrix(c(1, 0.5, 0.5, 1), 2)

## ---------------- Matrices of the paper ---------------------------------
XtXg <- ginv(crossprod(X))
P  <- X %*% XtXg %*% t(X); Q <- diag(n) - P
A  <- solve(t(Z) %*% Q %*% Z)                       # (Z'QZ)^{-1}
A22 <- A[(p + 1):(p + r), (p + 1):(p + r), drop = FALSE]
Mp <- t(Zp) %*% Q %*% Zp
D  <- solve(Mp, t(Zp) %*% Q %*% Zr)
H  <- Zp %*% D - Zr
W  <- Z %*% A %*% t(Z); Wp <- Zp %*% solve(Mp, t(Zp))
cat("D =", round(D, 4), "  A22 =", round(A22, 4), "\n")

# estimable contrast tau_1 - tau_2 (C' = (0, 1, -1, 0))
cvec <- c(0, 1, -1, 0)
G <- t(cvec) %*% XtXg %*% t(X)                       # 1 x n

# prediction: treatment 3, z0 = (1.5, 0.5, 0)
x0 <- c(1, 0, 0, 1); z0 <- c(1.5, 0.5, 0)
a  <- drop(X %*% XtXg %*% x0)
u  <- drop(z0 - t(Z) %*% a); up <- u[1:p]; ur <- u[(p + 1):(p + r)]
v  <- drop(t(D) %*% up - ur)

## ---------------- Exact traces of the MSE matrices ----------------------
b0 <- matrix(c(1, -0.5), 1)                          # direction of B_r
base <- sum(diag(solve(Sig, t(b0) %*% solve(A22, b0))))  # tr(Delta) for c = 1
trS <- sum(diag(Sig))
rows <- NULL
for (tD in c(0, 0.25, 0.5, 1, 2, 4, 8)) {
  Br <- sqrt(tD / base) * b0
  Bp_F <- sum(diag(kronecker(Sig, A[1:p, 1:p])))
  Bp_R <- sum(diag(kronecker(Sig, solve(Mp)))) + sum((D %*% Br)^2)
  ct_F <- trS * drop(G %*% (diag(n) + W) %*% t(G))
  ct_R <- trS * drop(G %*% (diag(n) + Wp) %*% t(G)) + sum((G %*% H %*% Br)^2)
  pr_F <- trS * (sum(a^2) + drop(t(u) %*% A %*% u) + 1)
  pr_R <- trS * (sum(a^2) + drop(t(up) %*% solve(Mp, up)) + 1) +
          sum((t(v) %*% Br)^2)
  ES11 <- (Sig + t(Br) %*% solve(A22, Br) / (n - s - p))[1, 1]
  rows <- rbind(rows, c(trDelta = tD, Bp_F = Bp_F, Bp_R = Bp_R,
                        ct_F = ct_F, ct_R = ct_R, pred_F = pr_F,
                        pred_R = pr_R, ESig11 = ES11))
}
cat("\n=== Exact traces ===\n")
print(round(rows, 4))

## ---------------- Monte Carlo check at tr(Delta) = 1 ---------------------
set.seed(20261026)
tD <- 1; Br <- sqrt(tD / base) * b0
Tm <- matrix(c(0, 1, 2, 3, 0, -1, 0, 1), 4, 2)       # any treatment matrix
Bp <- matrix(c(1, 0.5, -0.5, 1), 2)
EY <- X %*% Tm + Zp %*% Bp + Zr %*% Br
Rchol <- chol(Sig)
LBf <- (A %*% t(Z) %*% Q)[1:p, , drop = FALSE]       # full B_p estimator
LBr <- solve(Mp, t(Zp) %*% Q)                        # reduced B_p estimator
Lcf <- G %*% (diag(n) - W %*% Q)                     # full contrast
Lcr <- G %*% (diag(n) - Wp %*% Q)                    # reduced contrast
Lpf <- a + drop(t(u) %*% A %*% t(Z) %*% Q)           # full prediction weights
Lpr <- a + drop(t(up) %*% solve(Mp, t(Zp) %*% Q))    # reduced prediction
mu0 <- drop(x0 %*% Tm + z0 %*% rbind(Bp, Br))
ct0 <- drop(t(cvec) %*% Tm)
nsim <- 20000
acc <- numeric(6)
for (b in 1:nsim) {
  Y  <- EY + matrix(rnorm(n * d), n, d) %*% Rchol
  y0 <- mu0 + drop(rnorm(d) %*% Rchol)               # new observation
  acc <- acc + c(sum((LBf %*% Y - Bp)^2), sum((LBr %*% Y - Bp)^2),
                 sum((Lcf %*% Y - ct0)^2), sum((Lcr %*% Y - ct0)^2),
                 sum((drop(Lpf %*% Y) - y0)^2), sum((drop(Lpr %*% Y) - y0)^2))
}
cat("\n=== Monte Carlo (tr(Delta) = 1) vs exact ===\n")
print(round(rbind(MonteCarlo = acc / nsim,
                  Exact = rows[rows[, "trDelta"] == 1, 2:7]), 4))
