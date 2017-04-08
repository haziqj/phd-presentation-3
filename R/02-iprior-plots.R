## ---- prelim ----
library(iprior)
library(ggplot2)
library(gganimate)
library(reshape2)
library(directlabels)
gg_colour_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## ---- points ----
set.seed(123)
N <- 150
f <- function(x, truth = FALSE) {
  35 * dnorm(x, mean = 1, sd = 0.8) +
    65 * dnorm(x, mean = 4, sd = 1.5) +
    (x > 4.5) * (exp((1.25 * (x - 4.5))) - 1) +
    3 * dnorm(x, mean = 2.5, sd = 0.3)
}
x <- c(seq(0.2, 1.9, length = N * 5 / 8), seq(3.7, 4.6, length = N * 3 / 8))
x <- sample(x, size = N)
x <- x + rnorm(N, sd = 0.65)  # adding random fluctuation to the x
x <- sort(x)
y.err <- rt(N, df = 1)
y <- f(x) + sign(y.err) * pmin(abs(y.err), rnorm(N, mean = 4.1))  # adding random terms to the y

# True values
x.true <- seq(-2.1, 7, length = 1000)
y.true <- f(x.true, TRUE)

# Data for plot
dat <- data.frame(x, y)
dat.truth <- data.frame(x.true, y.true)

p1 <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y)) +
  scale_x_continuous(
    limits = c(min(x.true), max(x.true)),
    breaks = NULL, name = expression(italic(x))
  ) +
  scale_y_continuous(
    limits = c(min(y) - 5, max(y) + 5),
    breaks = NULL, name = expression(italic(y))
  ) +
  theme_bw()

## ---- plot.function1 ----
fnH4 <- function(x, y = NULL, l = 1) {
x <- scale(x, scale = FALSE)
if (is.vector(x))
  x <- matrix(x, ncol = 1)
n <- nrow(x)
A <- matrix(0, n, n)
index.mat <- upper.tri(A)
index <- which(index.mat, arr.ind = TRUE)
xcrossprod <- tcrossprod(x)
if (is.null(y)) {
  tmp1 <- diag(xcrossprod)[index[, 1]]
  tmp2 <- diag(xcrossprod)[index[, 2]]
  tmp3 <- xcrossprod[index]
  A[index.mat] <- tmp1 + tmp2 - 2 * tmp3
  A <- A + t(A)
  tmp <- exp(-A / (2 * l ^ 2))
} else {
  if (is.vector(y))
    y <- matrix(y, ncol = 1)
  else y <- as.matrix(y)
  y <- sweep(y, 2, attr(x, "scaled:center"), "-")
  m <- nrow(y)
  B <- matrix(0, m, n)
  indexy <- expand.grid(1:m, 1:n)
  ynorm <- apply(y, 1, function(z) sum(z ^ 2))
  xycrossprod <- tcrossprod(y, x)
  tmp1 <- ynorm[indexy[, 1]]
  tmp2 <- diag(xcrossprod)[indexy[, 2]]
  tmp3 <- as.numeric(xycrossprod)
  B[, ] <- tmp1 + tmp2 - 2 * tmp3
  tmp <- exp(-B / (2 * l ^ 2))
}
tmp
}

dev.SEkern <- function(theta, y = y) {
  alpha <- mean(y)
  lambda <- theta[1]
  psi <- theta[2]
  n <- length(y)
  H <- fnH4(x, l = theta[3])
  H2 <- H %*% H
  Vy <- psi * lambda ^ 2 * H2 + diag(1 / psi, n)
  tmp <- -(n / 2) * log(2 * pi) - (1 / 2) * determinant(Vy)$mod -
    (1 / 2) * (y - alpha) %*% solve(Vy, y - alpha)
  as.numeric(-2 * tmp)
}

plot1 <- function(kernel, no.of.draws = 100) {
  # Fit an I-prior model -------------------------------------------------------
  if (kernel == "SE") {
    mod.se <- optim(c(1, 1, 1), dev.SEkern, method = "L-BFGS", y = y,
                    lower = c(-Inf, 1e-9, 1e-9))
    n <- length(y)
    H <- fnH4(x, l = mod.se$par[3])
    H2 <- H %*% H
    alpha <- mean(y)
    lambda <- mod.se$par[1]
    psi <- mod.se$par[2]
    Vy <- psi * lambda ^ 2 * H2 + diag(1 / psi, n)
    w.hat <- psi * lambda * H %*% solve(Vy, y - alpha)
    VarY.inv <- solve(Vy)
    H.star <- fnH4(x = x, y = x.true, l = mod.se$par[3])
    y.fitted <- as.numeric(mean(y) + lambda * H.star %*% w.hat)
    y.fitted2 <- as.numeric(mean(y) + lambda * H %*% w.hat)
  } else {
    if (kernel == "Canonical") {
      mod <- iprior(kernL(y, x, model = list(kernel = "Canonical")))
      H.star <- fnH2(x = x, y = x.true)
    }
    if (kernel == "FBM") {
      mod <- fbmOptim(kernL(y, x, model = list(kernel = "FBM")))
      H.star <- fnH3(x = x, y = x.true, gamma = mod$ipriorKernel$model$Hurst[1])
    }
    # Estimated values  ----------------------------------------------------------
    y.fitted <- predict(mod, list(matrix(x.true, ncol = 1)))
    y.fitted2 <- fitted(mod)
    lambda <- mod$lambda
    psi <- mod$psi
    w.hat <- mod$w.hat
    H <- mod$ipriorKernel$Hl[[1]]
    H2 <- H %*% H
    Vy <- vary(mod)
    VarY.inv <- solve(Vy)
  }

  # Prepare random draws from prior and posterior ------------------------------
  draw.pri <- draw.pos <- matrix(NA, ncol = no.of.draws, nrow = nrow(H.star))
  L <- chol(VarY.inv)
  for (i in 1:no.of.draws) {
    w.draw <- w.hat + crossprod(L, rnorm(length(y)))
    draw.pos[, i] <- mean(y) + lambda * as.numeric(H.star %*% w.draw)
    draw.pri[, i] <- mean(y) +
      lambda * as.numeric(H.star %*% rnorm(length(y), sd = sqrt(psi)))
  }
  dat.f <- rbind(data.frame(x = x.true, y = y.fitted, type = "Posterior"),
                 data.frame(x = x.true, y = mean(y), type = "Prior"))
  melted.pos <- melt(data.frame(f = draw.pos, x = x.true), id.vars = "x")
  melted.pri <- melt(data.frame(f = draw.pri, x = x.true), id.vars = "x")
  melted <- rbind(cbind(melted.pri, type = "Prior"),
                  cbind(melted.pos, type = "Posterior"))

  # Posterior predictive covariance matrix -------------------------------------
  varystar <- psi * (lambda ^ 2) * H.star %*% t(H.star) +
    diag(1 / psi, nrow(H.star))
  covystary <- psi * (lambda ^ 2) * H.star %*% H
  VarY.stary <- varystar - covystary %*% VarY.inv %*% t(covystary)
  dat.fit <- data.frame(x.true, y.fitted, sdev = sqrt(diag(VarY.stary)),
                        type = "95% credible interval")

  # Prepare random draws for posterior predictive checks -----------------------
  VarY.hat <- Vy - (psi ^ 2) * (lambda ^ 4) * H2 %*% VarY.inv %*% H2
  ppc <- matrix(NA, ncol = no.of.draws, nrow = nrow(H))
  L <- chol(VarY.hat)
  for (i in 1:no.of.draws) {
    ppc[, i] <- y.fitted2 + crossprod(L, rnorm(n))
    # ppc[, i] <- y.fitted2 + rnorm(n, sd = sqrt(diag(VarY.hat)))
  }
  melted.ppc <- melt(data.frame(x = x, ppc = ppc), id.vars = "x")
  melted.ppc <- cbind(melted.ppc, type = "Posterior predictive check")

  # Random draws from prior and posterior function -----------------------------
  p2.tmp <- ggplot() +
    geom_point(data = dat, aes(x = x, y = y), col = "grey55", alpha = 0.5) +
    scale_x_continuous(
      limits = c(min(x.true), max(x.true)),
      breaks = NULL, name = expression(italic(x))
    ) +
    scale_y_continuous(
      limits = c(min(y) - 5, max(y) + 5),
      breaks = NULL, name = expression(italic(y))
    ) +
    theme_bw()
  p2 <- p2.tmp +
    geom_line(data = melted, aes(x = x, y = value, group = variable),
              col = "steelblue3", size = 0.19, alpha = 0.5) +
    facet_grid(type ~ .) +
    geom_line(data = dat.f, aes(x = x, y = y), size = 1, linetype = 2, col = "grey10")
  p2.prior <- p2.tmp +
    geom_line(data = subset(melted, type == "Prior"),
              aes(x = x, y = value, group = variable),
              col = "steelblue3", size = 0.19, alpha = 0.5) +
    facet_grid(type ~ .)
  p2.prior.line <- p2.prior +
    geom_line(data = subset(dat.f, type == "Prior"), aes(x = x, y = y),
              size = 1, linetype = 2, col = "grey10")
  p2.posterior <- p2.tmp +
    geom_line(data = subset(melted, type == "Posterior"),
              aes(x = x, y = value, group = variable),
              col = "steelblue3", size = 0.19, alpha = 0.5) +
    facet_grid(type ~ .)
  p2.posterior.line <- p2.posterior +
    geom_line(data = subset(dat.f, type == "Posterior"), aes(x = x, y = y),
              size = 1, linetype = 2, col = "grey10")

  # Confidence band for predicted values  --------------------------------------
  p3 <- p1 +
    geom_line(data = dat.fit, aes(x = x.true, y = y.fitted), col = "grey50",
              size = 0.9, linetype = 2) +
    geom_ribbon(data = dat.fit, fill = "grey70", alpha = 0.5,
                aes(x = x.true, ymin = y.fitted - 1.96 * sdev,
                    ymax = y.fitted + 1.96 * sdev)) +
    facet_grid(type ~ .)

  p4 <- p2 +
    geom_line(data = dat.truth, aes(x = x.true, y = y.true, col = "Fitted"),
              size = 1, alpha = 0.75) + theme(legend.position = "none")

  # Posterior predictive checks ------------------------------------------------
  p5 <- ggplot() +
  scale_x_continuous(breaks = NULL, name = expression(italic(y))) +
  scale_y_continuous(breaks = NULL) +
  geom_line(data = melted.ppc,
            aes(x = value, group = variable, col = "yrep", size = "yrep"),
            stat = "density", alpha = 0.5) +
  geom_line(data = dat, aes(x = y, col = "y", size = "y"), stat = "density") +
    theme(legend.position = "bottom") +
  scale_colour_manual(
    name = NULL, labels = c("Observed", "Replications"),
    values = c("grey10", "steelblue3")
  ) +
  scale_size_manual(
    name = NULL, labels = c("Observed", "Replications"),
    values = c(1.1, 0.19)
  ) +
  facet_grid(type ~ .) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.5))

  list(p2 = p2, p2.prior = p2.prior, p2.posterior = p2.posterior,
       p2.prior.line = p2.prior.line,
       p2.posterior.line = p2.posterior.line, p3 = p3, p4 = p4, p5 = p5)
}

## ---- canonical.kernel ----
plot.can <- plot1("Canonical")

## ---- fbm.kernel ----
plot.fbm <- plot1("FBM")

## ---- se.kernel.mle ----
plot.se <- plot1("SE")

## ---- variational.comparison ----
# The pdf to approximate
rx <- function(x) {
  exp(-(x ^ 2) / 2) / (1 + exp(-(20 * x + 4)))
}
const <- 1 / integrate(rx, -Inf, Inf)$value
px <- function(x) {
  rx(x) * const
}
dev <- function(x) -2 * log(px(x))
EX <- integrate(function(x) x * px(x), -Inf, Inf)$value
EX2 <- integrate(function(x) (x ^ 2) * px(x), -Inf, Inf)$value
VarX <- EX2 - EX ^ 2
x.hat <- optim(0, dev, method = "BFGS")$par

# Laplace approximation
tmp <- optim(0, function(x) -log(rx(x)), method = "BFGS", hessian = TRUE)
mode.p <- tmp$par
var.p <- 1 / tmp$hessian
lap.approx <- function(x) dnorm(x, mean = mode.p, sd = sqrt(var.p))
lap.approx.dev <- function(x) -2 * log(lap.approx(x))
x.lap <- optim(0, lap.approx.dev, method = "BFGS")$par

# Variational approximation
var.approx <- function(x) dnorm(x, mean = EX, sd = sqrt(VarX))
var.approx.dev <- function(x) -2 * log(var.approx(x))
x.var <- optim(0, var.approx.dev, method = "BFGS")$par

# Density plots
p6 <- function(mode = FALSE, laplace = FALSE, mean = FALSE, variational = FALSE) {
  x <- seq(-2, 3, length = 1000)
  den.df.full <- data.frame(x = x, Truth = px(x), Laplace = lap.approx(x),
                            Variational = var.approx(x))
  dl.list <- list("top.bumptwice",
                  dl.move("Truth", 0.3, 0.69),
                  dl.move("Laplace", -0.54, 0.55),
                  dl.move("Variational", 1.2, 0.60))
  thecol <- gg_colour_hue(3)
  if (isTRUE(laplace)) {
    if (isTRUE(variational)) {
      den.df <- den.df.full
    } else {
      den.df <- den.df.full[, -4]
      thecol <- thecol[1:2]
    }
  } else {
    den.df <- den.df.full[, -(3:4)]
    thecol <- thecol[1]
  }
  den.df <- melt(den.df, id.vars = "x")

  p <- ggplot(data = den.df, aes(x = x, y = value, group = variable)) +
    geom_line(aes(col = variable)) +
    scale_colour_manual(values = thecol) +
    scale_fill_manual(values = thecol) +
    scale_x_continuous(
      breaks = NULL, name = expression(italic(z))
    ) +
    scale_y_continuous(
      limits = c(-0.03, px(x.hat) + 0.05),
      breaks = NULL, name = "Density"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  if (isTRUE(mode)) {
    p <- p +
      geom_vline(xintercept = x.hat, linetype = 2, col = thecol[1]) +
      annotate("text", label = "mode", col = thecol[1], x = x.hat - 0.2, y = -0.03)
  }
  if (isTRUE(mean)) {
    p <- p +
      geom_vline(xintercept = EX, linetype = 2, col = thecol[1]) +
      annotate("text", label = "mean", col = thecol[1], x = EX + 0.2, y = -0.03)
  }

  p <- p +
    geom_area(aes(fill = variable), alpha = 0.3, position = "identity") +
    geom_dl(aes(label = variable, col = variable), method = dl.list)

  p
}

# Deviance plots
p7 <- function() {
  x <- seq(-0.5, 1, length = 1000)
  dev.df.full <- data.frame(x = x, Truth = dev(x), Laplace = lap.approx.dev(x),
                            Variational = var.approx.dev(x))
  dev.df <- melt(dev.df.full, id.vars = "x")
  thecol <- gg_colour_hue(3)

  p <- ggplot(data = dev.df, aes(x = x, y = value, group = variable)) +
    geom_line(aes(col = variable), size = 0.9) +
    geom_segment(x = x.hat, y = dev(x.hat), xend = x.hat, yend = -2,
                 col = thecol[1], linetype = 2) +
    geom_segment(x = x.lap + 0.01, y = lap.approx.dev(x.lap + 0.01),
                 xend = x.lap + 0.01, yend = -2,
                 col = thecol[2], linetype = 2) +
    geom_segment(x = x.var, y = var.approx.dev(x.var), xend = x.var, yend = -2,
                 col = thecol[3], linetype = 2) +
    scale_colour_manual(values = thecol) +
    scale_x_continuous(
      limits = c(-0.5, 1.18),
      breaks = NULL, name = expression(italic(z))
    ) +
    scale_y_continuous(
      limits = c(-1, 7),
      breaks = NULL, name = expression(paste("Deviance (", -2, " x Log-density)"))
    ) +
    geom_dl(aes(label = variable, col = variable), method = "last.bumpup") +
    theme_bw() +
    theme(legend.position = "none")

  p
}


## ---- variational.example ----
# Simulate data from N(10, 4)
set.seed(123)
n <- 30
mu <- 0
psi <- 1
ydat <- rnorm(n, mean = 0, sd = sqrt(1 / psi))
s <- var(ydat) * (n - 1) / n

# pdf of normal gamma
dnormgamma <- function(x, mu = 0, lambda = 1, alpha = 1, beta = 1) {
  const <- alpha * log(beta) + 0.5 * log(lambda) - lgamma(alpha) - 0.5 * log(2 * pi)
  res <- const + (alpha - 0.5) * log(x[2]) - beta * x[2] -
    lambda * x[2] * (x[1] - mu) ^ 2 / 2
  exp(res)
}

# pdf of variational approximation
dvarapprox <- function(x, mean.mu, var.mu, shape.psi, rate.psi) {
  dnorm(x[1], mean = mean.mu, sd = sqrt(var.mu)) *
    dgamma(x[2], shape = shape.psi, rate = rate.psi)
}

# Exact parameters of the posterior
lambda0 <- 0.01
mu0 <- 0
alpha0 <- 0.01
beta0 <- 0.01
mu.post <- (lambda0 * mu0 + n * mean(ydat)) / (lambda0 + n)
lambda.post <- lambda0 + n
alpha.post <- alpha0 + n / 2
beta.post <- beta0 + 0.5 * (n * s + lambda0 * n * (mean(ydat) - mu0) ^ 2 / (lambda0 + n))
mu.post  # estimate of mean
alpha.post / beta.post  # estimate of precision

# Data for contour plot
plotxlim <- c(-0.4, 0.5)
plotylim <- c(0.5, 1.8)
x <- seq(plotxlim[1], plotxlim[2], length = 50)  # mean
y <- seq(plotylim[1], plotylim[2], length = 50)  # variance
xy <- expand.grid(x, y); names(xy) <- c("x", "y")
z <- apply(xy, 1, dnormgamma, mu = mu.post, lambda = lambda.post,
           alpha = alpha.post, beta = beta.post)
contour.df <- data.frame(xy, z)
p1 <- ggplot(contour.df, aes(x = x, y = y, z = z)) + #geom_contour(size = 1) +
  xlim(plotxlim[1], plotxlim[2]) + coord_cartesian(ylim = plotylim) +
  labs(x = expression(Mean~(mu)), y = expression(Precision~(psi))) +
  geom_raster(aes(fill = density)) +
  theme_bw()
p1


# Function for variational inference contour data
varContour <- function(niter = 10) {
  # Set up result list
  res <- NULL

  # Initialise
  q.psi.c <- 60; q.psi.d <- 40
  q.mu.a <- 0.4; q.mu.b <- 400
  q.mu.var <- q.psi.d / (q.psi.c * q.mu.b)
  res[1] <- apply(xy, 1, dvarapprox, mean.mu = q.mu.a,
                  var.mu = q.mu.var, shape.psi = q.psi.c, rate.psi = q.psi.d)
}


# Variational approximation


# Update mu
q.mu.a <- (lambda0 * mu0 + n * mean(ydat)) / (lambda0 + n)
q.mu.var <- q.psi.d / (q.psi.c *  (lambda0 + n))
zq10 <- apply(xy, 1, dvarapprox, mean.mu = q.mu.a,
              var.mu = q.mu.var, shape.psi = q.psi.c,
              rate.psi = q.psi.d)
contour.df <- data.frame(xy, z = zq10)

# Update psi
q.psi.c <- n / 2
q.psi.d <-  0.5 * (sum(ydat ^ 2) - 2 * mean(ydat) * q.mu.a + q.mu.var + q.mu.a ^ 2)
zq11 <- apply(xy, 1, dvarapprox, mean.mu = q.mu.a,
              var.mu = q.mu.var, shape.psi = q.psi.c,
              rate.psi = q.psi.d)
contour.df <- data.frame(xy, z = zq11)

# Update mu
q.mu.a <- (lambda0 * mu0 + n * mean(ydat)) / (lambda0 + n)
q.mu.var <- q.psi.d / (q.psi.c *  (lambda0 + n))
zq21 <- apply(xy, 1, dvarapprox, mean.mu = q.mu.a,
              var.mu = q.mu.var, shape.psi = q.psi.c,
              rate.psi = q.psi.d)
contour.df <- data.frame(xy, z = zq21)

# Update psi
q.psi.c <- n / 2
q.psi.d <-  0.5 * (sum(ydat ^ 2) - 2 * mean(ydat) * q.mu.a + q.mu.var + q.mu.a ^ 2)
zq22 <- apply(xy, 1, dvarapprox, mean.mu = q.mu.a,
              var.mu = q.mu.var, shape.psi = q.psi.c,
              rate.psi = q.psi.d)
contour.df <- data.frame(xy, z = zq22)

# Update mu
q.mu.a <- (lambda0 * mu0 + n * mean(ydat)) / (lambda0 + n)
q.mu.var <- q.psi.d / (q.psi.c *  (lambda0 + n))
zq32 <- apply(xy, 1, dvarapprox, mean.mu = q.mu.a,
              var.mu = q.mu.var, shape.psi = q.psi.c,
              rate.psi = q.psi.d)
contour.df <- data.frame(xy, z = zq32)

# Update psi
q.psi.c <- n / 2
q.psi.d <-  0.5 * (sum(ydat ^ 2) - 2 * mean(ydat) * q.mu.a + q.mu.var + q.mu.a ^ 2)
zq33 <- apply(xy, 1, dvarapprox, mean.mu = q.mu.a,
              var.mu = q.mu.var, shape.psi = q.psi.c,
              rate.psi = q.psi.d)
contour.df <- data.frame(xy, z = zq33)

# Contour plot
p1 + geom_contour(data = contour.df, aes(x = x, y = y, z = z), col = "red")

test.df <- rbind(
  data.frame(xy, z = zq00, iteration = "Iteration 0"),
  data.frame(xy, z = zq10, iteration = "Iteration 1 (updating mean)"),
  data.frame(xy, z = zq11, iteration = "Iteration 2 (updating precision)")
)

p <- ggplot(test.df, aes(x = x, y = y, z = z, group = iteration, frame = iteration)) + geom_contour() +
  xlim(plotxlim[1], plotxlim[2]) + coord_cartesian(ylim = plotylim) + theme_bw()
gganimate_save(gganimate(p), "var.gif")


## ---- save.plots.for.presentation ----
ggsave("figure/points.pdf", p1, width = 6.5, height = 6.5 / 2.25)
# ggsave("figure/can-prior.pdf", plot.can$p2.prior.line,
#        width = 6.5, height = 6.5 / 1.5)
# ggsave("figure/can-posterior.pdf", plot.can$p2.posterior.line,
#        width = 6.5, height = 6.5 / 1.5)
ggsave("figure/fbm-prior.pdf", plot.fbm$p2.prior.line,
       width = 6.5, height = 6.5 / 1.5)
ggsave("figure/fbm-posterior.pdf", plot.fbm$p2.posterior.line,
       width = 6.5, height = 6.5 / 1.5)
ggsave("figure/fbm-posterior-truth.pdf", {
  plot.fbm$p2.posterior +
    geom_line(data = dat.truth, aes(x = x.true, y = y.true), size = 1,
              alpha = 0.75, col = "red3") +
    annotate("text", label = "Truth", col = "red3", x = max(x.true),
             y = max(y.true) + 1)
  }, width = 6.5, height = 6.5 / 1.5)
ggsave("figure/credible-interval.pdf", plot.fbm$p3, width = 6.5, height = 6.5 / 1.5)
ggsave("figure/ppc.pdf", plot.fbm$p5, width = 6.5, height = 6.5 / 1.5)

ggsave("figure/compare1.pdf", p6(), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare2.pdf", p6(T), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare3.pdf", p6(T, T), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare4.pdf", p6(T, T, T), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare5.pdf", p6(T, T, T, T), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare6.pdf", p7(), width = 6.5, height = 6.5 / 1.5)









