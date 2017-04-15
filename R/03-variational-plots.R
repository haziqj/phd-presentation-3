source("01-prelim.R")

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

## ---- save.plots.for.presentation ----
ggsave("figure/compare1.pdf", p6(), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare2.pdf", p6(T), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare3.pdf", p6(T, T), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare4.pdf", p6(T, T, T), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare5.pdf", p6(T, T, T, T), width = 6.5, height = 6.5 / 1.5)
ggsave("figure/compare6.pdf", p7(), width = 6.5, height = 6.5 / 1.5)
