source("01-prelim.R")

## ---- data.smoke ----
smoking <- read.table(
  "data/gum.txt",
  sep = "\t", header = TRUE
)
smoking.small <- smoking[,]

makeSmoke <- function(x) {
  rbind(
    data.frame(y = rep(1, x[2]), Study = x[1], Group = "Treated"),
    data.frame(y = rep(0, x[3] - x[2]), Study = x[1], Group = "Treated"),
    data.frame(y = rep(1, x[4]), Study = x[1], Group = "Control"),
    data.frame(y = rep(0, x[5] - x[4]), Study = x[1], Group = "Control")
  )
}
dat.smoke <- NULL
for (i in 1:nrow(smoking.small)) {
  dat.smoke <- rbind(dat.smoke, makeSmoke(smoking.small[i, ]))
}
dat.smoke$y <- as.factor(dat.smoke$y)
levels(dat.smoke$y) <- c("Remain", "Quit")
str(dat.smoke)

## ---- sumtab.smoke ----
tab <- with(smoking, {as.data.frame(rbind(
  Control = c(min(n0), round(mean(n0), 0), max(n0), round(mean(d0 / n0), 3),
              round(mean(d0 / n0) / (1 - mean(d0 / n0)), 3)),
  Treated = c(min(n1), round(mean(n1), 0), max(n1), round(mean(d1 / n1), 3),
              round(mean(d1 / n1) / (1 - mean(d1 / n1)), 3))
))})
colnames(tab) <- c("Min.", "Avg.", "Max.", "Prop. quit", "Odds quit")
knitr::kable(tab, format = "latex", align = "r")

## ---- mod.smoke ----
# This takes quite a long time, but saved object mod in data/smokingmod
# mod <- iprobit2(y = dat.smoke$y, dat.smoke$Group, dat.smoke$Study)
# mod1 <- iprobit(y = dat.smoke$y, dat.smoke$Group)
# mod2 <- iprobit2(y = dat.smoke$y, dat.smoke$Group, dat.smoke$Study,
#                  interactions = FALSE, maxit = 1000)

## ---- fit.smoke ----
load("data/smokingmod"); l3 <- logLik(mod)
load("data/smokingmod1"); l1 <- logLik(mod1)
load("data/smokingmod2"); l2 <- logLik(mod2)
calcOdds <- function(mod) {
  probs <- fitted(mod)$prob
  probs.lo <- fitted(mod, "lower")$prob
  probs.up <- fitted(mod, "upper")$prob
  latent <- mod$ystar
  studies <- as.character(unique(dat.smoke$Study))
  res <- matrix(NA, nrow = length(studies), ncol = 2 * 6 + 1)
  for (i in 1:length(studies)) {
    probs.treated <- probs[dat.smoke$Study == studies[i] &
                             dat.smoke$Group == "Treated"]
    probs.control <- probs[dat.smoke$Study == studies[i] &
                             dat.smoke$Group == "Control"]
    probs.treated.lo <- probs.lo[dat.smoke$Study == studies[i] &
                                   dat.smoke$Group == "Treated"]
    probs.control.lo <- probs.lo[dat.smoke$Study == studies[i] &
                                   dat.smoke$Group == "Control"]
    probs.treated.up <- probs.up[dat.smoke$Study == studies[i] &
                                   dat.smoke$Group == "Treated"]
    probs.control.up <- probs.up[dat.smoke$Study == studies[i] &
                                   dat.smoke$Group == "Control"]
    latent.treated <- latent[dat.smoke$Study == studies[i] &
                               dat.smoke$Group == "Treated"]
    latent.control <- latent[dat.smoke$Study == studies[i] &
                               dat.smoke$Group == "Control"]
    current.prob <- c(mean(probs.treated), mean(probs.control))
    current.prob.lo <- c(mean(probs.treated.lo), mean(probs.control.lo))
    current.prob.up <- c(mean(probs.treated.up), mean(probs.control.up))
    current.odds <- current.prob / (1 - current.prob)
    current.odds.lo <- current.prob.lo / (1 - current.prob.lo)
    current.odds.up <- current.prob.up / (1 - current.prob.up)
    current.lodds <- log(current.odds)
    current.latent <- c(mean(latent.treated), mean(latent.control))
    res[i, ] <- c(
      current.prob, current.odds, current.odds[1] / current.odds[2],
      current.odds.lo, current.odds.up, current.lodds, current.latent
    )
  }
  res <- data.frame(studies, res)
  colnames(res) <- c(
    "Study", "Prob.T", "Prob.C", "Odds.T", "Odds.C", "Odds.ratio",
    "Odds.T.lo", "Odds.C.lo", "Odds.T.up", "Odds.C.up", "Log-Odds.T",
    "Log-Odds.C", "Latent.T", "Latent.C"
  )
  res

  list(odds = res[, c(1, 4:10)], avg.odds = mean(res[, "Odds.ratio"]))
}
tmp <- calcOdds(mod)  # f(x) + f(j) + f(x,j)
BrierFun <- function(x) mean((fitted(x)$prob - x$y) ^ 2 )
Brier.score <- c(
  BrierFun(mod1), BrierFun(mod2), BrierFun(mod)
)

## ---- mod.compare.smoke ----
tab.compare <- as.data.frame(cbind(
  c("$f_1$", "$f_1 + f_2$",
    "$f_1 + f_2 + f_{12}$"),
  decPlac(c(l1, l2, l3), 2),
  decPlac(Brier.score, 4),
  c(1, 2, 2)
))
colnames(tab.compare) <- c("Model", "Lower bound", "Brier score",
                           "No. of RKHS\\newline param.")
knitr::kable(tab.compare, format = "latex", align = c("l", "r", "r", "R{2.2cm}"),
             row.names = TRUE, escape = FALSE)

## ---- plot.smoke ----
set.seed(456)
odds.dat.orig <- tmp$odds[sample(1:nrow(tmp$odds), 10), ]
odds.dat <- odds.dat.orig[, -(5:8)]
plot.dat <- reshape2::melt(odds.dat, id.vars = c("Study", "Odds.ratio"),
                           measure.vars = c("Odds.C", "Odds.T"))
plot.dat <- cbind(plot.dat,
                  lo = c(odds.dat.orig[, 6], odds.dat.orig[, 5]),
                  up = c(odds.dat.orig[, 8], odds.dat.orig[, 7]))
lab.msg <- paste("avg. odds\nratio = ", decPlac(tmp$avg.odds, 3))
ggplot(plot.dat, aes(x = variable, y = value, col = Study, group = Study)) +
  geom_point() +
  geom_line() +
  # geom_errorbar(aes(x = variable, ymin = lo, ymax = up, width = 0.05)) +
  directlabels::geom_dl(aes(label = decPlac(Odds.ratio, 3)),
                        method = list("last.polygons")) + #,
                                      #directlabels::dl.trans(x = x + 0.2))) +
  directlabels::geom_dl(aes(label = Study),
                        method = list("first.bumpup",
                                      directlabels::dl.trans(x = x - 0.2))) +
  annotate("text", x = 2.4, y = 0.13, label = lab.msg, size = 4.3,
           col = "grey30", lineheight = 0.9) +
  scale_x_discrete(labels = c("Control", "Nicotine gum treatment")) +
  labs(x = NULL, y = "Model predicted odds") +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 12))

## ---- plot.smoke.all ----
odds.dat.orig <- tmp$odds
odds.dat <- odds.dat.orig[, -(5:8)]
plot.dat <- reshape2::melt(odds.dat, id.vars = c("Study", "Odds.ratio"),
                           measure.vars = c("Odds.C", "Odds.T"))
plot.dat <- cbind(plot.dat,
                  lo = c(odds.dat.orig[, 6], odds.dat.orig[, 5]),
                  up = c(odds.dat.orig[, 8], odds.dat.orig[, 7]))
lab.msg <- paste("avg. odds\nratio = ", decPlac(tmp$avg.odds, 3))
ggplot(plot.dat, aes(x = variable, y = value, col = Study, group = Study)) +
  geom_point() +
  geom_line() +
  directlabels::geom_dl(aes(label = decPlac(Odds.ratio, 3)),
                        method = list("last.polygons")) +
  directlabels::geom_dl(aes(label = Study),
                        method = list("first.bumpup",
                                      directlabels::dl.trans(x = x - 0.2))) +
  annotate("text", x = 2.2, y = 0.03, label = lab.msg, size = 3.9,
           col = "grey30", lineheight = 0.9) +
  scale_x_discrete(labels = c("Control", "Nicotine gum treatment")) +
  labs(x = NULL, y = "Model predicted odds") +
  theme_bw() +
  theme(legend.position = "none")

## ---- tab.smoke.all ----
odds.dat.orig <- tmp$odds
comb.odds <- function(x, k = 2) {
  xx <- c(x[1], decPlac(x[-1], k))
  odds.T <- paste0(xx[2], " (", xx[5], ", ", xx[7], ")")
  odds.C <- paste0(xx[3], " (", xx[6], ", ", xx[8], ")")
  odds.R.loup <- decPlac(sort(c(max(x[5], x[7]) / min(x[6], x[8]),
                                min(x[5], x[7]) / max(x[6], x[8]))), k)
  odds.R <- paste0(xx[4], " (", odds.R.loup[1], ", ", odds.R.loup[2], ")")
  c(odds.C, odds.T, odds.R)
}
tmp1 <- split(odds.dat.orig, seq(nrow(odds.dat.orig)))
tmp1 <- lapply(tmp1, comb.odds)
comb.odds.df <- data.frame(cbind(
  as.character(odds.dat.orig$Study),
  matrix(unlist(tmp1), ncol = 3, byrow = TRUE)
))
colnames(comb.odds.df) <- c("Study", "Control", "Treated", "Odds ratio")
levels(comb.odds.df[, 4])[27] <- "3.04 (1.83, 11.6)"  # manual fix signif.
knitr::kable(comb.odds.df, align = c("l", "r", "r", "r"))
