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
colnames(tab) <- c("Min", "Average", "Max", "Prop. quit", "Odds quit")
knitr::kable(tab, format = "latex")

## ---- mod.smoke ----
# This takes quite a long time, but saved object mod in data/smokingmod
mod <- iprobit2(y = dat.smoke$y, dat.smoke$Group, dat.smoke$Study)
mod1 <- iprobit(y = dat.smoke$y, dat.smoke$Group)
mod2 <- iprobit2(y = dat.smoke$y, dat.smoke$Group, dat.smoke$Study,
                 interactions = FALSE, maxit = 1000)

## ---- fit.smoke ----
load("data/smokingmod")
load("data/smokingmod1")
load("data/smokingmod2")
calcOdds <- function(mod) {
  probs <- unique(fitted(mod)$prob)
  # probs.normalised <- {
  #   tmp <- NULL
  #   ind <- seq(1, length(probs), by = 2)
  #   for (i in ind) {
  #     normalis <- probs[i] + probs[i + 1]
  #     tmp <- c(tmp, probs[i] / normalis, probs[i + 1] / normalis)
  #   }
  #   tmp
  # }
  # probs <- probs.normalised
  latent <- unique(mod$ystar)
  probs.dat <- data.frame(
    Study = levels(dat.smoke$Study),
    Treated = probs[seq(1, length(probs), by = 4)],
    Control = probs[seq(3, length(probs), by = 4)]
  )
  latent.dat <- data.frame(
    Study = levels(dat.smoke$Study),
    Treated = latent[seq(1, length(latent), by = 4)],
    Control = latent[seq(3, length(latent), by = 4)]
  )
  odds.dat <- data.frame(
    Study = levels(dat.smoke$Study),
    Treated = probs.dat$Treated / (1 - probs.dat$Treated),
    Control =  probs.dat$Control / (1 - probs.dat$Control)
  )
  log.odds.dat <- odds.dat <- data.frame(
    Study = levels(dat.smoke$Study),
    Treated = log(odds.dat$Treated),
    Control =  log(odds.dat$Control)
  )
  avg.odds <- apply(odds.dat[, -1], 2, mean)
  avg.odds["Treated"] / avg.odds["Control"]
  odds.dat <- cbind(odds.dat,
                    "Odds.ratio" = round(odds.dat$Treated / odds.dat$Control, 3))

  list(odds = odds.dat, avg.odds = avg.odds)
}
tmp <- calcOdds(mod)  # f(x) + f(j) + f(x,j)
odds.dat <- tmp$odds
tmp1 <- calcOdds(mod1)$avg.odds  # f(x)
tmp2 <- calcOdds(mod2)$avg.odds  # f(x) + f(j)

## ---- mod.compare.smoke ----
tab.compare <- as.data.frame(cbind(
  c("f(x)", "f(x) + f(j)", "f(x) + f(j) + f(x,j)"),
  round(c(logLik(mod), logLik(mod), logLik(mod)), 2),
  round(tmp$avg.odds[1] / tmp$avg.odds[2], 3)
))
colnames(tab.compare) <- c("Model", "Lower bound", "Avg. odds ratio")
tab.compare

## ---- plot.smoke ----
set.seed(456)
odds.dat <- odds.dat[sample(1:nrow(odds.dat), 10), ]
plot.dat <- reshape2::melt(odds.dat, id.vars = c("Study", "Odds.ratio"),
                           measure.vars = c("Control", "Treated"))
ggplot(plot.dat, aes(x = variable, y = value, col = Study, group = Study)) +
  geom_point() +
  geom_line() +
  directlabels::geom_dl(aes(label = Odds.ratio),
                        method = list("last.polygons")) + #,
                                      #directlabels::dl.trans(x = x + 0.2))) +
  directlabels::geom_dl(aes(label = Study),
                        method = list("first.bumpup",
                                      directlabels::dl.trans(x = x - 0.2))) +
  scale_x_discrete(labels = c("Control", "Nicotine gum treatment")) +
  labs(x = NULL, y = "Model predicted odds") +
  theme_bw() +
  theme(legend.position = "none")




