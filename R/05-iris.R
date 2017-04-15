source("01-prelim.R")

## ---- data.iris ----
data(iris)
y <- ifelse(iris$Species == "setosa", 1, 0)
setosa <- as.factor(y)
levels(setosa) <- c("Others", "Setosa")
y <- setosa
X <- iris[, -(3:5)]
n <- length(y)
H <- fnH2(X)  # canonical kernel
ggplot(data = cbind(X, Class = setosa),
       aes(x = Sepal.Length, y = Sepal.Width, col = Class)) +
  geom_point(size = 3) +
  theme_bw()

## ---- mod.iris ----
set.seed(123)
system.time(
  (mod <- iprobit(y, X, maxit = 10000))
)

## ---- fake.mod.iris ----
system.time(
  (mod <- iprobit(y, X))
)
cat(
  "\n|==================                              |  61%\n",
  "Converged after 6141 iterations.\n",
  "Training error rate: 0 %\n",
  "   user  system elapsed\n",
  " 67.857   6.396  74.277"
)

## ---- save.plot.iris ----
p1 <- iplot_decbound(mod)
p1 <- p1 +
  scale_x_continuous(
    breaks = NULL, name = NULL
  ) +
  scale_y_continuous(
    breaks = NULL, name = NULL
  ) + theme_classic() + theme(legend.position = "none")
p2 <- iplot_prob(mod, 2)
p2 <- p2 +
  scale_x_continuous(
    breaks = NULL, name = NULL
  ) +
  scale_y_continuous(
    breaks = NULL, name = NULL
  ) + theme_classic() + theme(legend.position = "none")
ggsave("../figure/iprior-probit1.pdf", p1, width = 4, height = 3)
ggsave("../figure/iprior-probit2.pdf", p2, width = 3, height = 2.25)
