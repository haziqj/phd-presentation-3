source("01-prelim.R")

data(iris)
y <- ifelse(iris$Species == "setosa", 1, 0)
X <- iris[, -5]
n <- length(y)
H <- fnH2(X)

mod <- "
data {
int<lower=0> n;
int<lower=0,upper=1> y[n];
matrix[n, n] H; // canonical RKHS
}
parameters {
real alpha;
real<lower=0> lambda;
vector[n] w;
}
transformed parameters {
vector<lower=0,upper=1>[n] pi;
pi = Phi_approx(alpha + lambda * H * w);
}
model {
target += bernoulli_lpmf(y | pi);
target += normal_lpdf(w | 0, 1);
target += cauchy_lpdf(lambda | 0, 2.5);
target += normal_lpdf(alpha | 0, 2.5);
}
"
m <- stan_model(model_code = mod)
m@model_name <- "iprior.probit"


fit.stan <- rstan::sampling(
  m, data = list(n = n, y = y, H = H), pars = c("alpha", "lambda", "w"),
  iter = 1000, chains = 8, thin = 5, control = list(adapt_delta = 0.8)
)

a <- stan2coda(fit.stan)
fit.ggs <- ggs(a)
p1 <- ggs_traceplot(fit.ggs, family = "lambda") +
  ggplot2::coord_cartesian(ylim = c(0, 20)) +
  ggplot2::theme_bw()

p2 <- ggs_density(fit.ggs, family = "lambda") +
  ggplot2::coord_cartesian(xlim = c(0, 20)) +
  ggplot2::theme_bw()

ggsave("figure/mcmc1.pdf", p1, width = 7.5, height = 2.5)
ggsave("figure/mcmc2.pdf", p2, width = 7.5, height = 2.5)
