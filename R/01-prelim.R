## ---- prelim ----
library(iprobit)
library(iprior)
library(ggplot2)
gg_colour_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# library(gganimate)
library(animation)
library(reshape2)
library(directlabels)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}
library(ggmcmc)
library(coda)

source("06-class-sim.R")

# Function to specify decimal places
decPlac <- function(x, k = 2) format(round(x, k), nsmall = k)

# Function to determine even numbers
isEven <- function(x) x %% 2 == 0
