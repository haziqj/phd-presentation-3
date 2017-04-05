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
isEven <- function(x) x %% 2 == 0
