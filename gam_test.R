library(ggplot2)
library(mgcv)
library(dHSIC)

# Random:
x <- rnorm(n=100)
y <- rnorm(n=100)
Sample_data <- data.frame(y,x)
ggplot(Sample_data, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x))
gam_y <- gam(y ~ s(x), method = "REML")
cat(summary(gam_y)$s.pv)


# y = sin(x) + err
x <- seq(0, pi * 4, 0.1)
sin_x <- sin(x)
y <- sin_x + rnorm(n = length(x), mean = 0, sd = sd(sin_x))
Sample_data <- data.frame(y,x)
ggplot(Sample_data, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x))
gam_y <- gam(y ~ s(x), method = "REML")
cat(summary(gam_y)$s.pv)

dh <- dhsic.test(x,y)
cat(dh$p.value)
