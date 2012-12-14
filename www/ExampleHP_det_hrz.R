library(highriskzone)

data(craterA)
spatstat.options(npixel=400)
## type: dist
hrzd1 <- det_hrz(craterA, type = "dist", criterion = "area", cutoff = 1000000, nxprob = 0.1)
hrzd2 <- det_hrz(craterA, type = "dist", criterion = "indirect", cutoff = 0.9, nxprob = 0.1)
hrzd3 <- det_hrz(craterA, type = "dist", criterion = "direct", cutoff = 100, nxprob = 0.1)

op <- par(mfrow = c(2, 2))
plot(craterA)
plot(hrzd1, zonecol = 2, win = craterA$window, plotwindow = TRUE)
plot(hrzd2, zonecol = 3,  win = craterA$window, plotwindow = TRUE)
plot(hrzd3, zonecol = 4,  win = craterA$window, plotwindow = TRUE)
par(op)

# or first calculate the distancemap and use it:
distm <- distmap(craterA)
hrzd <- det_hrz(craterA, type = "dist", criterion = "direct", cutoff = 100,
                distancemap = distm, nxprob = 0.1)
plot(hrzd, zonecol = 4,  win = craterA$window, plotwindow = TRUE)

## type: intens
hrzi1 <- det_hrz(craterA, type = "intens", criterion = "area", cutoff = 1000000, nxprob = 0.1)
hrzi2 <- det_hrz(craterA, type = "intens", criterion = "indirect", cutoff = 0.1, nxprob = 0.1)
hrzi3 <- det_hrz(craterA, type = "intens", criterion = "direct", cutoff = 0.0001, nxprob = 0.1)

op <- par(mfrow = c(2, 2))
plot(craterA)
plot(hrzi1, zonecol = 2, win = craterA$window, plotwindow = TRUE)
plot(hrzi2, zonecol = 3,  win = craterA$window, plotwindow = TRUE)
plot(hrzi3, zonecol = 4,  win = craterA$window, plotwindow = TRUE)
par(op)

# or first estimate the intensity and then use it:
intensity <- est_intens(craterA)
plot(intensity$intensest, main = "")

hrzii <- det_hrz(craterA, type = "intens", criterion = "area", cutoff = 1000000, 
                 intens = intensity$intensest, nxprob = 0.1)
plot(hrzii, zonecol = 4,  win = craterA$window, plotwindow = TRUE, main = "")

# or give covariance-matrix:
covmatrix <- Hscv(cbind(craterA$x, craterA$y))
hrzic <- det_hrz(craterA, type = "intens", criterion = "area", cutoff = 1000000, 
                 nxprob = 0.1, covmatrix = covmatrix)
plot(hrzic, zonecol = 4,  win = craterA$window, plotwindow = TRUE)
