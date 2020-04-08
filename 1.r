
# IRR Reporting in L2 Meta-Analyses (2014-2019):

i <- read.csv("https://raw.githubusercontent.com/hkil/m/master/L2.csv")$IRR

y <- table(i)/length(i)

plot(y, xlim = c(-.2, 4.2), ylim = c(0, .58), lend = 1, lwd = 10, 
     ylab = "Percentage of Use", xlab = "Type of IRR", font.lab = 2, 
     xaxt = "n", yaxt ="n", panel.f = abline(h = 0, col = 8), type = "h", 
     panel.l = axis(2, at = at <- axTicks(2), labels = paste0(at*1e2, "%"), las = 1))

text(0:4, y, paste0(round(y, 4)*1e2, "%"), font = 2, pos = 3, col = 4)
axis(1, at = 0:4, labels = c("NA", "%Agreement", "Kappa", "Mixed", "Unclear"))
