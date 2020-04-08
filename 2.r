
source("https://raw.githubusercontent.com/hkil/m/master/s.r")

dat <- data.frame(Coder1 = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0),
                  Coder2 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

lst <- lapply(0:8, function(x) {dat[0:x, ] <- 0; dat })
ESL <- sapply(lst, sum)
EFL <- sum(lengths(dat)) - ESL

inter <- sapply(lst, irr)

plot(1:9, inter[1,], type = "b", ylim = c(min(inter[1,]),  0.7), ylab = "Kappa", panel.f = abline(h = 0, col = 8), xaxt = "n", yaxt = "n",
     xlab = "ESL:EFL Frequency", font.lab = 2, pch = 19, las = 1, mgp = c(1.5, .55, 0), cex.axis = .7, cex.lab = .7)
axis(1, at = 1:9, labels = paste0("(", ESL, ":", EFL, ")"), cex.axis = .7, mgp = c(1, .3, 0))
axis(2, at = c(-.1, 0:6*.1), cex.axis = .7, mgp = c(1.5, .55, 0), las = 1, padj = .4)