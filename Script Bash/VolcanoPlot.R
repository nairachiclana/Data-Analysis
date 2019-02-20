options(warn=0)

library(calibrate)

load(file="savedData/DEA.RData")

args = commandArgs(trailingOnly=TRUE)
lfc = as.numeric(args[1])
pval = as.numeric(args[1])


tab = data.frame(logFC=diff.exp.df$logFC, negLogPval= -log10(diff.exp.df$adj.P.Val))
tab2 = data.frame(logFC=diff.exp.df$logFC, negLogPval= -log10(diff.exp.df$adj.P.Val), Gene=diff.exp.df$gene.name)

par(mar = c(5, 4, 4, 5))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
points(tab[(abs(tab$logFC) > lfc), ], pch = 16, cex = 0.8, col = "orange") 
points(tab[(tab$negLogPval > -log10(pval)), ], pch = 16, cex = 0.8, col = "green") 
points(tab[(abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval)), ], pch = 16, cex = 0.8, col = "red") 
abline(h = -log10(pval), col = "green3", lty = 2) 
#sobreexpresion lines (separate dysregulated genes)
abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
with(subset(tab2, negLogPval > -log10(pval) & abs(logFC)>lfc), textxy(logFC, negLogPval, labs=Gene, cex=.4))

cat("Se ha realizado el volcano plor para los datos lfc=", lfc, "y pval=", pval, "\n")
