#!/usr/bin/env Rscript
# genewise ml plot
# created 2018-03-02

args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]
plottitle = args[2]

#inputfile = "whelan-d16o_cat-gtr_genewise_logl.tab"

genelnldat = read.table(inputfile, header=TRUE, sep="\t")

colorset_clear = c("#fc8d62aa", "#66c2a5aa", "#377eb888")
colorset_dark = c("#fc8d62", "#66c2a5", "#377eb8")

t1lnl = genelnldat[,2]
t2lnl = genelnldat[,3]
t3lnl = genelnldat[,4]

numgenes = length(t1lnl)

rowmaxs = apply( genelnldat[,2:4], 1, which.max)

minlnl = apply( genelnldat[,2:4], 1, min)
maxsitelnl = apply( genelnldat[,2:4], 1, max)
medianln = apply( genelnldat[,2:4], 1, median )
first_v_second_lnl = maxsitelnl - medianln

xpositions = 1:length(genelnldat[,1])
ymax = max(first_v_second_lnl)

strongsites = (first_v_second_lnl >= (ymax/2))

if (!is.na(plottitle)) {
mainlab = plottitle
} else {
mainlab = inputfile
}

outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)
pdf(outputfile, width=7+numgenes/100, height=7)
par(mar=c(6,4.5,3,1))
plot(xpositions, first_v_second_lnl, type="h", col=colorset_clear[rowmaxs], lwd=250/numgenes, xlim=c(0.01,0.99)*numgenes, ylim=c(0,max(pretty(ymax))), ylab=expression(paste(Delta,"ln(L)",sep="")), xlab="", main=mainlab, axes=FALSE, cex.lab=1.4)
axis(2, cex.axis=1.3)
#axis(1, at=seq(0,1700,100), labels=seq(0,1700,100) , cex.axis=1.3)
abline(h=c(5), lwd=1, lty=2)
text(0,28,"Favors T1", cex=2, col="#fc8d62", pos=4)
text(0,25,"Favors T2", cex=2, col="#66c2a5", pos=4)
text(0,22,"Favors T3", cex=2, col="#377eb8", pos=4)
mtext(genelnldat[,1],side=1,at=c(1:numgenes),las=2, cex=0.7)
text( xpositions[strongsites], first_v_second_lnl[strongsites]+1, genelnldat[,1][strongsites], col=colorset_dark[rowmaxs][strongsites])

dev.off()

#