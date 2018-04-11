#!/usr/bin/env Rscript
# sitewise ml plot
# created 2018-04-11

args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]
plottitle = args[2]

#inputfile = "RAxML_perSiteLLs.Whelan_D16_Opisthokonta_site_lk.tab"

sitelnldat = read.table(inputfile, header=TRUE, sep="\t")

colorset_clear = c("#fc8d62aa", "#66c2a5aa", "#377eb888")
colorset_dark = c("#fc8d62", "#66c2a5", "#377eb8")

t1lnl = sitelnldat[,2]
t2lnl = sitelnldat[,3]
t3lnl = sitelnldat[,4]

numsites = length(t1lnl)

rowmaxs = apply( sitelnldat[,2:4], 1, which.max)

minlnl = apply( sitelnldat[,2:4], 1, min)
maxsitelnl = apply( sitelnldat[,2:4], 1, max)
medianln = apply( sitelnldat[,2:4], 1, median )
first_v_second_lnl = maxsitelnl - medianln

xpositions = 1:length(sitelnldat[,1])
ymax = max(first_v_second_lnl)

stronggenes = (first_v_second_lnl >= (ymax/2))

if (!is.na(plottitle)) {
mainlab = plottitle
} else {
mainlab = inputfile
}

outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)
pdf(outputfile, width=7+numsites/1000, height=7)
par(mar=c(1,4.5,3,1))
plot(xpositions, first_v_second_lnl, type="h", pch=16, col=colorset_clear[rowmaxs], lwd=150/numsites, xlim=c(0.03,0.97)*numsites, ylim=c(0,max(pretty(ymax))), ylab=expression(paste(Delta,"ln(L)",sep="")), xlab="", main=mainlab, axes=FALSE, cex.lab=1.4)
axis(2, cex.axis=1.3)
#axis(1, at=seq(0,1700,100), labels=seq(0,1700,100) , cex.axis=1.3)
abline(h=c(0.5), lwd=1, lty=2)
text(0,28,"Favors T1", cex=2, col="#fc8d62", pos=4)
text(0,25,"Favors T2", cex=2, col="#66c2a5", pos=4)
text(0,22,"Favors T3", cex=2, col="#377eb8", pos=4)
text( xpositions[stronggenes], first_v_second_lnl[stronggenes]+0.2, sitelnldat[,1][stronggenes], col=colorset_dark[rowmaxs][stronggenes])

dev.off()

#