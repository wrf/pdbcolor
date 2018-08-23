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
pretty_ymax = max(pretty(ymax))

is_strong = ( first_v_second_lnl >= 0.5 )

t1strong = first_v_second_lnl[is_strong][(rowmaxs==1)[is_strong]]
t2strong = first_v_second_lnl[is_strong][(rowmaxs==2)[is_strong]]
t3strong = first_v_second_lnl[is_strong][(rowmaxs==3)[is_strong]]

strongsites_labs = (first_v_second_lnl >= (ymax/2))

if (!is.na(plottitle)) {
mainlab = plottitle
} else {
mainlab = inputfile
}

outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)
pdf(outputfile, width=7+numsites/1000, height=7)
par(mar=c(1,4.5,3,1))
plot(xpositions, first_v_second_lnl, type="h", pch=16, col=colorset_clear[rowmaxs], lwd=150/numsites, xlim=c(0.03,0.97)*numsites, ylim=c(0,pretty_ymax), ylab=expression(paste(Delta,"ln(L)",sep="")), xlab="", main=mainlab, axes=FALSE, cex.lab=1.4)
axis(2, cex.axis=1.3)
#axis(1, at=seq(0,1700,100), labels=seq(0,1700,100) , cex.axis=1.3)
abline(h=c(0.5), lwd=1, lty=2)
text(0,0.9*pretty_ymax,paste("Favors T1 (",length(t1strong),")"), cex=2, col="#fc8d62", pos=4)
text(0,0.8*pretty_ymax,paste("Favors T2 (",length(t2strong),")"), cex=2, col="#66c2a5", pos=4)
text(0,0.7*pretty_ymax,paste("Favors T3 (",length(t3strong),")"), cex=2, col="#377eb8", pos=4)
text( xpositions[strongsites_labs], first_v_second_lnl[strongsites_labs]+0.2, sitelnldat[,1][strongsites_labs], col=colorset_dark[rowmaxs][strongsites_labs])

dev.off()

#