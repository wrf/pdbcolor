#!/usr/bin/env Rscript
# genewise ml plot
# created 2018-03-02
# last modified 2020-08-27

args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]
plottitle = args[2]

#inputfile = "whelan-d16o_cat-gtr_genewise_logl.tab"

genelnldat = read.table(inputfile, header=TRUE, sep="\t")

colorset_clear = c("#fc8d62aa", "#66c2a5aa", "#377eb888")
colorset_dark = c("#fc8d62", "#66c2a5", "#377eb8")

t1lnl = genelnldat[,2]
t2lnl = genelnldat[,3]

if (length(genelnldat) < 4) {
numtrees = 3
} else {
numtrees = 4
}

is_using_t3 = (length(genelnldat) == 4)


numgenes = length(t1lnl)

rowmaxs = apply( genelnldat[,2:numtrees], 1, which.max)

minlnl = apply( genelnldat[,2:numtrees], 1, min)
maxsitelnl = apply( genelnldat[,2:numtrees], 1, max)
medianln = apply( genelnldat[,2:numtrees], 1, median )
first_v_second_lnl = maxsitelnl - medianln

xpositions = 1:length(genelnldat[,1])
ymax = max(first_v_second_lnl)
pretty_ymax = max(pretty(ymax))

# 16 is based on control dataset
STRONGSCORE = 16
is_strong = ( first_v_second_lnl >= STRONGSCORE )

stronggenes_lab = (first_v_second_lnl >= (ymax/2))

t1strong = first_v_second_lnl[is_strong][(rowmaxs==1)[is_strong]]
t2strong = first_v_second_lnl[is_strong][(rowmaxs==2)[is_strong]]

if (is_using_t3) {
t3lnl = genelnldat[,4]
t3strong = first_v_second_lnl[is_strong][(rowmaxs==3)[is_strong]]
}


if (!is.na(plottitle)) {
mainlab = plottitle
} else {
mainlab = inputfile
}

outputfile = gsub("([\\w/]+)\\....$","\\1.genes.pdf",inputfile,perl=TRUE)
pdf(outputfile, width=7+numgenes/100, height=7)
par(mar=c(6,4.5,3,1))
plot(xpositions, first_v_second_lnl, type="h", col=colorset_clear[rowmaxs], lwd=250/numgenes, xlim=c(0.01,0.99)*numgenes, ylim=c(0,pretty_ymax), ylab=expression(paste(Delta,"ln(L)",sep="")), xlab="", main=mainlab, axes=FALSE, cex.lab=1.4)
axis(2, cex.axis=1.3)
#axis(1, at=seq(0,1700,100), labels=seq(0,1700,100) , cex.axis=1.3)
abline(h=c(STRONGSCORE), lwd=1, lty=2)
text(0, 0.9*pretty_ymax, paste("Favors T1 (",length(t1strong),"genes,",round(sum(t1strong),digits=1),")"), cex=2, col="#fc8d62", pos=4)
text(0, 0.8*pretty_ymax, paste("Favors T2 (",length(t2strong),"genes,",round(sum(t2strong),digits=1),")"), cex=2, col="#66c2a5", pos=4)
if (is_using_t3) {
text(0, 0.7*pretty_ymax, paste("Favors T3 (",length(t3strong),"genes,",round(sum(t3strong),digits=1),")"), cex=2, col="#377eb8", pos=4)
}
mtext(genelnldat[,1],side=1,at=c(1:numgenes),las=2, cex=0.7)
text( xpositions[stronggenes_lab], first_v_second_lnl[stronggenes_lab]+1, genelnldat[,1][stronggenes_lab], col=colorset_dark[rowmaxs][stronggenes_lab])

dev.off()



t1_hist = hist(first_v_second_lnl[(rowmaxs==1)], plot=FALSE)
t2_hist = hist(first_v_second_lnl[(rowmaxs==2)], plot=FALSE)
t3_hist = hist(first_v_second_lnl[(rowmaxs==3)], plot=FALSE)

outputfile = gsub("([\\w/]+)\\....$","\\1.hist.pdf",inputfile,perl=TRUE)
pdf(outputfile, width=8, height=7)
par(mar=c(4.5,4.5,2,2))
plot(0,0, type="n", xlim=c(0,max(first_v_second_lnl)), ylim=c(0,numgenes/3), xlab=expression(paste(Delta,"ln(L)",sep="")), ylab="Number of genes", main=mainlab, frame.plot=FALSE, cex.lab=1.4, cex.axis=1.4)
lines(t1_hist$mids, t1_hist$counts, lwd=4, col="#fc8d62")
lines(t2_hist$mids, t2_hist$counts, lwd=4, col="#66c2a5")
lines(t3_hist$mids, t3_hist$counts, lwd=4, col="#377eb8")
abline(v=c(STRONGSCORE), lwd=1, lty=2)
dev.off()


#