x = read.table("total.fst", na.strings = c("na"))
y = as.matrix(x)

#y = ifelse(abs(y)<=abs(min(y, na.rm = T)), NA, y)
y = ifelse(abs(y)<=0.01, NA, y)
hist(y, br=100)

filter = apply(y, 1, function(x) sum(x, na.rm=T)==0)
length(which(!filter))
y = y[which(!filter),]

#filter = apply(y, 2, function(x) sum(x, na.rm=T)==0)
#length(which(!filter))

z = as.data.frame(y, row.names=F)
z$pos = as.numeric(rownames(y))
long = reshape(z, varying=1:10, v.names="fst", timevar="group_pair", times=names(z)[1:10], direction="long")
long = long[,c(2,1,3)]
row.names(long)=NULL
long = long[!is.na(long$fst),]

library(lattice)
xyplot(fst ~ pos | group_pair, data=long, type="l", layout=c(2,5))



x = read.table("afr-other.fst")
colnames(x) = c("pos", "fst")
y=x[which(x$fst>0.01),]
#xyplot(y$fst~y$pos, type = 'l')

z = scan("afr-other-exon.txt", what = integer())

fst_afr = x
fst_afr.sig = y
fst_afr.sig.exon = z

plot(fst_afr.sig$pos,fst_afr.sig$fst, col = "gray", las =1, xlab="SNP position", ylab = expression(F[st]), main = "PVT1 SNPs with Fst > 0.01 for AFR vs non-AFR", cex.main = 0.8)
points(fst_afr.sig.exon,fst_afr.sig[which(fst_afr.sig$pos %in% fst_afr.sig.exon),2], col="magenta")

#pi by windows
x = read.table("pvt1.windowed.pi", head=T)
x$pos = x$BIN_START
x = x[,c(6,5)]
plot(x$pos, x$PI*100, type = 'l', col="gray", las=1, xlab = "SNP position", ylab = "PI*100", cex.main = 0.8)

points(fst_afr.sig$pos, fst_afr.sig$fst, col=2, cex=0.5)

# pi by stits:
x = read.table("pvt1.sites.pi", head=T)
y = x[which(x$POS %in% fst_afr.sig$pos),]

plot(y$POS, y$PI, col="gray", las=1, xlab = "SNP position", ylab=expression(pi), cex=0.5, ylim=c(-0.55,0.55), pch =16)
points(fst_afr.sig$pos, -1*fst_afr.sig$fst, col="blue", cex=0.5, pch =16)


exon = read.table("exon.txt")

plot(y$POS, y$PI, las=1, xlab = "SNP position", ylab=expression(pi), ylim=c(-0.55,0.55), type='n')

rect(exon[,1],-0.55,exon[,2],0.55, col=2, border=NA)

points(fst_afr.sig$pos, -1*fst_afr.sig$fst, col="blue", cex=0.5, pch =16)

points(y$POS, y$PI, col="gray", cex=0.5, pch =16)



x= read.table("fst-pi-count", sep="\t")
plot(x[,3], x[,2], xlab= expression(pi), ylab=expression(F[st]), col="magenta")
points(x[,3], x[,4], col = "green")
abline(v=0.1, lty = 2, col=2)
abline(h=0.1, lty = 2, col=2)



x = read.table("myc.percent")
max = max(x)
plot(rownames(x), x$AFR, cex=0.5, col='red', ylim = c(0,max))
points(rownames(x), x$AMR, cex=0.5, col='green')
points(rownames(x), x$EAS, cex=0.5, col='blue')
points(rownames(x), x$EUR, cex=0.5, col='cyan')
points(rownames(x), x$SAS, cex=0.5, col='magenta')


x=read.table("myc-dot")
x$pi=x$PI*100
x$fst = x$WEIR_AND_COCKERHAM_FST *500
x=x[,c(1:6, 9,10)]

x=x[which(x$fst>30),]
max = max(x)
plot(rownames(x), x$fst, cex=0.5, ylim = c(0,max), type ="b")
points(rownames(x), x$AFR, cex=0.5, col='red')
points(rownames(x), x$AMR, cex=0.5, col='green')
points(rownames(x), x$EAS, cex=0.5, col='blue')
points(rownames(x), x$EUR, cex=0.5, col='cyan')
points(rownames(x), x$SAS, cex=0.5, col='magenta')

x=read.table("myc-dot")
y=read.table("myc-afr.eur.fst", row.names=1)
z=intersect(rownames(x), rownames(y))

xx = data.frame(overall=x[z,8], afr_eur=y[z,1], id=y[z,2], row.names=z)
range = range(xx[,c(1,2)])
plot(rownames(xx), xx$overall, cex=0.5, ylim = range)
points(rownames(xx), xx$afr_eur, cex=0.5, col='red')

yy = xx[order(xx$afr_eur),]


#hardy in archaic
x = read.table("arc.hwe")
y = as.matrix(x)
y = -log10(y)
range = range(y)
x = as.data.frame(y)

plot(rownames(x), x$AFR, cex=0.5, ylim = range, type ="b")
points(rownames(x), x$AMR, cex=0.5, type="b", col=2)
points(rownames(y), x$EAS, cex=0.5, type="b", col=3)
points(rownames(y), x$EUR, cex=0.5, type="b", col=4)
points(rownames(y), x$SAS, cex=0.5, type="b", col=5)


#sample
sample = read.table("sample.dat", header = T, row.names=1)
pop_group = read.table("pop-group.txt")


#heatmap for snps:
library(gplots)
frq_pop = read.table("archaic/frq-pop.txt", header=T, row.names=1, check.names=F)
pop_id = read.table("archaic/pos-id", row.names=1)

x = as.data.frame(t(as.matrix(frq_pop)))
rownames(x) = pop_id[rownames(x),1]
x = t(as.matrix(x))

#hc.pop=hclust(as.dist(1-cor(x, use="pairwise.complete.obs")))
#hc.seq=hclust(as.dist(1-cor(t(x), use="pairwise.complete.obs")))
rlab=as.character(pop$color)

heatmap.2(x, scale="none", Colv=F, col=colorpanel(1000,"white","darkcyan"), RowSideColors=rlab, cexCol=0.6, cexRow=0.6, trace="none")

#hc.rows = hclust(dist(x))
#plot(hc.rows)


# hardy
hardy = read.table("archaic_26pop.hwe")


group99$mg = ifelse(group99$group %in% tail(t$tip.label,19), "archaic", "cosmo")
write.table(group99, "group99.txt", col.names = F, row.names=F, sep="\t", quote = F)

hardy_mg = read.table("count-archaic-cosmo.txt", header = T, row.names=1)

hardy <- function(x) {
  n_ind = sum(x);
  freq_ind = x/n_ind;
  freq_gene = c((x[1]*2+x[2]), (x[2]+x[3]*2))/(2*n_ind)
  freq_exp = c(freq_gene[1]^2, freq_gene[1]*freq_gene[2]*2, freq_gene[2]^2)
  count_exp = freq_exp * n_ind
  diff = (x-count_exp)^2/count_exp
  p.value.under = pchisq(sum(diff), df = 2)
  p.value.over = pchisq(sum(diff), df = 2, lower.tail = F)
}
apply(hardy_mg, 1, function(x) hardy(x))


## fst by site using vcftools, 26 populations
fst = read.table("pvt1.fst.site.vs.others.vcftools.26pop", header = T, row.names=1)

yrange = c(0.02,range(fst)[2])

plot(rownames(fst), fst[,1], las=1, ylab = expression(F[st]), xlab = 'SNP Position', ylim =yrange , type="n")

for (i in 1:dim(exon)[1]) {
  rect(exon[i,1], yrange[1], exon[i,2], yrange[2], col="gray", border="gray")
  text((exon[i,1]+exon[i,2])/2, yrange[2]-0.005, rownames(exon)[i], pos=3, cex=0.7, col="gray")
}

legend("topleft", rownames(super), pch=1, col=super[,2], cex=0.8, ncol=5)

for (i in 1:dim(fst)[2]){points(rownames(fst), fst[,i], cex=0.3, col=pop[colnames(fst)[i],2])}

# export 6x14

fst_site_vs_all_vcftools_26pop = fst
