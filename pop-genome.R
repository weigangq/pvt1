# Use PopGenome
library(PopGenome)
g = readVCF("pvt1.recode.vcf.gz", 1000, "8", 127890000, 128102000)
pops = split(sample[,1], sample[,2]) # create a list of populations
g = set.populations(g, pops, diploid = T) # set population names

# by windows
slide = sliding.window.transform(g, width = 100, jump = 20) # nsnps, not actual length
slide = F_ST.stats(slide, mode = "nucleotide")
snp.pos = slide@region.data@biallelic.sites # SNP positions
win.num = length(slide@region.names)
win.start = numeric()
for (i in 1:win.num) {win.start[i] = snp.pos[[i]][1]}
fst = slide@nuc.F_ST.vs.all 
pop.names = names(slide@populations) # population names
plot(win.start, fst[,1], type ="n", las = 1, ylab = expression(F[st]), xlab = "SNP Position", ylim = c(0,0.4))
for (i in 1:length(slide@populations)) {
  lines(win.start, fst[,i], type = "l", col = pop.group[pop.names[i],4])
}
arch.coords=c(127982050, 127992931)
abline(v = arch.coords, col = "orange")
#rect(xleft = arch.coords[1], ybottom = -1, xright = arch.coords[2], ytop = 0.5, border = "transparent", col = 2)



