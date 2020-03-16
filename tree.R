library(tidyverse)

# count haplotypes in each supergroup
# 3/13/2020
sample <- read_tsv("sample.dat", col_names = T) # sample	pop	super_pop	gender		
group99 <- read_tsv("group.txt2", col_names=c("id","allale", "group")) # HG01767	1	g_1

group99 <- group99 %>% left_join(sample, c("id" = "sample")) 

hap.ct <- group99 %>% group_by(group, super_pop) %>% count()

hap.sum <- hap.ct %>% group_by(group) %>% summarise(sum=sum(n))
hap.50 <- hap.sum %>% filter(sum>=50)
hap.50 %>% ggplot(aes(x=reorder(group, sum), y=sum)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + xlab("haplotype") + ylab("total counts (out of 5008 chroms)") # "hap50"

hap.ct2 <- hap.ct %>% left_join(hap.sum, "group")
hap.ct2 <- hap.ct2 %>% mutate(freq = n/sum)
hap.ct2 %>% ggplot(aes(x=sum)) + geom_histogram(bins =10) + scale_x_log10()

# pick haps having at least 50 individuals:
# hap.ct3 <- hap.ct2 %>% filter(sum >=100)
# hap.ct.wide <- hap.ct3 %>% select(1,2,5) %>% spread(key = "super_pop", value = "freq", fill = 0)

hap.ct3<- hap.ct2 %>% filter(group %in% hap.50$group) 
hap.afr <- hap.ct3 %>% filter(super_pop == 'AFR') %>% arrange(freq)
hap.ct3$group <- factor(hap.ct3$group, levels = hap.afr$group)
hap.ct3 %>% ggplot(aes(x=group, y=freq, fill=super_pop)) +  geom_bar(stat = "identity", position = "stack") + coord_flip() + xlab("haplotype") + ylab("frequency") + theme_bw()

# HWE test for g_62
g62 <- c(rep("A/A", 19), rep("A/G", 180), rep("G/G", 462))
g62.geno <- genotype(g62)
HWE.test(g62.geno)

# HWE test for g_30
g30 <- c(rep("A/A", 12), rep("A/G", 571), rep("G/G", 1260))
g30.geno <- genotype(g62)
HWE.test(g30.geno)

################
# END 3/13/2020
################
#group99$super_pop = sample$super_pop[group99$id]
#group99$pop = sample$pop[group99$id]

count.group = table(group99$group, group99$super_pop)
count.group = as.data.frame.matrix(count.group)


library(ape)
t = read.tree("tree0.dnd")
plot(t)
library("phangorn", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
tm = midpoint(t)
#plot(t)
#plot(t, show.tip.label = F)
write.tree(tm, "tree.dnd")

t = read.tree("tree.dnd")

plot(t, show.tip.label=F, x.lim=2)

for (i in rownames(count.group)){
  x1 = getphylo_x(t, i)
  y1 = getphylo_y(t, i)
  segments(x1,y1, 1.5,y1, col="lightgray", lty=2)
}

for (i in 1:5){
  text(rep(0.8+0.15*(i-1), 99), 1:99, count.group[t$tip.label,i], cex=0.5, pos=2, col=i+1)
}

text(rep(1.55,99), 1:99, apply(count.group,1, sum)[t$tip.label], cex=0.5, pos=2)

text(c(0.8,0.95,1.1,1.25,1.4,1.55),rep(100.5,6), c(colnames(count.group),"total"), cex=0.5, pos=2, col="gray")

text(c(0.8,0.95,1.1,1.25,1.4,1.55),rep(-0.5,6), c(apply(count.group,2,sum), sum(apply(count.group,2,sum))), cex=0.5, pos=2, col="gray")

getphylo_x <- function(tree, node) {
  if(is.character(node)) {
    node <- which(c(tree$tip.label, tree$node.label)==node)
  }
  pi <- tree$edge[tree$edge[,2]==node, 1]
  if (length(pi)) {
    ei<-which(tree$edge[,1]==pi & tree$edge[,2]==node)
    tree$edge.length[ei] + Recall(tree, pi)
  } else {
    if(!is.null(tree$root.edge)) {
      tree$root.edge
    } else {
      0
    }
  }
}

getphylo_y <- function(tree, node) {
  if(is.character(node)) {
    node <- which(c(tree$tip.label, tree$node.label)==node)
  }
  ci <- tree$edge[tree$edge[,1]==node, 2]
  if (length(ci)==2) {
    mean(c(Recall(tree, ci[1]), Recall(tree, ci[2])))
  } else if (length(ci)==0) {
    Ntip <- length(tree$tip.label)
    which(tree$edge[tree$edge[, 2] <= Ntip, 2] == node)
  } else {
    stop(paste("error", length(ci)))
  }
}
