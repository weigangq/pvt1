sample <- read.table("sample.dat", header = T) # sample	pop	super_pop	gender		
group99 = read.table("group.txt2", col.names=c("id","allale", "group")) # HG01767	1	g_1
group99$super_pop = sample$super_pop[group99$id]
group99$pop = sample$pop[group99$id]

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
