library(ggplot2)
library(UpSetR)
library(progress)
require(compiler)
library(gtools)

load("~/global.nt.RData")          #  This dataset is downloadable in the repository.
load("~/list.matrix.org.RData")    #  These files are generated in SupplementaryPhaseSync_Step1.R
load("~/list.matrix.org.wl.RData") #  These files are generated in SupplementaryPhaseSync_Step1.R
load("~/list.matrix.RData")        #  These files are generated in SupplementaryPhaseSync_Step1.R
global.nt <- global.nt.RData
#######################################################################################################################

change.name.obj <- function(objname, newname)
{
  
  assign(newname, get(objname), envir = .GlobalEnv)
  remove(objname, envir = .GlobalEnv)
  
}

saveit <- function(..., string, file) {
  
  x <- list(...)
  names(x) <- string
  save(list=names(x), file=file, envir=list2env(x))
  
}

compute.phase.sync <- function(p.on, v.r)
{
  n.osc.net <- colnames(p.on)
  p.on <-  p.on[rowSums(p.on==0, na.rm=TRUE)<ncol(p.on), ]
  v.r  <-  v.r[rowSums(v.r==0, na.rm=TRUE)<ncol(v.r), ]
  p.on1 <-  p.on[,-(which(colSums(p.on ) == 0))] 
  v.r  <-  v.r[, -(which(colSums(p.on ) == 0))] 
  p.on <- p.on1
  #View(v.r)
  #View(p.on)
  l <-  dim(p.on)[1]
  v.m.m      <- matrix(0, l, 3)
  v.m.m.tot  <- matrix(0, l, 3)
  v.m.md     <- matrix(0, l, 3)
  v.m.md.tot <- matrix(0, l, 3)
  v.m.s      <- matrix(0, l, 3)
  v.m.s.tot  <- matrix(0, l, 3)
  
  
  for(i in 1:l)
  {
    v.m.s[i,1]   <- i
    v.m.m[i,1]   <- i
    v.m.md[i,1]  <- i
    v.m.s[i,3]   <- as.numeric(n.osc.net[i])
    v.m.m[i,3]   <- as.numeric(n.osc.net[i])
    v.m.md[i,3]  <- as.numeric(n.osc.net[i])
    
    
    v.m.m[i,2]   <- mean(v.r[i,  p.on[i,] > 0])
    v.m.md[i,2]  <- median(v.r[i,  p.on[i,] > 0])
    v.m.s[i,2]   <- sd(v.r[i,  p.on[i,] > 0])
    
    v.m.s.tot[i,1] <- i
    v.m.m.tot[i,1] <- i
    v.m.md.tot[i,1] <- i
    v.m.s.tot[i,3]   <- as.numeric(n.osc.net[i])
    v.m.m.tot[i,3]   <- as.numeric(n.osc.net[i])
    v.m.md.tot[i,3]  <- as.numeric(n.osc.net[i])
    
    v.m.m.tot[i,2] <- mean(v.r[i, p.on[i,] <= 0])
    v.m.md.tot[i,2] <- median(v.r[i, p.on[i,] <= 0])
    v.m.s.tot[i,2] <- sd(v.r[i, p.on[i,] <= 0])
  }
  
  v.m.s.tot   <-  v.m.s.tot[!is.na(v.m.s[,2]),]
  v.m.m.tot   <-  v.m.m.tot[!is.na(v.m.s[,2]),]
  v.m.md.tot  <-  v.m.md.tot[!is.na(v.m.s[,2]),]
  v.m.m       <-  v.m.m[!is.na(v.m.s[,2]),]
  v.m.md      <-  v.m.md[!is.na(v.m.s[,2]),]
  v.m.s       <-  v.m.s[!is.na(v.m.s[,2]),]

  
  v.m.s.tot   <-  v.m.s.tot[!is.na(v.m.s[,3]),]
  v.m.m.tot   <-  v.m.m.tot[!is.na(v.m.s[,3]),]
  v.m.md.tot  <-  v.m.md.tot[!is.na(v.m.s[,3]),]
  v.m.m       <-  v.m.m[!is.na(v.m.s[,3]),]
  v.m.md      <-  v.m.md[!is.na(v.m.s[,3]),]
  v.m.s       <-  v.m.s[!is.na(v.m.s[,3]),]
  
  return(list(v.m.m     = v.m.m,  v.m.m.tot  = v.m.m.tot,
              v.m.md    = v.m.md, v.m.md.tot = v.m.md.tot,
              v.m.s     = v.m.s,
              v.m.s.tot = v.m.s.tot))
}



compute.phase.sync.compiled <- cmpfun(compute.phase.sync)


l.org.n <- c("bce:", "bsu:", "bth:", "cac:", "cje:", "eco:", "hpy:", "mtu:", "pae:", "sme:", "stm:")

#Example plot for  a within organisms velocity of propagation boxplot of the i-th organism.

phase.v1 <- compute.phase.sync.compiled(list.matrix.org[[1]][[1]][[1]], list.matrix.org.wl[[1]][[1]][[1]])
stat.freq.get <- function(A,B){
  stat.freq.osc <- c(mean(as.numeric(A[,2])), median(as.numeric(A[,2])), sd(as.numeric(A[,2])),
                     mean(as.numeric(B[,2])), median(as.numeric(B[,2])), sd(as.numeric(B[,2])))
  names(stat.freq.osc) <- c("m", "md", "sd", "m.ps", "md.ps", "sd.ps")
  
  return(stat.freq.osc)
}

#plot the box plot
A = data.frame(x = phase.v1$v.m.m.tot[,1], y= phase.v1$v.m.m.tot[,2], KEGG.id = phase.v1$v.m.m.tot[,3])
B = data.frame(x = phase.v1$v.m.m[,1], y= phase.v1$v.m.m[,2], KEGG.id = phase.v1$v.m.m[,3] )
v1 <- stat.freq.get(A,B)
A["ps"] = 0
B["ps"] = 1
df0 = rbind(A,B)
df0["Label"] = gsub(":","",l.org.n)[[1]]

fun_mean <- function(x){ return(round(data.frame(y=mean(x),label=mean(x,na.rm=T)), digits=3))}

p1 = ggplot(df0, aes(x=ps, y=y, group=ps)) + 
  geom_boxplot(aes(fill=Label)) +
  #facet_wrap( ~ Label, scales="free") +
  coord_flip() +
  stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.6, hjust=1) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position = "none") #+ geom_text(aes(label = KEGG.id), na.rm = TRUE, vjust=+6, hjust=1) #outliers KEGG.id
p1

