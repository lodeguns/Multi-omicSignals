####################################################################
# Multi-omic normalization and rescaling
require(igraph)
require(KernSmooth)
require(grDevices)
require(moments)
require(progress)
require(classInt)
require(nortest) # https://www.r-bloggers.com/normality-tests-for-continuous-data/

work.dir <- "~/Data Comb/data_extracted/"
work.dir2 <- "~/Data Extraction/"
work.dir3 <- ""
load("~/Data Extraction/lorg.RData")
source(paste(work.dir3, "FunzioniAccessorie.R", sep = ""))
selected.org <- c(2, 3, 4, 7, 8, 10, 11, 12, 13, 14, 15, 17)



l.org[selected.org, ]

change.name.obj <- function(objname, newname)
{
  
  assign(newname, get(objname), envir = .GlobalEnv)
  remove(objname, envir = .GlobalEnv)
  
}

doanes <- function(data)
{  require(moments)
   n <- length(data)
   Doane <- 1 + log(n) + log(1 + abs(sqrt((6*(n-2)) /((n+1)*(n+3))) / skewness(data) ))
   Doane <- signif(x = Doane, digits = 1)
   return(Doane)
}


list.orgs <- list()

for (i in 1:length(selected.org))
{
  print(selected.org[[i]])
  gene_ass <- paste(work.dir, l.org[selected.org[[i]], 5], "_gene_assoc.RData", sep = "")
  paths <- paste(work.dir, l.org[selected.org[[i]], 5], "_pathways.RData", sep = "")
  operons <- paste(work.dir2, l.org[selected.org[[i]], 3], ".opr.Rdata", sep = "")
  expm <- paste(work.dir2, l.org[selected.org[[i]], 3], ".RData", sep = "")
  load.c <- c(gene_ass, paths, operons, expm)
  print(load.c)
  file.c <- c()
  
  for (h in 1:length(load.c))
  {
    file.c <- c(file.c, load(load.c[[h]]))
  }
  
  # df0 esperimenti df1 CAI, gene order, ids df2 operoni
  change.name.obj(file.c[[3]], "df2")
  change.name.obj(file.c[[4]], "df0")
  
  omic1 <- as.vector(scale(df1$cai))
  omic2 <- as.vector(scale(df1$mw))
  oon <- gsub(paste(l.org[selected.org[[i]], 5], ":", sep = ""), "", df1$gene_loc)
  names(omic1) <- oon
  names(omic2) <- oon
  
  omic3 <- scale(df0$exprdata)  ## questi sono più dati quindi scalo CAI e mw
  omic3r <- omic3[complete.cases(omic3), ]
  omic3r <- omic3r[rownames(omic3r) %in% names(omic1), ]
  
  omic1r <- omic1[names(omic1) %in% rownames(omic3r)]
  omic2r <- omic2[names(omic2) %in% rownames(omic3r)]
  oon <- names(omic1r)
  
  rr1 <- range(omic3, na.rm = TRUE, finite = FALSE)
  library(oce)
  omic1r <- rescale(omic1r, rlow = rr1[[1]], rhigh = rr1[[2]])
  omic2r <- rescale(omic2r, rlow = rr1[[1]], rhigh = rr1[[2]])
  names(omic1r) <- oon
  names(omic2r) <- oon
  
  list.intr <- list(CAI_sc = omic1r, MW_sc = omic2r, exp_sc = omic3r, gene_assoc = gene_ass, paths = paths, operons = operons, expm = expm)
  
  list.orgs[length(list.orgs) + 1] <- list(list.intr)
}

save(list.orgs, file = "list.orgs.RData")




#################################################################
# consensus tests for the levels of discretization



cai.test.result     <- c()
mw.test.result      <- c()
exp.test.result     <- c()
nclass.norm.cai     <- c()
nclass.norm.mw      <- c()
nclass.norm.exp     <- c()
nclass.not.norm.cai <- c()
nclass.not.norm.mw  <- c()
nclass.not.norm.exp <- c()
length.org <- c()

for (j in 1:length(list.orgs))
{ length.org <- c(length.org, length(list.orgs[[j]]$CAI_sc))
  ad.t1  <- ad.test(list.orgs[[j]]$CAI_sc)
  nc.t1  <- doanes(list.orgs[[j]]$CAI_sc)
  nc.t11 <- nclass.scott(list.orgs[[j]]$CAI_sc)
  ad.t2  <- ad.test(list.orgs[[j]]$MW_sc)
  nc.t2  <- doanes(list.orgs[[j]]$MW_sc)
  nc.t22 <- nclass.scott(list.orgs[[j]]$MW_sc)
  cai.test.result <- c(cai.test.result, ad.t1$p.value)
  mw.test.result  <- c(mw.test.result, ad.t2$p.value)
  nclass.not.norm.cai <- c(nclass.not.norm.cai, nc.t1)
  nclass.not.norm.mw  <- c(nclass.not.norm.mw, nc.t2)
  nclass.norm.cai <- c(nclass.norm.cai, nc.t11)
  nclass.norm.mw  <- c(nclass.norm.mw, nc.t22)
  
  for (k in 1:length(list.orgs[[j]]$exp_sc[1, ]))
  {
    ad.t3  <- ad.test(list.orgs[[j]]$exp_sc[, k])
    nc.t3  <- doanes(list.orgs[[j]]$exp_sc[, k])
    nc.t33 <- nclass.scott(list.orgs[[j]]$exp_sc[, k])
    exp.test.result <- c(exp.test.result, ad.t3$p.value)
    
    nclass.norm.exp     <- c(nclass.norm.exp, nc.t33)
    nclass.not.norm.exp <- c(nclass.not.norm.exp, nc.t3)
    
    
    
    
  }
  
  
}



alpha = 0.05
cai.test.result < 0.05
mw.test.result  < 0.05
exp.test.result < 0.05

prop.table(table(round(cai.test.result))) * 100
prop.table(table(round(mw.test.result))) * 100
prop.table(table(round(exp.test.result))) * 100

prop.table(table(nclass.not.norm.exp)) * 100
prop.table(table(nclass.not.norm.cai)) * 100
prop.table(table(nclass.not.norm.mw)) * 100


















