require(compiler)
library(gtools)
library(ggplot2)
library(UpSetR)
library(progress)

load("~/global.nt.RData")          #download this dataset from the repository
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


get.circuits.pathways <- function(score, org, j = 3, n.path=NULL, static.omics=TRUE, global.n)
{
  global.n <- global.nt
  if(static.omics==TRUE){
    print(paste("Org: ", org,  sep=" " ))
    code   <- c("n1", "n2", "n3", "n4", "o1", "o2", "o3", "o4")
    nomi   <- c("cai.mw", "cai.exp", "mw.exp", "cai.mw.exp", "opr.cai.mw", "opr.cai.exp", "opr.mw.exp", "opr.cai.mw.exp")
    print(paste("Multi-omic oscillation score threshold: ", score,  sep=" " ))
    print(paste("Dynamic/static  multi-omic combination: ", nomi[j], " : ", code[j],  sep=" " ))
    
  } else
  {
    print(paste("Org: ", org,  sep=" " ))
    code   <- c("n2", "n3", "n4", "o2", "o3", "o4")
    nomi   <- c("cai.exp", "mw.exp", "cai.mw.exp", "opr.cai.exp", "opr.mw.exp", "opr.cai.mw.exp")
    print(paste("Multi-omic oscillation score threshold: ", score,  sep=" " ))
    print(paste("Dynamic/static  multi-omic combination: ", nomi[j], " : ", code[j],  sep=" " ))
    
  }
  
  mask_j     <- global.n$score >= score & global.n$code == code[j] & global.n$kegg.id == org
  global.n1  <- global.n[mask_j, ]
  # dd         <- duplicated(global.n1) 
  # global.n1  <- global.n1[!dd,  ]
  
  paths.g    <- table(global.n$n.path)
  print(paths.g)
  paths.p    <- global.n1$n.path
  paths.t    <- table(paths.p) #all the paths >= score
  exps.p     <- paste(global.n1$exp.ref, global.n1$exp.ctr, sep = "|<---->|")
  class.p    <- global.n1$class
  func.p     <- global.n1$func
  map.p      <- global.n1$pathway_map
  
  exps.t     <- table(exps.p)  # all the exps            >= score
  class.t    <- table(class.p) # all the pathws class    >= score
  func.t     <- table(func.p)  # all the pathws funcs    >= score
  map.t      <- table(map.p)   # all the pathws maps     >= score
  
  paths.sum  <- sum(paths.t)
  
  
  
  mask_u     <- global.n$kegg.id == org & global.n$code == code[j]
  global.nu  <- global.n[mask_u, ]
  dd         <- duplicated(global.nu) 
  global.nu  <- global.nu[!dd,  ]
  exps.g     <- paste(global.nu$exp.ref, global.nu$exp.ctr, sep = "|<---->|")
  exps.g     <- table(exps.g)
  class.g    <- table(global.nu$class)
  func.g     <- table(global.nu$func)
  map.g      <- table(global.nu$pathway_map)
  
  
  
  
  
  
  exps.g.sum  <- sum(exps.g)
  class.g.sum <- sum(class.g)
  func.g.sum  <- sum(func.g)
  map.g.sum   <- sum(map.g)
  
  
  cov.exp.perc    <- length(names(exps.t))/length(names(exps.g)) * 100
  cov.exp.tot     <- sum(exps.t)/sum(exps.g) * 100
  
  cov.class.perc    <- length(names(class.t))/length(names(class.g)) * 100
  cov.class.tot     <- sum(class.t)/sum(class.g) * 100
  
  cov.func.perc    <- length(names(func.t))/length(names(func.g)) * 100
  cov.func.tot     <- sum(func.t)/sum(func.g) * 100
  
  cov.map.perc    <- length(names(map.t))/length(names(map.g)) * 100
  cov.map.tot     <- sum(map.t)/sum(map.g) * 100
  
  print("Oscillations %% of exp coverage")
  print(paste(" experiments:", cov.exp.tot))
  print(paste(" experiments (for tipology):", cov.exp.perc)) 
  print("Oscillations %% of pathways coverage for class, function and map")
  print(paste("pathways classes involved:", cov.class.tot))
  print(paste("pathways classes (for tipology):", cov.class.perc)) 
  print(paste("pathways functions involved:", cov.func.tot))
  print(paste("pathways functions (for tipology):", cov.func.perc)) 
  print(paste("pathways maps involved:", cov.map.tot))
  print(paste("pathways maps (for tipology):", cov.map.perc)) 
  
  # median value of the length of decay/window for the signal
  circuit.w    <- matrix(0, length(exps.p), length(paths.g))
  # median score absolute for each window = shifted half a period 
  circuit.s    <- matrix(0, length(exps.p), length(paths.g)) 
  # number of shifted half a period divided by the length of the step signal
  circuit.v    <- matrix(0, length(exps.p), length(paths.g))
  rownames(circuit.w ) <- exps.p
  colnames(circuit.w ) <- names(paths.g)
  
  rownames(circuit.s) <- exps.p
  colnames(circuit.s) <- names(paths.g)
  
  rownames( circuit.v ) <- exps.p
  colnames( circuit.v ) <- names(paths.g)
  
  for(k in 1:length(paths.g))
  {
    #pb$tick()
    
    mask_k     <- #global.n$score   >= score   & 
      global.n$code    == code[j] & 
      global.n$kegg.id == org     &
      global.n$n.path  == names(paths.g)[k] 
    
    #here
    
    global.nk  <- global.n[mask_k, ]
    dd         <- duplicated(global.nk) 
    global.nk  <- global.nk[!dd,  ]
    
    exp.c1     <- paste(global.nk$exp.ref, global.nk$exp.ctr, sep = "|<---->|")
    
    sdiff1 <- setdiff(exp.c1, rownames(circuit.w))
    w.out <- which(exp.c1 %in% sdiff1)
    
    if(length(sdiff1) == 0)
    {
      circuit.w[exp.c1, names(paths.g)[k]]      <- global.nk$med.w
    } else
    {
      exp.c <- exp.c1[-w.out]
      global.nl <- global.nk[-w.out, ]
      
      circuit.w[exp.c, names(paths.g)[k]]      <- global.nl$med.w
    }
    
    
    sdiff1 <- setdiff(exp.c1, rownames(circuit.v))
    w.out <- which(exp.c1 %in% sdiff1)
    
    if(length(sdiff1) == 0)
    {
      circuit.v[exp.c1, names(paths.g)[k]]      <- global.nk$v.change.w
    } else
    {
      exp.c <- exp.c1[-w.out]
      global.nl <- global.nk[-w.out, ]
      
      circuit.v[exp.c, names(paths.g)[k]]      <- global.nl$v.change.w
    }
    
    
    
    sdiff1 <- setdiff(exp.c1, rownames(circuit.s))
    w.out <- which(exp.c1 %in% sdiff1)
    
    if(length(sdiff1) == 0)
    {
      circuit.s[exp.c1, names(paths.g)[k]]      <- global.nk$med.s
    } else
    {
      exp.c <- exp.c1[-w.out]
      global.nl <- global.nk[-w.out, ]
      
      circuit.s[exp.c, names(paths.g)[k]]      <- global.nl$med.s
    }
    
    
    
  }
  
  
  
  
  
  
  
  
  
  
  paths.not.common <-  setdiff(names(paths.g), names(table(global.nu$n.path)))
  
  circuit.p <- matrix(0, length(exps.p), length(paths.g))
  rownames(circuit.p) <- exps.p
  colnames(circuit.p) <- names(paths.g)
  circuit.p[, paths.not.common] <- 0   #paths not comprised in the org path set. better -1? In the second function (phase sync) they are removed. 
  
  # pb <- progress_bar$new(
  #   format = "Processing circuit associations [:bar] :percent in :elapsed ",
  #   total = length(paths.t), clear =TRUE, width= 100)
  #  # Circuit associations between pathways and experiments
  
  
  
  for(k in 1:length(paths.t))
  {
    # pb$tick()
    # 
    mask_k     <- global.n$score   >= score   & 
      global.n$code    == code[j] & 
      global.n$kegg.id == org     &
      global.n$n.path  == names(paths.t)[k] 
    
    global.nk  <- global.n[mask_k, ]
    dd         <- duplicated(global.nk) 
    global.nk  <- global.nk[!dd,  ]
    
    exp.cc     <- paste(global.nk$exp.ref, global.nk$exp.ctr, sep = "|<---->|")
    circuit.p[exp.cc, names(paths.t)[k]]      <- circuit.p[exp.cc, names(paths.t)[k]] + 1
    
    
    
  }
  
  
  # Circuit associations between pathways and pathway classes
  circuit.c    <- matrix(0, length(exps.p), length(class.g))
  
  
  rownames(circuit.c) <- exps.p
  colnames(circuit.c) <- names(class.g)
  
  # pb <- progress_bar$new(
  #   format = "Processing circuit associations [:bar] :percent in :elapsed ",
  #   total = length(class.t), clear =TRUE, width= 100)
  # Circuit associations between pathways and experiments
  
  for(k in 1:length(class.t))
  {
    #pb$tick()
    
    mask_k     <- global.n$score   >= score   & 
      global.n$code    == code[j] & 
      global.n$kegg.id == org     &
      global.n$class  == names(class.t)[k] 
    
    global.nk  <- global.n[mask_k, ]
    dd         <- duplicated(global.nk) 
    global.nk  <- global.nk[!dd,  ]
    
    exp.cc     <- paste(global.nk$exp.ref, global.nk$exp.ctr, sep = "|<---->|")
    circuit.c[exp.cc, names(class.t)[k]]      <- circuit.c[exp.cc, names(class.t)[k]] + 1
    
    
  }
  
  
  ##pathways for the i-esim experiment.
  ##circuit.p[rownames(circuit.p)[[1]],which(c(circuit.p[rownames(circuit.p)[[1]], ]) != 0)]
  
  # Circuit associations between pathways and pathway functions
  circuit.f <- matrix(0, length(exps.p), length(func.g))
  rownames(circuit.f) <- exps.p
  colnames(circuit.f) <- names(func.g)
  # pb <- progress_bar$new(
  #   format = "Processing circuit associations [:bar] :percent in :elapsed ",
  #   total = length(func.t), clear =TRUE, width= 100)
  
  for(k in 1:length(func.t))
  {
    # pb$tick()
    # 
    mask_k     <- global.n$score   >= score   & 
      global.n$code    == code[j] & 
      global.n$kegg.id == org     &
      global.n$func  == names(func.t)[k] 
    
    global.nk  <- global.n[mask_k, ]
    dd         <- duplicated(global.nk) 
    global.nk  <- global.nk[!dd,  ]
    
    exp.cc     <- paste(global.nk$exp.ref, global.nk$exp.ctr, sep = "|<---->|")
    circuit.f[exp.cc, names(func.t)[k]]      <- circuit.f[exp.cc, names(func.t)[k]] + 1
    
    
    
  }
  
  
  
  # Circuit associations between pathways and pathway maps
  circuit.m <- matrix(0, length(exps.p), length(map.g))
  rownames(circuit.m) <- exps.p
  colnames(circuit.m) <- names(map.g)
  # pb <- progress_bar$new(
  #   format = "Processing circuit associations [:bar] :percent in :elapsed ",
  #   total = length(map.t), clear =TRUE, width= 100)
  for(k in 1:length(map.t))
  {
    # pb$tick()
    # 
    mask_k     <- global.n$score   >= score   & 
      global.n$code    == code[j] & 
      global.n$kegg.id == org     &
      global.n$pathway_map == names(map.t)[k] 
    
    global.nk  <- global.n[mask_k, ]
    dd         <- duplicated(global.nk) 
    global.nk  <- global.nk[!dd,  ]
    
    exp.cc     <- paste(global.nk$exp.ref, global.nk$exp.ctr, sep = "|<---->|")
    circuit.m[exp.cc, names(map.t)[k]]      <- circuit.m[exp.cc, names(map.t)[k]] + 1
    
    
    
  }
  
  
  
  return(list(pathways.on = circuit.p, class.on = circuit.c,
              function.on = circuit.f, 
              maps.on = circuit.m,
              med.path.window = circuit.w, 
              med.score.window = circuit.s,
              ratio.window.length = circuit.v))
  
  
  
  
}



get.circuits.pathways.compiled <- cmpfun(get.circuits.pathways)


################################################################################################
#### Here are computed and collected the intersections of active pahtways between organisms
################################################################################################
l.org.n <- c("bce:", "bsu:", "bth:", "cac:", "cje:", "eco:", "hpy:", "mtu:", "pae:", "sme:", "stm:")
list.matrix.org <- vector("list", 6)
list.circ       <- vector("list", 11)

list.ko.level.org <- vector("list", 4)


list.matrix.org.wl   <- vector("list", 6)
list.circ.wl         <- vector("list", 11)
list.ko.level.org.wl <- vector("list", 4)


for(moc in 1:6) #1:6
{ for(KO.LVL in 1:4) #1:3
  {
    score_setted = 0.8
    set.org <-  l.org.n[[1]]
    invisible(capture.output(circ.capt0
                             <- get.circuits.pathways.compiled(score  = score_setted,   
                                                               org    = set.org,  
                                                               moc, 
                                                               n.path = NULL, 
                                                               static.omics = FALSE, 
                                                               global.nt)))
    list.circ[1] <- list(circ.capt0[[KO.LVL]])
    list.circ.wl[1] <- list(circ.capt0$ratio.window.length)

    for(org in 2:11)
    {
      set.org <-  l.org.n[[org]]
      print(paste("moc:", moc, " org:", set.org, sep = " " ))
      
      
      invisible(capture.output(circ.capt
                               <- get.circuits.pathways.compiled(score  = score_setted,   
                                                                 org    = set.org,  
                                                                 moc, 
                                                                 n.path = NULL, 
                                                                 static.omics = FALSE, 
                                                                 global.nt)))
      
      
      list.circ[org] <- list(circ.capt[[KO.LVL]])
      list.circ.wl[org] <- list(circ.capt$ratio.window.length)

      
      
    }

    list.ko.level.org[KO.LVL] <- list(list.circ)
    list.ko.level.org.wl[KO.LVL] <- list(list.circ.wl)
  }
  
  list.matrix.org[moc] <- list(list.ko.level.org)
  list.matrix.org.wl[moc] <- list(list.ko.level.org.wl)
  
}


saveit(list.matrix.org = list.matrix.org, string = "list.matrix.org", file ="list.matrix.org.RData")
saveit(list.matrix.org.wl = list.matrix.org.wl, string = "list.matrix.org.wl", file ="list.matrix.org.wl.RData")
###################################################################################################################
###################################################################################################################
View(list.matrix.org[[1]][[1]])
View(list.matrix.org[[1]][[1]][[1]])

####
View(list.matrix.org[[1]][[1]])
View(list.matrix.org[[1]][[1]][[1]])




###################################################################################################################
###################################################################################################################
###################################################################################################################
# Here are computed and collected the intersections of active pahtways within organisms
l.org.n <- c("bce:", "bsu:", "bth:", "cac:", "cje:", "eco:", "hpy:", "mtu:", "pae:", "sme:", "stm:")


list.matrix     <- vector("list", 6) 
list.circ       <- vector("list", 11)
list.ko.level   <- vector("list", 3)


for(moc in 1:6) #1:6
{ 
  for(KO.LVL in 1:4) #1:3
  {
    score_setted = 0.8
    set.org <-  l.org.n[[1]]
    invisible(capture.output(circ.capt0
                             <- get.circuits.pathways.compiled(score  = score_setted,   
                                                               org    = set.org,  
                                                               moc, 
                                                               n.path = NULL, 
                                                               static.omics = FALSE, 
                                                               global.nt)))
    list.circ[1] <- list(circ.capt0[[KO.LVL]])
    list.circ.wl[1] <- list(circ.capt0$ratio.window.length)
    FLAG.capt = TRUE
    for(org in 2:11)
    {
      set.org <-  l.org.n[[org]]
      print(paste("moc:", moc, " org:", set.org, sep = " " ))
      
      
      invisible(capture.output(circ.capt
                               <- get.circuits.pathways.compiled(score  = score_setted,   
                                                                 org    = set.org,  
                                                                 moc, 
                                                                 n.path = NULL, 
                                                                 static.omics = FALSE, 
                                                                 global.nt)))
      
      # name.c.ko.lvl <- colnames(bind.circ.cai.mw.exp[[KO.LVL]])
      
      # mat1 <- matrix(ncol=length(name.c.ko.lvl), nrow=length(l.org.n ))
      
      list.circ[org] <- list(circ.capt[[KO.LVL]])
      
      if(FLAG.capt)
      {
        bind.circ1  <- smartbind(circ.capt0[[KO.LVL]], circ.capt[[KO.LVL]], fill=0)
        FLAG.capt = FALSE
      }
      else
      {
        bind.circ1  <- smartbind(bind.circ1 , circ.capt[[KO.LVL]], fill=0)
      }
      
      
    }
    bind.circ1 <- bind.circ1[apply(bind.circ1[,-1], 1, function(x) !all(x==0)),]
    list.ko.level[KO.LVL] <- list(bind.circ1)
  }
  
  list.matrix[moc] <- list(list.ko.level)
  
  
}

###################################################################################################################

saveit(list.matrix = list.matrix, string = "list.matrix", file ="list.matrix.RData")

# ####
# View(list.matrix.org[[1]][[1]])
# View(list.matrix.org[[1]][[1]][[1]])
# 
# ####
# View(list.matrix.org[[1]][[1]])
# View(list.matrix.org[[1]][[1]][[1]])



###
### Plot the UpSets
library(progress)
load("~/list.matrix.org.RData")
load("~/list.matrix.RData")

DIM <- function( ... ){
  args <- list(...)
  lapply( args , function(x) { if( is.null( dim(x) ) )
    return( length(x) )
    dim(x) } )
}


#View(list.matrix[[2]][[3]][[12]])
#View(list.matrix.org[[6]][[3]][[11]])
for(i in 1:6){
  for (j in 1:4){
    for (k in 1:11){
      set.test <- list.matrix.org[[i]][[j]][[k]][apply(list.matrix.org[[i]][[j]][[k]][,-1], 1, function(x) !all(x==0)),]
      pdf(file=paste("~/upset/upsets/upset_", i, "_", j, "_", k, ".pdf", sep=""), onefile=FALSE)
      ust <- upset(as.data.frame(set.test), nsets = 12 , main.bar.color = "black", order.by = "freq")
      dev.off()
    }
  }
}

for (i in 1:6){
  for (j in 1:4){
    set.test <- list.matrix[[i]][[j]][apply(list.matrix[[i]][[j]][,-1], 1, function(x) !all(x==0)),]
    pdf(file=paste("~/upset/upsets_1/upset_", i, "_", j, ".pdf", sep=""), onefile=FALSE)
    ust <- upset(as.data.frame(set.test), nsets = 12 , main.bar.color = "black", order.by = "freq")
    dev.off()
  }
}






set.test <- list.matrix[[1]][[1]][apply(list.matrix[[1]][[1]][,-1], 1, function(x) !all(x==0)),]
ust <- upset(as.data.frame(set.test), nsets = 12 , main.bar.color = "black", order.by = "freq")
