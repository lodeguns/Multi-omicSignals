library(progress)

##################### ########################################################
#COLOMBOS parser compendium
##############################################################################
parseCompendium <- function(destfile){
  require(Rcolombos)
  out_dir <- strsplit(destfile, "\\.")[[1]][1]
  unzip(destfile, exdir=out_dir) # unzip the files in the proper directory 
  files <- dir(path=out_dir, pattern="colombos")
  temp <- paste(out_dir, files[grep("exprdata", files)], sep="/")
  my_cols <- na.omit(scan(temp, nlines=1, sep="\t", what="c", na.strings="", quiet=TRUE))
  exprdata <- read.csv(temp, row.names=1, skip=7, stringsAsFactors=FALSE, sep="\t", header=FALSE)
  exprdata <- exprdata[,c(2:dim(exprdata)[[2]])] 
  colnames(exprdata) = my_cols; exprdata <- exprdata[,c(2:dim(exprdata)[[2]])]
  ## condition annotations 
  if(getOption("REST.version")=="http://rest.colombos.net/"){
    temp <- paste(out_dir, files[grep("refannot", files)], sep="/")
    refannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
    temp <- paste(out_dir, files[grep("testannot", files)], sep="/")
    testannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
    return( list(exprdata=exprdata, refannot=refannot, testannot=testannot) )
    
  } else if(getOption("REST.version")=="http://rest.legacyv2.colombos.net/") {
    temp <- paste(out_dir, files[grep("condannot", files)], sep="/")
    condannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
    return( list(exprdata=exprdata, condannot=condannot) )
  } else return(NULL)
}

parseCompendium_MetaAll <- function(destfile){
  require(Rcolombos)
  out_dir <- strsplit(destfile, "\\.")[[1]][1]
  unzip(destfile, exdir=out_dir) # unzip the files in the proper directory 
  files <- dir(path=out_dir)
  temp <- paste(out_dir, files[grep("exprdata", files)], sep="/")
  my_cols <- na.omit(scan(temp, nlines=1, sep="\t", what="c", na.strings="", quiet=TRUE))
  exprdata <- read.csv(temp, row.names=1, skip=7, stringsAsFactors=FALSE, sep="\t", header=FALSE)
  exprdata <- exprdata[,c(2:dim(exprdata)[[2]])] 
  colnames(exprdata) = my_cols; exprdata <- exprdata[,c(2:dim(exprdata)[[2]])]
  ## condition annotations 
  if(getOption("REST.version")=="http://rest.colombos.net/"){
    temp <- paste(out_dir, files[grep("refannot", files)], sep="/")
    refannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
    temp <- paste(out_dir, files[grep("testannot", files)], sep="/")
    testannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
    meta_allv2 <- list(exprdata=exprdata, refannot=refannot, testannot=testannot) 
  }
  
  return(meta_allv2)
 
}

# require(Rcolombos)
# out_dir  <- strsplit(destfile, "\\.")[[1]][1]
# unzip(destfile, exdir = out_dir) # unzip the files in the proper directory 
# files    <- dir(path=out_dir, pattern="colombos")
# temp     <- paste(out_dir, files[grep("colombos_[a-z]+_exprdata_[0-9]+.txt", files)], sep="/")
# my_cols  <- na.omit(scan(temp, nlines=1, sep="\t", what="c", na.strings="", quiet=TRUE))
# exprdata <- read.csv(temp, row.names=1, skip=7, stringsAsFactors=FALSE, sep="\t", header=FALSE)
# exprdata <- exprdata[,c(2:dim(exprdata)[[2]])] 
# colnames(exprdata) = my_cols; exprdata <- exprdata[,c(2:dim(exprdata)[[2]])]
# ## condition annotations 
# if(getOption("REST.version")=="http://rest.colombos.net/"){
#   temp <- paste(out_dir, files[grep("colombos_[a-z]+_refannot_[0-9]+.txt", files)], sep="/")
#   refannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
#   temp <- paste(out_dir, files[grep("colombos_[a-z]+_testannot_[0-9]+.txt", files)], sep="/")
#   testannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
#   return( list(exprdata=exprdata, refannot=refannot, testannot=testannot) )
#   
# } else if(getOption("REST.version")=="http://rest.legacyv2.colombos.net/") {
#   temp <- paste(out_dir, files[grep("colombos_[a-z]+_condannot_[0-9]+.txt", files)], sep="/")
#   condannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
#   return( list(exprdata=exprdata, condannot=condannot) )
# } else return(NULL)


#######################################################################################
# Pathways multi-organism data extraction KEGG
#######################################################################################
library(KEGGREST)
require(KEGGgraph)
require(igraph)
ext.path <- function(KEGG.id)
{
  p <- keggLink(KEGG.id, "pathway")
  return(unique(names(p)))
}


ext.path.kegg <- function(KEGG.id.path){
   ll <- list()
   for(i in 1:length(KEGG.id.path))
  { print(paste("i", KEGG.id.path[[i]], sep="-"))
     ll[length(ll)+1] <- list(ext.path(KEGG.id.path[[i]]))
  
  }

  list.pathways <- ll
  names(list.pathways) <- KEGG.id.path
  save(list.pathways , file="list.pathways.RData")

  list.common     <- list()
  list.KEGG.id    <- list()
  for(j in 1:length(list.pathways)){
     list.common[length(list.common) +1]     <- list(gsub("([:a-z:])", "", list.pathways[[j]]))
     list.KEGG.id[length(list.KEGG.id) + 1]  <- unique(list(gsub("path", "", gsub("([:0-9:])", "", list.pathways[[j]])))[[1]])
   }

return(list(list.common = list.common, list.KEGG.id = list.KEGG.id))}

kegg.get.path <- function(pathid, path.org="eco"){
  g<-NULL
  #print(pathid)
  gKGML <- getKGMLurl(pathid,  organism=path.org )
  tryCatch( g <- parseKGML2Graph(gKGML, genesOnly = TRUE, expandGenes = TRUE),
            error = function(c) {
              c$message <- paste(c$message, " pathway:",gKGML, ")", sep=" ")
              print(c)
              return(NULL)
            })
  
  return(g)
}

get.igraph.kegg   <- function(pathid.l,  path.org="eco")
{ 
  graph1<- NULL
  if(!is_igraph(pathid.l)){
    #print(paste(pathid.l))
    tryCatch( graph1<-
      igraph.from.graphNEL(kegg.get.path(pathid.l, path.org)
                           , name = TRUE, weight = FALSE, unlist.attrs = TRUE),
      error = function(c) {
        c$message <- paste(c$message, " pathway:",pathid.l, ")", sep=" ")
        print(c)
        return(NULL)
      })
}
  
  if(is_igraph(graph1))
    return(graph1)
  return(NULL)
}

get.list.paths.org <- function(i, list.paths.kegg)
{ require(progress)
  list.paths00 <- list()
  pb <- progress_bar$new(
    format = "Pathways downloading [:bar] :percent in :elapsed",
    total = length(list.paths.kegg$list.common[[i]]), clear = FALSE, width= 60)
  
  for( j in 1:length(list.paths.kegg$list.common[[i]]))
  {
    list.paths00[length(list.paths00) +1] <- list(get.igraph.kegg(list.paths.kegg$list.common[[i]][[j]],  list.paths.kegg$list.KEGG.id[[i]]))
    
    pb$tick()
    
  }
  
  names(list.paths00) <- list.paths.kegg$list.common[[i]]
  return(list.paths00)
}


get.kegg.ncbi.id <- function(ll1)
{
  require(KEGGREST)
  kegg.ncbi.list <- list()
  
  for(i in 1:length(ll1))
  { ll<- NULL
    flag = FALSE
    tryCatch( ll <- keggConv(ll1[[i]], "ncbi-geneid"),
              error = function(c) {
                c$message <- paste(c$message, " :",ll1[[i]], ")", sep=" ")
                print(c)
                flag = TRUE
              })
    if(flag)
    { kegg.ncbi.list[length(kegg.ncbi.list)+1] <- list() }
    else{
    kegg.ncbi.list[length(kegg.ncbi.list)+1] <- list(ll) }
    
  }
  
  names(kegg.ncbi.list) <- ll1
  
  return(kegg.ncbi.list) 
}


########################################################################
## Estrai tutti i dati relativi 
##

get.nt.seq <- function(str){
  tryCatch({lstr1 <- keggGet(str, "ntseq")[[1]]
  return(paste(lstr1))},
  warning = function(w) {print(paste("nt invalid parameter:: ", str)); return(0);},
  error = function(e) {print(paste("nt invalid identifier:: ", str)); return(0);} ) 
}


get.cai.trn <- function(osc.ribs.eco.static)
{
  gentest <- data.frame(locus<-character(0), cds<-character(0), desc<-character(0))
  pb <- progress_bar$new(
    format = "CAI training set [:bar] :percent in :elapsed",
    total = length(osc.ribs.eco.static), clear = FALSE, width= 60)
  count<-0
  for(i in 1:length(osc.ribs.eco.static))
  {
    row<-osc.ribs.eco.static[i]
    kggett<-get.nt.seq(row)
    if(kggett!=0){
      gentest<-rbind(gentest, data.frame(locus = as.character(row),
                                         cds = as.character(tolower(kggett[[1]])),
                                         desc = as.character(paste(count))
      ))
      count<-count+1}
    pb$tick()
  }
  
  return(gentest)
  
}






get.gene.info <- function(str, ids.org1, cai.tr.set, emw){
  require(XML)
  require(reutils)
  require(rentrez)
  require(CAIpackage)

  
  fetch  <- efetch(ids.org1, db = "gene")
  xmlp   <- xmlParse(fetch$content)
  a      <- getNodeSet(xmlp, "//Entrezgene_gene/Gene-ref")
  to.ret <- xpathApply(a[[1]], "//Gene-ref_locus", xmlValue)
  if(length(to.ret) == 0)
  {
    to.ret<-xpathApply(a[[1]], "//Gene-ref_locus-tag", xmlValue)
  }
  
  from    <- as.numeric(xmlValue(getNodeSet(xmlp, "//Seq-interval_from/text()[1]")[[1]]))
  to      <- as.numeric(xmlValue(getNodeSet(xmlp, "//Seq-interval_to/text()[1]")[[1]]))
  dir1    <- getNodeSet(xmlp, "//Na-strand") 
  to.ret2 <- xpathSApply(dir1[[1]], "//Na-strand/@value")[[1]]
  
  dir2    <- getNodeSet(xmlp, "//Entrezgene_type") 
  to.ret3 <- xpathSApply(dir1[[1]], "//Entrezgene_type/@value")[[1]]
  
  to.ret4 <- NULL
  tryCatch({ a1      <- getNodeSet(xmlp, "//Prot-ref_name_E")[[1]]
  to.ret4 <- xpathApply(a1[[1]], "//Prot-ref_name_E", xmlValue)[[1]]
  },
  warning = function(w) {  to.ret4 <- to.ret3 },
  error = function(e) {    to.ret4 <- to.ret3  }) 
  
  lstr1<- NULL
  tryCatch({lstr1 <- keggGet(str, "aaseq")[[1]]
  },
  warning = function(w) { lstr1 <- NULL},
  error = function(e)   { lstr1 <- NULL} ) 
  
  if(!is.null(lstr1)){
  amw <- nchar(paste(lstr1))*110
  } else {amw=0}
  
  tryCatch({lstr2 <- keggGet(str, "ntseq")[[1]]
  },
  warning = function(w) {print(paste("nt invalid parameter:: ", str)); lstr2 <- NULL},
  error = function(e) {print(paste("nt invalid identifier:: ", str)); lstr2 <- NULL} ) 
  
  if(!is.null(lstr2)){
  train.set <- cai.tr.set
  gene.test <- data.frame(locus<-character(0), cds<-character(0), desc<-character(0))
  count<-0
  gene.test <- rbind(gene.test, data.frame( locus = as.character(paste(to.ret)),
                                            cds = as.character(tolower(lstr2)),
                                           desc = as.character(paste(ids.org1))))
                                     
                                     osc.CAI <-CAIObject()
                                     RefSet(osc.CAI)<- train.set
                                     CAISet(osc.CAI)<- gene.test
                                     gene_assoc_cai <- CAIcalc(osc.CAI)
                                     cai = gene_assoc_cai@caitable[['CAI']]
  } else {cai <- NULL}
  
  if(!is.null(lstr1)){
  calc.mw  <- calc_molecular_weight(paste(lstr1), emw)
  } else
  {calc.mw  <- NULL}
  
  return(list(  gene_sym = paste(to.ret),  from   = from, to = to, 
                strand = to.ret2, aa.seq = lstr1, nt.seq = lstr2, 
                cai = cai, mw = calc.mw,
                amw = amw, gene_type = to.ret3, gene_dsc = to.ret4))
}






get.aa.seq <- function(str, org){
  tryCatch({lstr1 <- keggGet(paste(org, ":",str, sep=""), "aaseq")[[1]]
  return(paste(lstr1))},
  warning = function(w) {print(paste("aa invalid parameter:: ", str)); return(0);},
  error = function(e) {print(paste("aa invalid identifier:: ", str)); return(0);} ) 
  
}



#exact molecular weight
emw <- read.table("/home/lodeguns/Dropbox/Lavoro FRANPA 2017/aa_molecular_weight", sep="\t", header = T)
calc_molecular_weight <- function(str, emw)
{ sum_w <- 0
vect_aa_names <- as.vector(emw[,3])
vect_aa_mw    <- as.numeric(emw[,4])
for(i in 1:length(vect_aa_mw))
{
  found_aa <- sapply(regmatches(str, gregexpr(vect_aa_names[i], str)), length)
  
  sum_w <- sum_w + (found_aa * vect_aa_mw[i])
  
}

return(sum_w)
}









## test the integrity of these KEGG nodes against gene.col of EcoCyc
# 
# require(RCurl)
# require(KEGGREST)
# require(KEGGgraph)
# require(XML)
# require(reutils)
# require(rentrez)
# diff.genes <-  setdiff( v.kegg.conn, v.ecocy )
# diff.genes.service <- c()
# for(j in 1:length(diff.genes)){
#   cc<-keggConv("ncbi-geneid", diff.genes[[j]])
#   print(cc)
#   if(length(cc)!= 0){
#     st1 <- strsplit(cc, ":")[[1]][2]
#     
#     fetch<-efetch(st1, db = "gene")
#     
#     xmlp<-xmlParse(fetch$content)
#     
#     a = getNodeSet(xmlp, "//Entrezgene_gene/Gene-ref")
#     
#     xp  <- xpathApply(a[[1]], "//Gene-ref_locus", xmlValue)
#     xp1 <- xpathApply(a[[1]], "//Gene-ref_syn/*", xmlValue)
#   } else
#   {xp =""}
#   
#   gg1 <- grep(xp, genes.col$GENE.NAME)
#   gg2 <- grep(xp, genes.col$y.name)
#   
#   for(kk0 in 1:length(xp))
#   {
#     for(dfi in 1:ncol(genes.col))
#     {   vc1 <- as.vector(genes.col[,dfi])
#     gg1<-grep(xp[[kk0]], vc1)
#     if(length(gg1)!=0)
#     { break;}
#     }
#   }
#   if(length(gg1)==0)
#     for(kk1 in 1:length(xp1))
#     {
#       for(dfi in 1:ncol(genes.col))
#       {   vc1 <- as.vector(genes.col[,dfi])
#       gg1<-grep(xp1[[kk1]], vc1)
#       if(length(gg1)!=0)
#       { break;}
#       }
#     }
#   
#   if(length(gg1)!=0)
#   { if(genes.col[gg1[1],]$b == gsub("eco:", "", diff.genes[[j]]))
#   {
#     diff.genes.service <- c(diff.genes.service, 0)
#   } 
#     else
#     {diff.genes.service <- c(diff.genes.service, 1)}
#   } else
#   {
#     diff.genes.service <- c(diff.genes.service, 1)
#   }
#   
# }
# 
# #possible incongruency nodes   
# incongruency.nodes <- diff.genes[diff.genes.service==1]   
# # Ok, there isn't incongruency, are nodes not accounted.
# # check also on http://genexpdb.ou.edu/databases/genexpdb/ 
# 
# 
# 
# fname <- "./list.kegg.path.igraph.Rdata"
# list.kegg.path.igraph<-m.l
# save(list.kegg.path.igraph, file=fname)
# 
# fname <- "./kegg.path.merged.igraph.Rdata"
# kegg.path.merged.igraph<-m.u.connected
# save(kegg.path.merged.igraph, file=fname)








