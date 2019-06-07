##
##    Median periodicity search window estimator
##    = theta's estimator
##


library(compiler)

get.window.median.prev <- function(altern, l.sc, N){
  window <- rep(0, l.sc)
  
  for(k in seq_len(l.sc-1))
  {  l <- k
  c.w          <- 1
  alter        <- abs(altern[k] - altern[k+1])
  median.prev  <-  alter
  if(k+1 < l.sc){
    while(alter <= median(median.prev))
    {    j = k+1
    if(j+1 < l.sc)
    {
      alter <- abs(alter   - altern[j+1])
      median.prev <- c(median.prev, alter)
      k = j
    } else { break;}
    c.w <- c.w +1
    }

    median.prev <- c(0)
    window[l] <- c.w+1
  }
  
  
  }
  return(window)
}

get.window = cmpfun(get.window.median.prev) 

multi.omic.signal = c(1,5,8,1,5,8,1,5,7,1,4,6,1,4,2,2,1,2,1,8,1,6,7)
IOAC.bins         = 9
mo.windows <- get.window(multi.omic.signal, length(multi.omic.signal), IOAC.bins)
#mean(mo.windows)
#sd(mo.windows)
algo1.period.estimation = theta = max(mo.windows)   #output 

