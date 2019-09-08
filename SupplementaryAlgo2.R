library(compiler)
library(ggplot2)
 
oscillation.waveform.scores  <- function(altern, l.sc, t0, N, n.path){
  score.test = TRUE
  while(score.test){
    l.sc <- length(altern)
    s.q <- rep(0, l.sc-1)
    s.a <- rep(0, l.sc-1)
    flag <- 0
    count <- 1
    
    for(i in seq_len(l.sc-1))
    {
      
      if(altern[i] > altern[i+1])
      {
        s.q[i] <- -1
        s.a[i] <- abs(altern[i] - altern[i+1])
        
      }
      
      if(altern[i] < altern[i+1])
      {
        s.q[i] <- +1
        s.a[i] <- abs(altern[i] - altern[i+1])
      }
      
      if(altern[i] == altern[i+1])
      {
        s.q[i] <- 0
        s.a[i] <- 0
      }
    }
    
    
    buff.a <-  s.a[1]
    buff.q <-  s.q[1]
    
    s.s  <-  rep(NA, l.sc-1)
    s.f  <-  rep(NA, l.sc-1)
    
    l.s.q <- length(s.q)
    
    i=2
    j=1
    while(!is.na(s.q[i]))
    { 
      
      skip = FALSE
      
      l.b.q <- length(buff.q)
      l.b.a <- length(buff.a)
      
      if((s.q[i] < 0 & buff.q[l.b.q] > 0) || (s.q[i] > 0 & buff.q[l.b.q] < 0))
      {
        s.s[j] <- sum(buff.a)
        s.f[j] <- l.b.a
        j = j + 1
        buff.a <- s.a[i]
        buff.q <- s.q[i]
        i = i + 1
        skip = TRUE
      } else 
        if(s.q[i]  <= 0 & buff.q[l.b.q] <= 0 & skip == FALSE)
        {
          if(l.b.q <= t0)
          {
            buff.a <- c(buff.a, s.a[i])
            buff.q <- c(buff.q, s.q[i])
            i = i+1
            skip = TRUE
          } else
          {
            s.s[j] <- sum(buff.a)
            s.f[j] <- l.b.a
            j = j + 1
            buff.a <- s.a[i]
            buff.q <- s.q[i]
            i = i + 1
          }
          
        } else if(s.q[i] >= 0 & buff.q[l.b.q] >= 0 & skip == FALSE)
        {
          if(l.b.q <= t0)
          {
            buff.a <- c(buff.a, s.a[i])
            buff.q <- c(buff.q, s.q[i])
            i = i+1
            
          } else
          {
            
            s.s[j] <-  sum(buff.a)
            s.f[j] <-  l.b.a
            j =j +1
            buff.a <- s.a[i]
            buff.q <- s.q[i]
            i = i+ 1
          }
          
        } else {print("Something goes wrong. ")}
      
      skip = FALSE
      
    } #while
    
    
    if(length(buff.a) > 0 )
    { s.s[j] <-   sum(buff.a) 
    s.f[j] <-   length(buff.a)
    
    }
    
    s.s <- s.s[!is.na(s.s)]
    s.f <- s.f[!is.na(s.f)]
    p.w <- s.s * s.f
    p.l <- sum(s.f)
    l.s <- length(s.s)
    p.s <- 0
    for(h in seq_len(l.s))
    {
      p.s <- p.s +p.w[h]
    }
    s <- p.s/(p.l * (N-1))
    c.v <- l.s/p.l
    # print(j)
    # print(s.s)
    # print(s.f)
    score.test = FALSE
    if( s > 1 && t0 > 1)
    { 
      t0 <- t0 - 1 
      score.test = TRUE
    }
    
    
  }
  
  ret <- data.frame(score=s, m.s=mean(s.s), med.s=median(s.s), sd.s=sd(s.s),
                    m.w=mean(s.f), med.w=median(s.f),  sd.w=sd(s.f), 
                    change.w=l.s,
                    v.change.w= c.v, path.l= l.sc,
                    n.path=n.path,
                    stringsAsFactors=FALSE)
  
  colnames(ret) <- c("osc_s", "mean(osc_s)", "median(osc_s)", "sd(osc_s)", "multi-omic values mean(w)", "median(w)", "sd(w)", "number of changing windows", "osc_k", "signal length", "name of the pathway")

  
  return(ret)
}


get.o.w.s.cmp = cmpfun(oscillation.waveform.scores)

# Ex1
name.pathway  = "pathway number X"
multi.omic.signal = c(1,5,7,1,5,8,1,5,7,1,5,7,1,5,8,1,5,7,1,5,7,1,5,8)
IOAC.bins         = 8
algo1.period.estimation = theta = 3
get.o.w.s.cmp(multi.omic.signal, length(multi.omic.signal), algo1.period.estimation, IOAC.bins, name.pathway)
multi.omic.signal.df <- as.data.frame(multi.omic.signal)
multi.omic.signal.df["gene.order"] <- 1:length(multi.omic.signal)
ggplot(multi.omic.signal.df, aes(x=gene.order, y=multi.omic.signal )) +
  geom_line()+
  geom_point()+
  scale_linetype_manual(values=c("twodash", "dotted"))

# Ex2
name.pathway  = "pathway number Y"
multi.omic.signal = c(1,5,7,3,5,8,2,3,5,1,4,6,1,4,3,2,1,2,1,8,1,6,7)
IOAC.bins         = 8
algo1.period.estimation = theta = 3
get.o.w.s.cmp(multi.omic.signal, length(multi.omic.signal), algo1.period.estimation, IOAC.bins, name.pathway)

multi.omic.signal.df <- as.data.frame(multi.omic.signal)
multi.omic.signal.df["gene.order"] <- 1:length(multi.omic.signal)
ggplot(multi.omic.signal.df, aes(x=gene.order, y=multi.omic.signal )) +
  geom_line()+
  geom_point()+
  scale_linetype_manual(values=c("twodash", "dotted"))

