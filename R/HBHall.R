globalVariables(c("..count..","..density..","ind", "values",
                  "Success Probability", "Expected Data Transmissions",
                  "Expected ACK Transmissions", "Expected Total Transmissions",
                  "Expected Data Receptions", "Expected ACK Receptions",
                  "Expected Total Receptions"))

####################################################################################
#THEORETICAL
####################################################################################

HBH <- function(p1,p2,L,N)
{
  if(p1 >= 1 | p1 <= 0) stop("p1 must be a real number in (0,1)")
  if(p2 >= 1 | p2 <= 0) stop("p2 must be a real number in (0,1)")
  if(L != Inf && (L%%1!=0 | L<0)) stop("L must be a positive integer")
  if(N%%1!=0 | N<1) stop("N must be a positive integer")

  cat("   ","\n")
  cat("###########################################","\n")
  cat("HOP BY HOP - THEORETICAL RESULTS","\n")
  cat("###########################################","\n")
  cat("   ","\n")

  cat(paste("Data success probability         p1 = ",p1),"\n")
  cat(paste("ACK success probability          p2 = ",p2),"\n")
  cat(paste("Maximum number of transmissions  L  = ",L),"\n")
  cat(paste("Number of Hops                   N  = ",N),"\n")
  cat("   ","\n")

  P = function(y,p1,p2,L)
  {
    pp = p1*p2
    return(pp*(1-pp)^(y-1) + ifelse(y==L,(1-pp)^L,0))
  }

  if(L==Inf)
  {
    expect1 = 1/(p1*p2)
    expect2 = p1*expect1
    ET1 = (1+p1)/(p1*p2)

    REC_expect1 = 1/p2
    REC_expect2 = 1
    REC_ET1 = (1+p2)/p2
  }else
  {
    y = seq(1,L)
    expect1 = sum(y*sapply(X = y,FUN = P,p1,p2,L))
    expect2 = p1*expect1
    ET1  = expect1 + expect2

    REC_expect0 = sum(y*sapply(X = y,FUN = P,p1,p2,L))
    REC_expect1 = p1*REC_expect0
    REC_expect2 = p2*REC_expect1
    REC_ET1  = REC_expect1 + REC_expect2
  }

  Pr1    = 1-(1-p1)^L
  PrS    = Pr1^N
  PrSV   = Pr1^seq(1,N)
  w      = Pr1^seq(0,N-1)
  geo    = sum(w)

  ETData = expect1*geo
  ETACK  = expect2*geo
  ETS    = ET1*geo

  ETDataV = expect1*w
  ETACKV  = expect2*w
  ETSV    = ET1*w

  REC_ETData = REC_expect1*geo
  REC_ETACK  = REC_expect2*geo
  REC_ETS    = REC_ET1*geo

  REC_ETDataV = REC_expect1*w
  REC_ETACKV  = REC_expect2*w
  REC_ETSV    = REC_ET1*w

  res    = round(matrix(data = c(c(PrSV,PrS),c(ETDataV,ETData),c(ETACKV,ETACK),c(ETSV,ETS),c(REC_ETDataV,REC_ETData),c(REC_ETACKV,REC_ETACK),c(REC_ETSV,REC_ETS)),nrow = 7,ncol = N+1,byrow = T),3)
  rownames(res)=c("Success Probability","Expected Data Transmissions","Expected ACK Transmissions","Expected Total Transmissions","Expected Data Receptions","Expected ACK Receptions","Expected Total Receptions")
  colnames(res)= c(paste("Hop ",seq(1,N),"/",N,sep = ""),"Total")
  return(res)
}

####################################################################################
#MONTE CARLO
####################################################################################

MCHBH = function(p1,p2,L,N,M=5000)
{
  if(p1 >= 1 | p1 <= 0) stop("p1 must be a real number in (0,1)")
  if(p2 >= 1 | p2 <= 0) stop("p2 must be a real number in (0,1)")
  if(L != Inf && (L%%1!=0 | L<0)) stop("L must be a positive integer")
  if(N%%1!=0 | N<1) stop("N must be a positive integer")
  if(M%%1!=0 | M<1) stop("N must be a positive integer")

  cat("   ","\n")
  cat("###########################################","\n")
  cat("HOP BY HOP - MONTE CARLO SIMULATION RESULTS","\n")
  cat("###########################################","\n")
  cat("   ","\n")

  cat(paste("Data success probability         p1 = ",p1),"\n")
  cat(paste("ACK success probability          p2 = ",p2),"\n")
  cat(paste("Maximum number of transmissions  L  = ",L),"\n")
  cat(paste("Number of Hops                   N  = ",N),"\n")
  cat(paste("Monte Carlo Simulations          M  = ",M),"\n")
  cat("   ","\n")

  prog2 = function(p1,p2,L,N)
  {
    pos=1
    transD = 0
    transA = 0
    RtransD = 0
    RtransA = 0
    paso = TRUE
    TX = matrix(data = 0,nrow = N,ncol = 3)
    RX = matrix(data = 0,nrow = N,ncol = 3)

    while(pos<=N && paso == TRUE)
    {
      trial = 1
      paso = FALSE
      resp = FALSE
      while(trial<=L && resp == FALSE)
      {
        transD = transD +1
        TX[pos,1] = TX[pos,1] + 1
        u1 = runif(1)
        if(u1<p1)
        {
          paso = TRUE
          RtransD = RtransD +1
          RX[pos,1] = RX[pos,1] + 1
          transA = transA +1
          TX[pos,2] = TX[pos,2] + 1
          u2 = runif(1)
          if(u2<p2)
          {
            RtransA = RtransA +1
            RX[pos,2] = RX[pos,2] + 1
            resp=TRUE
          }
        }
        trial=trial+1
      }
      if(paso==TRUE){pos = pos + 1}
    }
    TX[,3] = TX[,1] + TX[,2]
    RX[,3] = RX[,1] + RX[,2]
    return(list(pos=pos,tDATA = transD,tACK = transA,trans = transD+transA,TX = TX,RDATA = RtransD,RACK = RtransA,Rtrans = RtransD+RtransA,RX = RX))
    cat("   ","\n")
  }

  resul = matrix(data = NA,nrow = M,ncol = 7)
  resulA = resulB = array(data = NA,dim = c(M,N,3))
  resulB = array(data = NA,dim = c(M,N,3))
  pb <- txtProgressBar(min = 0, max = M, style = 3)

  for(k in 1:M)
  {
    run = prog2(p1,p2,L,N)
    resul[k,1]=run$pos
    resul[k,2]=run$tDATA
    resul[k,3]=run$tACK
    resul[k,4]=run$trans
    resul[k,5]=run$RDATA
    resul[k,6]=run$RACK
    resul[k,7]=run$Rtrans
    resulA[k,,] = run$TX
    resulB[k,,] = run$RX
    setTxtProgressBar(pb, k)
  }
  close(pb)

  counts = matrix(data = 0,nrow = N+1,ncol = 1)
  counts[sort(unique(resul[,1]))] = c(as.data.frame(table(resul[,1]))$Freq)

  success = rev(cumsum(rev(counts)))[-1]/M
  mean2 = apply(X = resulA,MARGIN = 2:3,FUN = mean)
  mean = apply(X = resul[,2:4],MARGIN = 2,FUN = mean)
  Rmean2 = apply(X = resulB,MARGIN = 2:3,FUN = mean)
  Rmean = apply(X = resul[,5:7],MARGIN = 2,FUN = mean)
  success[is.na(success)] = 1

  table = round(rbind(c(success,success[N]),cbind(t(mean2),mean),cbind(t(Rmean2),Rmean)),3)
  rownames(table)=c("MC Success Probability","MC Mean Data Transmissions","MC Mean ACK Transmissions","MC Mean Total Transmissions","MC Mean Data Receptions","MC Mean ACK Receptions","MC Mean Total Receptions")
  colnames(table)= c(paste("Hop ",seq(1,N),"/",N,sep = ""),"Total")
  return(table)
  cat("   ","\n")
}

# #RUN
#
# #Theoretical
# HBH(p1=0.65,p2=0.4,L=7,N=5)
# #MonteCarlo Simulations
# MCHBH(p1=0.65,p2=0.4,L=7,N=5,M=1000)

HBH0 <- function(p1,p2,L,N)
{
  P = function(y,p1,p2,L)
  {
    pp = p1*p2
    return(pp*(1-pp)^(y-1) + ifelse(y==L,(1-pp)^L,0))
  }

  if(L==Inf)
  {
    expect1 = 1/(p1*p2)
    expect2 = p1*expect1
    ET1 = (1+p1)/(p1*p2)

    REC_expect1 = 1/p2
    REC_expect2 = 1
    REC_ET1 = (1+p2)/p2
  }else
  {
    y = seq(1,L)
    expect1 = sum(y*sapply(X = y,FUN = P,p1,p2,L))
    expect2 = p1*expect1
    ET1  = expect1 + expect2

    REC_expect0 = sum(y*sapply(X = y,FUN = P,p1,p2,L))
    REC_expect1 = p1*REC_expect0
    REC_expect2 = p2*REC_expect1
    REC_ET1  = REC_expect1 + REC_expect2
  }

  Pr1    = 1-(1-p1)^L
  PrS    = Pr1^N
  w      = Pr1^seq(0,N-1)
  geo    = sum(w)

  ETData = expect1*geo
  ETACK  = expect2*geo
  ETS    = ET1*geo

  REC_ETData = REC_expect1*geo
  REC_ETACK  = REC_expect2*geo
  REC_ETS    = REC_ET1*geo

  res    = round(matrix(data = c(PrS,ETData,ETACK,ETS,REC_ETData,REC_ETACK,REC_ETS),nrow = 7,ncol = 1,byrow = T),3)
  #rownames(res)=c("Success Probability","Expected Data Transmissions","Expected ACK Transmissions","Expected Total Transmissions","Expected Data Receptions","Expected ACK Receptions","Expected Total Receptions")
  #colnames(res)= c(paste("Hop ",seq(1,N),"/",N,sep = ""),"Total")
  return(res)
}

#HBH0(p1 = 0.7,p2 = 0.6,L = 7,N = 5)

stochastic_HBH = function(dist1,p11,p12,dist2,p21,p22,L,N,M=10^5,printout=TRUE,plotspdf=TRUE){
  if(L != Inf && (L%%1!=0 | L<0)) stop("L must be a positive integer")
  if(N%%1!=0 | N<1) stop("N must be a positive integer")
  if(dist1 == "uniform"){
    p1 = runif(M,p11,p12)
  }else{
    if(dist1 == "beta"){
      p1 = rbeta(M,p11,p12)
    }else{
      stop("p1 distribution must be either 'uniform' or 'beta'.")
    }
  }
  if(dist2 == "uniform"){
    p2 = runif(M,p21,p22)
  }else{
    if(dist2 == "beta"){
      p2 = rbeta(M,p21,p22)
    }else{
      stop("p1 distribution must be either 'uniform' or 'beta'.")
    }
  }

  out = apply(X = cbind(p1,p2),MARGIN = 1, function(x) HBH0(x[1], x[2], L, N))
  outsum = matrix(apply(out,1,mean),7,1)
  rownames(outsum) = c("Success Probability", "Expected Data Transmissions",
                       "Expected ACK Transmissions", "Expected Total Transmissions",
                       "Expected Data Receptions", "Expected ACK Receptions",
                       "Expected Total Receptions")
  colnames(outsum) = c("Total")

  df = data.frame(p1,p2,t(out))
  colnames(df) = c("p1","p2",rownames(outsum))

  #print(summary(df))
  #install.packages("pastecs")
  #library(pastecs)

  stats = as.data.frame(t(stat.desc(df))[,-c(1:3)])

  p1 = ggplot(df,aes(p1)) +
    geom_histogram(aes(y=..density.., fill = "p1"),alpha = 0.4,color="gray40",breaks = seq(min(p1),max(p1),length.out = round(diff(range(p1))/0.0625)+1)) +
    geom_histogram(aes(x = p2, y = ..density..,fill = "p2"),alpha = 0.4,color="gray40",breaks = seq(min(p2),max(p2),length.out = round(diff(range(p2))/0.0625)+1))+
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(fill = "")+
    xlab("probability") +
    ylab("density") +
    ggtitle("Data and ACK Success Probabilities") +
    xlim(0,1) +
    theme_classic()

  p2 = ggplot(df,aes(x=`Success Probability`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 1,
                   breaks = seq(min(df$`Success Probability`),max(df$`Success Probability`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_vline(xintercept = outsum[1],color="gray40",lwd=1.2)  +
    ggtitle(paste("Success Probability (mean = ",round(outsum[1],3),")",sep = "")) +
    ylab("relative frequency")+ theme_classic()


  p3 = ggplot(df,aes(x=`Expected Data Transmissions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 2,
                   breaks = seq(min(df$`Expected Data Transmissions`),max(df$`Expected Data Transmissions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[2],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected Data Transmissions (mean = ",round(outsum[2],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  p4 = ggplot(df,aes(x=`Expected ACK Transmissions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 3,
                   breaks = seq(min(df$`Expected ACK Transmissions`),max(df$`Expected ACK Transmissions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[3],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected ACK Transmissions (mean = ",round(outsum[3],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  p5 = ggplot(df,aes(x=`Expected Total Transmissions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 4,
                   breaks = seq(min(df$`Expected Total Transmissions`),max(df$`Expected Total Transmissions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[4],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected Total Transmissions (mean = ",round(outsum[4],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  p6 = ggplot(df,aes(x=`Expected Data Receptions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 5,
                   breaks = seq(min(df$`Expected Data Receptions`),max(df$`Expected Data Receptions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[5],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected Data Receptions (mean = ",round(outsum[5],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  p7 = ggplot(df,aes(x=`Expected ACK Receptions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 6,
                   breaks = seq(min(df$`Expected ACK Receptions`),max(df$`Expected ACK Receptions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[6],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected ACK Receptions (mean = ",round(outsum[6],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()


  p8 = ggplot(df,aes(x=`Expected Total Receptions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 7,
                   breaks = seq(min(df$`Expected Total Receptions`),max(df$`Expected Total Receptions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[7],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected Total Receptions (mean = ",round(outsum[7],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  df2 = stack(df[,c(4:9)])

  p9 = ggplot(df2, aes(x = ind, y = values)) +
    geom_boxplot(fill = rev(2:7),alpha = 0.2,color="gray40") +
    coord_flip() +
    #scale_y_continuous(trans='sqrt') +
    xlab("") +
    scale_x_discrete(limits = rev(levels(df2$ind))) +
    theme_classic()

  if(isTRUE(printout)){
    cat(paste("Monte Carlo simulations          M  = ", M), "\n")
    cat(paste("Maximum number of transmissions  L  = ", L), "\n")
    cat(paste("Number of Hops                   N  = ", N), "\n")
    print(stats)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
    print(p8)
    print(p9)
  }

  if(isTRUE(plotspdf)){
    plotsPath = paste("HBH",format(Sys.time(),"%d%m%y_%H%M%S"),".pdf",sep="")
    pdf(file=plotsPath,width = 8.27,height = 5.83)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
    print(p8)
    print(p9)
    dev.off()
    print(paste("Plots file ",plotsPath," saved in working directory ",getwd(),".",sep = ""))
  }
  return(list(data=df,stats = stats))
}

# #We now consider p1 ~ Uniform(0.2,0.6)
# dist1 = "uniform"
# p11 = 0.2
# p12 = 0.6
#
# #and p2 ~ Beta(3,1)
# dist2 = "beta"
# p21 = 3
# p22 = 1
#
# library(pastecs)
# library(ggplot2)
# out = stochastic_HBH(dist1,p11,p12,dist2,p21,p22,L=7,N=5)
