#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

start_time <- Sys.time()

## import arguments from command line
compl <- as.character(args[1])
thresh <- as.numeric(args[2])
bootNum <- as.numeric(args[3])
montecarlo <- as.numeric(args[4])
strat <- as.numeric(args[5])
array_num <- as.numeric(args[6])
if(is.na(as.numeric(args[7]))){
  rand<-NULL
} else{
  rand<-as.numeric(args[7])
}
if(is.na(as.numeric(args[8]))){
  expo<-NULL
} else{
  expo<-as.numeric(args[8])
}
if(is.na(as.numeric(args[9]))){
  cens<-NULL
} else{
  cens<-as.numeric(args[9])
}
if(is.na(as.numeric(args[10]))){
  int<-NULL
} else{
  int<-as.numeric(args[10])
}

cat("EAGeR Compliance Variable \n")
compl
cat("Threshold Defining Compliance \n")
thresh
cat("Number of Bootstraps \n")
bootNum
cat("MC Sample \n")
montecarlo
cat("Array Number \n")
array_num
cat("Stratify by Eligibility? (1 = true) \n")
strat
cat("Randomization Assignment \n")
rand
cat("Compliance Assignment \n")
expo
cat("Censoring \n")
cens
cat("Interaction \n")
int

# load libraries/functions
setwd(".")
packages <- c("data.table","tidyverse","haven","lubridate","here","VGAM","splines","parallel","cmprsk")
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

cat("Check if library directory exists",'\n')
ifelse(!dir.exists(userLib), dir.create(userLib), FALSE)

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

dr_here()

a <- read_csv(here("data","2020_03_06-eager_weekly.csv"))
a

cat("The Data", '\n')
print(head(data.frame(a)))

# set seed (0 does no bootstrap resampling)
seed <- array_num
# function to bootstrap the g Computation algorithm
g_boot1<-function(seed){
  set.seed(seed)
  cat("Now Running SEED",seed,'\n')
  cat('\n')
  cat("Resampling Data",'\n')
  clusters <- as.numeric(names(table(a$id)))
  index <- sample(1:length(clusters), length(clusters), replace=TRUE)
  bb <- table(clusters[index])
  boot <- NULL
  if(seed==0){
    boot<-a
  } else{
    for(zzz in 0:max(bb)){ # this counter goes from zero to select empirical data (no resample)
      cc <- a[a$id %in% names(bb[bb %in% c(zzz:max(bb))]),]
      cc$bid<-paste0(cc$id,zzz)
      boot <- rbind(boot, cc)
    }}

  if(seed==0){
    print(head(boot))
  }

  boot$last_id <- as.numeric(!duplicated(boot$id,fromLast = T))
  head(data.frame(boot))
  sum(boot$last_id)

  print(table(boot[boot$last_id==1,]$R))
  print(table(boot[boot$last_id==1,]$E))

  table(boot$X)

  boot$X  <- as.factor(1*as.numeric(boot$X==0&boot$R==0) +
                         1*as.numeric(boot$X==0&boot$R==1) +
                         2*as.numeric(boot$X==1&boot$R==0) +
                         3*as.numeric(boot$X==1&boot$R==1))


  head(data.frame(boot))

  boot$Xl  <- as.factor(1*as.numeric(boot$Xl==0&boot$R==0) +
                          1*as.numeric(boot$Xl==0&boot$R==1) +
                          2*as.numeric(boot$Xl==1&boot$R==0) +
                          3*as.numeric(boot$Xl==1&boot$R==1))

  boot$Xl1 <- as.factor(1*as.numeric(boot$Xl1==0&boot$R==0) +
                          1*as.numeric(boot$Xl1==0&boot$R==1) +
                          2*as.numeric(boot$Xl1==1&boot$R==0) +
                          3*as.numeric(boot$Xl1==1&boot$R==1))

  print(table(boot$X))
  print(table(boot$Xl))
  print(table(boot$Xl1))

  cat("Are resampled data identical to original for seed=", seed,"?:",identical(boot,a2),'\n')

  cat("Centering and scaling the time-scale variable jj", '\n')
  boot$jj<-scale(boot$j)
  boot$jjsq<-boot$jj*boot$jj
  j_mean<-attributes(boot$jj)$`scaled:center`
  j_sd<-attributes(boot$jj)$`scaled:scale`

  print(head(data.frame(boot)))

  cat('\n')
  cat("Fitting parametric models",'\n')
  mY_0<-function(k){
    fitY<-glm(Y~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+Xl1+B+Bl+N+Nl
              +jj,
              data=boot,subset=Z==1&R==0,family=binomial)
    return(fitY)
  }
  mY_1<-function(k){
    fitY<-glm(Y~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+Xl1+B+Bl+N+Nl
              +ns(jj,df=3),
              data=boot,subset=Z==1&R==1,family=binomial)
    return(fitY)
  }
  mD_0<-function(k){
    fitD<-glm(D~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+B+Bl+N+Nl
              +bs(jj,df=3),
              data=boot,subset=Z==1&R==0,family=binomial(link="logit"))
    return(fitD)
  }
  mD_1<-function(k){
    fitD<-glm(D~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+B+Bl+N+Nl
              +bs(jj,df=3),
              data=boot,subset=Z==1&R==1,family=binomial(link="logit"))
    return(fitD)
  }
  mC_0<-function(k){
    fitC<-glm(C~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+B+Bl+N+Nl
              +bs(jj,df=3),
              data=boot,subset=R==0,family=binomial(link="logit"))
    return(fitC)
  }
  mC_1<-function(k){
    fitC<-glm(C~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+B+Bl+N+Nl
              +bs(jj,df=3),
              data=boot,subset=R==1,family=binomial(link="logit"))
    return(fitC)
  }
  mS_0<-function(k){
    fitS<-glm(S~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+B+Bl+N+Nl
              +as.numeric(jj)+as.numeric(jjsq),
              data=boot,subset=Z==0&R==0,family=binomial(link="logit"))
    return(fitS)
  }
  mS_1<-function(k){
    fitS<-glm(S~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+B+Bl+N+Nl
              +as.numeric(jj),
              data=boot,subset=Z==0&R==1,family=binomial(link="logit"))
    return(fitS)
  }
  mZ_0<-function(k){
    fitZ<-glm(Z~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+Xl1+B+Bl+N+Nl
              +as.numeric(jj)+as.numeric(jjsq),
              family=binomial(link="logit"),data=boot,subset=Zl==0&R==0)
    return(fitZ)
  }
  mZ_1<-function(k){
    fitZ<-glm(Z~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +X+Xl+Xl1+B+Bl+N+Nl
              +as.numeric(jj)+as.numeric(jjsq),
              family=binomial(link="logit"),data=boot,subset=Zl==0&R==1)
    return(fitZ)
  }

  mX_0<-function(k){
    fitX<-vglm(X~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
               +B+Bl+N+Nl
               +as.numeric(jj)+as.numeric(jjsq),
               family=multinomial(refLevel = 1),data=boot)
    return(fitX)
  }
  mX_1<-function(k){ # same as mX_0
    fitX<-vglm(X~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
               +B+Bl+N+Nl
               +as.numeric(jj)+as.numeric(jjsq),
               family=multinomial(refLevel = 1),data=boot)
    return(fitX)
  }

  mB_0<-function(k){
    fitB<-glm(B~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +Xl+Xl1+Bl+N+Nl
              +as.numeric(jj)+as.numeric(jjsq),
              family=binomial,data=subset(boot,R==0))
    return(fitB)
  }
  mB_1<-function(k){
    fitB<-glm(B~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +Xl+Xl1+Bl+N+Nl
              +as.numeric(jj)+as.numeric(jjsq),
              family=binomial,data=subset(boot,R==1))
    return(fitB)
  }
  mN_0<-function(k){
    fitN<-glm(N~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +Xl+Xl1+Bl+Nl
              +as.numeric(jj)+as.numeric(jjsq),
              family=binomial,data=subset(boot,R==0))
    return(fitN)
  }
  mN_1<-function(k){
    fitN<-glm(N~V1+V2+V3+V4+V5+V6+bs(V7,df=3)+bs(V8,df=3)
              +Xl+Xl1+Bl+Nl
              +as.numeric(jj)+as.numeric(jjsq),
              family=binomial,data=subset(boot,R==1))
    return(fitN)
  }
  mL<-c(mN_0,mB_0,mX_0,mZ_0,mC_0,mS_0,mD_0,mY_0,
        mN_1,mB_1,mX_1,mZ_1,mC_1,mS_1,mD_1,mY_1)

  fitR<-lapply(1:16,function(x) mL[[x]](k))

  m1<-list(fitR[[1]],fitR[[2]],fitR[[3]],fitR[[4]],fitR[[5]],fitR[[6]],fitR[[7]],fitR[[8]])
  m2<-list(fitR[[9]],fitR[[10]],fitR[[11]],fitR[[12]],fitR[[13]],fitR[[14]],fitR[[15]],fitR[[16]])
  fitR<-list(m1,m2)

  cat('\n')
  cat("Creating monte carlo dataset",'\n')
  ## create monte carlo dataset
  # select first obs for each person to obtain joint empirical distribution of baseline covariates
  MC0<-boot[boot$j==1,]

  cat('\n')
  cat("Sampling (with replacement) Monte Carlo Sample from each treatment arm",'\n')
  # sample with replacement from each treatment arm
  spl <- split(MC0, list(MC0$R))
  samples <- lapply(spl, function(x) x[sample(1:nrow(x), montecarlo/length(spl), replace=T),])
  MC <- rbindlist(samples)
  MC$id<-1:montecarlo

  if(seed==0){
    cat('\n')
    cat("The Monte Carlo Dataset",'\n')
    print(head(MC))

    cat('\n')
    cat("The Monte Carlo Dataset: Randomization Table",'\n')
    print(table(MC$R))
  }

  cat('\n')
  cat("Simulating follow-up using parametric fits and baseline resampled data",'\n')
  # PREDICT FOLLOW-UP BASED ON G FORMULA USING PGF FUNCTION
  pgf<-function(ii, mc_data, lngth, randomization = NULL, exposure = NULL, censoring=NULL, interaction=NULL, stratum=0){
    pFunc<-function(mod,ndat){as.numeric(predict(mod,newdata=ndat,type="response")>runif(1))}
    d <- mc_data
    d <- d[d$id==ii,]
    cat("Observation",ii,"from Monte Carlo Data",'\n')
    print(d)
    lngth <- lngth
    Vp <- Rp <- Ep <- Bp <- Np <- Zp <- Xp <- Cp <- Sp <- Dp <- Yp <- Y2p <- mm <- cJ <- GA <- numeric()
    Xp <- factor(levels=c("1","2","3"))
    mm[1] <- j <- 1
    id <- d$id
    Vp <- d[,paste("V",1:8,sep="")]
    cJ<-99

    if(stratum==0){
      Ep <- d$E
      if (is.null(randomization)) {
        Rp <- d$R
      } else {
        Rp <- randomization
      }
    } else {
      # for stratum == 1 R and E are switched (this is good)
      Rp <- d$R

      if (is.null(randomization)) {
        Ep <- d$E
      } else {
        Ep <- randomization
      }
    }

    Bp[1] <- d$B
    Np[1] <- d$N
    Zp[1] <- d$Z

    if(Zp[1]==1){cJ<-1}

    Xp[1] <- d$X

    # withdrawal
    if(is.null(censoring)){
      Cp[1] <- d$C
    } else{
      Cp[1] <- censoring
    }

    # efuwp
    if(Zp[1]==0&Cp[1]==0){
      Sp[1] <- d$S
    } else{
      Sp[1] <- 0
    }

    Dp[1] <- Yp[1] <- 0

    for (j in 2:lngth) {
      cat("Iteration",j,"for observation",ii,"from Monte Carlo Data",'\n')
      if (Yp[j - 1]==0&Dp[j - 1]==0&Cp[j - 1]==0&Sp[j - 1]==0) {
        cat("Conditions",'\n')
        if (j==2){
          Bl1<-Nl1<-Zl1<-0
          Xl1<-factor(3,levels=c("1","2","3"))
        } else{
          Xl1 <- Xp[j - 2]
          Bl1 <- Bp[j - 2]
          Nl1 <- Np[j - 2]
          Zl1 <- Zp[j - 2]
        }
        Xl=Xp[j-1];Zl=Zp[j-1];Bl=Bp[j-1];Nl=Np[j-1]

        cat("Generating Nausea",'\n')
        dNp <- data.table(Vp,E=Ep,Xl,Xl1,Zl,Zl1,Bl,Bl1,Nl,Nl1,R=Rp,jj=as.numeric((j-j_mean)/j_sd),jjsq=(as.numeric((j-j_mean)/j_sd))^2)
        Np[j] <- pFunc(fitR[[Rp+1]][[1]], dNp)

        cat("Generating Bleeding",'\n')
        dBp <- data.table(Vp,E=Ep,Xl,Xl1,Zl,Zl1,Bl,Bl1,N=Np[j],Nl,Nl1,R=Rp,jj=as.numeric((j-j_mean)/j_sd),jjsq=(as.numeric((j-j_mean)/j_sd))^2)
        Bp[j] <- pFunc(fitR[[Rp+1]][[2]], dBp)

        dXp <- data.table(Vp,E=Ep,Xl,Xl1,Zl,Zl1,B=Bp[j],Bl,Bl1,N=Np[j],Nl,Nl1,R=Rp,jj=as.numeric((j-j_mean)/j_sd),jjsq=(as.numeric((j-j_mean)/j_sd))^2)
        cat("Generating Compliance",'\n')
        full_pred <- predict(fitR[[Rp+1]][[3]], newdata=dXp, type="response")
        if(Rp==1){
          tmp <- full_pred[,1]/sum(full_pred[,c(1,3)])
          Xp[j] <- if_else(as.numeric(runif(1)>tmp)*3==3,3,1)
        } else{
          tmp <- full_pred[,1]/sum(full_pred[,c(1,2)])
          Xp[j] <- if_else(as.numeric(runif(1)>tmp)*2==2,2,1)
        }

        cat("Generating Conception",'\n')
        dZp <- data.table(Vp,E=Ep,X=Xp[j],Xl,Xl1,B=Bp[j],Bl,Bl1,N=Np[j],Nl,Nl1,R=Rp,jj=as.numeric((j-j_mean)/j_sd),jjsq=(as.numeric((j-j_mean)/j_sd))^2)
        if (Zp[j - 1] == 0){
          Zp[j] <- pFunc(fitR[[Rp+1]][[4]], dZp)
        } else {
          Zp[j] <- 1
        }

        if(Zp[j-1] == 0 & Zp[j] == 1){
          cJ <- j
        }

        # withdrawal
        cat("Generating Withdrawal",'\n')
        dCp <- data.table(Vp,E=Ep,X=Xp[j],Xl=Xp[j-1],Xl1,B=Bp[j],Bl=Bp[j-1],Bl1,N=Np[j],Nl=Np[j-1],Nl1,R=Rp,jj=as.numeric((j-j_mean)/j_sd),jjsq=(as.numeric((j-j_mean)/j_sd))^2)
        if(is.null(censoring)){
          Cp[j] <- pFunc(fitR[[Rp+1]][[5]], dCp)
        } else{
          Cp[j] <- 0
        }

        # efuwp
        cat("Generating EFUWP",'\n')
        dSp <- data.table(Vp,E=Ep,X=Xp[j],Xl=Xp[j-1],Xl1,B=Bp[j],Bl=Bp[j-1],Bl1,N=Np[j],Nl=Np[j-1],Nl1,R=Rp,jj=as.numeric((j-j_mean)/j_sd),jjsq=(as.numeric((j-j_mean)/j_sd))^2)
        if(Zp[j]==0&Cp[j]==0&j<50){
          Sp[j] <- pFunc(fitR[[Rp+1]][[6]], dSp)
        } else{
          Sp[j] <- 0
        }

        # pregnancy loss and LB
        dYp <- data.table(Vp,E=Ep,X=Xp[j],Xl=Xp[j-1],Xl1,B=Bp[j],Bl=Bp[j-1],Bl1,N=Np[j],Nl=Np[j-1],Nl1,R=Rp,j,jj=(j-j_mean)/j_sd,jjsq=(as.numeric((j-j_mean)/j_sd))^2)
        Dp[j] <- 0
        if(Zp[j]==1&Cp[j]==0){
          Dp[j] <- pFunc(fitR[[Rp+1]][[7]],dYp)
        }

        Yp[j] <- 0
        if(Zp[j]==1&Cp[j]==0&Dp[j]==0){
          if(j>=40&j-cJ>20){ #
            Yp[j] <- pFunc(fitR[[Rp+1]][[8]],dYp)
          }
        }

        if(Yp[j]==1){
          GA <- j - cJ
        } else{
          GA <- 99
        }
      } else {
        break
      }
      mm[j] <- j
    }
    boot_num<-seed
    gdat <- data.table(boot_num,id,mm,Vp,Rp,Ep,Np,Bp,Zp,Xp,Cp,Sp,Yp,Dp)
    gdat$last<-as.numeric(!(gdat$Cp==0)|!(gdat$Sp==0)|!(gdat$Yp==0)|!(gdat$Dp==0)|gdat$mm==lngth)
    return(gdat)
  } # end pgf function
  r<-mclapply(1:montecarlo, function(i) pgf(i,MC,60,randomization = rand,exposure = expo,censoring=cens,interaction=int,stratum=strat),mc.cores=cores)
  results<-do.call(rbind,r)
  cat('\n')
  cat("End of run for SEED",seed,'\n')
  return(results)
} # end g_boot1 function

for(ii in 0:bootNum){ #bootNum

  cat('\n')
  cat("Now running bootstrap",ii,'\n')

  results<-g_boot1(ii)

  cat('\n')
  cat("Save simulated follow-up data",'\n')
  if(is.null(rand)&is.null(expo)){
    filename<-"gS-NaturalCourse.txt"
  } else if(is.null(int)&expo==3) {
    filename<-"gS-compliant_placebo-noInt.txt"
  } else if(is.null(int)&expo==4) {
    filename<-"gS-compliant_treatment-noInt.txt"

  }else if(is.null(int)&rand==1&expo==1) {
    filename<-"gS-exposed-noInt.txt"
  } else if(int==1){
    filename<-"gS-PreConceptionExposure.txt"
  } else if(int==2){
    filename<-"gS-PostConceptionExposure.txt"
  } else if(int==6){
    filename<-"gS-PostConception6Exposure.txt"
  } else if(int==8){
    filename<-"gS-PostConception8Exposure.txt"
  } else if(int==12){
    filename<-"gS-PostConception12Exposure.txt"
  } else if(int==20){
    filename<-"gS-PostConception20Exposure.txt"
  }

  print(filename)

  if(bootNum==0){
    write.table(results,filename,sep="\t",row.names=F)
  } else{
    write.table(results,filename,sep="\t",row.names=F,col.names=F,append=T)
  }
}

elapsed_time<-Sys.time()-start_time
print(elapsed_time)
