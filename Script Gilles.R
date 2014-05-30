
setwd("C:/WORKS/Logiciels/msms/bin")
library(vegan)
source("ms.r")
source("FunctionNielsen.r")

#### model Parameters        
ndraws <- 1   # Number of replicates
theta_com <- 10 # theta=2 J L v      where J is the size of the pooled localities, L is the number of localities and v is the ??????
SamplePop <- c(150,150,150,150,150,150)   ### Gives the number of sampled individuals in each localities
nsam <- sum(SamplePop)  ## total number of sampled individuals
L <- length(SamplePop)  ## number of localities sampled
I <- 5 # I=2 J m     where m is the migration rate
Ratio <- 100   #### Ratio mutation rate/speciation rate

### msms paramaters
theta <- (1+Ratio)*theta_com/L
Npop <- L
LimitSp <-1-1/Ratio
Mig<-I

###### Coalescence simulation using 'ms'
if (Npop>1) SIM <- paste("msms",nsam,ndraws,"-t",theta,"-I",Npop,apply(as.matrix(SamplePop),2,paste,collapse=" "), Mig,   " > out2.txt")
if (Npop==1)SIM <-  paste("msms",nsam,ndraws,"-t",theta," > out2.txt")

shell(SIM)

OUTPUT <- ms.inp.multi(nsam=nsam,ndraws=ndraws,ms.output.file="out2.txt")

OUTPUT_SP <- list()
OUTPUT_G <- list()

j=1
# for (j in 1:ndraws) {
  
  ###### ##############################
  ##### Partitionning sample data into species and haplotypes according to the threshold given
  MatSp <- OUTPUT[[j]]$sample[,which(OUTPUT[[j]]$mutations>LimitSp)]
  MatG <- OUTPUT[[j]]$sample[,which(OUTPUT[[j]]$mutations<=LimitSp)]
  
  Output_Sp=list()
  Output_G=list()
  Output_Sp[[1]] <- length(which(OUTPUT[[j]]$mutations>LimitSp)); Output_G[[1]] <- length(which(OUTPUT[[j]]$mutations<=LimitSp)) ### Number of mutations/speciation events
  Output_Sp[[2]] <- OUTPUT[[j]]$mutations[which(OUTPUT[[j]]$mutations>LimitSp)]; Output_G[[2]] <- OUTPUT[[j]]$mutations[which(OUTPUT[[j]]$mutations<=LimitSp)]  ### Position of speciation/mutation events in the array
  Output_Sp[[3]] <- MatSp; Output_G[[3]] <- MatG
  names(Output_Sp) <- names(OUTPUT[[j]]) ; names(Output_G) <- names(OUTPUT[[j]])
  
  OUTPUT_SP[[j]] <- Output_Sp; OUTPUT_G[[j]] <- Output_G   ####### Datasets arranged into format expected by function makeStats 

  ##############################
  ##### Species level stats
  ##############################
  SpeciesRich <- matrix(NA,1,Npop+1)
  NielsenSpDiv <- matrix(NA,1,Npop+1)
  
  Stats_Sp=list()
  ListeSp <- matrix(0,length(table(apply(MatSp,1,paste,collapse=" "))),Npop+2)  #### ListeSp gives the list of species with abundance and 'sequence' composition
  ListeSp[,1] <- table(apply(MatSp,1,paste,collapse=" "))
  ListeSp[,2] <- names(table(apply(MatSp,1,paste,collapse=" ")))
  NAMES <- NULL; for (i in 1:nrow(ListeSp)) NAMES <- c(NAMES,paste("Sp",i));rownames(ListeSp) <- NAMES
  #rownames(MatSp)=seq(1, nrow(MatSp),1); for (i in 1:length(ListeSp[,1])) rownames(MatSp)[which(apply(MatSp,1,paste,collapse=" ")==ListeSp[i,2])]=rownames(ListeSp)[i]
  
  if (Npop>1) {nsamBrk<-SamplePop*0   ##### Defining interval for extracting individual by population
  nsamBrk[1]<- 1
  for(i in 2:length(SamplePop))
    nsamBrk[i]<- nsamBrk[i-1]+SamplePop[i-1]
  nsamBrk<- c(nsamBrk, sum(SamplePop)+1)}else{nsamBrk<-c(1,SamplePop+1)}
  
  for (k in 1:Npop){    ##### Counting individuals by populations
    SubMat <- table(apply(MatSp[nsamBrk[k]:(nsamBrk[k+1]-1),],1,paste,collapse=" "))
    for (i in 1:length(SubMat)) ListeSp[which(ListeSp[,2]==names(SubMat)[i]),k+2] <- as.numeric(SubMat[i])
  }
  
  Stats_Sp[[1]] <- as.matrix(ListeSp[,3:(3+Npop-1)]); NAMES=NULL; for (i in 1:Npop) NAMES <- c(NAMES,paste("pop",i)); colnames(Stats_Sp[[1]]) <- NAMES; mode(Stats_Sp[[1]]) <- "numeric"
  Stats_Sp[[2]] <- ListeSp[,2]; names(Stats_Sp[[2]]) <- "Species identity"; names(Stats_Sp[[2]]) <- rownames(ListeSp)
  Stats_Sp[[3]] <- dim(ListeSp)[1]; names(Stats_Sp[[3]]) <- " Species richness"
  
  SpeciesRich[1,1] <- Stats_Sp[[3]]
  SpeciesRich[1,2:(Npop+1)] <- as.vector(specnumber(t(Stats_Sp[[1]])))
  
  NielsenSpDiv[1,1] <- Nielsen(t(as.matrix(rowSums(Stats_Sp[[1]]))))
  NielsenSpDiv[1,2:(Npop+1)] <- Nielsen(t(Stats_Sp[[1]]))
  
  #### Matching species names with genetic sequences
  rownames(Output_Sp[[3]]) <- seq(1,nrow(Output_Sp[[3]])); for (i in 1:length(ListeSp[,1])) rownames(Output_Sp[[3]])[which(apply(Output_Sp[[3]],1,paste,collapse=" ")==ListeSp[i,2])]=rownames(ListeSp)[i]
  rownames(Output_G[[3]]) <- rownames(Output_Sp[[3]])
  
  ##############################
  ##### Genetic level stats
  ##############################
  Stats_G_Full=list()
  Stats_G=list()
  ListeG <- matrix(0,length(table(apply(MatG,1,paste,collapse=" "))),Npop+2)
  ListeG[,1] <- table(apply(MatG,1,paste,collapse=" "))
  ListeG[,2] <- names(table(apply(MatG,1,paste,collapse=" ")))
  
  for (k in 1:Npop){  ##### Counting individuals by populations
#     MatG[nsamBrk[k]:(nsamBrk[k+1]-1),]
    SubMat <- table(apply(MatG[nsamBrk[k]:(nsamBrk[k+1]-1),],1,paste,collapse=" "))
    for (i in 1:length(SubMat)) ListeG[which(ListeG[,2]==names(SubMat)[i]),k+2] <- as.numeric(SubMat[i])
  }

  Stats_G[[1]]<- as.matrix(ListeG[,3:(3+Npop-1)]); colnames(Stats_G[[1]])=NAMES; mode(Stats_G[[1]]) <- "numeric"
  Stats_G[[2]]<- ListeG[,2]; names(Stats_G[[2]])="Haplotype identity"
  Stats_G[[3]]<- dim(ListeG)[1]; names(Stats_G[[3]])=" Haplotype richness"
  Stats_G_Full[[1]]<- Stats_G
  names(Stats_G_Full)[[1]]<- "All species"

  HapRich <- matrix(NA,Stats_Sp[[3]],Npop+1)  ## will store haplotype richness for each species global and within each locality
  HapNielsen <- matrix(NA,Stats_Sp[[3]],Npop+1) ### will store Nielsen estimate of each species global and for each locality
  Abundances <- matrix(NA,Stats_Sp[[3]],Npop+1) ### will store the total number of individuals and the nbe of individuals in each locality

  for (h in 1:Stats_Sp[[3]]) {  #### Extracting haplotype abundance and composition for all species
    Stats_G <- list()
    SubMatG <- Output_G[[3]][which(rownames(Output_G[[3]])==rownames(Stats_Sp[[1]])[h]),]
    if(is.matrix(SubMatG)) {SubMatG=SubMatG}else{SubMatG=t(as.matrix(SubMatG)); rownames(SubMatG) <-  rownames(Stats_Sp[[1]])[h]}     
    ListeG <- matrix(0,length(table(apply(SubMatG,1,paste,collapse=" "))),Npop+2)
    ListeG[,1] <- table(apply(SubMatG,1,paste,collapse=" "))
    ListeG[,2] <- names(table(apply(SubMatG,1,paste,collapse=" ")))

    if (Npop==1){   ### if there is only one locality
      ListeG[,3]=ListeG[,1]
    }else{
      for (k in 1:Npop){  ##### Counting individuals by locality
        #     MatG[nsamBrk[k]:(nsamBrk[k+1]-1),]
        SubMatG <- Output_G[[3]][nsamBrk[k]:(nsamBrk[k+1]-1),]
        SubMatG <- SubMatG[which(rownames(SubMatG)==rownames(Stats_Sp[[1]])[h]),]
        if(!is.matrix(SubMatG)) {SubMatG=t(as.matrix(SubMatG))}  #; rownames(SubMatG) <-  rownames(Stats_Sp[[1]])[h]
        if (nrow(SubMatG)>0){         #### If Species 'h' is present in locality 'k'
          SubMat <- table(apply(SubMatG,1,paste,collapse=" "))
          for (i in 1:length(SubMat)) ListeG[which(ListeG[,2]==names(SubMat)[i]),k+2]=as.numeric(SubMat[i])
        }
      }   # end k 
    }

    Stats_G[[1]] <- as.matrix(ListeG[,3:(3+Npop-1)]); if(dim(ListeG)[1]==1) {Stats_G[[1]]=t(Stats_G[[1]])}; colnames(Stats_G[[1]])=NAMES ;mode(Stats_G[[1]]) <- "numeric"
    Stats_G[[2]] <- ListeG[,2]; names(Stats_G[[2]])="Haplotype identity"
    Stats_G[[3]] <- dim(ListeG)[1]; names(Stats_G[[3]])=" Haplotype richness"
    
    ## Haplotype richness
    HapRich[h,1] <- Stats_G[[3]]
    HapRich[h,2:(Npop+1)] <- as.vector(specnumber(t(Stats_G[[1]])))
    
    ## Nielsen estimate of diversity
    if(dim(Stats_G[[1]])[1]>1)  HapNielsen[h,1] <- round(Nielsen(t(as.matrix(rowSums(Stats_G[[1]])))),1)
    HapNielsen[h,2:(Npop+1)] <- round(Nielsen(t(Stats_G[[1]])),1)
    
      ## Haplotype abundances
    Abundances[h,1] <- sum(Stats_G[[1]])
    Abundances[h,2:(Npop+1)] <- as.vector(colSums(Stats_G[[1]]))
    
    Stats_G_Full[[h+1]] <- Stats_G
    names(Stats_G_Full)[[h+1]] <- paste("sp",h)
  } # end h


# } # end j loop --> number of draws

NAMES=c(); NAMES=c(NAMES,"SpRichTot"); for (i in 1:Npop) NAMES=c(NAMES,paste("SpRichLoc",i)); colnames(SpeciesRich) <- NAMES
NAMES=c(); NAMES=c(NAMES,"SpNielsenTot"); for (i in 1:Npop) NAMES=c(NAMES,paste("SpNielsenLoc",i)); colnames(NielsenSpDiv) <- NAMES
NAMES=c(); NAMES=c(NAMES,"HapRichTot"); for (i in 1:Npop) NAMES=c(NAMES,paste("HapRichLoc",i)); colnames(HapRich) <- NAMES
NAMES=c(); NAMES=c(NAMES,"HapNielsenTot"); for (i in 1:Npop) NAMES=c(NAMES,paste("HapNielsenLoc",i)); colnames(HapNielsen) <- NAMES
NAMES=c(); NAMES=c(NAMES,"AbTot"); for (i in 1:Npop) NAMES=c(NAMES,paste("AbLoc",i)); colnames(Abundances) <- NAMES


#### Stats at community between-species level
SpeciesRich
NielsenSpDiv

### Stats at population intra-species level
HapRich
HapNielsen
Abundances




colMeans(HapNielsen, na.rm=T)

NielAbSp=c()
for (h in 1:Npop) NielAbSp=c(NielAbSp,mean(HapNielsen[which(Abundances[,2]/sum(Abundances[,2])>0.05),h+1])) 




