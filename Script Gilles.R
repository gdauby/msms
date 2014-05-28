
setwd("C:/WORKS/Logiciels/msms/bin")
source("ms.r")


#### Input Parameter for msms simulation

ndraws=1
theta=10
SamplePop <- c(5,5)
nsam=sum(SamplePop)
Npop <- length(SamplePop)
Mig=3

##### Input parameter for paritionning species and genetic diversity
LimitSp <- 5   #### Ratio mutation rate/speciation rate
LimitSp <-1-1/LimitSp


###### Coalescence simulation using 'ms'
if (Npop>1) SIM <- paste("msms",nsam,ndraws,"-t",theta,"-I",Npop,apply(as.matrix(SamplePop),2,paste,collapse=" "), Mig,   " > out2.txt")
if (Npop==1)SIM <-  paste("msms",nsam,ndraws,"-t",theta," > out2.txt")

shell(SIM)

OUTPUT <- ms.inp.multi(nsam=nsam,ndraws=ndraws,ms.output.file="out2.txt")
#read.ms.output(file="out2.txt")

GeneticDiversity=list()
OUTPUT_SP=list()
OUTPUT_G=list()
if(Npop>1) FST=matrix(NA,ndraws,2)
for (j in 1:ndraws) {
  ###### ##############################
  ##### Partitionning sample data into species and haplotypes according to the threshold given
  MatSp <- OUTPUT[[j]]$sample[,which(OUTPUT[[j]]$mutations>LimitSp)]
  MatG <- OUTPUT[[j]]$sample[,which(OUTPUT[[j]]$mutations<=LimitSp)]
  
  Output_Sp=list()
  Output_G=list()
  Output_Sp[[1]]=length(which(OUTPUT[[j]]$mutations>LimitSp)); Output_G[[1]]=length(which(OUTPUT[[j]]$mutations<=LimitSp)) ### Number of mutations/speciation events
  Output_Sp[[2]]=OUTPUT[[j]]$mutations[which(OUTPUT[[j]]$mutations>LimitSp)]; Output_G[[2]]=OUTPUT[[j]]$mutations[which(OUTPUT[[j]]$mutations<=LimitSp)]  ### Position of speciation/mutation events in the array
  Output_Sp[[3]]=MatSp; Output_G[[3]]=MatG
  names(Output_Sp)=names(OUTPUT[[j]]) ; names(Output_G)=names(OUTPUT[[j]])
  
  OUTPUT_SP[[j]]=Output_Sp; OUTPUT_G[[j]]=Output_G   ####### Datasets arranged into format expected by function makeStats 
  
  ##### Species level stats
  Stats_Sp=list()
  ListeSp <- matrix(0,length(table(apply(MatSp,1,paste,collapse=" "))),Npop+2)  #### ListeSp gives the list of species with abundance and 'sequence' composition
  ListeSp[,1]=table(apply(MatSp,1,paste,collapse=" "))
  ListeSp[,2]=names(table(apply(MatSp,1,paste,collapse=" ")))
  NAMES =NULL; for (i in 1:nrow(ListeSp)) NAMES=c(NAMES,paste("Sp",i));rownames(ListeSp)=NAMES
  #rownames(MatSp)=seq(1, nrow(MatSp),1); for (i in 1:length(ListeSp[,1])) rownames(MatSp)[which(apply(MatSp,1,paste,collapse=" ")==ListeSp[i,2])]=rownames(ListeSp)[i]
  
  if (Npop>1) {nsamBrk<-SamplePop*0   ##### Defining interval for extracting individual by population
  nsamBrk[1]<-1
  for(i in 2:length(SamplePop))
    nsamBrk[i]<-nsamBrk[i-1]+SamplePop[i-1]
  nsamBrk<-c(nsamBrk, sum(SamplePop)+1)}else{nsamBrk<-c(1,SamplePop+1)}
  
  for (k in 1:Npop){    ##### Counting individuals by populations
#     MatSp[nsamBrk[k]:(nsamBrk[k+1]-1),]
    SubMat<-table(apply(MatSp[nsamBrk[k]:(nsamBrk[k+1]-1),],1,paste,collapse=" "))  
    for (i in 1:length(SubMat)) ListeSp[which(ListeSp[,2]==names(SubMat)[i]),k+2]=as.numeric(SubMat[i])
  }
  
  Stats_Sp[[1]]=as.matrix(ListeSp[,3:(3+Npop-1)]); NAMES=NULL; for (i in 1:Npop) NAMES=c(NAMES,paste("pop",i)); colnames(Stats_Sp[[1]])=NAMES; mode(Stats_Sp[[1]]) <- "numeric"
  Stats_Sp[[2]]=ListeSp[,2]; names(Stats_Sp[[2]])="Species identity"; names(Stats_Sp[[2]])=rownames(ListeSp)
  Stats_Sp[[3]]=dim(ListeSp)[1]; names(Stats_Sp[[3]])=" Species richness"
  
  #### Matching species names with genetic sequences
  rownames(Output_Sp[[3]])=seq(1,nrow(Output_Sp[[3]])); for (i in 1:length(ListeSp[,1])) rownames(Output_Sp[[3]])[which(apply(Output_Sp[[3]],1,paste,collapse=" ")==ListeSp[i,2])]=rownames(ListeSp)[i]
  rownames(Output_G[[3]])=rownames(Output_Sp[[3]])
  
  ##### Genetic level stats
  Stats_G_Full=list()
  Stats_G=list()
  ListeG <- matrix(0,length(table(apply(MatG,1,paste,collapse=" "))),Npop+2)
  ListeG[,1]=table(apply(MatG,1,paste,collapse=" "))
  ListeG[,2]=names(table(apply(MatG,1,paste,collapse=" ")))
  
  for (k in 1:Npop){  ##### Counting individuals by populations
    MatG[nsamBrk[k]:(nsamBrk[k+1]-1),]
    SubMat<-table(apply(MatG[nsamBrk[k]:(nsamBrk[k+1]-1),],1,paste,collapse=" "))
    for (i in 1:length(SubMat)) ListeG[which(ListeG[,2]==names(SubMat)[i]),k+2]=as.numeric(SubMat[i])
  }
  
  Stats_G[[1]]=as.matrix(ListeG[,3:(3+Npop-1)]); colnames(Stats_G[[1]])=NAMES; mode(Stats_G[[1]]) <- "numeric"
  Stats_G[[2]]=ListeG[,2]; names(Stats_G[[2]])="Haplotype identity"
  Stats_G[[3]]=dim(ListeG)[1]; names(Stats_G[[3]])=" Haplotype richness"
  Stats_G_Full[[1]]=Stats_G
  names(Stats_G_Full)[[1]]="All species"
  
  HapRich=c()
  for (h in 1:Stats_Sp[[3]]) {  #### Extracting haplotype abundance and composition for all species
    Stats_G=list()
    SubMatG=Output_G[[3]][which(rownames(Output_G[[3]])==rownames(Stats_Sp[[1]])[h]),]
    if(is.matrix(SubMatG)) {SubMatG=SubMatG}else{SubMatG=t(as.matrix(SubMatG))}      
        
    Stats_G[[1]]=table(apply(SubMatG,1,paste,collapse=" ")); for (i in 1:length(Stats_G[[1]])) names(Stats_G[[1]])[i]=paste("Hap",i)
    Stats_G[[2]]=names(table(apply(SubMatG,1,paste,collapse=" "))); names(Stats_G[[2]])=names(Stats_G[[1]])
    Stats_G[[3]]=length(table(apply(SubMatG,1,paste,collapse=" "))); names(Stats_G[[3]])=" Haplotype richness"
    HapRich=c(HapRich,Stats_G[[3]])
    Stats_G_Full[[h+1]]=Stats_G
    names(Stats_G_Full)[[h+1]]=paste("sp",h)
  }

  GeneticDiversity[[j]]=cbind(HapRich,as.vector(rowSums(Stats_Sp[[1]])))
  names(GeneticDiversity)[[j]]=j
  rownames(GeneticDiversity[[j]])= rownames(Stats_Sp[[1]])
  colnames(GeneticDiversity[[j]])= c("Haplotype richness","Species abundance")


}

STATS_Sp <- makeStats(OUTPUT_SP, SamplePop)
STATS_G <- makeStats(OUTPUT_G, SamplePop)
FST[,1]=STATS_Sp$fst
FST[,2]=STATS_G$fst



