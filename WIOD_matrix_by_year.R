

## 2015-07-04
rm(list = ls())    # Remove all
setwd("D:/Copy/GVC/WIOD/R")
rawd <- "D:/Copy/GVC/WIOD/wiot/"    # Raw data Directory

library(foreign)
library(matlab)
library(xlsx)


## WIOT includes 41 countries (including RoW) with 35 industries over 1995-2011
S <- 35
N <- 41
SN <- S*N
period <- c(1995:2011)


## Import data from WIOD and assign IDs
wiot  <- read.csv(paste0(rawd,"wiot1995_meta.csv"), header=FALSE, skip=6)
cid  <- toupper(unique(wiot$V3[1:SN]))
iid  <- c(1:35)
ciid <- toupper(paste(wiot$V3[1:SN], as.character(iid), sep=""))
industry <- as.character(wiot$V2[1:SN])
save(cid, iid, ciid, industry, file="wiot_meta.Rdata")
#writeMat(wiot_meta.mat, cid=cid, iid=iid, ciid=ciid, industry=industry)


## Compute the values for defined matrices

for(yr in period) {
      wiot <- read.csv(paste0(rawd,"wiot",yr,".csv"), header=FALSE, as.is=TRUE)
      wiot <- data.matrix(wiot)
      
      y <- wiot[nrow(wiot), 1:SN]                           # y = World Output Vector
      names(y) <- ciid
      y.nzero <- y[y>0]
      ciid.nzero <- names(y.nzero)
            
      va <- y - colSums(wiot[1:SN, 1:SN])                   # va = World VA Vector
      names(va) <- ciid
      va.nzero <- va[va>0]
      
      M <- wiot[1:SN, 1:SN]                                 # M = Intermediate IO Matrix
      dimnames(M) <- list(ciid, ciid)
      
      A <- zeros(SN, SN)                                    # A = Input Coefficient Matrix
      dimnames(A) <- list(ciid, ciid)
      A[, ciid.nzero] <- M[, ciid.nzero] / (ones(SN,1) %*% y.nzero)
      
      r  <- 1 - colSums(A)                                  # r = Ratio of Value-added to Total Output
      names(r) <- ciid
      
      Leon <- diag(length(ciid)) - A
      LeonInv <- solve(Leon)                                # LeonInv = Leontief Inverse

      FD <- zeros(SN, N)                                    # FD = Final Demand Matrix for N countries
      dimnames(FD) <- list(ciid, cid)
      for(i in 1:N) {
            FD[, i] <- wiot[1:SN, (SN+1+5*(i-1)):(SN+5*i)] %*% ones(5,1)
      }
      
      FDD <- zeros(SN, SN)                                  # FDD = partial diagonalization of Final Demand
      dimnames(FDD) <- list(ciid, ciid)
      for (j in 1:N) {                                      # j indicates column
            for (i in 1:N) {                                # i indicates row
                  FDD[(1+S*(i-1)):(S*i), (1+S*(j-1)):(S*j)] <- diag(FD[(1+S*(i-1)):(S*i), j])
            }
      }
      
      Y.alloc <- LeonInv %*% FDD                            # Y.alloc = Output Allocation Matrix
      dimnames(Y.alloc) <- list(ciid, ciid)
      
      VA.alloc <- diag(r) %*% Y.alloc                       # VA.alloc = Value-added Allocation Matrix
      dimnames(VA.alloc) <- list(ciid, ciid)
      
      yInv <- zeros(SN,1)
      names(yInv) <- ciid
      yInv[ciid.nzero] <- 1/y.nzero
      OS <- diag(yInv) %*% Y.alloc                          # OS = Output Share Matrix
      dimnames(OS) <- list(ciid, ciid)
      
      vaInv <- zeros(SN,1)
      names(vaInv) <- ciid
      vaInv[ciid.nzero] <- 1/va.nzero
      VAS <- diag(vaInv) %*% VA.alloc                       # VAS = Value-added Share Matrix
      dimnames(VAS) <- list(ciid, ciid)

      
      # Save objects
      save(wiot,S,N,SN,cid,iid,ciid,industry,M,A,y,va,r,LeonInv,FD,FDD,Y.alloc,VA.alloc,OS,VAS, 
           file=paste0("WIOD_matrix_",yr,".RData"))
}


