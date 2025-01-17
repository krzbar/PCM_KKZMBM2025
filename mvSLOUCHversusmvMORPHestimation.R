## These are the R scripts and numerical results accompanying 
## CHECK FINAL AUTHORS
## Bartoszek, Brahmantio, Munoz-Duran, Fuentes-Gonzalez, Chi Kiang, Pienaar, and Polly
## TITLE?
## "TITLE?"

## The R setup for the manuscript was as follows: R version 4.4.1 (2024-06-14) 
## Platform: x86_64-pc-linux-gnu (64-bit) Running under: openSUSE Leap 15.6

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## with a singular tree due to a single short tip branch numerical mvMORPH will return estimates that are very
## different from those when this branch is removed. mvMORPH in the near singular case might also return 
## comment that convergence attained. mvSLOUCH returns -Inf for likelihood prompting the user to hopefully
## investigate the tree

## TODO: set random number generator

library(ape)
library(mvSLOUCH)
library(TreeSim)
library(mvMORPH)
source("add_cherry_tip.R")

rexp(1)
n<-30
rm(.Random.seed)
my_random_seed<-1.890968652176264530596
#my_random_seed<-rexp(1)
set.seed(my_random_seed)
phyltreen<-TreeSim::sim.bd.taxa(n,1,1,0)[[1]]
#plot(phyltreen)

#tip_id<-which(phyltreen$tip.label=="t40") # 1000
tip_label<-"t15"
#tip_label<-paste0("t",sample(1:n,1))
tip_id<-which(phyltreen$tip.label==tip_label) #54321
phyltreen1<-add_new_tip_split(phyltreen,tip_id,0.000001)
BMdatan1<-mvSLOUCH::simulBMProcPhylTree(phyltreen1, X0=matrix(c(10,0),ncol=1), Sigma=rbind(c(0.1,0),c(-10,10)))
BMdatan<-BMdatan1[-1,]
resmvsln1<-mvSLOUCH::BrownianMotionModel(phyltreen1,BMdatan1);resmvsln1$ParamSummary$LogLik;resmvsln1$ParamsInModel$vX0;resmvsln1$ParamSummary$StS
resmvsln<-mvSLOUCH::BrownianMotionModel(phyltreen,BMdatan);resmvsln$ParamSummary$LogLik;resmvsln$ParamsInModel$vX0;resmvsln$ParamSummary$StS
resmvmorBM1n1<-mvMORPH::mvBM(phyltreen1,BMdatan1, model="BM1", scale.height = FALSE) 
resmvmorBM1n<-mvMORPH::mvBM(phyltreen,BMdatan, model="BM1", scale.height = FALSE) 


print("Original tree and data, mvSLOUCH raises error.")
print(resmvsln1)
print("==========================================")

print("Tree and data with short tip branch species removed, mvSLOUCH estimates without error.")
print(resmvsln)
print("==========================================")

print("Tree and data with short tip branch species removed, mvMORPH estimates same as mvSLOUCH.")
print(resmvmorBM1n)
print("==========================================")

print("Original tree and data, mvMORPH returns different estimates.")
print(resmvmorBM1n1)
print("==========================================")
