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

library(ape)
library(mvSLOUCH)
library(TreeSim)
library(mvMORPH)
source("add_cherry_tip.R")

vN<-c(50,100,500) ## perhaps other tree sizes? 1000 is the limit for mvMORPH rpf
num_repeats<-100 ## see what can be done
BMparams<-
f_singularbranchlen<-0.000001
f_normalbranchlen<-0.01

sapply(vN,function(n,){
sapply(1:numrepeats,function(i,){
at each iteration save the seed
randomseed<-.Random.seed;set.seed(randomseed, "L'Ecuyer-CMRG");RNG_kind<-RNGkind();RNG_version<-getRversion()

    phyltreen<-TreeSim::sim.bd.taxa(n,1,1,0)[[1]]	
    phyltreen1singular<-add_new_tip_split(phyltreen,tip_id,f_singularbranchlen)
    phyltreen1normal<-add_new_tip_split(phyltreen,tip_id,f_normalbranchlen) ## will need investigating which tip_id can be split, i.e., has height>f_normalbranchlen if none-resimulate tree
    
    BMdatanInternal<-mvSLOUCH::simulBMProcPhylTree(phyltreen1, X0=BMparams$X0, Sigma=BMparams$Sigma,dropInternal = FALSE)
    ## identify split node of new additional nodes, get their values and simulate the new added node by YUIMA NEED to set step size for singular branch length
    BNdatan1singular<-
    BMdatan1normal<-
    BMdatan <- drop internal from BMdatanInternal

    ## we can additionally remove a random tip
    ## instead of keeping a single f_normalbranchlen we can have a vector of them-and check how the liklihood changes with shortening and shortening of this branch
    
    resmvsln<-mvSLOUCH::BrownianMotionModel(phyltreen,BMdatan) ## if -Inf rerun but keep counter of errors
    resmvsln1singular<-mvSLOUCH::BrownianMotionModel(phyltreen1singular,BMdatan1singular) ## should be -Inf
    resmvsln1normal<-mvSLOUCH::BrownianMotionModel(phyltreen1normal,BMdatan1normal) ## if -Inf rerun but keep counter of errors
    compare resmvsln$ParamSummary$LogLik and resmvsln1normal$ParamSummary$LogLik and estimated parameters

    resmvmorBM1n1singular_pic<-mvMORPH::mvBM(phyltreen1singular,BMdatan1singular, model="BM1", scale.height = FALSE,method="pic") 
    resmvmorBM1n1normal_pic<-mvMORPH::mvBM(phyltreen1normal,BMdatan1normal, model="BM1", scale.height = FALSE,method="pic") 
    resmvmorBM1n_pic<-mvMORPH::mvBM(phyltreen,BMdatan, model="BM1", scale.height = FALSE,method="pic") 

    resmvmorBM1n1singular_rpf<-mvMORPH::mvBM(phyltreen1singular,BMdatan1singular, model="BM1", scale.height = FALSE,method="rpf") 
    resmvmorBM1n1normal_rpf<-mvMORPH::mvBM(phyltreen1normal,BMdatan1normal, model="BM1", scale.height = FALSE,method="rpf") 
    resmvmorBM1n_rpf<-mvMORPH::mvBM(phyltreen,BMdatan, model="BM1", scale.height = FALSE,method="rpf") 

    compare loglikelihood between resmvmorBM1n1singular_pic, resmvmorBM1n1normal_pic, resmvmorBM1n_pic, resmvmorBM1n1singular_rpf, resmvmorBM1n1normal_rpf, resmvmorBM1n_rpf, 
    illustrate sensitvity, robustness and cf with mvSLOUCH
    
    save in a RData file also random seed things randomseed=randomseed,RNG_kind=RNG_kind;RNG_version=RNG_version,i=i
}
}
