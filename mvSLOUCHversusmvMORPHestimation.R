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

phyltree50<-sim.bd.taxa(50,1,1,0)[[1]]
plot(phyltree50)
which(phyltree50$tip.label=="t50") ## here we perhaps need to change the node dependen on the random seed
phyltree51<-add_new_tip_split(phyltree,4,0.00001)
BMdata51<-simulBMProcPhylTree(phyltree51, X0=matrix(c(0,0),ncol=1), Sigma=diag(1,2,2))
BMdata50<-BMdata[-1,]
resmvsl50<-BrownianMotionModel(phyltree50,BMdata50);resmvsl50$ParamSummary$LogLik;resmvsl50$ParamsInModel$vX0;resmvsl50$ParamSummary$StS
resmvsl51<-BrownianMotionModel(phyltree51,BMdata51);resmvsl51$ParamSummary$LogLik;resmvsl51$ParamsInModel$vX0;resmvsl51$ParamSummary$StS
resmvmorBM150<-mvMORPH::mvBM(phyltree50,BMdata50, model="BM1", scale.height = FALSE) ## notice the big difference betweem the two estimates
resmvmorBM151<-mvMORPH::mvBM(phyltree51,BMdata51, model="BM1", scale.height = FALSE) ## perhaps saying something very different

