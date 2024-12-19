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

## TODO:
## 1) decide on dataset
## 2) tree as .tre file
## 3) data as csv file

library(PCMBase)
library(PCMBaseCpp)
library(mvSLOUCH)
library(mvMORPH)
library(ape)

## mvData<-read.table("",header=TRUE,sep=)
## phyltree<-ape::read.tree("")
set.seed(5)

## below we observe output with errorsome likelihood
## this is not an error on mvSLOUCH's side but a 
## defense mechanism against situations with potentially
## high numerical errors in the calculation
## with the aim of prompting the user to investigate what is the root
## of the problem
## it should be pointed out that mvSLOUCH has analytical estimators in the BM case
mvslBMres_error<-mvSLOUCH::BrownianMotionModel(phyltree, mvData)

## mvMORPH returns without errors 
## but mvMORPH (from investigation of source code) finds the BM's model's
## parameters through numerical optimization of the likelihood function
## furthermore, mvMORPH uses a different way to calculate the likelihood
## which will be robust to short internal branch lenghts, but has
## far greater computational complexity
mvmorBM1result_1<-mvMORPH::mvBM(phyltree,mvData, model="BM1", scale.height = FALSE)

## There are a number of ways that we can rectify the situation
## so that mvSLOUCH produces a valid likelihood.
## The first uses PCMBase's mechanism, avialable through mvSLOUCH, to treat short branches
## as polytomies (or in other words skipped). The parameter min_bl defines what is meant by a short branch.
## It can be set to the smallest value that returns a valid likelihood.
## Notice that we are changing the tree, but given, that the original
## branch was short anyway we should not suspect that this would have
## a significant impact on the value of the likelihood. This method
## works for all models in mvSLOUCH.
mvslBMres_validlik_1<-mvSLOUCH::BrownianMotionModel(phyltree, mvData, min_bl=0.0105)

## We can also equivalenty set the PCMBase option prior to calling mvSLOUCH
## saving options to be adjusted
PCMBase_orig_threshold<-PCMBase::PCMOptions()$PCMBase.Threshold.Skip.Singular
PCMBase_orig_skip<-PCMBase::PCMOptions()$PCMBase.Skip.Singular
## ============================================================================
options(PCMBase.Threshold.Skip.Singular=0.0105)
options(PCMBase.Skip.Singular=TRUE) ## logical statement if such branches should be treated as polytomies/skipped
mvslBMres_validlik_2<-mvSLOUCH::BrownianMotionModel(phyltree, mvData)
## restoring original options
options(PCMBase.Threshold.Skip.Singular=PCMBase_orig_threshold)
options(PCMBase.Skip.Singular=PCMBase_orig_skip) 
## ============================================================================

## We can achieve the same by manually adjusting the offending branch lengths in the phylogeny.
## Notice that this makes the phylogeny non-ultrametric, however, mvSLOUCH is agnostic
## to this.
## Again we are changing the tree, but given, that the original
## branch was short anyway we should not suspect that this would have
## a significant impact on the value of the likelihood. This method
## works for all models in mvSLOUCH.
phyltree2<-phyltree
phyltree2$edge.length<-0.0105
mvslBMres_validlik_3<-mvSLOUCH::BrownianMotionModel(phyltree2, mvData)

## We can also collapse the branch into a polytomy.
phyltree3<-ape::di2multi(phyltree,tol=0.01)
mvslBMres_validlik_4<-mvSLOUCH::BrownianMotionModel(phyltree3, mvData)

## An alternaitive approach is to try to linearly rescale the data
## so that the variance along a branch will kill the short branch effect.
## This however, might not work for OU models, as there branch lengths
## are exponentially transfored in variance calculations.
## However we have to remember to rescale model parameters afterwords.
rescale_factor<-0.0105/0.00999495
mvData2<-rescale_factor*mvData2
mvslBMres_validlik_5<-mvSLOUCH::BrownianMotionModel(phyltree,mvData2)
mvmorBM1result_2<-mvMORPH::mvBM(phyltree,mvData2, model="BM1", scale.height = FALSE)
## Correcting mvSLOUCH's estimates for original data scale, this is much more complex for OU models.
## The likelihood remains uneffected by an affine transformation (as spelled out in our Evolution ms).
mvslBMres_validlik_5_originalscale<-mvslBMres_validlik_5
mvslBMres_validlik_5_originalscale$ParamsInModel$vX0<-mvslBMres_validlik_5_originalscale$ParamsInModel$vX0/rescale_factor
mvslBMres_validlik_5_originalscale$ParamsInModel$Sxx<-mvslBMres_validlik_5_originalscale$ParamsInModel$Sxx/rescale_factor
mvslBMres_validlik_5_originalscale$ParamSummary$StS<-mvslBMres_validlik_5_originalscale$ParamSummary$StS/(rescale_factor^2)
mvslBMres_validlik_5_originalscale$ParamSummary$confidence.interval$regression.summary$X0.regression.confidence.interval<-mvslBMres_validlik_5_originalscale$ParamSummary$confidence.interval$regression.summary$X0.regression.confidence.interval/rescale_factor
mvslBMres_validlik_5_originalscale$ParamSummary$confidence.interval$regression.summary$regression.covariance.matrix<-mvslBMres_validlik_5_originalscale$ParamSummary$confidence.interval$regression.summary$regression.covariance.matrix/(rescale_factor^2)



## We can end by comparing the results of mvSLOUCH and mvMORPH
print("mvSLOUCH invalid likelihood")
print(mvslBMres_error$ParamsInModel);print(mvslBMres_error$ParamSummary$LogLik)
print("=====================================================================")

print("mvSLOUCH valid likelihood: collapse internally short branches")
print(mvslBMres_validlik_1$ParamsInModel);print(mvslBMres_validlik_1$ParamSummary$LogLik)
print("=====================================================================")

print("mvSLOUCH valid likelihood: set PCMBase's option globally ")
print(mvslBMres_validlik_2$ParamsInModel);print(mvslBMres_validlik_2$ParamSummary$LogLik)
print("=====================================================================")

print("mvSLOUCH valid likelihood: extend branch lengths in phylogeny ")
print(mvslBMres_validlik_3$ParamsInModel);print(mvslBMres_validlik_3$ParamSummary$LogLik)
print("=====================================================================")

print("mvSLOUCH valid likelihood: collapse internal node in phylogeny to polytomy ")
print(mvslBMres_validlik_4$ParamsInModel);print(mvslBMres_validlik_4$ParamSummary$LogLik)
print("=====================================================================")

print("mvSLOUCH valid likelihood: rescale data ")
print("mvSLOUCH results for scaled data ")
print(mvslBMres_validlik_5$ParamsInModel);print(mvslBMres_validlik_5$ParamSummary$LogLik)
print("mvSLOUCH results for data on original scale")
print(mvslBMres_validlik_5_originalscale$ParamsInModel);print(mvslBMres_validlik_5_originalscale$ParamSummary$LogLik)
print("mvMORPH results for scaled data ")
print("=====================================================================")

print("mvMORPH results ")
print("mvMORPH results for original data ")
print(mvmorBM1result_1)
print("mvMORPH results for scaled data ")
print(mvmorBM1result_2)
print("=====================================================================")



