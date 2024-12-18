# Set up

library(PCMBaseCpp)
library(mvSLOUCH)
library(mvMORPH)

set.seed(5)

# Load workspace (Workspace.RData) with load() function

# Global regime
OU1Results<-estimate.evolutionary.model(mvStree,mvData,regimes = NULL,repeats=5,
                                        model.setups=model_setups,doPrint=TRUE,
                                        pESS=NULL,maxiter=c(10,50,100))

# Sort models according to AICc
OU1AIC<-rep(NA,length(OU1Results$testedModels))
for (i in 1:length(OU1AIC)) {
  OU1AIC[i]<-OU1Results$testedModels[[i]]$aic.c
}
# Note that many models failed to run ("Inf")
sort(OU1AIC)
# BM between them
OU1Results$testedModels[1]
# BM works under mvMORPH though
mvBM(Tree,mvData, model="BM1", scale.height = T)





# Global regime with rescaled data
OU1Rescale<-estimate.evolutionary.model(mvStree,mvData*10,regimes = NULL,repeats=5,
                                        model.setups=model_setups,doPrint=TRUE,
                                        pESS=NULL,maxiter=c(10,50,100))

# Sort models according to AICc
OUrAIC<-rep(NA,length(OU1Rescale$testedModels))
for (i in 1:length(OUrAIC)) {
  OUrAIC[i]<-OU1Rescale$testedModels[[i]]$aic.c
}
# Note that all models ran this time
sort(OUrAIC)
# Including BM
OU1Rescale$testedModels[1]
# mvSLOUCH & mvMORPH report matching results under rescaled data
mvBM(Tree,mvData*10, model="BM1", scale.height = T)



















































