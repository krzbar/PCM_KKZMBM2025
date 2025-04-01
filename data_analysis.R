library(ape)
library(mvSLOUCH)
library(mvMORPH)

# load data and tree
load("mvData.RData")
Tree <- read.tree("Tree.tre")

# ================================ fit mvSLOUCH ================================
ML_mvs <- BrownianMotionModel(phyltree = Tree, mData = mvData)

# ----------------- mvSLOUCH with min_bl
ML_mvs_mbl <- BrownianMotionModel(phyltree = Tree, mData = mvData, min_bl = 0.0105)

# ----------------- polytomy
# create a new tree based on this polytomy
l_min <- Tree$edge.length[Tree$edge.length < 0.0105]
Tree_pol <- di2multi(Tree, tol = 0.0105)

# fix the tip branch to make it ultrametric
tips_id <- which(Tree_pol$tip.label %in% c("LGY", "LCU"))

Tree_pol$edge.length[which(Tree_pol$edge[,2] %in% tips_id)] <- 
  Tree_pol$edge.length[which(Tree_pol$edge[,2] %in% tips_id)] + l_min
plot(Tree_pol)

ML_mvs_pol <- BrownianMotionModel(phyltree = Tree_pol, mData = mvData)


# ============== check the log-likelihood by hand ==============
llik_byhand <- function(Tree, X0_hat, Sigma_hat, mData){
  Tmat <- ape::vcv(Tree)
  Tmat_inv <- solve(Tmat)
  n <- Ntip(Tree)
  
  # maximum log-likelihood
  V <- kronecker(Sigma_hat,Tmat)
  D <- kronecker(diag(2), rep(1,n))
  v1 <- c(mData) - D%*%X0_hat
  mloglik <- -(n*2)/2*log(2*pi) - log(det(V))/2 - (t(v1)%*%solve(V)%*%v1)/2
  
  c(mloglik)
}



# ----------------- Try multiplying data by a constant
# Multiplier constant c
c <- 0.0105/min(Tree$edge.length)

ML_mvs_mod <- BrownianMotionModel(phyltree = Tree, mData = c*mvData)
ML_mvs_mod_res <- list()

# don't forget to rescale back
ML_mvs_mod_res[[1]] <- ML_mvs_mod$ParamsInModel$vX0/c
ML_mvs_mod_res[[2]] <- ML_mvs_mod$ParamSummary$StS/c^2
ML_mvs_mod_res[[3]] <- llik_byhand(Tree, ML_mvs_mod_res[[1]], 
                                   ML_mvs_mod_res[[2]], mvData)



# ================================ fit mvMORPH ================================
# pic
ML_pic <- mvBM(tree = Tree, data = mvData, method = "pic")

# rpf
ML_rpf <- mvBM(tree = Tree, data = mvData, method = "rpf")



# ========================= create a table of results ==========================
rnames <- c("mvS", "mvS_mbl", "poly", "scaled", "pic", "rpf")
cnames <- c("X0_1", "X0_2", "StS_11", "StS_12", "StS_22", "loglik")

Mres <- matrix(nrow = length(rnames), ncol = length(cnames), 
               dimnames = list(rnames, cnames))

# retrieve results for mvSLOUCH
mvs_res <- function(obj){
  c(obj$ParamsInModel$vX0, 
    c(obj$ParamSummary$StS[c(1,2,4)]),
    obj$ParamSummary$LogLik)
}

# retrieve results for mvMORPH
mvm_res <- function(obj){
  c(obj$theta, c(obj$sigma)[c(1,2,4)], obj$LogLik)
}


# assign results row by row
# mvSLOUCH (default)
Mres[1,] <- mvs_res(ML_mvs)

# mvSLOUCH with min_bl option
Mres[2,] <- mvs_res(ML_mvs_mbl)

# mvSLOUCH on a tree with polytomy
Mres[3,] <- mvs_res(ML_mvs_pol)

# mvSLOUCH where the data is multiplied by a constant
Mres[4,] <- unlist(ML_mvs_mod_res)[-4]

# mvMORPH (pic)
Mres[5,] <- mvm_res(ML_pic)

# mvMORPH (rpf)
Mres[6,] <- mvm_res(ML_rpf)

# print results
round(Mres, 6)

