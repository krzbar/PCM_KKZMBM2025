library(ape)
library(TreeSim)
library(mvSLOUCH)
library(mvMORPH)

# random seed
random_seeds <- 123 + 0:2

# ==================== create a tree with 4 and 100 tips =======================
tree_4 <- read.tree(text = "(y1, ((y2, y3), y4));") 
tree_4$edge.length <- c(1, 0.5, 0.25, 0.25, 0.25, 0.5)

# ============= create a tree with 100 tips and simulate data ==================
set.seed(random_seeds[1])
tree_100 <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]] 
# normalize
tree_100$edge.length <- tree_100$edge.length/max(nodeHeights(tree_100))

# ====================== simulate traits at tips, 10 times =====================

# true values
X0 <- matrix(c(0, 0), ncol = 1) # ancestral trait value at the root
Sigma <- matrix(c(1,2,2,1), nrow = 2)
StS <- Sigma %*% t(Sigma)

# vector of true parameter values
truevalues <- c(X0, c(StS)[c(1,2,4)])

# simulate, n = 4
Xsim_4 <- list()
set.seed(random_seeds[2])
for (i in 1:10){
  Xsim_4[[i]] <- mvSLOUCH::simulBMProcPhylTree(phyltree = tree_4, X0 = X0, 
                                               Sigma = Sigma)
}

# simulate, n = 100
set.seed(random_seeds[3])
Xsim_100 <- list()
for (i in 1:10){
  Xsim_100[[i]] <- mvSLOUCH::simulBMProcPhylTree(phyltree = tree_100, X0 = X0, 
                                                 Sigma = Sigma)
}


# ======================= get the short tips branches ids ======================
# find the shortest tip branches in the tree
get_short_tip <- function(tree){
  # find which tip branch is the shortest
  tip_id <- which(tree$edge[,2] %in% 1:Ntip(tree))
  edge_id <- tip_id[which.min(tree$edge.length[tip_id])]
  parent_node <- tree$edge[edge_id,1]
  ids_desc <- which(tree$edge[,1] == parent_node) # ids of tip branch
  id_anc <- which(tree$edge[,2] == parent_node) # id of ancestor
  return(c(id_anc, ids_desc))
}

# n = 4
st_4 <- get_short_tip(tree_4)

# n = 5
st_100 <- get_short_tip(tree_100)


# ======================= sequence of tip branch length ========================
n_seq <- 100

# for tree_4
l_seq_4 <- sapply(1:n_seq, function(n){0.25*((2/3)^(n-1))})

# for tree_100, subset of l_seq_4 which is lower than l_stip
l_stip <- tree_100$edge.length[get_short_tip(tree_100)[2]]
l_seq_100 <- c(l_stip, l_seq_4[l_seq_4 < l_stip])



# =========================== create a list of trees ===========================

TREES_4 <- list()
TREES_100 <- list()

# for n = 4
tree_new <- tree_4
for (i in 1:length(l_seq_4)){
  tree_new$edge.length[st_4[2:3]] <- l_seq_4[i]
  tree_new$edge.length[st_4[1]] <- 0.5 - l_seq_4[i]
  TREES_4[[i]] <- tree_new
}

# for n = 100
tree_new <- tree_100
l_anc_edge <- tree_new$edge.length[st_100[1]]
for (i in 1:length(l_seq_100)){
  tree_new$edge.length[st_100[2:3]] <- l_seq_100[i]
  tree_new$edge.length[st_100[1]] <- l_anc_edge - l_seq_100[i]
  TREES_100[[i]] <- tree_new
}




# ============ ML function without and with known measurement error ============

ml_params <- function(tree, Xsim, method, M.error = NULL){
  # return a vector: c(x01, x02, StS11, StS12, StS22)
  
  ntips <- Ntip(tree)
  if (ntips == 4){
    st <- st_4
  }else{
    st <- st_100
  }
  
  if (!is.null(M.error)){
    if (M.error == "single"){
      Merror <- sapply(1:ntips, function(x){matrix(0,2,2)}, simplify=FALSE)
      Merror[[tree$edge[st[3],2]]] <- tree$edge.length[st[1]]*StS
      
      merror_matrix <- matrix(0, nrow = ntips, ncol = 2)
      merror_matrix[tree$edge[st[3],2],] <- diag(tree$edge.length[st[1]]*StS)
    }
    
    if (M.error == "both"){
      Merror <- sapply(1:ntips, function(x){matrix(0,2,2)}, simplify=FALSE)
      Merror[[tree$edge[st[2],2]]] <- tree$edge.length[st[1]]*StS*0.5
      Merror[[tree$edge[st[3],2]]] <- tree$edge.length[st[1]]*StS*0.5
      
      merror_matrix <- matrix(0, nrow = ntips, ncol = 2)
      merror_matrix[tree$edge[st[2:3],2],] <- diag(tree$edge.length[st[1]]*StS*0.5)
    }
  }else{
    Merror <- M.error
    merror_matrix <- M.error
  }
  
  if (method == "mvSLOUCH"){
    res <- mvSLOUCH::BrownianMotionModel(phyltree = tree, mData = Xsim,
                                         M.error = Merror)
    if (is.infinite(res$ParamSummary$LogLik) | (res$ParamSummary$LogLik <= -1e5)){
      return(rep(NA, 6))
    }else{
      return(c(c(res$ParamsInModel$vX0), c(res$ParamSummary$StS)[c(1,2,4)],
               res$ParamSummary$LogLik))
    }
    
  }else{
    res <- try(mvMORPH::mvBM(tree = tree, data = Xsim, method = method, 
                             echo = FALSE, diagnostic = FALSE,
                             error = merror_matrix), 
               silent = TRUE)
    if ("try-error" %in% class(res)){
      return(rep(NA, 6))
    }else{
      return(c(res$theta, res$sigma[c(1,2,4)], res$LogLik))
    }
  }
}



# ------------------------ run maximum likelihood fits -------------------------

# ----------------------- no measurement errors -----------------------
# --------------------- n = 4 --------------------- 
# for measuring time
times_4_1 <- c()

start <- Sys.time()
ml1_mvs_4 <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "mvSLOUCH")})
  }
)
times_4_1[1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml1_mvs_4, file = "ml1_mvs_4.RData")
rm(ml1_mvs_4)

start <- Sys.time()
ml1_rpf_4 <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "rpf")})
}
)
times_4_1[2] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml1_rpf_4, file = "ml1_rpf_4.RData")
rm(ml1_rpf_4)

start <- Sys.time()
ml1_pic_4 <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "pic")})
}
)
times_4_1[3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml1_pic_4, file = "ml1_pic_4.RData")
rm(ml1_pic_4)

save(times_4_1, file = "times_4_1.RData")

# --------------------- n = 100 ---------------------
# for measuring time
times_100_1 <- c()

start <- Sys.time()
ml1_mvs_100 <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "mvSLOUCH")})
}
)
times_100_1[1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml1_mvs_100, file = "ml1_mvs_100.RData")
rm(ml1_mvs_100)

start <- Sys.time()
ml1_rpf_100 <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "rpf")})
}
)
times_100_1[2] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml1_rpf_100, file = "ml1_rpf_100.RData")
rm(ml1_rpf_100)

start <- Sys.time()
ml1_pic_100 <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "pic")})
}
)
times_100_1[3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml1_pic_100, file = "ml1_pic_100.RData")
rm(ml1_pic_100)

save(times_100_1, file = "times_100_1.RData")



# ------ measurement errors prop to StS and lengthened branch, single ----------

# --------------------- n = 4 --------------------- 
# for measuring time
times_4_2_single <- c()

start <- Sys.time()
ml2_mvs_4_single <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "mvSLOUCH",
              M.error = "single")})
}
)

times_4_2_single[1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_mvs_4_single, file = "ml2_mvs_4_single.RData")
rm(ml2_mvs_4_single)

start <- Sys.time()
ml2_rpf_4_single <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "rpf",
              M.error = "single")})
}
)

times_4_2_single[2] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_rpf_4_single, file = "ml2_rpf_4_single.RData")
rm(ml2_rpf_4_single)

start <- Sys.time()
ml2_pic_4_single <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "pic",
              M.error = "single")})
}
)

times_4_2_single[3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_pic_4_single, file = "ml2_pic_4_single.RData")
rm(ml2_pic_4_single)

save(times_4_2_single, file = "times_4_2_single.RData")

# --------------------- n = 100 ---------------------
# for measuring time
times_100_2_single <- c()

start <- Sys.time()
ml2_mvs_100_single <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "mvSLOUCH",
              M.error = "single")})
}
)

times_100_2_single[1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_mvs_100_single, file = "ml2_mvs_100_single.RData")
rm(ml2_mvs_100_single)

start <- Sys.time()
ml2_rpf_100_single <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "rpf",
              M.error = "single")})
}
)

times_100_2_single[2] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_rpf_100_single, file = "ml2_rpf_100_single.RData")
rm(ml2_rpf_100_single)

start <- Sys.time()
ml2_pic_100_single <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "pic",
              M.error = "single")})
}
)

times_100_2_single[3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_pic_100_single, file = "ml2_pic_100_single.RData")
rm(ml2_pic_100_single)

save(times_100_2_single, file = "times_100_2_single.RData")




# -------- measurement errors prop to StS and lengthened branch, both ----------
# --------------------- n = 4 --------------------- 
# for measuring time
times_4_2_both <- c()

start <- Sys.time()
ml2_mvs_4_both <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "mvSLOUCH",
              M.error = "both")})
}
)

times_4_2_both[1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_mvs_4_both, file = "ml2_mvs_4_both.RData")
rm(ml2_mvs_4_both)

start <- Sys.time()
ml2_rpf_4_both <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "rpf",
              M.error = "both")})
}
)

times_4_2_both[2] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_rpf_4_both, file = "ml2_rpf_4_both.RData")
rm(ml2_rpf_4_both)

start <- Sys.time()
ml2_pic_4_both <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[k]], method = "pic",
              M.error = "both")})
}
)

times_4_2_both[3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_pic_4_both, file = "ml2_pic_4_both.RData")
rm(ml2_pic_4_both)

save(times_4_2_both, file = "times_4_2_both.RData")

# --------------------- n = 100 ---------------------
# for measuring time
times_100_2_both <- c()

start <- Sys.time()
ml2_mvs_100_both <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "mvSLOUCH",
              M.error = "both")})
}
)

times_100_2_both[1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_mvs_100_both, file = "ml2_mvs_100_both.RData")
rm(ml2_mvs_100_both)

start <- Sys.time()
ml2_rpf_100_both <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "rpf",
              M.error = "both")})
}
)

times_100_2_both[2] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_rpf_100_both, file = "ml2_rpf_100_both.RData")
rm(ml2_rpf_100_both)

start <- Sys.time()
ml2_pic_100_both <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[k]], method = "pic",
              M.error = "both")})
}
)

times_100_2_both[3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
save(ml2_pic_100_both, file = "ml2_pic_100_both.RData")
rm(ml2_pic_100_both)

save(times_100_2_both, file = "times_100_2_both.RData")





# ================= ML function with unknown measurement error =================


# objective function for mvSLOUCH: given v, return loglik
# fopt_mvs <- function(v, ntips, ids, return_obj = FALSE){
#   # convert v to M.error
#   
#   # create SPD matrix
#   L <- matrix(c(v[1], v[2], 0, v[3]), nrow = 2)
#   M <- L %*% t(L)
# 
#   # create an empty list of error matrices
#   Merror <- sapply(1:ntips, function(x){matrix(0,2,2)}, simplify=FALSE)
# 
#   if (length(ids) == 1){
#     Merror[[ids]] <- M
#   }else{
#     Merror[[ids[1]]] <- M*0.5
#     Merror[[ids[2]]] <- M*0.5
#   }
# 
#   res <- BrownianMotionModel(phyltree = tree, mData = Xsim,
#                              M.error = Merror)
#   
#   if (return_obj){
#     return(res)
#   }
# 
#   if (is.infinite(res$ParamSummary$LogLik)){
#     return(-1e5)
#   }else{
#     return(res$ParamSummary$LogLik)
#   }
# }
# 
# tree <- tree_4
# Xsim <- Xsim_4[[1]]
# res <- fopt_mvs(v = c(0,1,0), ntips = 4, ids = c(2,3))
# res_obj <- fopt_mvs(v = c(0,1,0), ntips = 4, ids = c(2,3), return_obj = TRUE)
# 
# 
# res1 <-  mvBM(tree = tree, data = Xsim, echo = FALSE, diagnostic = FALSE)
# init_rpf <- c(c(t(chol(res1$sigma)))[c(1,2,4)], c(res1$theta), rep(0.1,2))



# # design matrix
# D <- cbind(rep(c(1,0), each = ntips), rep(c(0,1), each = ntips))
# 
# # variance covariance matrix
# V <- kronecker(res1$ParamSummary$StS, vcv(tree))
# 
# mvLL(V, data = Xsim, method = "rpf",
#      param = list(estim = FALSE, D = D, mu = as.vector(pars$vX0),
#                   sigma = res1$ParamSummary$StS),
#      error = merror4)$logl
# 
# 
# 
# 
# 
# 
# # objective function for rpf: given v, return loglik
# fopt_rpf <- function(v, ntips, ids){
#   
#   # Sigma (StS) from v
#   Sxx <- matrix(c(v[1], v[2], 0, v[3]), nrow = 2)
#   Sigma <- Sxx %*% t(Sxx)
#   
#   # mean ancestral
#   vX0 <- v[4:5]
#   
#   # matrix of errors
#   Merror <- matrix(0, nrow = ntips, ncol = 2)
#   
#   if (length(ids) == 1){
#     Merror[ids,] <- c(v[6], v[7])
#   }else{
#     Merror[ids[1],] <- c(v[6], v[7])*0.5
#     Merror[ids[2],] <- c(v[6], v[7])*0.5
#   }
#   
#   # design matrix
#   D <- cbind(rep(c(1,0), each = ntips), rep(c(0,1), each = ntips))
#   
#   # variance covariance matrix
#   V <- kronecker(Sigma, vcv(tree))
#   
#   res <- try(mvLL(V, data = Xsim, method = "rpf",
#               param = list(estim = FALSE, D = D, mu = vX0,
#                            sigma = Sigma),
#               error = Merror)$logl, silent = TRUE)
#   
#   if ("try-error" %in% class(res)){
#     return(-1e6)
#   }else{
#     return(res)
#   }
# }
# 
# 
# lb <- c(0, -100, 0, -100, -100, 0, 0)
# ub <- c(100, 100, 100, 100, 100, 100, 100)
# optim(par = init_rpf, fn = fopt_rpf, method = "L-BFGS-B",
#       lower = lb, upper = ub, control = list(fnscale = -1),
#       ntips = 4, ids = 3)
# 
# fopt_rpf(v=init_rpf+rnorm(7, sd=0.01), ntips=4, ids=3)
# 
# 
# 
# main function
# ml_me <- function(tree, Xsim, method, single = TRUE){
# 
#   # number of tips
#   ntips <- Ntip(tree)
# 
#   # short tip indices
#   st <- if(ntips==4) st_4 else st_100
# 
#   # branches indices to tips indices
#   if (single == TRUE){
#     ids <- tree$edge[st[3],2]
#   }else{
#     ids <- tree$edge[st[2:3],2]
#   }
# 
#   if (method == "mvSLOUCH"){
#     # initial ML estimates with zero measurement errors
#     res1 <- mvSLOUCH::BrownianMotionModel(phyltree = tree, mData = Xsim)
#     
#     if (res1$ParamSummary$LogLik < -1e6){
#       return(rep(NA,9))
#     }
#     
#     # set lower and upper bound of optimization
#     lb <- c(0, -100, 0)
#     ub <- c(100, 100, 100)
# 
#     # run optim
#     opt_res <- optim(par = c(0,0,0), fn = fopt_mvs, method = "L-BFGS-B",
#                      lower = lb, upper = ub, control = list(fnscale = -1),
#                      ntips = ntips, ids = ids)
#     
#     
# 
#     # return optim results if loglik is better
#     if (opt_res$value > res1$ParamSummary$LogLik){
#       res2 <- fopt_mvs(v = opt_res$par, ntips = ntips, ids = ids, 
#                        return_obj = TRUE)
#       return(c(c(res2$ParamsInModel$vX0), c(res2$ParamSummary$StS)[c(1,2,4)],
#                opt_res$par, res2$ParamSummary$LogLik))
#     }else{
#       return(c(c(res1$ParamsInModel$vX0), c(res1$ParamSummary$StS)[c(1,2,4)],
#                rep(0,3), res1$ParamSummary$LogLik))
#     }
#   }
# 
#   if (method == "rpf"){
#     # initial ML estimates
#     res1 <- try(mvBM(tree = tree, data = Xsim,
#                      echo = FALSE, diagnostic = FALSE), silent = TRUE)
# 
#     if ("try-error" %in% class(res1)){
#       return(rep(NA,7))
#     }
# 
#     # if not, extract initial param
#     init_rpf <- c(c(t(chol(res1$sigma)))[c(1,2,4)], c(res1$theta), rep(0.1,2))
# 
#     # set lower and upper bound of optimization
#     lb <- c(0, -100, 0, -100, -100, 0, 0)
#     ub <- c(100, 100, 100, 100, 100, 100, 100)
# 
#     # run optim
#     opt_res <- optim(par = init_rpf, fn = fopt_rpf, method = "L-BFGS-B",
#                      lower = lb, upper = ub, control = list(fnscale = -1),
#                      ntips = ntips, ids = ids)
# 
#     # return optim results if loglik is better
#     if (opt_res$value > res1$LogLik){
#       
#       return(c(opt_res$par, opt_res$value))
#     }else{
#       return(c(c(res1$theta), c(res1$sigma)[c(1,2,4)],
#                rep(0,2), res1$LogLik))
#     }
#   }
# }
# 
# 
# # initial ML optimization
# tree <- TREES_4[[2]]
# Xsim <- Xsim_4[[2]]
#  
# times_4_3_single <- c()
# 
# timestamp()
# start <- Sys.time()
# test_mvs <- ml_me(tree = tree, Xsim = Xsim, method = "mvSLOUCH",
#                   single = FALSE)
# times_4_3_single[1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
# timestamp()
# 
# 
# 
# tree <- TREES_100[[2]]
# Xsim <- Xsim_100[[2]]
# times_100_3_single <- c()
# timestamp()
# start <- Sys.time()
# test_mvs <- ml_me(tree = tree, Xsim = Xsim, method = "mvSLOUCH",
#                   single = FALSE)
# times_100_3_single[1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
# timestamp()


# test_rpf <- ml_me(tree, Xsim, method = "rpf")








