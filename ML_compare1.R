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
