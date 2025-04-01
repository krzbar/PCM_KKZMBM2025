# ----- compare the ML results from different new trait assignments ------------
# ------ on the pruned tree ----------------------------------------------------
# experiment parameters:
# - Number of tips (n): 5, 25, 50, 100
# - number of simulations: 100
# - Trait assignment: average, left, right

library(ape)
library(TreeSim)
library(phytools)
library(mvSLOUCH)
library(mvMORPH)

# random seed 
random_seeds <- 1234 + 0:4

# ========================= generate a list of trees =========================
n_tips <- c(5, 25, 50, 100)
n_sim <- 100
l <- 1e-5 # close to the threshold of mvSLOUCH

set.seed(random_seeds[1])
TREES <- lapply(n_tips, function(n){
  TreeSim::sim.bd.taxa(n, numbsim = n_sim, lambda = 1, mu = 0)})


# =========== convert trees into trees with short tips =========================
short_tip_tree <- function(tree, l = 1e-5){
  # normalize the height of the tree into unit height
  tree$edge.length <- tree$edge.length/max(phytools::nodeHeights(tree))
  
  # find which tip branch is the shortest
  tip_id <- which(tree$edge[,2] %in% 1:Ntip(tree))
  edge_id <- tip_id[which.min(tree$edge.length[tip_id])]
  parent_node <- tree$edge[edge_id,1]
  ids_desc <- which(tree$edge[,1] == parent_node)
  id_anc <- which(tree$edge[,2] == parent_node)
  
  # change the values: shorten tip branches, lengthen parent branch
  l_orig <- tree$edge.length[ids_desc][1]
  tree$edge.length[ids_desc] <- l
  tree$edge.length[id_anc] <- tree$edge.length[id_anc] + l_orig - l
  
  return(tree)
}

# ------------ trees with short tip branch -----> use this trees
TREES_st <- list()
for (i in 1:length(n_tips)){
  TREES_st[[i]] <- lapply(TREES[[i]], short_tip_tree)
}


# -------------- retrieve labels of the short tips
get_st_labs <- function(tree, l = 1e-5){
  st_id <- tree$edge[which(tree$edge.length == l),2]
  st_labels <- tree$tip.label[st_id]
  
  return(setNames(c(st_labels, st_id),
                  c("tip1", "tip2", "id1", "id2")))
}

# retrieve labels
st_labs <- list()
for (i in 1:length(n_tips)){
  st_labs[[i]] <- t(sapply(TREES_st[[i]], get_st_labs))
}



# =========== prune one of the two short tip branches in the tree ==============
# run ML fit on pruned trees
# prune the short tip branch
prTREES <- list()

# ---> the lower index of the two tips is removed
for (i in 1:length(n_tips)){
  prTREES[[i]] <- lapply(1:n_sim, function(j){drop.tip(TREES_st[[i]][[j]], 
                                                       as.numeric(st_labs[[i]][j,3]))})
}




# ============= simulate data for each tree (with short tips) ==================
# true parameters:
X0 <- matrix(c(0, 0), ncol = 1) # ancestral trait value at the root
Sigma <- matrix(c(1,2,2,1), nrow = 2) # BM diffusion
StS <- Sigma %*% t(Sigma)
truevalues <- c(X0, c(StS)[c(1,2,4)])

# simulate data at the tips
X_SIM <- list()
for (i in 1:length(n_tips)){
  set.seed(random_seeds[i+1])
  X_SIM[[i]] <- lapply(TREES_st[[i]], 
                       function(tree){simulBMProcPhylTree(phyltree = tree, 
                                                          X0 = X0, 
                                                          Sigma = Sigma)})
}


# ------------ modify the simulated data to fit the pruned tree
mod_X <- function(X, tip_labs, method = c("avg", "left", "right")){
  traits_tips <- X[tip_labs[1:2],]
  
  if (method == "avg"){
    new_traits <- colMeans(traits_tips)
  } else if (method == "left"){
    new_traits <- traits_tips[1,]
  } else if (method == "right"){
    new_traits <- traits_tips[2,]
  }
  
  # remove the trait of left tip in X, 
  # replace by new_traits
  X_new <- X[-which(rownames(X) == tip_labs[1]),]
  X_new[tip_labs[2],] <- new_traits
  
  return(X_new)
}


# case 1: average of the two tips
X_SIM_avg <- list() 

# case 2: left of the two tips
X_SIM_left <- list() 

# case 3: right of the two tips
X_SIM_right <- list() 


for (i in 1:length(n_tips)){
  X_SIM_avg[[i]] <- lapply(1:n_sim, 
                           function(j){mod_X(X = X_SIM[[i]][[j]], 
                                             tip_labs = st_labs[[i]][j,],
                                             method = "avg")}) 
  
  X_SIM_left[[i]] <- lapply(1:n_sim, 
                            function(j){mod_X(X = X_SIM[[i]][[j]], 
                                              tip_labs = st_labs[[i]][j,],
                                              method = "left")}) 
  
  X_SIM_right[[i]] <- lapply(1:n_sim, 
                             function(j){mod_X(X = X_SIM[[i]][[j]], 
                                               tip_labs = st_labs[[i]][j,],
                                               method = "right")}) 
}



# ===== find ML parameters for the pruned trees based on different cases =======

# function to run ML fit given a tree and simulated data
# returns a vector of X0's, StS values, and max log-lik
MLparams <- function(tree, Xsim){
  res <- mvSLOUCH::BrownianMotionModel(phyltree = tree, mData = Xsim)
  return(c(c(res$ParamsInModel$vX0), c(res$ParamSummary$StS)[c(1,2,4)],
           res$ParamSummary$LogLik))
}


# list for the pruned tree
MLavg <- list()
MLleft <- list()
MLright <- list()

for (i in 1:length(n_tips)){
  # average of the two short branch tips
  MLavg[[i]] <- t(sapply(1:n_sim, function(j){MLparams(prTREES[[i]][[j]], 
                                                       X_SIM_avg[[i]][[j]])}))
  
  # left of the two short branch tips
  MLleft[[i]] <- t(sapply(1:n_sim, function(j){MLparams(prTREES[[i]][[j]], 
                                                        X_SIM_left[[i]][[j]])}))
  
  # right of the two short branch tips
  MLright[[i]] <- t(sapply(1:n_sim, function(j){MLparams(prTREES[[i]][[j]], 
                                                         X_SIM_right[[i]][[j]])}))
  
  # print the progress
  print(paste0("n_tips: ", n_tips[i]))
}


# save the ML results
# save(MLavg, file = "MLavg.RData")
# save(MLleft, file = "MLleft.RData")
# save(MLright, file = "MLright.RData")



# ============================ plot the results ================================

xlims <- rbind(c(-5,5),
               c(-5,5),
               c(0,10),
               c(-10,0),
               c(0,10))

cols <- c(rgb(1,0.1,0.1,0.4),
          rgb(0.1,1,0.1,0.4),
          rgb(0.1,0.1,1,0.4))

mains <- c(expression(hat(X[0])^{(1)}),
           expression(hat(X[0])^{(2)}),
           expression(paste(hat(Sigma), "[1,1]")),
           expression(paste(hat(Sigma), "[1,2]")),
           expression(paste(hat(Sigma), "[2,2]")))

Mplots <- rbind(matrix(1:20, nrow = 4), rep(21, 5))
Mrow <- c(0, max(Mplots) + 1:4, 0)
Mcol <- max(Mrow) + 1:5


{

pdf(file = "removed_tip.pdf", width = 16, height = 10)
layout(cbind(Mrow, rbind(Mcol, Mplots)), 
       heights = c(1,rep(5,4),1), widths = c(2, rep(5,5)))

# -------------------------- boxplots --------------------------
par(mar = c(1.1, 2.1, 1.1, 2.1))
for (i in 1:length(truevalues)){
  for (j in 1:length(n_tips)){
    df1 <- cbind(MLavg[[j]][,i], MLleft[[j]][,i], MLright[[j]][,i])
    boxplot(df1, col = rgb(0,0,0,0), border = NA, xaxt = "n",
            cex.axis = 1.5)
    grid()
    boxplot(df1, col = c("dodgerblue4", "skyblue2", "lightcyan"), add = TRUE,
            boxwex = 0.5,  xaxt = "n", yaxt = "n")
    abline(h = truevalues[i], lwd = 3, col = "red", lty = 2)
  }
}

# -------------------------- legend --------------------------

par(mar = c(0,0,0,0))
plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
legend("center", legend = c("true value", "average", "left tip", "right tip"),
       lty = c(2,rep(NA, 3)), pch = c(NA, rep(15, 3)), col = c("red",rep(NA,3)), 
       lwd = rep(2,4), horiz = TRUE, x.intersp = 1, text.width = 1.5, cex = 2, bty = "n",
       border = c(NA, rep("black", 3)), 
       fill = c(NA, "dodgerblue4", "skyblue2", "lightcyan"))


# -------------------------- row names (n_tips) --------------------------
for (i in 1:length(n_tips)){
  plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
  text(0,0, paste0("n = ", n_tips[i]), font = 2, cex  = 1.5)
}


# -------------------------- row names (param. names) --------------------------
for (i in 1:length(truevalues)){
  plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
  text(0,0, mains[i], cex  = 1.5)
}

dev.off()

}
