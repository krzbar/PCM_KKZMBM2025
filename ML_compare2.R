# ------ Calculate maximum likelihood given trees with short branch tips -------
# two cases:
# - small tree: n = 4
# - large tree: n = 100
library(ape)
library(TreeSim)
library(mvSLOUCH)
library(mvMORPH)

# random seed
random_seeds <- 123 + 0:2

# ========================== create a tree with 4 tips =========================
tree_4 <- read.tree(text = "(y1, ((y2, y3), y4));") 
tree_4$edge.length <- c(1, 0.5, 0.25, 0.25, 0.25, 0.5)

# ===================== simulate a tree with 100 tips ==========================
set.seed(random_seeds[1])
tree_100 <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]] 
# normalize
tree_100$edge.length <- tree_100$edge.length/max(nodeHeights(tree_100))


# ========= create list of trees with decreasing short branch lengths ==========

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

# ----------- create the sequence of short branch lengths
# length of the sequence
n_seq <- 100

# for tree_4
l_seq_4 <- sapply(1:n_seq, function(n){0.25*((2/3)^(n-1))})

# for tree_100, subset of l_seq_4 which is lower than l_stip
l_stip <- tree_100$edge.length[get_short_tip(tree_100)[2]]
l_seq_100 <- c(l_stip, l_seq_4[l_seq_4 < l_stip])


# ----------- create list of trees

TREES_4 <- list()
TREES_100 <- list()

# for n = 4
tree_new <- tree_4
for (i in 1:length(l_seq_4)){
  tree_new$edge.length[c(4,5)] <- l_seq_4[i]
  tree_new$edge.length[3] <- 0.5 - l_seq_4[i]
  TREES_4[[i]] <- tree_new
}

# for n = 100
tree_new <- tree_100
ids_edges <- get_short_tip(tree_new)
l_anc_edge <- tree_new$edge.length[ids_edges[1]]
for (i in 1:length(l_seq_100)){
  tree_new$edge.length[ids_edges[2:3]] <- l_seq_100[i]
  tree_new$edge.length[ids_edges[1]] <- l_anc_edge - l_seq_100[i]
  TREES_100[[i]] <- tree_new
}



# =================== simulate traits at tips for each tree ====================
# --------------------- define the true values ---------------------
X0 <- matrix(c(0, 0), ncol = 1) # ancestral trait value at the root
Sigma <- matrix(c(1,2,2,1), nrow = 2) # BM fluctuations
StS <- Sigma %*% t(Sigma)

# vector of true parameter values
truevalues <- c(X0, c(StS)[c(1,2,4)])

# --------------------- simulate traits at tips ---------------------
# prepare matrix of results
Xsim_4 <- list()
Xsim_100 <- list()

set.seed(random_seeds[2])

# simulate for n = 4
for (i in 1:length(l_seq_4)){
  X4 <- list()
  for (j in 1:10){
    X4[[j]] <- mvSLOUCH::simulBMProcPhylTree(phyltree = TREES_4[[i]], 
                                             X0 = X0, Sigma = Sigma)
  }
  Xsim_4[[i]] <- X4
}

for (i in 1:length(l_seq_100)){
  X100 <- list()
  for (j in 1:10){
    X100[[j]] <- mvSLOUCH::simulBMProcPhylTree(phyltree = TREES_100[[i]], 
                                               X0 = X0, Sigma = Sigma)
  }
  Xsim_100[[i]] <- X100
}


# ================= ML functions given tree and simulated data =================
# retrieve maximum likelihood parameters given
# branch length and simulated data
ml_params <- function(tree, Xsim, method = c("mvSLOUCH", "rpf", "pic")){
  if (length(method) == 2){
    method <- "mvSLOUCH"
    cat("mvSLOUCH is used \n")
  }
  
  # return a vector: c(x01, x02, StS11, StS12, StS22)
  if (method == "mvSLOUCH"){
    res <- mvSLOUCH::BrownianMotionModel(phyltree = tree, mData = Xsim)
    if (is.infinite(res$ParamSummary$LogLik)){
      return(rep(NA, 6))
    }else{
      return(c(c(res$ParamsInModel$vX0), c(res$ParamSummary$StS)[c(1,2,4)],
               res$ParamSummary$LogLik))
    }
    
  }else{
    if (method %in% c("rpf", "pic")){
      res <- try(mvMORPH::mvBM(tree = tree, data = Xsim, method = method, 
                               echo = FALSE), silent = TRUE)
      if ("try-error" %in% class(res)){
        return(rep(NA, 6))
      }else{
        return(c(res$theta, res$sigma[c(1,2,4)], res$LogLik))
      }
    }
  }
}


# ================== run ML inference for tree with n = 4 ======================

# mvSLOUCH
mlparams_mvs_4 <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[i]][[k]],
              method = "mvSLOUCH")
  })
})

# rpf (mvMORPH)
mlparams_rpf_4 <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[i]][[k]],
              method = "rpf")
  })
})

# pic (mvMORPH)
mlparams_pic_4 <- lapply(1:length(l_seq_4), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_4[[i]], Xsim = Xsim_4[[i]][[k]],
              method = "pic")
  })
})

# save results
save(mlparams_mvs_4, file = "mlparams_mvs_4.RData")
save(mlparams_rpf_4, file = "mlparams_rpf_4.RData")
save(mlparams_pic_4, file = "mlparams_pic_4.RData")


# ================= run ML inference for tree with n = 100 =====================

# mvSLOUCH
mlparams_mvs_100 <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[i]][[k]],
              method = "mvSLOUCH")
  })
})

# rpf (mvMORPH)
mlparams_rpf_100 <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[i]][[k]],
              method = "rpf")
  })
})

# pic (mvMORPH)
mlparams_pic_100 <- lapply(1:length(l_seq_100), function(i){
  sapply(1:10, function(k){
    ml_params(tree = TREES_100[[i]], Xsim = Xsim_100[[i]][[k]],
              method = "pic")
  })
})

# save results
save(mlparams_mvs_100, file = "mlparams_mvs_100.RData")
save(mlparams_rpf_100, file = "mlparams_rpf_100.RData")
save(mlparams_pic_100, file = "mlparams_pic_100.RData")



# ============================== process results ===============================

# load results
load(file = "mlparams_mvs_4.RData")
load(file = "mlparams_rpf_4.RData")
load(file = "mlparams_pic_4.RData")

load(file = "mlparams_mvs_100.RData")
load(file = "mlparams_rpf_100.RData")
load(file = "mlparams_pic_100.RData")


# n = 4
res_mlparams_mvs_4 <- sapply(mlparams_mvs_4, rowMeans)
res_mlparams_rpf_4 <- sapply(mlparams_rpf_4, rowMeans)
res_mlparams_pic_4 <- sapply(mlparams_pic_4, rowMeans)

# combine into one variable
# combine the mean into one variable for each n
res_4 <- list(res_mlparams_mvs_4, 
              res_mlparams_rpf_4, 
              res_mlparams_pic_4)

# n = 100
res_mlparams_mvs_100 <- sapply(mlparams_mvs_100, rowMeans)
res_mlparams_rpf_100 <- sapply(mlparams_rpf_100, rowMeans)
res_mlparams_pic_100 <- sapply(mlparams_pic_100, rowMeans)

res_100 <- list(res_mlparams_mvs_100, 
                res_mlparams_rpf_100, 
                res_mlparams_pic_100)



# ========================== start plotting routine ============================

# prepare common variables and functions

# colors
cols <- sapply(c("skyblue1","midnightblue", "orange1"), col2rgb)/255
col_mvs <- rgb(cols[1,3], cols[2,3], cols[3,3], 0.7)
col_rpf <- rgb(cols[1,2], cols[2,2], cols[3,2], 0.7)
col_pic <- rgb(cols[1,1], cols[2,1], cols[3,1], 0.7)

# title for the subplots
mains <- c(expression(hat(X[0])^{(1)}),
           expression(hat(X[0])^{(2)}),
           expression(paste(hat(Sigma), "[1,1]")),
           expression(paste(hat(Sigma), "[1,2]")),
           expression(paste(hat(Sigma), "[2,2]")),
           "max. log-lik.")

# prepare limits
xid <- c(1, seq(20, 100, by = 20))
i_zoomed <- 15
xid_zoomed <- c(1, seq(5, i_zoomed, by = 5))


# function to plot one column of results
plot_column <- function(res, ylims, full){
  
  n <- dim(res[[1]])[2]
  
  l_seq <- if(n==100) l_seq_4 else l_seq_100
  ids <- if(full) xid else xid_zoomed
  max_id <- ifelse(full, length(l_seq), max(ids))
  
  for (i in 1:6){
    plot(res[[3]][i,1:max_id], pch = 16, col = rgb(0,0,0,0), axes = FALSE, 
         xlab = "", ylab = "", type = "o", ylim = ylims[i,])
    grid()
    abline(h = truevalues[i], lty = 2, lwd = 2)
    lines(res[[3]][i,1:max_id], pch = 16, col = col_pic, type = "o", cex = 2,
          lwd = 2)
    lines(res[[2]][i,1:max_id], pch = 5, col = col_rpf, type = "o", cex = 2,
          lty = 2, lwd = 2)
    lines(res[[1]][i,1:max_id], pch = 17, col = col_mvs, type = "o", cex = 2)
    
    axis(1, at = ids, labels = formatC(l_seq[ids], format = "e", digits = 2), 
         cex.axis = 2, line = 1.5)
    axis(2, cex.axis = 2, line = 1.5)
    mtext(bquote("\U2113"[i]), side = 1, line = 5, cex = 2)
  }
}

# prepare the layout matrix
M4 <- matrix(1:(2*6), ncol = 2)
M100 <- M4 + max(M4)
Mplots <- rbind(M4, c(0,0), M100)
Mlegend <- rep(max(Mplots)+1, 3)
Mtitle4 <- rep(max(Mlegend)+1, 2)
Mtitle100 <- rep(max(Mlegend)+2, 2)
Mrow <- c(0, max(Mtitle100)+(1:6), c(0,0), max(Mtitle100)+(1:6)+6)
M <- rbind(cbind(Mrow, rbind(Mtitle4, M4, c(0,0), Mtitle100, M100)), Mlegend)

heights <- c(1, rep(5, 6), 2, 1, rep(5, 6), 3)
widths <- c(2, 5, 5)




{
  svg(filename = "sc2.svg", width = 30, height = 40)
  
  # -------------------------- results plots --------------------------
  layout(mat = M, heights = heights, widths = widths)
  par(mar = c(5.1, 4.1, 3.1, 4.1))
  
  # -------------------------- n = 4 --------------------------
  ylims_part <- matrix(c(-2,2,
                         -2,2,
                         0,10,
                         0,10,
                         0,10,
                         -15, 0), 
                       ncol = 2, byrow = TRUE)
  ylims_full <- matrix(c(-2,2,
                         -2,2,
                         0,10,
                         0,10,
                         0,10,
                         -15, 30), 
                       ncol = 2, byrow = TRUE)
  
  plot_column(res = res_4, ylims = ylims_part, full = FALSE)
  plot_column(res = res_4, ylims = ylims_full, full = TRUE)
  
  
  # -------------------------- n = 100 --------------------------
  ylims_part <- matrix(c(-2,2,
                         -2,2,
                         0,10,
                         0,10,
                         0,10,
                         -250, -230), 
                       ncol = 2, byrow = TRUE)
  ylims_full <- matrix(c(-2,2,
                         -2,2,
                         0,10,
                         0,10,
                         0,10,
                         -250, -200), 
                       ncol = 2, byrow = TRUE)
  
  plot_column(res = res_100, ylims = ylims_part, full = FALSE)
  plot_column(res = res_100, ylims = ylims_full, full = TRUE)
  
  
  # -------------------------- legend --------------------------
  par(mar = c(0,0,0,0))
  plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
  legend("center", legend = c("true value", "mvSLOUCH", "mvMORPH (rpf)", "mvMORPH (pic)"),
         lty = c(2,1,2,1), pch = c(NA, 17, 5, 16), 
         col = c("black", "orange1", "midnightblue", "skyblue1"), lwd = c(2,1,1,1),
         horiz = TRUE, x.intersp = 3, text.width = 1.5, cex = 3, bty = "n")
  
  
  # -------------------------- titles --------------------------
  plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
  text(0,0, paste0("n = 4"), font = 2, cex  = 5)
  
  plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
  text(0,0, paste0("n = 100"), font = 2, cex  = 5)
  
  
  # -------------------------- parameter names --------------------------
  for (k in 1:2){
    for (i in 1:6){
      plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
      text(0,0, mains[i], font = 2, cex  = 4)
    }
  }
  
  dev.off()
}
