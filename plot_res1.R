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



# ------------------------ Case 1: no measurement errors -----------------------

# load results
load("ml1_mvs_4.RData")
load("ml1_rpf_4.RData")
load("ml1_pic_4.RData")
load("ml1_mvs_100.RData")
load("ml1_rpf_100.RData")
load("ml1_pic_100.RData")

# extract the means from results
res_ml1_mvs_4 <- sapply(ml1_mvs_4, rowMeans)
res_ml1_rpf_4 <- sapply(ml1_rpf_4, rowMeans)
res_ml1_pic_4 <- sapply(ml1_pic_4, rowMeans)
res_ml1_mvs_100 <- sapply(ml1_mvs_100, rowMeans)
res_ml1_rpf_100 <- sapply(ml1_rpf_100, rowMeans)
res_ml1_pic_100 <- sapply(ml1_pic_100, rowMeans)

# combine the mean into one variable for each n
res_4 <- list(res_ml1_mvs_4, 
              res_ml1_rpf_4, 
              res_ml1_pic_4)

res_100 <- list(res_ml1_mvs_100, 
                res_ml1_rpf_100, 
                res_ml1_pic_100)

{
  svg(filename = "sc1_case1.svg", width = 30, height = 40)
  
  # -------------------------- results plots --------------------------
  layout(mat = M, heights = heights, widths = widths)
  par(mar = c(5.1, 4.1, 3.1, 4.1))
  
  # -------------------------- n = 4 --------------------------
  ylims_part <- matrix(c(-2,2,
                      -2,2,
                      0,80,
                      0,80,
                      0,80,
                      -17, -8), 
                    ncol = 2, byrow = TRUE)
  ylims_full <- matrix(c(-2,2,
                        -2,2,
                        0,max(res_4[[3]][3,], na.rm=TRUE),
                        0,max(res_4[[3]][3,], na.rm=TRUE),
                        0,max(res_4[[3]][3,], na.rm=TRUE),
                        min(res_4[[3]][6,], na.rm=TRUE), max(res_4[[3]][6,], na.rm=TRUE)), 
                      ncol = 2, byrow = TRUE)
  
  plot_column(res = res_4, ylims = ylims_part, full = FALSE)
  plot_column(res = res_4, ylims = ylims_full, full = TRUE)
  
  
  # -------------------------- n = 100 --------------------------
  ylims_part <- matrix(c(-2,2,
                         -2,2,
                         0,15,
                         0,15,
                         0,15,
                         range(res_100[[3]][6,1:15])), 
                       ncol = 2, byrow = TRUE)
  ylims_full <- matrix(c(-2,2,
                         -2,2,
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         range(res_100[[3]][6,], na.rm = TRUE)), 
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


rm("ml1_mvs_4")
rm("ml1_rpf_4")
rm("ml1_pic_4")
rm("ml1_mvs_100")
rm("ml1_rpf_100")
rm("ml1_pic_100")

rm("res_ml1_mvs_4")
rm("res_ml1_rpf_4")
rm("res_ml1_pic_4")
rm("res_ml1_mvs_100")
rm("res_ml1_rpf_100")
rm("res_ml1_pic_100")

rm("res_4")
rm("res_100")




# ------------- Case 2: with measurement errors (on single tip) ----------------

# load results
load("ml2_mvs_4_single.RData")
load("ml2_rpf_4_single.RData")
load("ml2_pic_4_single.RData")
load("ml2_mvs_100_single.RData")
load("ml2_rpf_100_single.RData")
load("ml2_pic_100_single.RData")

# extract the means from results
res_ml2_mvs_4_single <- sapply(ml2_mvs_4_single, rowMeans)
res_ml2_rpf_4_single <- sapply(ml2_rpf_4_single, rowMeans)
res_ml2_pic_4_single <- sapply(ml2_pic_4_single, rowMeans)
res_ml2_mvs_100_single <- sapply(ml2_mvs_100_single, rowMeans)
res_ml2_rpf_100_single <- sapply(ml2_rpf_100_single, rowMeans)
res_ml2_pic_100_single <- sapply(ml2_pic_100_single, rowMeans)

# combine the mean into one variable for each n
res_4 <- list(res_ml2_mvs_4_single, 
              res_ml2_rpf_4_single, 
              res_ml2_pic_4_single)

res_100 <- list(res_ml2_mvs_100_single, 
                res_ml2_rpf_100_single, 
                res_ml2_pic_100_single)

{
  svg(filename = "sc1_case2.svg", width = 30, height = 40)
  
  # -------------------------- results plots --------------------------
  layout(mat = M, heights = heights, widths = widths)
  par(mar = c(5.1, 4.1, 3.1, 4.1))
  
  # -------------------------- n = 4 --------------------------
  ylims_part <- matrix(c(-2,2,
                         -2,2,
                         0,80,
                         0,80,
                         0,80,
                         -17, -8), 
                       ncol = 2, byrow = TRUE)
  ylims_full <- matrix(c(-2,2,
                         -2,2,
                         0,max(res_4[[3]][3,], na.rm=TRUE),
                         0,max(res_4[[3]][3,], na.rm=TRUE),
                         0,max(res_4[[3]][3,], na.rm=TRUE),
                         min(res_4[[3]][6,], na.rm=TRUE), max(res_4[[3]][6,], na.rm=TRUE)), 
                       ncol = 2, byrow = TRUE)
  
  plot_column(res = res_4, ylims = ylims_part, full = FALSE)
  plot_column(res = res_4, ylims = ylims_full, full = TRUE)
  
  
  # -------------------------- n = 100 --------------------------
  ylims_part <- matrix(c(-2,2,
                         -2,2,
                         0,15,
                         0,15,
                         0,15,
                         range(res_100[[3]][6,1:15])), 
                       ncol = 2, byrow = TRUE)
  ylims_full <- matrix(c(-2,2,
                         -2,2,
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         range(res_100[[3]][6,], na.rm = TRUE)), 
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


rm("ml2_mvs_4_single")
rm("ml2_rpf_4_single")
rm("ml2_pic_4_single")
rm("ml2_mvs_100_single")
rm("ml2_rpf_100_single")
rm("ml2_pic_100_single")

rm("res_ml2_mvs_4_single")
rm("res_ml2_rpf_4_single")
rm("res_ml2_pic_4_single")
rm("res_ml2_mvs_100_single")
rm("res_ml2_rpf_100_single")
rm("res_ml2_pic_100_single")

rm("res_4")
rm("res_100")






# -------------- Case 2: with measurement errors (on both tips) ----------------

# load results
load("ml2_mvs_4_both.RData")
load("ml2_rpf_4_both.RData")
load("ml2_pic_4_both.RData")
load("ml2_mvs_100_both.RData")
load("ml2_rpf_100_both.RData")
load("ml2_pic_100_both.RData")

# extract the means from results
res_ml2_mvs_4_both <- sapply(ml2_mvs_4_both, rowMeans)
res_ml2_rpf_4_both <- sapply(ml2_rpf_4_both, rowMeans)
res_ml2_pic_4_both <- sapply(ml2_pic_4_both, rowMeans)
res_ml2_mvs_100_both <- sapply(ml2_mvs_100_both, rowMeans)
res_ml2_rpf_100_both <- sapply(ml2_rpf_100_both, rowMeans)
res_ml2_pic_100_both <- sapply(ml2_pic_100_both, rowMeans)

# combine the mean into one variable for each n
res_4 <- list(res_ml2_mvs_4_both, 
              res_ml2_rpf_4_both, 
              res_ml2_pic_4_both)

res_100 <- list(res_ml2_mvs_100_both, 
                res_ml2_rpf_100_both, 
                res_ml2_pic_100_both)

{
  svg(filename = "sc1_case3.svg", width = 30, height = 40)
  
  # -------------------------- results plots --------------------------
  layout(mat = M, heights = heights, widths = widths)
  par(mar = c(5.1, 4.1, 3.1, 4.1))
  
  # -------------------------- n = 4 --------------------------
  ylims_part <- matrix(c(-2,2,
                         -2,2,
                         0,80,
                         0,80,
                         0,80,
                         -20, 5), 
                       ncol = 2, byrow = TRUE)
  ylims_full <- matrix(c(-2,2,
                         -2,2,
                         0,max(res_4[[3]][3,], na.rm=TRUE),
                         0,max(res_4[[3]][3,], na.rm=TRUE),
                         0,max(res_4[[3]][3,], na.rm=TRUE),
                         min(res_4[[2]][6,], na.rm=TRUE), 5), 
                       ncol = 2, byrow = TRUE)
  
  plot_column(res = res_4, ylims = ylims_part, full = FALSE)
  plot_column(res = res_4, ylims = ylims_full, full = TRUE)
  
  
  # -------------------------- n = 100 --------------------------
  ylims_part <- matrix(c(-2,2,
                         -2,2,
                         0,15,
                         0,15,
                         0,15,
                         range(res_100[[3]][6,1:15])), 
                       ncol = 2, byrow = TRUE)
  ylims_full <- matrix(c(-2,2,
                         -2,2,
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         0,max(res_100[[3]][3,], na.rm=TRUE),
                         range(res_100[[3]][6,], na.rm = TRUE)), 
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


rm("ml2_mvs_4_both")
rm("ml2_rpf_4_both")
rm("ml2_pic_4_both")
rm("ml2_mvs_100_both")
rm("ml2_rpf_100_both")
rm("ml2_pic_100_both")

rm("res_ml2_mvs_4_both")
rm("res_ml2_rpf_4_both")
rm("res_ml2_pic_4_both")
rm("res_ml2_mvs_100_both")
rm("res_ml2_rpf_100_both")
rm("res_ml2_pic_100_both")

rm("res_4")
rm("res_100")


