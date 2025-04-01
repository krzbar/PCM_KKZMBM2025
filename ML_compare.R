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

# =============== create a tree with 4 tips and simulate data ==================
tree_4 <- read.tree(text = "(y1, ((y2, y3), y4));") 
tree_4$edge.length <- c(1, 0.5, 0.25, 0.25, 0.25, 0.5)

# ----------- new tree (for plotting)
# make the tip branch shorter
tree_new <- tree_4
l_new <- 0.245
tree_new$edge.length <- c(1, 0.5, 0.25+l_new, 0.25-l_new, 0.25-l_new, 0.5)

pdf("four_tips.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
plot.phylo(tree_4, main = "A", label.offset = 0.1, edge.width = 2)
axis(1, at = seq(0, 1, by = 0.2))
plot.phylo(tree_new,  main = "B", label.offset = 0.1, edge.width = 2)
axis(1, at = seq(0, 1, by = 0.2))
dev.off()

# --------------------- simulate traits at tips
set.seed(random_seeds[1])
X0 <- matrix(c(0, 0), ncol = 1) # ancestral trait value at the root
Sigma <- matrix(c(1,2,2,1), nrow = 2)
StS <- Sigma %*% t(Sigma)

# vector of true parameter values
truevalues <- c(X0, c(StS)[c(1,2,4)])

# simulate
Xsim_4 <- mvSLOUCH::simulBMProcPhylTree(phyltree = tree_4, X0 = X0, Sigma = Sigma)

# ============= create a tree with 100 tips and simulate data ==================
set.seed(random_seeds[2])
tree_100 <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]] 
# normalize
tree_100$edge.length <- tree_100$edge.length/max(nodeHeights(tree_100))

# --------------------- simulate traits at tips
set.seed(random_seeds[3])
# simulate
Xsim_100 <- mvSLOUCH::simulBMProcPhylTree(phyltree = tree_100, X0 = X0, Sigma = Sigma)

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


# ======================= sequence of tip branch length ========================
n_seq <- 100

# for tree_4
l_seq_4 <- sapply(1:n_seq, function(n){0.25*((2/3)^(n-1))})

# for tree_100, subset of l_seq_4 which is lower than l_stip
l_stip <- tree_100$edge.length[get_short_tip(tree_100)[2]]
l_seq_100 <- c(l_stip, l_seq_4[l_seq_4 < l_stip])


# =============== ML functions given tree and branch length  ===================
# retrieve maximum likelihood parameters given
# branch length and simulated data
ml_params <- function(l, tree, Xsim, method = c("mvSLOUCH", "rpf", "pic")){
  if (length(method) == 2){
    method <- "mvSLOUCH"
    cat("mvSLOUCH is used \n")
  }
  
  # number of tips
  ntips <- Ntip(tree)
  
  # change the tip branch length
  if (ntips == 4){
    tree_new <- tree
    tree_new$edge.length[c(4,5)] <- l
    tree_new$edge.length[3] <- 0.5 - l
  }else{
    tree_new <- tree
    ids_edges <- get_short_tip(tree_new)
    l_anc_edge <- tree_new$edge.length[ids_edges[1]]
    tree_new$edge.length[ids_edges[2:3]] <- l
    tree_new$edge.length[ids_edges[1]] <- l_anc_edge - l
  }
  
  
  # return a vector: c(x01, x02, StS11, StS12, StS22)
  if (method == "mvSLOUCH"){
    res <- mvSLOUCH::BrownianMotionModel(phyltree = tree_new, mData = Xsim)
    if (is.infinite(res$ParamSummary$LogLik)){
      return(rep(NA, 6))
    }else{
      return(c(c(res$ParamsInModel$vX0), c(res$ParamSummary$StS)[c(1,2,4)],
               res$ParamSummary$LogLik))
    }
    
  }else{
    if (method %in% c("rpf", "pic")){
      res <- try(mvMORPH::mvBM(tree = tree_new, data = Xsim, method = method, 
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
mlparams_mvs_4 <- sapply(l_seq_4, ml_params, tree = tree_4, Xsim = Xsim_4, 
                         method = "mvSLOUCH")

# rpf (mvMORPH)
mlparams_rpf_4 <- sapply(l_seq_4, ml_params, tree = tree_4, Xsim = Xsim_4, 
                         method = "rpf")

# pic (mvMORPH)
mlparams_pic_4 <- sapply(l_seq_4, ml_params, tree = tree_4, Xsim = Xsim_4, 
                         method = "pic")


# ================= run ML inference for tree with n = 100 =====================

# mvSLOUCH
mlparams_mvs_100 <- sapply(l_seq_100, ml_params, tree = tree_100, Xsim = Xsim_100, 
                           method = "mvSLOUCH")

# rpf (mvMORPH)
mlparams_rpf_100 <- sapply(l_seq_100, ml_params, tree = tree_100, Xsim = Xsim_100, 
                           method = "rpf")

# pic (mvMORPH)
mlparams_pic_100 <- sapply(l_seq_100, ml_params, tree = tree_100, Xsim = Xsim_100, 
                           method = "pic")




# ============================= plot the results ===============================

# prepare the colors
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

# prepare the layout
Mplots_ <- matrix(1:(4*6), nrow = 6, byrow = TRUE)
Mplots <- cbind(Mplots_[,1:2], rep(0, length(truevalues)+1), Mplots_[,3:4]) 
Mlegend <- rep(max(Mplots)+1, 5)
Mcol <- c(rep(max(Mlegend)+1, 2), 0, rep(max(Mlegend)+2, 2))
Mrow <- c(0, (max(Mcol)+1):(max(Mcol)+6), 0)

M <- cbind(Mrow, rbind(Mcol, Mplots, Mlegend))
heights <- c(1, rep(5, 6), 3)
widths <- c(3, 5,5, 1, 5,5)


# ----------> plot in svg to preserve the \ell symbol
# ----------> later converted to pdf
{
# ---------------------------- start plot ---------------------------

svg("ML_compare.svg", width = 20, height = 25) 
 
layout(mat = M, heights = heights, widths = widths)

xid <- c(1, seq(20, 100, by = 20))
i_zoomed <- 15
xid_zoomed <- c(1, seq(5, i_zoomed, by = 5))


ylims_zoomed_4 <- rbind(c(-2, 2),
                        c(-2, 2),
                        c(4, 85),
                        c(-24, 5),
                        c(1, 11),
                        c(-17, -10))

ylims_zoomed_100 <- rbind(c(-2, 2),
                        c(-2, 2),
                        c(0, 10),
                        c(0, 10),
                        c(0, 20),
                        c(-310, -220))


# ---------------------------- plot the results ---------------------------
par(mar = c(4.1, 4.1, 3.1, 2.1))
for (i in 1:6){
  # ------------------------------ plot for n = 4 ------------------------------
  # ---------------------- zoomed in scale ----------------------
  plot(mlparams_pic_4[i,1:i_zoomed], pch = 16, col = rgb(0,0,0,0), axes = FALSE, 
       xlab = "", ylab = "", type = "o", ylim = ylims_zoomed_4[i,])
  grid()
  lines(mlparams_pic_4[i,1:i_zoomed], pch = 16, col = col_pic, type = "o", cex = 1.5,
        lwd = 2)
  lines(mlparams_rpf_4[i,1:i_zoomed], pch = 5, col = col_rpf, type = "o", cex = 1.5,
        lty = 2, lwd = 1.5)
  lines(mlparams_mvs_4[i,1:i_zoomed], pch = 17, col = col_mvs, type = "o", cex = 1.5)
  
  axis(1, at = xid_zoomed, labels = formatC(l_seq_4[xid_zoomed], format = "e", digits = 2), 
       cex.axis = 1.75, line = 1.5)
  axis(2, cex.axis = 2, line = 1.5)
  mtext(bquote("\U2113"[i]), side = 1, line = 5, cex = 1.5)
  abline(h = truevalues[i], lty = 2, lwd = 2)
  
  # ---------------------- full scale ----------------------
  plot(mlparams_pic_4[i,], pch = 16, col = rgb(0,0,0,0), axes = FALSE, 
       xlab = "", ylab = "", type = "o")
  grid()
  lines(mlparams_pic_4[i,], pch = 16, col = col_pic, type = "o", cex = 1.5,
        lwd = 2)
  lines(mlparams_rpf_4[i,], pch = 5, col = col_rpf, type = "o", cex = 1.5,
        lty = 2, lwd = 1.5)
  lines(mlparams_mvs_4[i,], pch = 17, col = col_mvs, type = "o", cex = 1.5)
  
  axis(1, at = xid, labels = formatC(l_seq_4[xid], format = "e", digits = 2), 
       cex.axis = 1.75, line = 1.5)
  axis(2, cex.axis = 2, line = 1.5)
  mtext(bquote("\U2113"[i]), side = 1, line = 5, cex = 1.5)
  
  
  
  # ----------------------------- plot for n = 100 -----------------------------
  # ---------------------- zoomed in scale ----------------------
  plot(mlparams_pic_100[i,1:i_zoomed], pch = 16, col = rgb(0,0,0,0), axes = FALSE, 
       xlab = "", ylab = "", type = "o", ylim = ylims_zoomed_100[i,])
  grid()
  abline(h = truevalues[i], lty = 2, lwd = 2)
  lines(mlparams_pic_100[i,1:i_zoomed], pch = 16, col = col_pic, type = "o", cex = 1.5,
        lwd = 2)
  lines(mlparams_rpf_100[i,1:i_zoomed], pch = 5, col = col_rpf, type = "o", cex = 1.5,
        lty = 2, lwd = 1.5)
  lines(mlparams_mvs_100[i,1:i_zoomed], pch = 17, col = col_mvs, type = "o", cex = 1.5)
  
  axis(1, at = xid_zoomed, labels = formatC(l_seq_100[xid_zoomed], format = "e", digits = 2), 
       cex.axis = 1.75, line = 1.5)
  axis(2, cex.axis = 2, line = 1.5)
  mtext(bquote("\U2113"[i]), side = 1, line = 5, cex = 1.5)
  
  
  # ---------------------- full scale ----------------------
  plot(mlparams_pic_100[i,], pch = 16, col = rgb(0,0,0,0), axes = FALSE, 
       xlab = "", ylab = "", type = "o")
  grid()
  lines(mlparams_pic_100[i,], pch = 16, col = col_pic, type = "o", cex = 1.5,
        lwd = 2)
  lines(mlparams_rpf_100[i,], pch = 5, col = col_rpf, type = "o", cex = 1.5,
        lty = 2, lwd = 1.5)
  lines(mlparams_mvs_100[i,], pch = 17, col = col_mvs, type = "o", cex = 1.5)
  
  axis(1, at = xid, labels = formatC(l_seq_100[xid], format = "e", digits = 2), 
       cex.axis = 1.75, line = 1.5)
  axis(2, cex.axis = 2, line = 1.5)
  mtext(bquote("\U2113"[i]), side = 1, line = 5, cex = 1.5)
  
}

# ---------------------------- plot legend ---------------------------
par(mar = c(0,0,0,0))
plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
legend("center", legend = c("true value", "mvSLOUCH", "mvMORPH (rpf)", "mvMORPH (pic)"),
       lty = c(2,1,2,1), pch = c(NA, 17, 5, 16), 
       col = c("black", "orange1", "midnightblue", "skyblue1"), lwd = c(2,1,1,1),
       horiz = TRUE, x.intersp = 2, text.width = 1.5, cex = 3, bty = "n")

# ---------------------------- plot column names ---------------------------
plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
text(0,0, paste0("n = 4"), font = 2, cex  = 3)

plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
text(0,0, paste0("n = 100"), font = 2, cex  = 3)


# ---------------------------- plot row names ---------------------------
for (i in 1:6){
  plot(NA, xlim = c(-5,5), ylim = c(-5,5), axes = FALSE, xlab = "", ylab = "")
  text(-1,0, mains[i], font = 1, cex  = 3)
}


dev.off()

}


# =================== check ML estimates ===================
# minimum length of l for mvSLOUCH:
l_seq_4[max(which(!is.na(mlparams_mvs_4[6,])))]
l_seq_100[max(which(!is.na(mlparams_mvs_100[6,])))]

# parameter values
mlparams_mvs_4[,max(which(!is.na(mlparams_mvs_4[6,])))]
mlparams_mvs_100[,max(which(!is.na(mlparams_mvs_100[6,])))]

# # minimum length of l for rpf:
l_seq_4[tail(which(!is.na(mlparams_rpf_4[6,])), 2)]
mlparams_rpf_4[,max(which(!is.na(mlparams_rpf_4[6,])))]

# # minimum length of l for pic:
l_seq_100[max(which(!is.na(mlparams_pic_100[6,])))]
mlparams_pic_4[,100]
mlparams_pic_100[,dim(mlparams_pic_100)[2]]





