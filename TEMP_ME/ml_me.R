## These are the R scripts and numerical results accompanying 
## Bartoszek, Brahmantio, Munoz-Duran, Fuentes-Gonzalez,Pienaar, and Polly
## Short branch singularities in phylogenetic comparative methods

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


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

# pdf("four_tips.pdf", width = 8, height = 4)
# par(mfrow = c(1,2))
# plot.phylo(tree_4, main = "A", label.offset = 0.1, edge.width = 2)
# axis(1, at = seq(0, 1, by = 0.2))
# plot.phylo(tree_new,  main = "B", label.offset = 0.1, edge.width = 2)
# axis(1, at = seq(0, 1, by = 0.2))
# dev.off()



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


# --------------------- simulate traits at tips
set.seed(random_seeds[3])
# simulate
Xsim_100 <- mvSLOUCH::simulBMProcPhylTree(phyltree = tree_100, X0 = X0, Sigma = Sigma)




# ======================= sequence of tip branch length ========================
n_seq <- 100

# for tree_4
l_seq_4 <- sapply(1:n_seq, function(n){0.25*((2/3)^(n-1))})

# for tree_100, subset of l_seq_4 which is lower than l_stip
l_stip <- tree_100$edge.length[get_short_tip(tree_100)[2]]
l_seq_100 <- c(l_stip, l_seq_4[l_seq_4 < l_stip])



# helper function, convert vector of length 3 into a measurement errors list
v2Me <- function(v, ntips, id){
  L <- matrix(c(v[1], v[2], 0, v[3]), ncol = 2, nrow = 2)
  Ve <- L %*% t(L)
  Me <- sapply(1:ntips,function(x){matrix(0,2,2)}, simplify=FALSE)
  Me[[2]] <- Ve
  return(Me)
}

# objective function
fopt <- function(v){
  # create a list of params from the first 6 indices
  pars <- list(Sxx = matrix(v[1:4], 2, 2), vX0 = matrix(v[5:6], 2, 1))
  
  # create a list of measurement errors from the last 3 indices
  Me <- v2Me(v[7:9])
  
  # calculate log-likelihood
  SummarizeBM(tree, Xsim, res$ParamsInModel, M.error = Me)$t_1$LogLik
}


# function to find ML with measurement errors
ml_fit_me <- function(tree, Xsim){
  
  # number of tips
  ntips <- Ntip(tree)
  
  # find the id for the tip with short branch
  if (ntips == 4){
    id <- 3
  }else{
    ids_edges <- get_short_tip(tree)
    id <- tree$edge[ids_edges[3],2]
  }
  
  # initial fit for obtaining a good seed for the next optim
  res <- mvSLOUCH::BrownianMotionModel(phyltree = tree, mData = Xsim)
  
  # extract results to vector for initial optim value
  v_init <- c(c(res$ParamsInModel$Sxx)[-3], res$ParamsInModel$vX0, 0, 0, 0)
  
  fopt <- function(v){
    # create a list of params from the first 5 indices
    pars <- list(Sxx = matrix(c(v[1], v[2], 0, v[3]), 2, 2), 
                 vX0 = matrix(v[4:5], 2, 1))
    
    # create a list of measurement errors from the last 3 indices
    Me <- v2Me(v[6:8], ntips = ntips, id = id)
    
    # calculate log-likelihood
    SummarizeBM(tree, Xsim, pars, M.error = Me)$t_1$LogLik
  }
  
  lower_lim <- c(0, -5*v_init[2], 0, -abs(res$ParamsInModel$vX0), 
                 0, -5*v_init[2], 0)
  upper_lim <- c(5*v_init[1:3], 5*abs(res$ParamsInModel$vX0), 5*v_init[1:3])
  
  optim(par = v_init, fn = fopt, method = "L-BFGS-B",
        lower = lower_lim, upper = upper_lim,
        control = list(fnscale = -1))
}


start <- Sys.time()
res_opt_4 <- ml_fit_me(TREES_4[[10]], Xsim_4)
as.numeric(difftime(Sys.time(), start, units = "secs"))


start <- Sys.time()
res_opt_100 <- ml_fit_me(TREES_100[[10]], Xsim_100)
as.numeric(difftime(Sys.time(), start, units = "secs"))



res <- mvSLOUCH::BrownianMotionModel(phyltree = TREES_100[[30]], mData = Xsim_100)
res$ParamsInModel
res$ParamSummary$LogLik

SummarizeBM(TREES_100[[30]], Xsim_100, res$ParamsInModel)
