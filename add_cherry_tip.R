add_new_tip_split<-function(phyltree,sister_tip,new_branch_len){
    new_bl<-new_branch_len
    n<-max(phyltree$edge)-phyltree$Nnode
    if (sister_tip>n){stop(paste0("sister_tip is not a tip node, maximum possible is: ",n))}
    b_id<-which(phyltree$edge[,2]==sister_tip)
    if(new_bl<phyltree$edge.length[b_id]){
	new_id<-max(phyltree$edge)+2	
	phyltree$edge<-phyltree$edge+1
	anc_id<-phyltree$edge[b_id,1]
	phyltree$edge[b_id,]<-c(new_id,sister_tip+1)
	phyltree$edge<-rbind(phyltree$edge,rbind(c(new_id,1),c(anc_id,new_id)))
	phyltree$edge.length<-c(phyltree$edge.length,new_bl,phyltree$edge.length[b_id]-new_bl)
	phyltree$edge.length[b_id]<-new_bl
    }else{stop(paste0("New branch length is too long, has to be shorter than: ",new_bl))}
    phyltree
}

Need to find error in below
> phyltree<-sim.bd.taxa(3,1,1,0)[[1]]
> newtree<-add_new_tip_split(phyltree,2,0.01)
Error in add_new_tip_split(phyltree, 2, 0.01) : 
  object 'phytree' not found
> source("add_cherry_tip.R")
> newtree<-add_new_tip_split(phyltree,2,0.01)
Error in phyltree$edge[b_id, 1] <- c(new_id, sister_tip) : 
  number of items to replace is not a multiple of replacement length
> source("add_cherry_tip.R")
> newtree<-add_new_tip_split(phyltree,2,0.01)
> plot(newtree)
Error in plot.phylo(newtree) : 
  tree badly conformed; cannot plot. Check the edge matrix.
> phyltree$edge
     [,1] [,2]
[1,]    4    5
[2,]    5    2
[3,]    5    1
[4,]    4    3
> newtree$edge
     [,1] [,2]
[1,]    5    6
[2,]    7    2
[3,]    6    2
[4,]    5    4
[5,]    7    1
[6,]    6    7
> source("add_cherry_tip.R")
> newtree<-add_new_tip_split(phyltree,2,0.01)
> plot(newtree)
Error in plot.phylo(newtree) : 
  tree badly conformed; cannot plot. Check the edge matrix.
> newtree$edge
     [,1] [,2]
[1,]    5    6
[2,]    7    3
[3,]    6    2
[4,]    5    4
[5,]    7    1
[6,]    6    7
