## These are the R scripts and numerical results accompanying 
## CHECK FINAL AUTHORS
## Bartoszek, Brahmantio, Munoz-Duran, Fuentes-Gonzalez, Chi Kiang, Pienaar, and Polly
## TITLE?
## "TITLE?"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .
add_new_tip_split<-function(phyltree,sister_tip,new_branch_len){
    new_bl<-new_branch_len
    n<-max(phyltree$edge)-phyltree$Nnode ## number of tips
    if (sister_tip>n){stop(paste0("sister_tip is not a tip node, maximum possible is: ",n))}
    b_id<-which(phyltree$edge[,2]==sister_tip) ## id branch to split
    if(new_bl<phyltree$edge.length[b_id]){
	new_id<-max(phyltree$edge)+2 ## new tip node id
	phyltree$edge<-phyltree$edge+1 ## increase all current node ids by 1 as new tip will be added
	anc_id<-phyltree$edge[b_id,1] ## node id of direct ancestor of sister_tip
	phyltree$edge[b_id,]<-c(new_id,sister_tip+1) ## tip branch of sister_tip, entries of edge update sister_tip not
	phyltree$edge<-rbind(phyltree$edge,rbind(c(new_id,1),c(anc_id,new_id))) ## create edge leading to new tip and edge leading to its direct ancestor, i.e., split tip edge of sister_tip
	phyltree$edge.length<-c(phyltree$edge.length,new_bl,phyltree$edge.length[b_id]-new_bl) ## add lengths of new edges
	phyltree$edge.length[b_id]<-new_bl ## update length of tip edge of sister_tip
	phyltree$tip.label<-c(paste0("t",n+1),phyltree$tip.label) ## update tip labels
	phyltree$Nnode<-phyltree$Nnode+1 ## update number of internal nodes 
    }else{stop(paste0("New branch length is too long, has to be shorter than: ",new_bl))}
    phyltree
}

