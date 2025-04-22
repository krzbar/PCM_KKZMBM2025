## These are the R scripts and numerical results accompanying 
## Bartoszek, Brahmantio, Munoz-Duran, Fuentes-Gonzalez,Pienaar, and Polly
## Short branch singularities in phylogenetic comparative methods

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(ape)
library(geiger)
library(mvSLOUCH)



# Phylogenetic tree
tree = read.tree("TreeNonUltra.tre")
# Chronogram
Tree<-chronos(tree,lambda=0.1)

# Shape variables
dat<-read.csv("RW.csv",header = T,row.names = 1)
name.check(Tree,dat)

# Alignment
row.names(dat) == Tree$tip.label
Data <- dat[match(Tree$tip.label, row.names(dat)), ]
row.names(Data) == Tree$tip.label


# mvSLOUCH setup
# Data
mvData<-data.matrix(Data)
save(mvData, file="mvData.RData")

ape::write.tree(Tree,file = "Tree.tre")
