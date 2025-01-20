## These are the R scripts and numerical results accompanying 
## CHECK FINAL AUTHORS
## Bartoszek, Brahmantio, Munoz-Duran, Fuentes-Gonzalez, Chi Kiang, Pienaar, and Polly
## TITLE?
## "TITLE?"

## The R setup for the manuscript was as follows: R version 4.4.1 (2024-06-14) 
## Platform: x86_64-pc-linux-gnu (64-bit) Running under: openSUSE Leap 15.6

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(ape)
library(mvSLOUCH)
library(TreeSim)
library(mvMORPH)

vN<-c(5,10,30,50,100,200,1000) ##2000,5000,10000,20000,50000)
numreps<-30
BMparams<-list(X0=matrix(0,ncol=1,nrow=4),Sigma=rbind(c(1,0,0,0),c(0.5,1,0,0),c(0.5,0.5,1,0),c(0.5,0.5,0.5,1)))
filepath<-"."
fileprefix<-"mvSLOUCHmvMORPH_BM_timings"
vcolours<-c("red","blue","orange")
b_dosimulation<-TRUE
c_plotlogscale<-TRUE

if (b_dosimulation){
    rexp(1)

    ## simulate, reestimate, get timings
    simulres<-sapply(vN,function(n,numreps,BMparams,filepath,fileprefix){
	sink(paste0(filepath,"/",fileprefix,"_progress.txt"),append=TRUE)
        cat(paste0("Doing n=",n,"\n"))
	num_errors<-0
        res_n<-sapply(1:numreps,function(i,BMparams,filepath,fileprefix){
    	    cat(paste0("Doing iteration i=",i," ",Sys.time(),"\n"))
            randomseed<-.Random.seed;set.seed(randomseed, "L'Ecuyer-CMRG");RNG_kind<-RNGkind();RNG_version<-getRversion()
	    phyltree<-TreeSim::sim.bd.taxa(n,1,1,0)[[1]]
    	    BMdata<-mvSLOUCH::simulBMProcPhylTree(phyltree, X0=BMparams$X0, Sigma=BMparams$Sigma)
	    start_mvslouch<-Sys.time()
	    resmvslouch<-mvSLOUCH::BrownianMotionModel(phyltree,BMdata)	
	    end_mvslouch<-Sys.time()
	    time_mvslouch<-end_mvslouch-start_mvslouch
	    cat("time mvSLOUCH: ",time_mvslouch)
	    cat("\n")
	    if (is.infinite(resmvslouch$ParamSummary$LogLik)){cat(paste0("Error in iteration: ",i,"\n"));time_mvslouch<-start_mvslouch<-end_mvslouch<-NA}
    		start_mvmorph_pic<-Sys.time()
		resmvmorph_pic<-mvMORPH::mvBM(phyltree,BMdata, model="BM1",method="pic") 
		end_mvmorph_pic<-Sys.time()
	        time_mvmorph_pic<-end_mvmorph_pic-start_mvmorph_pic
		cat("time mvMORPH pic: ",time_mvmorph_pic)
		cat("\n")
		start_mvmorph_rpf<-Sys.time()
		resmvmorph_rpf<-mvMORPH::mvBM(phyltree,BMdata, model="BM1",method="rpf") 
    		end_mvmorph_rpf<-Sys.time()
	        time_mvmorph_rpf<-end_mvmorph_rpf-start_mvmorph_rpf
		cat("time mvMORPH rpf: ",time_mvmorph_rpf)
        	cat("\n")
		list(i=i,n=n,phyltree=phyltree,BMdata=BMdata,start_mvslouch=start_mvslouch,end_mvslouch=end_mvslouch,start_mvmorph_pic=start_mvmorph_pic,end_mvmorph_pic=end_mvmorph_pic,start_mvmorph_rpf=start_mvmorph_rpf,end_mvmorph_rpf=end_mvmorph_rpf,time_mvslouch=time_mvslouch,time_mvmorph_pic=time_mvmorph_pic,time_mvmorph_rpf=time_mvmorph_rpf,resmvslouch=resmvslouch,resmvmorph_pic=resmvmorph_pic,resmvmorph_rpf=resmvmorph_rpf,randomseed=randomseed,RNG_kind=RNG_kind,RNG_version=RNG_version)
	},BMparams=BMparams,filepath=filepath,fileprefix=fileprefix,simplify=FALSE)
        mTimings<-t(sapply(res_n,function(x){c(x$time_mvslouch,x$time_mvmorph_pic,x$time_mvmorph_rpf)},simplify=TRUE))
	colnames(mTimings)<-c("mvSLOUCH","mvMORPH_pic","mvMORPH_rpf")
        num_errors<-length(which(is.na(mTimings[,1])))
	meantime<-apply(mTimings,2,mean,na.rm=TRUE)
    	vartime<-apply(mTimings,2,var,na.rm=TRUE)
	mediantime<-apply(mTimings,2,median,na.rm=TRUE)
    	interquant_time<-cbind(quantile(mTimings[,1],na.rm=TRUE)[c(2,4)],quantile(mTimings[,2],na.rm=TRUE)[c(2,4)],quantile(mTimings[,3],na.rm=TRUE)[c(2,4)])
	colnames(interquant_time)<-c("mvSLOUCH","mvMORPH_pic","mvMORPH_rpf")
	save(res_n,mTimings,meantime,vartime,mediantime,interquant_time,file=paste0(filepath,"/",fileprefix,"_n_",n,"_results.RData"))
	cat("========================================================= \n")
	sink()
    },numreps=numreps,BMparams=BMparams,filepath=filepath,fileprefix=fileprefix,simplify=FALSE)
}


## make plot and table


l_timings<-vector("list",5)
names(l_timings)<-c("mMeans","mVars","mMedians","lInterQuantiles","NAsmvSLOUCH")
l_timings$lInterQuantiles<-vector("list",length(vN))
names(l_timings$lInterQuantiles)<-paste0("n_",vN)
i<-1
for (n in vN){
    load(paste0(filepath,"/",fileprefix,"_n_",n,"_results.RData"))
    ##n_i<-which(names(l_timings)==paste0("n_",n))
    l_timings$mMeans<-rbind(l_timings$mMeans,meantime)
    l_timings$mVars<-rbind(l_timings$mVars,vartime)
    l_timings$mMedians<-rbind(l_timings$mMedians,mediantime)
    l_timings$lInterQuantiles[[i]]<-interquant_time
    l_timings$NAsmvSLOUCH<-c(l_timings$NAsmvSLOUCH,sum(is.na(mTimings[,1])))
    i<-i+1
}

png(paste0(filepath,"/",fileprefix,".png"))
sink(paste0(filepath,"/",fileprefix,".txt"))
v_ylim<-c(0.01,(max(l_timings$mMedians[i-1,1]+l_timings$lInterQuantiles[[i-1]][2,1],l_timings$mMedians[i-1,2]+l_timings$lInterQuantiles[[i-1]][2,2])+0.05))
c_ylab<-"time[s]"
if (c_plotlogscale){v_ylim<-log(v_ylim);c_ylab<-"log(time[s])"}
plot(NA,xlim=c(0,max(vN)+1),ylim=v_ylim,main="",xlab="n",ylab=c_ylab)

    for (j in 1:3){ ## mvSLOUCH, mvMORPH
	i<-1
	v_q1<-c()
	v_q2<-c()
	v_med<-c()

	for (n in vN){
	    len_bar<-0.025
	    #https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
    	    median_time<-l_timings$mMedians[i,j]
    	    interquant_time[1]<-median_time-1.5*(median_time-l_timings$lInterQuantiles[[i]][1,j])
	    interquant_time[2]<-median_time+1.5*(l_timings$lInterQuantiles[[i]][2,j]-median_time)
	    if (c_plotlogscale){interquant_time[1]<-log(interquant_time[1]);interquant_time[2]<-log(interquant_time[2])}
	    segments(n,interquant_time[1], n, interquant_time[2], col=vcolours[j],lwd=2)
	    segments(n-len_bar,interquant_time[1], n+len_bar, interquant_time[1], col=vcolours[j],lwd=2)
	    segments(n-len_bar,interquant_time[2], n+len_bar, interquant_time[2], col=vcolours[j],lwd=2)
	    i<-i+1
	    v_q1<-c(v_q1,interquant_time[1])
	    v_q2<-c(v_q2,interquant_time[2])
	    v_med<-c(v_med,median_time)
	}	
	cat("n                            ")
	for (n in vN){
	    cat(paste0(n,"                  "))
	}	
	cat("\n")
	if (j==1){
	    cat("short branch errors mvSLOUCH ")
	    for (i in 1:length(vN)){cat(paste0(l_timings$NAsmvSLOUCH[i],"               "))}	
	    cat("\n")
	}
	switch(j,
		    
	    cat("lower quantile mvSLOUCH    "),
	    cat("lower quantile mvMORPH pic "),
	    cat("lower quantile mvMORPH rpf ")	    
	)
	for (i in 1:length(vN)){cat(paste0(v_q1[i]," "))}	
	cat("\n")
	switch(j,
	    cat("median mvSLOUCH            "),
	    cat("median mvMORPH pic         "),
	    cat("median mvMORPH rpf         "),
	)
	for (i in 1:length(vN)){cat(paste0(v_med[i]," "))}	
	cat("\n")
	switch(j,
	    cat("upper quantile mvSLOUCH    "),
	    cat("upper quantile mvMORPH pic "),
	    cat("upper quantile mvMORPH rpf ")
	)
	for (i in 1:length(vN)){cat(paste0(v_q2[i]," "))}	
	cat("\n")
    }    
    legend("topleft",legend=c("mvSLOUCH 2.7.6","mvMORPH 1.2.1 pic","mvMORPH 1.2.1 rpf"),col=vcolours,pch=19,bty="n")
dev.off()
