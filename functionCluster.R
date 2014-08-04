extendBegin<-function(begin,sensGraph){
	
	test<-FALSE
	sens<-sensGraph[begin]
	while(test == FALSE){
		indiceValle<-begin -1
		nbegin<-max(which(sens != sensGraph[1:indiceValle]))
#		print(nbegin)
		if(begin == 51){ 
			test<-TRUE
		}else{
		
			tmpbegin<-nbegin-50
			test<- sum(sensGraph[tmpbegin:nbegin]) < -10
			if(is.na(test)){
				test<-TRUE
				nbegin<-0
				
			}
			begin<-nbegin-1
		}
		message<-paste(test," ",begin,sep="")
#		print(message)
	}
	nbegin<-nbegin+1
	return(nbegin)



}

extendEnd<-function(end,sensGraph){

	test<-FALSE
	sens<-sensGraph[end]


	while(test == FALSE){
		indiceDescente<-end+1
		nend<-min(which(sensGraph[end] != sensGraph[indiceDescente:length(sensGraph)]))	
		nend<-length(sensGraph[1:end]) +nend
		tmpend<- nend + 50
#######
		if(is.infinite(tmpend) ){
			nend<-length(sensGraph)	
			test<-TRUE
		}else{
			test<-sum(sensGraph[nend:tmpend])  > 10
		}
#Si on atteind le dernier indice
		if(is.na(test) ){
			nend<-length(sensGraph)	
			test<-TRUE
		}
		end<-nend+1
	}	
	if(nend > length(sensGraph)){
		nend<-length(sensGraph)
	}
	return(nend)


}

fusionCluster<-function(begin,end){

	interXY<-0
	clusterTotal<-length(begin)
	finaleBegin<-NULL
	finaleEnd<-NULL				

	while(length(begin) != 1){		
		b1<- begin[1]
		e1<- end[1]
		b2<- begin[2]
		e2<- end[2]

		x<-seq(b1:e1) + b1 - 1
		y<-seq(b2:e2) + b2 - 1
		interXY<-intersect(x,y)  
		if(b2 - e1 > 100 ){
			finaleBegin<-c(finaleBegin,b1)
			finaleEnd<-c(finaleEnd,e1)
			begin<-begin[-1]
			end<-end[-1]
		}else{
			rXY<-range(union(x,y))
			begin<-begin[-1]
			end<-end[-1]
			begin[1]<-rXY[1]
			end[1]<-rXY[2]

		} 
	}

	finaleBegin<-c(finaleBegin,begin)	
	finaleEnd<-c(finaleEnd,end)	

	coordClust<-cbind(finaleBegin,finaleEnd)
	print(finaleBegin)
	print(finaleEnd)
	return(coordClust)
}

extractCluster<-function(coordClust,pval, pref="clust",dir="./"){

	fileList<-paste(dir,"/",pref,"/list.txt",sep="")
	outNameList<-vector()
	for(i in 1:nrow(coordClust)){
		clust<- rownames(pval)[coordClust[i,1] : coordClust[i,2]]
		pathName<-paste(dir,"/",pref,sep="")
		if(!file.exists(pathName)){
			dir.create(pathName)
		}
		outname<-paste(pathName,"/",pref,"_",i,".txt",sep="")
		cat(outname,"\n",file=fileList,append=TRUE,sep="")
		write.table(clust,outname,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
		outNameList<-c(outNameList,outname)
	}
	return(outNameList)

}



sensClust<-function(coord,value){
	
	for(i in 1:nrow(coord)){
		clustSign<-sign(value[which(abs(value) == max(abs(value[coord$finaleBegin[i] : coord$finaleEnd[i]])))])
		if(clustSign == -1){
			coord$UP[i]<- -1
		}else{
			coord$UP[i]<- 1
			
		}
	}
	return(coord)
}


`graphMmobile` <-
function(filename,value,seuil=NULL,pas="",title="",clust=NULL,ylim=NULL){

if(pas==""){
	pas <- length(value)*0.7/100
	pas<-round(pas,0)
}
curveMobile<-rep(1/pas,pas)
fil<-filter(value,curveMobile)

png(filename=filename,width = 1300, height = 900, bg = "white", res = NA )
par(mar=c(0,0,0,0),oma=c(0,0,0,0),mgp=c(0,0,0))
yPos<-ylim-2
if(!is.null(ylim)){
	ylim=c(-ylim,ylim)
}else{
	ylim= range(value,na.rm=TRUE)
}


plot(value,pch=19,cex=0.8,cex.lab=1.5,cex.axis=3,col="azure3",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(1,length(value)),ylim=ylim,main=title)
lines(fil,lwd=5,col="orange2")
axis(side=2,tcl=0.5,labels = TRUE,pos=c(0,0),cex.axis=3)

if (! is.null(seuil)){
	abline(h=seuil,lwd=2,pch=19,cex=0.8,cex.lab=1.5,cex.axis=1.5,col="darkgreen",lty=2)
}

if( ! is.null(clust)){
	
	abline(h=-seuil,lwd=2,pch=19,cex=0.8,cex.lab=1.5,cex.axis=1.5,col="darkgreen",lty=2)
	abline(h=0,lwd=2,pch=19,cex=0.8,cex.lab=1.5,cex.axis=1.5,col="darkred")
	if(nrow(clust) !=0){
		col<-gsub(x=clust$UP,pattern=-1,replacement=3)
		col<-gsub(x=col,pattern=1,replacement=2)
		abline(v=clust[,1],lwd=2,pch=19,cex=0.8,cex.lab=1.5,cex.axis=1.5,col="red")
		abline(v=clust[,2],lwd=2,pch=19,cex=0.8,cex.lab=1.5,cex.axis=1.5,col="red")

		xPos<-clust$nProbes/2 + clust[,1] 
		yPos<-rep(yPos,nrow(clust))
		yPos<--clust$UP * yPos

		
		labels<-paste(rep("C ",nrow(clust)), clust$IdValid, rep(" : P : ",nrow(clust)), clust$nProbes,sep="")
		

		par(srt=90)
		text(xPos,yPos,labels=labels)
	}
}
dev.off()
}

options( show.error.messages=TRUE,echo=TRUE, verbose = FALSE, warn = -1  )
args <- commandArgs(trailingOnly =TRUE)
filePval <- args[1]
fileCluster <- args[2]
out_img <- args[3]

seuil <- as.numeric(args[4])
seuilCluster <- as.numeric(args[5]
pas <- as.numeric(args[6])

tabpval <- read.delim(filePval, row.names = 1, header =TRUE)
pval <- tabpval$P.Value

value<- -log10(pval)

curveMobile<-rep(1/pas,pas)
fil<-filter(value,curveMobile)
select<- fil > seuil
selectm1<-select[-1]
selectm1[length(select)]<-NA
	
	#indice de debut et de fin, +1 pour le debut car le TRUE et sur le "brin" de longeur n-1

begin <- which(select == FALSE & selectm1 == TRUE)
end <- which(select == TRUE & selectm1 == FALSE)
print(begin)
print(end)
if(length(begin) != 0){
	
	#sens graphe est de longeur n-1 donne le signe de l'indice n+1 (indice 1 de sens graphe donne le sens de l indice 2 de value
	sensGraph<-sign(diff(fil))
		
		#Preparation a l'extension des pics : Recherche de la prochaine valle pour le debut et la fin
	newBegin<-sapply(begin,extendBegin,sensGraph)
	newEnd<-sapply(end,extendEnd,sensGraph)
	
	tailleCluster<-newEnd - newBegin

	
	selectBegin<-newBegin[which(tailleCluster >= seuilCluster)]
	selectEnd<-newEnd[which(tailleCluster >= seuilCluster)]

	if( length(selectBegin) != 0  || length(selectEnd) !=0){

		coordClust<-fusionCluster(selectBegin,selectEnd)
		coordClust<-as.data.frame(coordClust)
		coordClust$UP<-rep(TRUE,nrow(coordClust))
		coordClust$nProbes<-coordClust[,2] -coordClust[,1]
		coordClust$IdInit<-seq(1:nrow(coordClust))
		coordClust$IdValid<-rep(0,nrow(coordClust))

	}else{
			
		coordClust<-data.frame()	
	}
			
		

	for(i in 1:nrow(coordClust)){
	clustselect <- rownames(tabpval)[coordClust[i,1]:coordClust[i,2]]
	  write(clustselect, paste("clust_",i,".txt",sep=""))
	}

	graphMmobile(out_img, value,seuil= seuil,pas="",title="",clust=coordClust)


}else{
	noCluster<-data.frame()
	print(noCluster)	
}



