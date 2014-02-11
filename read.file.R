read.Agilent.PseudoColor <- function(pSetup, sampleSheet){

#' @title read agilent files pseudo
#' @description read agilent files and return a Expression Set object

#' @param files : names of array files
#' @param sampleSheet : data frame with files and sample names association

#' @export

	RG<-read.maimages(unique(pData$nomFichiers), source = source, columns = cols, other.columns = other.col, path = get, verbose=FALSE)

}
read_Array <-
function(pData,namesGenes,type="AG",pathDir="",Flag=FALSE, rmCtrl = TRUE){
#Description : La fonction renvoie une matrice d'expression en lisant les fichiers
			# issus de logiciel d'extraction d'image avec dans chaque colonne
			# le signal median du canal 
#pData : data frame avec 2 colonnes obligatoire nomFichier Dye
#type : type de puce, agilent, nimbgene, gpr ou custum : (AG, NG, GPR, CUST) vecteur(variable) character
#namesGenes : nom des identifiants de g?nes vecteur character


#Verification des parametres d'entree :

	require(limma)
	if(missing(type))
		stop(" Type de puce manquantes ")
	if(! type == "AG" & ! type == "NG" & ! type == "GPR")
		stop("Type de puce non correcte : seul les types AG (agilent) NG (NimbelGene) et GPR (issue de GenePix) sont disponibles")
	
		if(missing(namesGenes))
		warning("Le nom pour les genes n'est pas donne. 
              Les colonnes GeneID et FeatureNum seront utilisees par defaut")
  
	dye<-length(unique(tolower(pData$Dye)))

	#Colones des flag a considérer sans le prefix r ou g
	flagCol<-c("IsFeatNonUnifOL","IsSaturated","IsPosandSignif","IsFeatPopnOL")
#	print(flagCol)

	#Setup des arguments pour la creation de la matrice.
	# type de puce et noms des colones a extraire.
	other.col<-c()
	source<-c()
	if( type == "AG"){
		other.col<-c("FeatureNum")
		source="agilent"

		if(Flag==TRUE){
			other.col<-c(other.col,paste("g",flagCol,sep=""))	
#				print(other.col)
		}
		if(dye == 1){
			cols<-list(R="gMedianSignal",G="gMedianSignal",Rb="gBGMedianSignal",Gb="gBGMedianSignal")

		}else{
			cols<-list(R="rMedianSignal",G="gMedianSignal",Rb="rBGMedianSignal",Gb="gBGMedianSignal")

			if(Flag==TRUE){
				other.col<-c(other.col,paste("r",flagCol,sep=""))	
#				print(other.col)
			}
		}
		
	} else if (type == "NG"){
		source="generic"
		dye<-1
		cols<-list(R="PM",G="PM",Rb="MM",Gb="MM") #Un seul canal par fichier pair mono couleur
		other.col=c("PROBE_ID")
	
	} else if (type == "GPR"){
		source="genepix.median"
		col<-list(R = "F635 Median", G = "F532 Median", Rb = "B635 Median", Gb = "B532 Median")
		other.col<- c("")
	}

	#parametre du nom des echantillons
	namesEch<-c()
	if ( dye == 2 && (type == "AG" || type == "GPR" ) ){
		namesEch<-c(rownames(pData)[which(tolower(pData$Dye) == "cy5")],
        rownames(pData)[which(tolower(pData$Dye) == "cy3")])
    
	}else {
		namesEch<-rownames(pData)
	}


	#Lecture des fichiers pour creer un objet expressionArray. (limma)
	RG<-read.maimages(unique(pData$nomFichiers),source=source,columns=cols,other.columns=other.col,path=pathDir,verbose=FALSE)

	if(dye == 2){
      pDataCy5<-pData[which(tolower(pData$Dye)== "cy5"),]
      pDataCy3<-pData[which(tolower(pData$Dye)== "cy3"),]
      R<-RG$R
      G<-RG$G
      R<-R[,intersect(removeExt(pDataCy5$nomFichiers),colnames(R))]
      G<-G[,intersect(removeExt(pDataCy3$nomFichiers),colnames(G))]

	data<-cbind(R,G)


			
	}else {
		data<-RG$G
		pDataCy3<-pData[which(tolower(pData$Dye)== "cy3"),]
	}
	#nom des echantillons
	
	colnames(data)<-namesEch

	#noms des genes
	geneNames <- paste(RG$genes$GeneName,"//",RG$genes$ProbeName,"//", RG$other$FeatureNum[1:ncol(data)],"//",RG$genes$ControlType, sep="")
#	rownames(data) <- geneNames
	ctrl <- which(RG$genes$ControlType != 0)
	data <-data[-ctrl, ]

	seuilFM<-round(0.5*length(namesEch))
#	seuilFM<-2

if(Flag){
##Elimination des Flags
	filesummary<-"summaryFlag.txt"
	rmProbesAll<-c()
	control<-which(RG$gene$ControlType != 0)
	PosControl<-which(RG$gene$ControlType == 1)
	NegControl<-which(RG$gene$ControlType == -1)
	data<-data[-control,]
	namesGenes<-namesGenes[-control]
	write(paste("Controle Positif : ", length(PosControl),sep=""),file=filesummary,append=TRUE)
	write(paste("Controle Negatif: ", length(NegControl),sep=""),file=filesummary,append=TRUE)

	if(!is.null(RG$other$gIsFeatNonUnifOL)){
	#Création de la matrice IsFeatNonUnifOL
		gIsFeatNonUnifOL<-RG$other$gIsFeatNonUnifOL
		gIsFeatNonUnifOL<-gIsFeatNonUnifOL[,intersect(removeExt(pDataCy3$nomFichiers),colnames(gIsFeatNonUnifOL))]
		IsFeatNonUnifOL<-cbind(gIsFeatNonUnifOL)
		if(dye==2){
			
			rIsFeatNonUnifOL<-RG$other$rIsFeatNonUnifOL
			rIsFeatNonUnifOL<-rIsFeatNonUnifOL[,intersect(removeExt(pDataCy5$nomFichiers),colnames(rIsFeatNonUnifOL))]
			IsFeatNonUnifOL<-cbind(rIsFeatNonUnifOL,gIsFeatNonUnifOL)
		}
		IsFeatNonUnifOL<-IsFeatNonUnifOL[-control,]
		rownames(IsFeatNonUnifOL)<-namesGenes
		colnames(IsFeatNonUnifOL)<-namesEch
		
		#nombre flag par sonde
		sumIsFeatNonUnifOL<-apply(IsFeatNonUnifOL,1,sum)
		
		#synthese
		write("IsFeatNonUnifOL",file=filesummary,append=FALSE)
		write.table(summary(IsFeatNonUnifOL == 1),file=filesummary,sep="\t",append=TRUE,quote=FALSE)
		summaryLaTeX<-summary(IsFeatNonUnifOL == 1)

		#Sonde à retirer si le dépasse le seuil
		rmProbesIsFeatNonUnifOL<-which(sumIsFeatNonUnifOL >= seuilFM)
		rmProbesAll<-c(rmProbesAll,rmProbesIsFeatNonUnifOL)
		
		write(paste("IsFeatNonUnifOL  ",length(rmProbesIsFeatNonUnifOL)," sondes éliminés",sep=""),file=filesummary,append=TRUE)
	}

	#Création de la matrice IsSaturated
	if(!is.null(RG$other$gIsSaturated)){
		gIsSaturated<-RG$other$gIsSaturated
		gIsSaturated<-gIsSaturated[,intersect(removeExt(pDataCy3$nomFichiers),colnames(gIsSaturated))]
		IsSaturated<-cbind(gIsSaturated)
		if(dye==2){
			rIsSaturated<-RG$other$rIsSaturated
			rIsSaturated<-rIsSaturated[,intersect(removeExt(pDataCy5$nomFichiers),colnames(rIsSaturated))]
			IsSaturated<-cbind(rIsSaturated,gIsSaturated)
		}
		IsSaturated<-IsSaturated[-control,]
		rownames(IsSaturated)<-namesGenes
		colnames(IsSaturated)<-namesEch

		sumIsSaturated<-apply(IsSaturated,1,sum)

		write("IsSaturated",file=filesummary,append=TRUE)
		write.table(summary(IsSaturated == 1),file=filesummary,sep="\t",quote=FALSE,append=TRUE)
		
		rmProbesIsSaturated<-which(sumIsSaturated >= seuilFM)
	#	print(rmProbesIsSaturated)
		if(length(rmProbesIsSaturated) != 0){
			rmProbesAll<-c(rmProbesAll,rmProbesIsSaturated)
		}
	#	print(length(rmProbesAll))
		write(paste("IsSaturated",length(rmProbesIsSaturated)," sondes éliminés",sep=""),file=filesummary,append=TRUE)
	}

	#Création de la matrice IsPosandSignif
	if(!is.null(RG$other$gIsPosandSignif)){
		gIsPosandSignif<-RG$other$gIsPosandSignif
		gIsPosandSignif<-gIsPosandSignif[,intersect(removeExt(pDataCy3$nomFichiers),colnames(gIsPosandSignif))]
		IsPosandSignif<-cbind(gIsPosandSignif)
		if(dye==2){
			rIsPosandSignif<-RG$other$rIsPosandSignif
			rIsPosandSignif<-rIsPosandSignif[,intersect(removeExt(pDataCy5$nomFichiers),colnames(rIsPosandSignif))]
			IsPosandSignif<-cbind(rIsPosandSignif,gIsPosandSignif)
		}
		IsPosandSignif<-IsPosandSignif[-control,]
		rownames(IsPosandSignif)<-namesGenes
		colnames(IsPosandSignif)<-namesEch
		sumIsPosandSignif<-sum(IsPosandSignif,1,sum)

	#	print("IsPosandSignif")
	#	print(summary(IsPosandSignif == 1))
		summaryLaTeX<-summary(IsPosandSignif == 1)

		#Sonde à retirer si le dépasse le seuil
		rmProbesIsPosandSignif<-which(sumIsPosandSignif >= seuilFM)
		rmProbesAll<-c(rmProbesAll,rmProbesIsPosandSignif)
	#	print(length(rmProbesIsPosandSignif))
	#	print(length(rmProbesAll))

	}

	#Création de la matrice IsFeatPopnOL
	if(!is.null(RG$other$gIsFeatPopnOL)){
		gIsFeatPopnOL<-RG$other$gIsFeatPopnOL
		gIsFeatPopnOL<-gIsFeatPopnOL[,intersect(removeExt(pDataCy3$nomFichiers),colnames(gIsFeatPopnOL))]
		IsFeatPopnOL<-cbind(gIsFeatPopnOL)
		if(dye==2){
			rIsFeatPopnOL<-RG$other$rIsFeatPopnOL
			rIsFeatPopnOL<-rIsFeatPopnOL[,intersect(removeExt(pDataCy5$nomFichiers),colnames(rIsFeatPopnOL))]
			IsFeatPopnOL<-cbind(rIsFeatPopnOL,gIsFeatPopnOL)
		}
		IsFeatPopnOL<-IsFeatPopnOL[-control,]
		rownames(IsFeatPopnOL)<-namesGenes
		colnames(IsFeatPopnOL)<-namesEch
		sumIsFeatPopnOL<-apply(IsFeatPopnOL,1,sum)
		write("IsFeatPopnOL",file=filesummary,append=TRUE)
		write.table(summary(IsFeatPopnOL== 1),file=filesummary,sep="\t",quote=FALSE,append=TRUE)

		summaryLaTeX<-summary(IsFeatPopnOL == 1)

		#Sonde à retirer si le dépasse le seuil
		rmProbesIsFeatPopnOL<-which(sumIsFeatPopnOL >= seuilFM)
		rmProbesAll<-c(rmProbesAll,rmProbesIsFeatPopnOL)
		write(paste("IsFeatPopnOL ",length(rmProbesIsFeatPopnOL)," sondes éliminés",sep=""),file=filesummary,append=TRUE)

	}

#Bilan nombre de flag dans chaque échantillons et part type de flag


#Elimination des controle

#rmProbesAll<-c(rmProbesAll,control)
write(paste("Sondes éliminées au total : ",length(unique(rmProbesAll)), " avec au moins ", seuilFM, " échantilons avec un flag",sep=""),file=filesummary,append=TRUE)

#Bilan nombre de controle
		if(is.null(rmProbesAll)){
			data<-data[-rmProbesAll,]
		}
	}
	return(data)
}

`read.cdt` <-
function(file){
	cdt<-read.delim(file, sep="\t")
	cdt<-cdt[,-1]
	cdt<-cdt[,-1]
	cdt<-cdt[,-2]
	cdt<-cdt[-1,]
	cdt<-cdt[-1,]
	row.names(cdt)<-cdt[,1]
	cdt<-cdt[,-1]
	cdt<-as.data.frame(cdt)
	return(cdt)
}

`read.matrixArray` <-
function(file){
	mat<-read.delim(file, sep="\t")
	row.names(mat)<-mat[,1]
	mat<-mat[,-1]
	return(mat)

}

