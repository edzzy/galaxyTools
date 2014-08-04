############FUNCTION################
my.lowess <-
function(matriceW,vecRef,f)
{
	index = which(is.na(matriceW))
	indexVal = which(!is.na(matriceW))

	if (length(index)>1) {
		vecRef = vecRef[-index]
		matriceW2 = matriceW[-index]

		result = lowess(vecRef,matriceW2,f=f)
		temp = result$y

		sortie = rep(NA,length(matriceW))
		sortie[indexVal]=temp
	}

	else{
		result = lowess(vecRef,matriceW,f=f)
		sortie = result$y
	}
	sortie
}


calculeNorm <-
function(increm, echantillon, lowessCurve, ref)
{
	result = echantillon[,increm]/lowessCurve[,increm]*ref
	return(result)
}


##Script########
library(getopt)
spec <- matrix(c(
								 	'i_matrix', 'm', 1, "character",
								 	'i_norm', 'n', 1, "character",
								 	'log', 'l', 0, "logical",
									'o_file','o',  1,  "character",
									'o_path', 'p', 1, "character",
									'verbose', 'v', 2, "integer",
									'help', 'h', 0, "logical"
									), byrow = TRUE, ncol = 4)
	
opt <- getopt(spec)
if (!is.null(opt$help)){
	cat(getopt(spec, usage = TRUE))
	q(status=1)
}

if(is.null(opt$i_matrix		)) 	{stop("Need input matrix")} 		
if(is.null(opt$i_norm)) 	{opt$i_norm = "L"} 
if(is.null(opt$log 						)) 	{opt$log = FALSE}
if(is.null(opt$o_path 						)) 	{opt$o_path = "./"}
if(is.null(opt$o_file			)) 	{opt$o_file = paste("matrix_norm_",format(Sys.Date(), "%y%m%d"), sep="")} 			

if (is.null(opt$verbose)){
	opt$verbose = FALSE
}else{
	opt$verbose = TRUE 
}

options(echo = opt$verbose)


i_matrix <- opt$i_matrix
i_norm <- opt$i_norm
o_file <- opt$o_file
o_path <- opt$o_path
log <- opt$log

if(!file.exists(o_path)){
	dir.create(o_path)
}

data <- read.delim(i_matrix, header=TRUE, row.names = 1)
if(log == TRUE){
	data <- log(data)
}

if(i_norm == "L"){
	
	nomGenes = rownames(data)
	nomEchan = colnames(data)
	
# Calcul du profil median
	profil.median = apply(data, 1, median, na.rm = TRUE)
	
######### Tri par ordre croissant des profils medians  #######
	
#On ordonne la matrice par ordre croissant des profils median
	ordre.profil = order(profil.median)
	
#On range la matrice de depart
	matOrdonne = data[ordre.profil,]
	
#On range le vecteur de profil median
	profilOrdonne = profil.median[ordre.profil]
	
#On range le vecteur des noms de genes
	nomGenes = as.character(nomGenes)
	nomGenes = nomGenes[ordre.profil]	
########################################################
	
############# Lowess #######################	
#Matrice lowess
# On passe les donnees en log
#pngDir<-"../,"02-lowess/"
	lowessCurve = NULL
	#	lowessCurve = cbind(lowessCurve,apply(log(matOrdonne),2,my.lowess,log(profilOrdonne),f=0.01))
	lowessCurve = cbind(lowessCurve,apply(matOrdonne,2,my.lowess,profilOrdonne,f=0.01))
	# Calcul de la matrice normalisee
	increm=c(1:dim(data)[2])
#	if(log==FALSE){
#		matNorm = sapply(increm, calculeNorm, matOrdonne, exp(lowessCurve), profilOrdonne)
#	}else{
	matNorm = sapply(increm, calculeNorm, matOrdonne, lowessCurve, profilOrdonne)
#	}
	profilMedNorm = apply(matNorm, 1, median, na.rm = TRUE)
	
	
#  if(graph==1){	
#  	########## Plots ###########################
#	require(multicore)
#  	diagonal=c(min(profilOrdonne,na.rm=T), max(profilOrdonne,na.rm=T))
#  	
#  	# Graphs avant normalisation pour tous les echantillons
#  	graph=mclapply(increm, traceGraph, matOrdonne, profilOrdonne, diagonal, lowessCurve,nomEchan,"1-Avant", pngDir,mc.cores=10)
#		
#  	# Graphs apres normalisation pour tous les echantillons
#  	graph=mclapply(increm, traceGraph, matNorm, profilMedNorm, diagonal, NULL,nomEchan,"2-Apres", pngDir,mc.cores=10)
#	
#  	#Graph valeurs brutes/valeurs normalisees
#  	graph=mclapply(increm, traceAvantApres, matOrdonne, matNorm, nomEchan, pngDir,mc.cores=10)
#	} ## If graph
	##############################
	
	### Exportation des resultats	
	rownames(matNorm)=nomGenes
	colnames(matNorm) = nomEchan

}


if(i_norm == "Q"){
	library(limma)

	matNorm <- normalizeBetweenArrays(as.matrix(data), method = "quantile")

}


matNorm<-round(matNorm,2)
write.table(matNorm, file = o_file, col.names = NA, row.names=TRUE, quote = FALSE, sep = "\t")

