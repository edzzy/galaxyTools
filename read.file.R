library(limma, quietly= TRUE )
options(echo = FALSE)
args <- commandArgs(trailingOnly =TRUE)
samplesData <- args[1]
colSignal <- args[2]
remove_control <- args[3]
out_file <- args[4]
input_files <- args[5:length(args)]

if(remove_control == "true"){
	remove_control = TRUE
}
pData <- read.delim(samplesData, header = TRUE, row.names = 1, sep = "\t")
print(pData)

files <- input_files[seq(from = 1, to = length(input_files), by = 2)]
print(files)
names(files) <- removeExt(input_files[seq(from = 2, to = length(input_files), by = 2)])




dye<-length(unique(tolower(pData$Dye)))
other.col<-c()
source<-c()
other.col<-c("FeatureNum")
source="agilent"
namesEch<-c()

if(dye == 2){
	cols<-list(R = paste("r", colSignal, sep = ""), G = paste("g",colSignal, sep = ""), Rb = paste("rBG", colSignal, sep = ""), Gb = paste("gBG", colSignal, sep = ""))
}else{
	cols<-list(R = paste("g", colSignal, sep = ""), G = paste("g",colSignal, sep = ""), Rb = paste("gBG", colSignal, sep = ""), Gb = paste("gBG", colSignal, sep = ""))
}

if ( dye == 2 ){
	namesEch<-c(rownames(pData)[which(tolower(pData$Dye) == "cy5")],
       		 rownames(pData)[which(tolower(pData$Dye) == "cy3")])
}else {
	namesEch<-rownames(pData)
}

RG<-read.maimages(files, source = source, columns = cols, other.columns = other.col,  verbose = FALSE)

if(dye == 2){
	pDataCy5<-pData[which(tolower(pData$Dye)== "cy5"),]
print(head(pDataCy5))
      	pDataCy3<-pData[which(tolower(pData$Dye)== "cy3"),]
print(head(pDataCy3))
      	R<-RG$R
      	G<-RG$G
	colnames(R) <- names(files)
	colnames(G) <- names(files)


      	R<-R[,intersect(removeExt(pDataCy5$nameFile),colnames(R))]
      	G<-G[,intersect(removeExt(pDataCy3$nameFile),colnames(G))]
	print(names)
	
	namesEch <- c(rownames(pDataCy5), rownames(pDataCy3))
	data<-cbind(R,G)


			
}else {
		data<-RG$G
		colnames(data) <- names(files)
		pDataCy3<-pData[which(tolower(pData$Dye)== "cy3"),]
      		data<-data[,intersect(removeExt(pDataCy3$nameFile),colnames(data))]
}
print(namesEch)
print(head(data))
colnames(data) <- namesEch
geneNames <- paste(RG$genes$GeneName,"|",RG$genes$ProbeName,"|", RG$other$FeatureNum[1:nrow(data)], sep="")
rownames(data) <- geneNames

if(remove_control){
	ctrl <- which(RG$genes$ControlType != 0)
	data <-data[-ctrl, ]
}

write.table(data, out_file , sep="\t", row.names = TRUE, col.names = NA, quote = FALSE)
