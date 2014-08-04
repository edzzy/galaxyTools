library(getopt)
suppressMessages(library(Biobase,verbose=FALSE))
spec <- matrix(c(
								 	'i_file', 'f', 1, "character",
								 	'i_sample', 's', 1, "character",
									'o_file','o',  1,  "character",
									'o_selected','p',  1,  "character",
									'o_img','i',  1,  "character",
									'factor1', '1', 1, "character", 
									'factor2', '2', 1, "character", 
									'attribut', 'a', 1, "character",
									'methodadjust', 'm', 1, "character",
									'log', 'l', 2, "logical",
									'verbose', 'v', 2, "integer",
									'threshold', 't', 1, "double",
									'help', 'h', 0, "logical"
									), byrow = TRUE, ncol = 4)
	
opt <- getopt(spec)
if (!is.null(opt$help)){
	cat(getopt(spec, usage = TRUE))
	q(status=1)
}

if(is.null(opt$i_file		)) 	{stop("Need input file")} 		
if(is.null(opt$i_sample		)) 	{stop("Need sample file")} 		
if(is.null(opt$factor1		)) 	{stop("Need factor 1")} 		
if(is.null(opt$factor2		)) 	{stop("Need factor 2")} 		
if(is.null(opt$attribut		)) 	{stop("Need attribut name (column in sample file")} 		
if(is.null(opt$methodadjust))	{opt$methodadjust = "BH"}
if(is.null(opt$threshold))	{opt$threshold = 0.05}
if(is.null(opt$o_file			)) 	{opt$o_file = paste("limmaAll-",format(Sys.Date(), "%y%m%d"), sep="")} 			
if(is.null(opt$o_selected			)) 	{opt$o_file = paste("limmaSemected-",format(Sys.Date(), "%y%m%d"), sep="")} 			
if(is.null(opt$o_img			)) 	{opt$o_file = paste("limma-",format(Sys.Date(), "%y%m%d"),".png", sep="")} 			

if (is.null(opt$verbose)){
	opt$verbose = FALSE
}else{
	opt$verbose = TRUE 
}

options(echo = opt$verbose)

i_file <- opt$i_file
i_sample <- opt$i_sample
attribut <- opt$attribut
threshold <- opt$threshold
factor1 <- opt$factor1
factor2 <- opt$factor2
methodadjust <- opt$methodadjust
o_file <- opt$o_file
o_selected <- opt$o_selected
o_img <- opt$o_img
if(is.null(opt$log)){
	log = FALSE
}else{
	log =TRUE
}
library(limma, verbose=FALSE)
mat <- read.delim(i_file, row.names = 1, header = TRUE, stringsAsFactor = FALSE)
if(log){
	mat =log(mat)
}
samples <- read.delim(i_sample, row.names = 1, header = TRUE, stringsAsFactor = FALSE)
samples <- samples[colnames(mat), ]
if(!all(colnames(mat) == rownames(samples))){
	stop("Samples in samples files and matrix files are not in same order")
	q(status=1)
}

subNames <- row.names(samples[which(samples==factor1 | samples == factor2, arr.ind = TRUE)[,2]])
subSamples <- samples[subNames, ]
subMatrix <- mat[, subNames]

eset <- new("ExpressionSet", exprs = as.matrix(subMatrix))
adf <- new("AnnotatedDataFrame", data = subSamples)
phenoData(eset) <- adf

design <- model.matrix(~0+as.factor(pData(phenoData(eset))[, attribut]))
levelsParameters <-  as.character(levels(as.factor(pData(phenoData(eset))[, attribut])))

colnames(design) <- levelsParameters

comparaison <- paste(factor1 ,factor2, sep="-") 
contrasts.matrix <- makeContrasts(contrasts = comparaison, levels = design)

comparaison_name <-paste(factor1, factor2,sep="-VS-")

fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts.matrix)
fit2 <- eBayes(fit2)
top <- topTable(fit2, coef = 1, adjust = methodadjust, p.value = threshold, number = nrow(mat))
topAll <- topTable(fit2, coef = 1, adjust = methodadjust, p.value = 1, number = nrow(mat))
nameFile=o_file
write.table(topAll, nameFile, row.names = TRUE, col.names =NA, sep="\t", quote = FALSE)
write.table(top, o_selected, row.names = TRUE, col.names =NA, sep="\t", quote = FALSE)
png(o_img)
volcanoplot(fit2)
dev.off()
