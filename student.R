library(getopt)
spec <- matrix(c(
								 	'i_file', 'f', 1, "character",
								 	'i_sample', 's', 1, "character",
									'factor1', '1', 1, "character",
									'factor2', '2', 1, "character",
									'adjust', 'a', 1, "character",
									'threshold', 't', 1, "double",
									'log', 'l', 2, "logical",
									'o_file','o',  1,  "character",
									'o_selected', 'p', 1, "character",
									'verbose', 'v', 2, "integer",
									'help', 'h', 0, "logical"
									), byrow = TRUE, ncol = 4)
	
opt <- getopt(spec)
if (!is.null(opt$help)){
	cat(getopt(spec, usage = TRUE))
	q(status=1)
}

if(is.null(opt$i_file	)) 	{stop("Need input file")} 		
if(is.null(opt$i_sample)) {stop("Need sample file")}
if(is.null(opt$factor1)) 	{stop("Need factor name 1")}
if(is.null(opt$factor2)) 	{stop("Need factor name 2")}
if(is.null(opt$adjust))		{opt$adjust = "BH"}
if(is.null(opt$threshold)){opt$threshold = "0.05"}
if(is.null(opt$log)) 			{opt$log = TRUE}
if(is.null(opt$o_selected	)) 	{opt$o_selected =paste("selected-pvalue",opt$threshold,"-",format(Sys.Date(), "%y%m%d"), sep="")}
if(is.null(opt$o_file	)) 	{opt$o_file = paste("allpvalue",format(Sys.Date(), "%y%m%d"), sep="")} 			

if (is.null(opt$verbose)){
	opt$verbose = FALSE
}else{
	opt$verbose = TRUE 
}

options(echo = opt$verbose)

i_file <- opt$i_file
i_sample <- opt$i_sample
o_file <- opt$o_file
o_selected <- opt$o_selected
log <- opt$log
padj.method <- opt$adjust
threshold <- opt$threshold
print(threshold)
factor1 <- opt$factor1
factor2 <- opt$factor2



library(genefilter, verbose=FALSE)
mat <- read.delim(i_file, header = TRUE, sep= "\t", stringsAsFactor = FALSE, row.names = 1) 
mat <- as.matrix(mat)
samples <- read.delim(i_sample, header = TRUE, sep = "\t", stringsAsFactor = FALSE, row.names = 1)

samples <- samples[colnames(mat), ]

if(!all(colnames(mat) == rownames(samples))){
	stop("Samples in samples files and matrix files are not in same order")
	q(status=1)
}

if(log){
	mat = log(mat)
}
	

comparaison_name <-paste(factor1, factor2,sep="-VS-")
finale<-data.frame(row.names(mat))
	
#samplesSubset <- samples[which(samples == factor1 | samples == factor2)]
    
subNames <- row.names(samples[which(samples==factor1 | samples == factor2, arr.ind = TRUE)[,2]])
subSamples <- samples[subNames, ]
subMatrix <- mat[, subNames]

	
  
factorsSamples <- as.factor(samples[which(samples == factor1 | samples == factor2, arr.ind = TRUE)])

tt <- rowFtests(subMatrix, factorsSamples, var.equal = FALSE)
tt$p.value.adjust <- p.adjust(tt$p.value, method=padj.method)
tt <- tt[order(tt$p.value.adjust), ] 
ttselected <- tt[which(tt$p.value.adjust <= threshold), ]
#nameFile=paste(o_path,"/",o_file,"-student-",comparaison_name,".txt",sep="")
write.table(tt, o_file, row.names = TRUE, col.names =NA, sep="\t", quote = FALSE)
write.table(ttselected, o_selected, row.names = TRUE, col.names =NA, sep="\t", quote = FALSE)

