
options( show.error.messages=FALSE,echo=FALSE, verbose = FALSE, warn = -1  )
##Script########
#library(parallel, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)
#library(methods, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)
#library(BiocGenerics, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)
suppressMessages(library(affy,  verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE))
#library(affy,verbose  = FALSE )

#options(echo = FALSE, show)
args <- commandArgs(trailingOnly =TRUE)



i_samples <- args[1]
normalize <- args[2]
bgcorrect <- args[3]
pmcorrect <- args[4]
summary <- args[5]
o_file <-args[6] 
files <- args[7:length(args)]
files <- unlist(strsplit(files, " "))

filesnames <- files[seq(1,length(files), 2)]
filesCEL <- files[seq(2,length(files), 2)]

ptm <- proc.time()
pData <- read.delim(i_samples, header=TRUE, stringsAsFactor= FALSE)
Data <- suppressMessages(ReadAffy(filenames=filesnames, sampleNames=as.character(pData$NameID)))
eset <- suppressMessages(expresso(Data, bgcorrect.method = bgcorrect, normalize.method = normalize, pmcorrect.method = pmcorrect, summary.method = summary, verbose=FALSE))
#eset <- suppressMessages(rma(Data,verbose = FALSE))
mat <- exprs(eset)
write.table(mat, file = o_file, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)
print(proc.time() - ptm)
