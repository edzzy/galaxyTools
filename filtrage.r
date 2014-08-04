options( echo=TRUE )
## open MadPro library
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 7){
	stop(paste(length(args), " parameters must be 7"),sep="")
}
sample_file <- 	args[1]
if(!file.exists(sample_file)){
	stop(paste(sample_file, "not found"), sep = "")
}
filterParameter <- args[2]
matrix_file <- args[3]
if(!file.exists(matrix_file)){
	stop(paste(matrix_file, "not found"), sep = "")
}
pourcent <- as.integer(args[4])

logged <- args[7]

output_file <- args[5]

output_img <- args[6]

library(limma,quietly=TRUE)
library(genefilter,quietly=TRUE)

pData = read.delim(sample_file, header = TRUE, row.names = 1, sep = "\t")
matrix = read.delim(matrix_file, header = TRUE, row.names = 1, sep = "\t")
if(logged){
	matrix <- log(matrix)
}
seuil = mean(apply(matrix, 2, median))
n = pourcent/100
classes.matrix = list()

if(filterParameter != "global"){
	factors = pData[,which(colnames(pData) == filterParameter)]
	factors = as.factor(unlist(factors))
	for(i in 1:nlevels(factors)){
		classes.matrix = c(classes.matrix, list(matrix[, which(factors == levels(factors)[i])]))
	}
}else{
	factors <- as.factor(unlist(rep("1", nrow(pData))))
	classes.matrix = list(matrix)
}



names(classes.matrix) = levels(factors)
current_filter = rep(FALSE, dim(classes.matrix[[1]])[1])
	

for(i in 1:nlevels(factors)){
	filtre = kOverA(dim(classes.matrix[[i]])[2]*n, seuil)
	result_filter = genefilter(classes.matrix[[i]],filtre)
	current_filter= result_filter | current_filter
}

png(output_img)
plot(density(as.numeric(unlist(matrix))))
abline(v=seuil)
dev.off()


matrix = matrix[which(current_filter != FALSE), ]

write.table(matrix, output_file, row.names = TRUE, col.names = NA, quote = FALSE, sep="\t") 




