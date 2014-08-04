options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
input_matrix <- args[1]
logged <- args[2]
output_img <- args[3]

mat = read.delim( input_matrix, row.names = 1 )
if(logged == "True"){
	mat <- log(mat)
}


bottomarg = nchar(max(colnames(mat))) #nombre de ligne pour la marge du bas
jpeg(filename = output_img, width = 1300, height = 900, quality = 100, bg = "white", res = NA)
par(mar=c(bottomarg + 5,5,3,3))
boxplot(as.data.frame(mat), col="darkolivegreen2", las="2",cex.lab=1.5,cex.axis=1.5)
devname= dev.off()
