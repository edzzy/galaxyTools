suppressMessages(library(gplots, quietly = TRUE))
options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
input_matrix <- args[1]
input_sample <- args[2]
logged <- args[3]
print(logged)
output_name <- args[4]
parameter <- args[5]

data <- read.delim(input_matrix, header = TRUE, row.names = 1)
samples <-read.delim(input_sample, header = TRUE, row.names = 1) 
data <- data[, rownames(samples)]
if(logged == "True"){
	data <- log(data)
}


class <- as.factor(unlist(samples[, which(colnames(samples) == parameter)]))


colclass<-c()
nSamples <- ncol(data)
width<- 5*500
height<-5*500
res = 500
pointsize = 8
bg_color <-"white"
bottomarg = nchar(max(colnames(data))) 
uniqueCol<-rainbow(length(unique(class)))
names(uniqueCol) <- unique(class)
for (a in 1:length(uniqueCol)){
	for (b in 1:length(class)){
		if(class[b] == names(uniqueCol)[a] ){
			 colclass[b]<-uniqueCol[a]
		}
	}
}
png(filename = output_name, width = width, height = height,  bg = bg_color, res =res, pointsize = pointsize)
par(mar=c(bottomarg + 5,5,3,3))
mat_data <- cor(data)
heatmap.2(mat_data, notecol="black",margins=c(12,9),RowSideColors=colclass, density.info="none",trace="none",col=redgreen(75), dendrogram="row", key = TRUE, scale ="none" )
par(lend = 1)
legend("top",
legend = unique(class),
col = uniqueCol,
lty = 1,
lwd = 10
 )
devname= dev.off()


