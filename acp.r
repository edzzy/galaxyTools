#options( show.error.messages=FALSE,error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
options(echo=TRUE)
library(ade4,quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)
library( RSVGTipsDevice,quietly=TRUE)
library(R2HTML)
args <- commandArgs(trailingOnly = TRUE)
input_matrix <- args[1]
input_sample <- args[2]
input_parameter <- args[3]
logged <- args[4]
output_img <- args[5]
print(output_img)

mat = read.delim( input_matrix, row.names = 1 )
samples <- read.delim(input_sample, row.name = 1)

if(logged){
	mat = log(mat)
}
class = samples[,which(colnames(samples) == input_parameter)]
names(class)<-row.names(samples)
class <- class[colnames(mat)]
classLevels <- levels(as.factor(as.character(class)))
colclass<-c()
colclass <- class
nf <- 3 
print(nf)

acp = dudi.pca(t(mat),scannf=FALSE,nf=nf)
#bg_color ="white"
#png(filename = output_img, width = width, height = height, bg = bg_color, res = NA)
#plot(as.data.frame(acp[which(names(acp)=="li")])[,1], as.data.frame(acp[which(names(acp)=="li")])[,2], pch = 19)
#dev.off()
cex<-3
pch<-20


col <- rainbow(length(classLevels))  
for (i in 1:length(classLevels)){
   	colclass <- gsub(pattern = paste("^",classLevels[i],"$", sep=""), replacement = col[i], x = as.character(colclass), perl=TRUE)
}

directory = getwd()
#htmlfile <- HTMLInitFile(outdir = getwd(), filename = output_img, extension="html")
#devSVGTips(output_img, toolTipMode = 1, title= paste ("ACP",input_parameter, sep=""), width=width, height=height)
#for(i in 1:(nf-1)){
#  for(j in (i+1):nf){
png(output_img, width = 1500, height = 750)
par(mfrow = c(1, 2))

s.class(dfxy = acp$li, fac = as.factor(class), col = rainbow(n = nlevels(as.factor(class))), xax = 1, yax = 2)
s.label(acp$li, xax =1, yax = 2)
dev.off()
#plot(c(min(acp$li),max(acp$li)),c(min(acp$li),max(acp$li)),type="n",xlab=paste("Axe ", i,  sep=""),ylab = paste("Axe ",j, sep=""),main = paste ("ACP ",input_parameter, " AXE : ", i, "-", j,sep=""))
#invisible(sapply(1:nrow(acp$li),function(x) {
#		 setSVGShapeToolTip(title=rownames(acp$li)[x], desc=class[x] );
#		 points(acp$li[x,i],acp$li[x,j],cex=cex,pch=pch,col=colclass[x])

#}))
#}
#}
 # plot(c(min(acp$li),max(acp$li)),c(min(acp$li),max(acp$li)),type="n",xlab="",ylab="",main = paste ("Legende ",input_parameter,sep=""), xaxt="n", yaxt="n")
 # legend("center", legend = as.character(classLevels), col=col, pch=20, cex=3)
#dev.off()
#HTMLInsertGraph("acp.svg", file=htmlfile, caption="ACP Axe 1")
#HTMLEndFile()

#save(acp, file="ACP.RData")


