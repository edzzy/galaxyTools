  options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
     ## open MadPro library
	source("/home/ehirchaud/MadProDev/R/graph_png.R")
	mat = read.delim( "testNorm.txt",row.names=1 )
	create_boxplot(mat,"boxplot.jpg")	
