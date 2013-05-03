  options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
     ## open MadPro library
    library(MadProDev)
      mat = read.delim( "matRaw.txt",row.names=1 )
  	type="L"
  	if(type=="Q"){
  		mat=normQuantile(mat=mat,pngDir="./",img=FALSE)	
  	}
  	if(type=="L"){
  		mat = LOWESS(data=mat,graph=0,projet="",pngDir="./")
  	}
    write.table(mat,"testNorm.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)
