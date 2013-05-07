  options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
     ## open MadPro library
	suppressPackageStartupMessages(library(MadPro,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE))	
