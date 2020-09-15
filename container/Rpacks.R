# Function taken at https://stackoverflow.com/a/19870272/12519542
install_if_absent <- function(x){
	for(i in x){
		#  require returns TRUE invisibly if it was able to load package
		if(! require(i , character.only = TRUE)){
			#  If package was not able to be loaded then re-install
			install.packages(i , repos="https://cran.ma.imperial.ac.uk/")
			#  Load package after installing
			require(i , character.only = TRUE)
		}
	}
}

required <- c("argparser", "UpSetR", "tibble", "dplyr", "readr", "purrr", "tidyr", "ggplot2", "rjson")
install_if_absent(required)
