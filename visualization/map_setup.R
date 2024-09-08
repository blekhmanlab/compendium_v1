##############################################################################
# SETUP FOR MAP on OSX
### Install udunits from command line
# brew install udunits
### Then find where it is:
# brew --prefix udunits
### Then include the path to that dir in these two places:
UDUNITS_PARAMS <- '--with-udunits2-include=/opt/homebrew/opt/udunits/include --with-udunits2-lib=/opt/homebrew/opt/udunits/lib'

install.packages('udunits2', type = 'source', repo = 'cran.rstudio.com',
                 configure.args=UDUNITS_PARAMS)
install.packages('devtools')
library(devtools)
devtools::install_github('edzer/units', type = 'source', configure.args=UDUNITS_PARAMS)
# brew install gdal
install.packages('sf')
install.packages('rnaturalearth')
install.packages("rnaturalearthdata")
install.packages('rgeos')
##############################################################################
