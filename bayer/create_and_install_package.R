# The first time devtools are installed, some programs may be missing
# the following command will install the missing programs:
# sudo apt-get install libcurl4-openssl-dev
#
# before installing the package, clear environment and restart R,
# otherwise there may be unexpected error messages
library(devtools)
load_all()
document()
devtools::check(document=F)
install(build_vignettes = FALSE)
library(crmPack)
#stats::app
#install.packages("rlang", dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))
update_packages()