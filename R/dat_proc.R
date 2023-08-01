library(tidyverse)
library(here)

# apasumdat ---------------------------------------------------------------

fl <- here('data/apasumdat.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apasumdat.RData', destfile = fl)

# misdat, missing data in apacp 2021 ----------------------------------------------------------

fl <- here('data/misdat.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/misdat.RData', destfile = fl)
