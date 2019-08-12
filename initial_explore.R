# supplied
# data files in "../xhpi_data/" are in three directories, one for each data-set. 
# "Neu_7_reports" neutron diffraction
# "Res_7_reports" high resolution
# "Temp_7_reports" high temperature 
#  .report files are the plain text files containing information about the XH/pi interactions. 
#  Each protein has 6 files
#  "7a3h_H.pdb"
#  "7a3h_H.pdb.backbone" 
#  "7a3h_H.pdb.backbone.pml"
#  "7a3h_H.pdb.pml"
#  "7a3h_H.pdb.report"  
#  "7a3h_H.pdb_A.pml"
#  
# Geometric constraints for this data was d<7 and theta<50, although cut-offs of
# d<4.3 angstroms and theta<25 deg are often used to define XH/pi interactions, 
# so many of the interactions reported are not actually XH/pi interactions.

# pkgs
library(tidyverse)

# "Res_7_reports" has 
dir("../xhpi_data/Res_7_reports/") %>% length() 
# 994 files and 7 of these accessory

cols <- c("X_res", "X_atom", "PDB FILE", "Resolution", "Aromatic", "XDist", "XTheta")

highres <- read.table("../xhpi_data/Res_7_reports/", 
                      header = FALSE, 
                      skip = 3,
                      sep = "\t",
                      col.names = cols,
                      stringsAsFactors = FALSE)
