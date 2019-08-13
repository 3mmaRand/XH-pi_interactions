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
library(data.table)

# "Res_7_reports" has 
dir("../xhpi_data/Res_7_reports/") %>% length() 
# 994 files and 7 of these accessory

cols <- c("X_res_type",
          "X_res_pos",
          "X_atom_type",
          "X_atom_number",
          "PDB FILE",
          "Resolution",
          "Aromatic_res",
          "Aromatic_number",
          "XDist",
          "XTheta")
# reading in one file to check
highres <- read.table("../xhpi_data/Res_7_reports/1eb6_H.pdb.report", 
                      header = FALSE, 
                      skip = 3,
                      col.names = cols,
                      stringsAsFactors = FALSE)

#################################
# IMPORT HIGH RES
#################################
directory <- "../xhpi_data/Res_7_reports"
reports <- list.files(path = directory, pattern = "*.pdb.report")
number <- length(reports)

highres <- do.call("rbind", 
                   lapply(paste0(directory,"/", reports[1:number]),
                          function(x) read.table(x, header = FALSE,
                                                 skip = 3,
                                                 col.names = cols,
                                                 stringsAsFactors = FALSE)))

# this results in 5142 interactions
# # 
# # a subset with the more stringent criteria of  d < 4.3 and theta < 25
# 
# highres2 <- highres %>% 
#   subset(XDist =< 4.3 ) %>%
#   subset(XTheta =< 25)
# # gives 5104 interactions


#################################
# IMPORT NEUTRON DIFFRACTION
##############################
directory <- "../xhpi_data/Neu_7_reports/"
reports <- list.files(path = directory, pattern = "*.pdb.report")
number <- length(reports)

neutron <- do.call("rbind", 
                   lapply(paste0(directory,"/", reports[1:number]),
                          function(x) read.table(x, header = FALSE,
                                                 skip = 3,
                                                 col.names = cols,
                                                 stringsAsFactors = FALSE)))

# this results in 791 interactions
# 
# a subset with the more stringent criteria of  d < 4.3 and theta < 25
# 
# neutron2 <- neutron %>% subset(XDist =< 4.3 & XTheta =< 25)
# # gives 786 interactions

#################################
# IMPORT HIGH TEMPERATURE
##############################
directory <- "../xhpi_data/Temp_7_reports/"
reports <- list.files(path = directory, pattern = "*.pdb.report")
number <- length(reports)

hightemp <- do.call("rbind", 
                   lapply(paste0(directory,"/", reports[1:number]),
                          function(x) read.table(x, header = FALSE,
                                                 skip = 3,
                                                 col.names = cols,
                                                 stringsAsFactors = FALSE)))

# this results in 2108 interactions
# 
# # a subset with the more stringent criteria of  d < 4.3 and theta < 25
# 
# hightemp2 <- hightemp %>% subset(XDist < 4.3 & XTheta < 25)
# # gives 2098 interactions

###############################
# Interaction types
###############################

# crosstabulate donor type and acceptor types
# highres
highressum <- highres %>% 
  group_by(X_atom_type, Aromatic_res) %>% 
  count() 

ggplot(highressum, aes(x = X_atom_type, y = Aromatic_res, fill = n)) +
  geom_tile()+
  ggtitle("High res")

# highres
neutronsum <- neutron %>% 
  group_by(X_atom_type, Aromatic_res) %>% 
  count()

ggplot(neutronsum, aes(x = X_atom_type, y = Aromatic_res, fill = n)) +
  geom_tile()+
  ggtitle("Neutron")

# hightemp
hightempsum <- hightemp %>% 
  group_by(X_atom_type, Aromatic_res) %>% 
  count()

ggplot(hightempsum, aes(x = X_atom_type, y = Aromatic_res, fill = n)) +
  geom_tile()+
  ggtitle("High temp")




###################################
# DONOR CLASSIFYING
###################################



