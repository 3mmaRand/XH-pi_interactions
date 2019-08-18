# supplied
# data files in "../xhpi_data_new/" are in three directories, one for each data-set. 
# "Neu_43_adv" neutron diffraction
# "Res_43_adv" high resolution
# "Temp_43_adv" high temperature 
#  .report files are the plain text files containing information about the XH/pi interactions. 
#  Each protein has 3 files
#  "7a3h.pdb"
#  "7a3h.hpml" 
#  "7a3h.pdb.report"  
  
# Geometric constraints for this data was d<7 and theta<50, although cut-offs of
# d<4.3 angstroms and theta<25 deg are often used to define XH/pi interactions, 
# so many of the interactions reported are not actually XH/pi interactions.

# pkgs
library(tidyverse)
library(data.table)

# "Res_43_adv" has 
dir("../xhpi_data_new/Res_43_adv/") %>% length() 
# 493 files and 1 of these accessory, i.e., 164 proteins


# reading in one file to check
neutron <- read.table("../xhpi_data_new/Neu_43_adv/6eym.pdb.report",
                      header = TRUE,
                      stringsAsFactors = FALSE)

#################################
# IMPORT HIGH RES
#################################
directory <- "../xhpi_data_new/Res_43_adv/"
reports <- list.files(path = directory, pattern = "*.pdb.report")
number <- length(reports)

highres <- do.call("rbind", 
                   lapply(paste0(directory,"/", reports[1:number]),
                          function(x) read.table(x, 
                                                 header = TRUE,
                                                 stringsAsFactors = FALSE)))

# this results in 4305 interactions


#################################
# IMPORT NEUTRON DIFFRACTION
##############################
directory <- "../xhpi_data_new/Neu_43_adv/"
reports <- list.files(path = directory, pattern = "*.pdb.report")
number <- length(reports)

neutron <- do.call("rbind", 
                   lapply(paste0(directory,"/", reports[1:number]),
                          function(x) read.table(x, 
                                                 header = TRUE,
                                                 stringsAsFactors = FALSE)))

# this results in 772 interactions

#################################
# IMPORT HIGH TEMPERATURE
##############################
directory <- "../xhpi_data_new/Temp_43_adv/"
reports <- list.files(path = directory, pattern = "*.pdb.report")
number <- length(reports)

hightemp <- do.call("rbind", 
                   lapply(paste0(directory,"/", reports[1:number]),
                          function(x) read.table(x, 
                                                 header = TRUE,
                                                 stringsAsFactors = FALSE)))

# this results in 1788 interactions

###############################
# Interaction types
###############################

# crosstabulate donor type and acceptor types
# highres
highressum <- highres %>% 
  group_by(x_atom_id, pi_res_id) %>% 
  count() 

ggplot(highressum, aes(x = x_atom_id, y = pi_res_id, fill = n)) +
  geom_tile()+
  ggtitle("High res")

# highres
neutronsum <- neutron %>% 
  group_by(x_atom_id, pi_res_id) %>% 
  count()

ggplot(neutronsum, aes(x = x_atom_id, y = pi_res_id, fill = n)) +
  geom_tile()+
  ggtitle("Neutron")

# hightemp
hightempsum <- hightemp %>% 
  group_by(x_atom_id, pi_res_id) %>% 
  count()

ggplot(hightempsum, aes(x = x_atom_id, y = pi_res_id, fill = n)) +
  geom_tile()+
  ggtitle("High temp")




###################################
# DONOR CLASSIFYING
###################################

###################################
# CLUSTERING MIXED DATA TYPES
###################################
# https://dpmartin42.github.io/posts/r/cluster-mixed-types
library(cluster)
library(Rtsne)
# distance calculation - gower distance
# clustering algorithm - partitioning aroung medoids
# selecting the number of clusters - silhouette width

glimpse(neutron)
# testing on neutron because it has resolution data (others have xxxx)

# some variables will not used to cluster
neutron2 <- neutron %>% 
  select(-x_res_num, 
         -x_atom_num, 
         -pdb, 
         -pi_res_num,
         -resolution)

# char need to be converted to factors
neutron2 <- neutron2 %>% mutate_if(is.character,as.factor)
glimpse(neutron2)

# calc distances without any transformations in first instance
gower_dist <- daisy(neutron2,
                    metric = "gower")
summary(gower_dist)

# which observations are least and most similar
# most similar i.e., min disimilarity
gower_mat <- as.matrix(gower_dist)
neutron2[which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
               arr.ind = TRUE)[1, ], ]

# least similar i.e., max disimilarity
neutron2[which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
               arr.ind = TRUE)[1, ], ]

# cluster with PAM partitioning round medoids

# selecting the number of clusters with silhouette width
# calculating silhouette width for clusters ranging from 2 to cnum 
# for the PAM algorithm, we see that 3 clusters yields the highest value.
sil_width <- c(NA)
cnum <- 40

for(i in 2:cnum){
  
  pam_fit <- pam(gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

# Plot sihouette width (higher is better)

plot(1:cnum, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:cnum, sil_width)

# vis
tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering))

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))

# consider nature of the variables before clustering

glimpse(neutron2)

hist(neutron2$x_width)
# log the bfactors
gower_dist <- daisy(neutron2,
                    metric = "gower",
                    type = list(logratio = 4,
                                logratio = 8))

# which observations are least and most similar
# most similar i.e., min disimilarity
gower_mat <- as.matrix(gower_dist)
neutron2[which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
               arr.ind = TRUE)[1, ], ]

# least similar i.e., max disimilarity
neutron2[which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
               arr.ind = TRUE)[1, ], ]

# cluster with PAM partitioning round medoids

# selecting the number of clusters with silhouette width
# calculating silhouette width for clusters ranging from 2 to cnum 
# for the PAM algorithm, we see that 3 clusters yields the highest value.
sil_width <- c(NA)
cnum <- 40

for(i in 2:cnum){
  
  pam_fit <- pam(gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

# Plot sihouette width (higher is better)

plot(1:cnum, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:cnum, sil_width)

