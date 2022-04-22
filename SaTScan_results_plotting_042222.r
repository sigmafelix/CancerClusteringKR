## SaTScan results plotting
library(pacman)
p_load(tidyverse, sf, nngeo)

source("base_functions.R")

# tmap_satscan = function(basemap, sats, threshold = 2, significance = 0.01, alpha = 0.4, return_ellipses = TRUE) {
#     rotate = function(a) {a = a*pi/180; matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)}

p1_res = read_rds("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_ASMR_uncontrolled.rds")
p1_labels = as.data.frame(str_split_fixed(p1_res$analysis_title, "_", 6)[,-6])
colnames(p1_labels) = c("cancertype", "period", "measure", "target", "sex")
p1_res = p1_res %>%
        bind_cols(p1_labels, .) %>%
        filter(pvalue <= 0.01 & number_locs >= 2)

p1_resl = p1_res %>% split(., .$analysis_title)
names(p1_resl)


tmap_satscan(covar_origin_00_fc, p1_resl[[7]])
tmap_satscan(covar_origin_00_fc, p1_resl[[13]])


covar_origin_10_fc = covar_origin_10_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430))
covar_origin_05_fc = covar_origin_05_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430))
covar_origin_00_fc = covar_origin_00_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430))


map_lisa(covar_origin_10_fc, "ragest_i_Lung_total_3")
map_lisa(covar_origin_10_fc, "ragest_i_Lung_male_3")
map_lisa(covar_origin_10_fc, "ragest_i_Lung_female_3")
map_lisa(covar_origin_10_fc, "ragest_i_Lung_total_3")
