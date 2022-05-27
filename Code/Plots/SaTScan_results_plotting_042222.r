## SaTScan results plotting
library(pacman)
p_load(tidyverse, sf, nngeo, spdep, tmap)

source("base_functions.R")
username = "sigma"
rdsdir = sprintf("/mnt/c/Users/%s/OneDrive/NCC_Project/CancerClustering/", username)

# tmap_satscan = function(basemap, sats, threshold = 2, significance = 0.01, alpha = 0.4, return_ellipses = TRUE) {
#     rotate = function(a) {a = a*pi/180; matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)}

clear_input = function(filepath, alpha = 0.01, locmin = 2) {
    # alpha: significance level
    # locmin: how many points should be included at least?
    res = read_rds(filepath)
    labels = as.data.frame(str_split_fixed(res$analysis_title, "_", 6))
    colnames(labels) = c("cancertype", "period", "measure", "target", "sex", "vset")
    res = res %>%
          bind_cols(labels, .) %>%
          filter(pvalue <= alpha & number_locs >= locmin)
    resl = res %>% split(., .$analysis_title)
    return(resl)
}

p_v0_resl = clear_input(filepath = str_c(rdsdir, "satscan_ASMR_uncontrolled.rds"))
p_v1_resl = clear_input(filepath = str_c(rdsdir, "satscan_ASMR_vset1.rds"))
p_v2_resl = clear_input(filepath = str_c(rdsdir, "satscan_ASMR_vset2.rds"))
p_v3_resl = clear_input(filepath = str_c(rdsdir, "satscan_ASMR_vset3.rds"))
p_v4_resl = clear_input(filepath = str_c(rdsdir, "satscan_ASMR_vset4.rds"))



p1_res = read_rds("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_ASMR_uncontrolled.rds")
p1_labels = as.data.frame(str_split_fixed(p1_res$analysis_title, "_", 6))
colnames(p1_labels) = c("cancertype", "period", "measure", "target", "sex", "vset")
p1_res = p1_res %>%
        bind_cols(p1_labels, .) %>%
        filter(pvalue <= 0.01 & number_locs >= 2)

p1_resl = p1_res %>% split(., .$analysis_title)

# names(p_v1_resl)


tmap_satscan(covar_origin_00_fc, p_v1_resl[[7]])
tmap_satscan(covar_origin_00_fc, p_v1_resl[[13]])


covar_origin_10_fcc = covar_origin_10_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430)) %>%
    rmapshaper::ms_simplify(keep = 0.2, keep_shapes = TRUE)
covar_origin_05_fcc = covar_origin_05_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430)) %>%
    rmapshaper::ms_simplify(keep = 0.2, keep_shapes = TRUE)
covar_origin_00_fcc = covar_origin_00_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430)) %>%
    rmapshaper::ms_simplify(keep = 0.2, keep_shapes = TRUE)


map_lisa(covar_origin_10_fcc, "ragest_i_Lung_total_3")
map_lisa(covar_origin_10_fcc, "ragest_i_Lung_male_3")
map_lisa(covar_origin_10_fcc, "ragest_i_Lung_female_3")

map_lisa(covar_origin_00_fcc, "ragest_i_Lung_total_1")
map_lisa(covar_origin_00_fcc, "ragest_i_Lung_male_1")
map_lisa(covar_origin_00_fcc, "ragest_i_Lung_female_1")

map_lisa(covar_origin_10_fcc, "ragest_i_Lung_total_3")
