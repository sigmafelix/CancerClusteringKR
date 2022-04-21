

username = 'sigma'

source('./base_functions.R')
source('./Cleaning_Population_011922.r')


## FIT MAIN ####
# Naming: smerc_[cancertype]_[period][incidence/mortality][sex][variable set/uncontrolled]
# minimal sociodemographic (v2): str_c('(^(p_65p|p_hbac)*.*_', sex_b, '$')
# intersection across all periods (v2): str_c(str_c('(^p_*.*_', sex_b, '$'), '^ap_', '^NDVI_)', sep = '|')
# period 2 and 3 (v3): str_c(str_c('(^p_*.*_', sex_b, '$'), '^r_', '^ap_', '^NDVI_)', sep = '|')
# period 3 only (v4): str_c(str_c('(^p_*.*_', sex_b, '$'), '^ap_', '^NDVI_', '^n_pw)', sep = '|')
load(str_c(drive, "Manuscript/Clustering_Base_sf_021722.RData"))


covar_origin_10_fc = covar_origin_10_fc %>%
    dplyr::select(-colnames(.)[grep(str_c("^(", str_c(colnames(morinc_clw)[-1], collapse = "|"), ")"), colnames(.))]) %>%
    left_join(morinc_clw %>% dplyr::select(1, ends_with("_3")))
covar_origin_05_fc = covar_origin_05_fc %>%
    dplyr::select(-colnames(.)[grep(str_c("^(", str_c(colnames(morinc_clw)[-1], collapse = "|"), ")"), colnames(.))]) %>%
    left_join(morinc_clw %>% dplyr::select(1, ends_with("_2")))
covar_origin_00_fc = covar_origin_00_fc %>%
    dplyr::select(-colnames(.)[grep(str_c("^(", str_c(colnames(morinc_clw)[-1], collapse = "|"), ")"), colnames(.))]) %>%
    left_join(morinc_clw %>% dplyr::select(1, ends_with("_1")))
save(list = ls()[grep('^covar_origin_', ls())], file = str_c(drive, "Manuscript/Clustering_Base_sf_042022.RData"), compress = 'xz', compression_level = 9)

