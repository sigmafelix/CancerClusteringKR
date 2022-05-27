### Run Clustering Analysis - normal
### 04/20/22
username = 'sigma'

source('./base_functions.R')
source('./Cleaning_Population_011922.r')


# knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, error = FALSE, fig.height = 7.5)

covar_origin_10 = get_basecovar(target_year = 2010, target_year_kcdc = 2010)
covar_origin_05 = get_basecovar(target_year = 2005, target_year_kcdc = 2008)
covar_origin_00 = get_basecovar(target_year = 2000, target_year_kcdc = 2000)

covar_origin_10_fc = clean_consolidated(cleaned_df = covar_origin_10) %>%
    mutate_at(.vars = vars(-geom),
              .funs = list(~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    left_join(morinc_clw, by = c('sgg_cd_c'))
covar_origin_05_fc = clean_consolidated(cleaned_df = covar_origin_05, target_year = 2005) %>%
    mutate_at(.vars = vars(-geom),
              .funs = list(~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    left_join(morinc_clw, by = c('sgg_cd_c'))
covar_origin_00_fc = clean_consolidated(cleaned_df = covar_origin_00, target_year = 2000) %>%
    mutate_at(.vars = vars(-geom),
              .funs = list(~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    left_join(morinc_clw, by = c('sgg_cd_c'))


## FIT MAIN ####
# Naming: smerc_[cancertype]_[period][incidence/mortality][sex][variable set/uncontrolled]
# minimal sociodemographic (v2): str_c('(^(p_65p|p_hbac)*.*_', sex_b, '$')
# intersection across all periods (v2): str_c(str_c('(^p_*.*_', sex_b, '$'), '^ap_', '^NDVI_)', sep = '|')
# period 2 and 3 (v3): str_c(str_c('(^p_*.*_', sex_b, '$'), '^r_', '^ap_', '^NDVI_)', sep = '|')
# period 3 only (v4): str_c(str_c('(^p_*.*_', sex_b, '$'), '^ap_', '^NDVI_', '^n_pw)', sep = '|')
load(str_c(drive, "Manuscript/Clustering_Base_sf_042022.RData"))
NTHREADS = 14L

## Excluding Ongjin 23320 Ulleung 37430
covar_origin_10_fc = covar_origin_10_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430))
covar_origin_05_fc = covar_origin_05_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430))
covar_origin_00_fc = covar_origin_00_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430))



# Available sociodemographic variables (2) for all periods ####
# Incidence (i)
sex_bb = 'total'
vset1 = str_c('^(p_65p|p_hbac)*.*_', sex_bb, '$')
smerc_lung_3it_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', type = 'normal', yvar = 'n_i_Lung_total_3', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_3it_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', type = 'normal', yvar = 'n_i_Stomach_total_3', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_2it_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', type = 'normal', yvar = 'n_i_Lung_total_2', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_2it_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', type = 'normal', yvar = 'n_i_Stomach_total_2', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_1it_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', type = 'normal', yvar = 'n_i_Lung_total_1', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_1it_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', type = 'normal', yvar = 'n_i_Stomach_total_1', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)

sex_bb = 'male'
vset1 = str_c('^(p_65p|p_hbac)*.*_', sex_bb, '$')
smerc_lung_3im_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', type = 'normal', yvar = 'n_i_Lung_male_3', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_3im_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', type = 'normal', yvar = 'n_i_Stomach_male_3', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_2im_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', type = 'normal', yvar = 'n_i_Lung_male_2', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_2im_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', type = 'normal', yvar = 'n_i_Stomach_male_2', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_1im_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', type = 'normal', yvar = 'n_i_Lung_male_1', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_1im_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', type = 'normal', yvar = 'n_i_Stomach_male_1', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)

sex_bb = 'female'
vset1 = str_c('^(p_65p|p_hbac)*.*_', sex_bb, '$')
smerc_lung_3if_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', type = 'normal', yvar = 'n_i_Lung_female_3', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_3if_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', type = 'normal', yvar = 'n_i_Stomach_female_3', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_2if_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', type = 'normal', yvar = 'n_i_Lung_female_2', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_2if_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', type = 'normal', yvar = 'n_i_Stomach_female_2', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_1if_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', type = 'normal', yvar = 'n_i_Lung_female_1', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_1if_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', type = 'normal', yvar = 'n_i_Stomach_female_1', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)

# unadjusted
smerc_lung_3itu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', type = 'normal', yvar = 'ragest_i_Lung_total_3', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_stom_3itu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', type = 'normal', yvar = 'ragest_i_Stomach_total_3', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_lung_2itu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', type = 'normal', yvar = 'ragest_i_Lung_total_2', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_stom_2itu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', type = 'normal', yvar = 'ragest_i_Stomach_total_2', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_lung_1itu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', type = 'normal', yvar = 'ragest_i_Lung_total_1', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_stom_1itu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', type = 'normal', yvar = 'ragest_i_Stomach_total_1', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)

smerc_lung_3imu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', type = 'normal', yvar = 'ragest_i_Lung_male_3', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_stom_3imu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', type = 'normal', yvar = 'ragest_i_Stomach_male_3', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_lung_2imu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', type = 'normal', yvar = 'ragest_i_Lung_male_2', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_stom_2imu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', type = 'normal', yvar = 'ragest_i_Stomach_male_2', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_lung_1imu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', type = 'normal', yvar = 'ragest_i_Lung_male_1', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_stom_1imu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', type = 'normal', yvar = 'ragest_i_Stomach_male_1', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)

smerc_lung_3ifu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', type = 'normal', yvar = 'ragest_i_Lung_female_3', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_stom_3ifu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', type = 'normal', yvar = 'ragest_i_Stomach_female_3', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_lung_2ifu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', type = 'normal', yvar = 'ragest_i_Lung_female_2', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_stom_2ifu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', type = 'normal', yvar = 'ragest_i_Stomach_female_2', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_lung_1ifu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', type = 'normal', yvar = 'ragest_i_Lung_female_1', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_stom_1ifu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', type = 'normal', yvar = 'ragest_i_Stomach_female_1', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)


# par(mfcol = c(1,2))
# plot(smerc_lung_3dt)
# plot(smerc_lung_3dtu)


# Mortality (d)
sex_bb = 'total'
vset1 = str_c('^(p_65p|p_hbac)*.*_', sex_bb, '$')
smerc_lung_3dt_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Lung_total_3', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_3dt_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Stomach_total_3', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_2dt_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_d_Lung_total_2', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_2dt_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_d_Stomach_total_2', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_1dt_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', yvar = 'n_d_Lung_total_1', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_1dt_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', yvar = 'n_d_Stomach_total_1', sex_b = 'total', adjust = TRUE, string_search = vset1, ncores = NTHREADS)

sex_bb = 'male'
vset1 = str_c('^(p_65p|p_hbac)*.*_', sex_bb, '$')
smerc_lung_3dm_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_3dm_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Stomach_male_3', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_2dm_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_d_Lung_male_2', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_2dm_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_d_Stomach_male_2', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_1dm_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', yvar = 'n_d_Lung_male_1', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_1dm_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', yvar = 'n_d_Stomach_male_1', sex_b = 'male', adjust = TRUE, string_search = vset1, ncores = NTHREADS)

sex_bb = 'female'
vset1 = str_c('^(p_65p|p_hbac)*.*_', sex_bb, '$')
smerc_lung_3df_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Lung_female_3', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_3df_v1 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Stomach_female_3', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_2df_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_d_Lung_female_2', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_2df_v1 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_d_Stomach_female_2', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_lung_1df_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', yvar = 'n_d_Lung_female_1', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)
smerc_stom_1df_v1 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', yvar = 'n_d_Stomach_female_1', sex_b = 'female', adjust = TRUE, string_search = vset1, ncores = NTHREADS)

# unadjusted
smerc_lung_3dtu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Lung_total_3', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_stom_3dtu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Stomach_total_3', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_lung_2dtu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_d_Lung_total_2', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_stom_2dtu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_d_Stomach_total_2', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_lung_1dtu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', yvar = 'n_d_Lung_total_1', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)
smerc_stom_1dtu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', yvar = 'n_d_Stomach_total_1', sex_b = 'total', adjust = FALSE, ncores = NTHREADS)

smerc_lung_3dmu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_stom_3dmu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Stomach_male_3', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_lung_2dmu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_d_Lung_male_2', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_stom_2dmu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_d_Stomach_male_2', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_lung_1dmu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', yvar = 'n_d_Lung_male_1', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)
smerc_stom_1dmu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', yvar = 'n_d_Stomach_male_1', sex_b = 'male', adjust = FALSE, ncores = NTHREADS)

smerc_lung_3dfu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Lung_female_3', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_stom_3dfu = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Stomach_female_3', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_lung_2dfu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_d_Lung_female_2', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_stom_2dfu = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_d_Stomach_female_2', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_lung_1dfu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', yvar = 'n_d_Lung_female_1', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)
smerc_stom_1dfu = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', yvar = 'n_d_Stomach_female_1', sex_b = 'female', adjust = FALSE, ncores = NTHREADS)





# Available sociodemographic and environmental variables for all periods ####
# Incidence (i)
NTHREADS = 28
sex_bb = 'total'
vset2 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^ap_', '^NDVI_', sep = '|')
smerc_lung_3it_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Lung_total_3', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_3it_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Stomach_total_3', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_2it_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_i_Lung_total_2', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_2it_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_i_Stomach_total_2', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_1it_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', yvar = 'n_i_Lung_total_1', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_1it_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', yvar = 'n_i_Stomach_total_1', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)

sex_bb = 'male'
vset2 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^ap_', '^NDVI_', sep = '|')
smerc_lung_3im_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Lung_male_3', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_3im_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Stomach_male_3', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_2im_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_i_Lung_male_2', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_2im_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_i_Stomach_male_2', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_1im_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', yvar = 'n_i_Lung_male_1', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_1im_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', yvar = 'n_i_Stomach_male_1', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)

sex_bb = 'female'
vset2 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^ap_', '^NDVI_', sep = '|')
smerc_lung_3if_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Lung_female_3', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_3if_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Stomach_female_3', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_2if_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_i_Lung_female_2', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_2if_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_i_Stomach_female_2', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_1if_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', yvar = 'n_i_Lung_female_1', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_1if_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', yvar = 'n_i_Stomach_female_1', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)

# par(mfcol = c(1,2))
# plot(smerc_lung_3dt)
# plot(smerc_lung_3dtu)


# Mortality (d)
sex_bb = 'total'
vset2 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^ap_', '^NDVI_', sep = '|')
smerc_lung_3dt_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Lung_total_3', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_3dt_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Stomach_total_3', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_2dt_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_d_Lung_total_2', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_2dt_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_d_Stomach_total_2', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_1dt_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', yvar = 'n_d_Lung_total_1', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_1dt_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_total_1', yvar = 'n_d_Stomach_total_1', sex_b = 'total', adjust = TRUE, string_search = vset2, ncores = NTHREADS)

sex_bb = 'male'
vset2 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^ap_', '^NDVI_', sep = '|')
smerc_lung_3dm_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_3dm_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Stomach_male_3', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_2dm_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_d_Lung_male_2', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_2dm_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_d_Stomach_male_2', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_1dm_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', yvar = 'n_d_Lung_male_1', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_1dm_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_male_1', yvar = 'n_d_Stomach_male_1', sex_b = 'male', adjust = TRUE, string_search = vset2, ncores = NTHREADS)

sex_bb = 'female'
vset2 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^ap_', '^NDVI_', sep = '|')
smerc_lung_3df_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Lung_female_3', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_3df_v2 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Stomach_female_3', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_2df_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_d_Lung_female_2', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_2df_v2 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_d_Stomach_female_2', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_lung_1df_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', yvar = 'n_d_Lung_female_1', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)
smerc_stom_1df_v2 = run_smerc_cancertype(data = covar_origin_00_fc, population = 'n_p_female_1', yvar = 'n_d_Stomach_female_1', sex_b = 'female', adjust = TRUE, string_search = vset2, ncores = NTHREADS)


#save(list = ls()[grep('^smerc_', ls())],
#     file = str_c(dbdir, 'Manuscript/Scan_SMERC_periods_1_3_vsets_1_2_Results_p005.RData'))
# save(list = ls()[grep('^dclust_*.*_[3]*.*', ls())],
#      file = 'Scan_DCLUST_2009_2013_Results.RData')
# save(list = ls()[grep('^(smerc|dclust|covar_origin_10_)', ls())],
#      file = 'Scan_2009_2013_Results.RData')





# period 2 only ####
# Incidence (i)
sex_bb = 'total'
vset3 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^r_(?!physmid)', '^ap_', '^NDVI_', sep = '|')
smerc_lung_3it_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Lung_total_3', sex_b = 'total', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_3it_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Stomach_total_3', sex_b = 'total', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_lung_2it_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_i_Lung_total_2', sex_b = 'total', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_2it_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_i_Stomach_total_2', sex_b = 'total', adjust = TRUE, string_search = vset3, ncores = NTHREADS)

sex_bb = 'male'
smerc_lung_3im_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Lung_male_3', sex_b = 'male', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_3im_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Stomach_male_3', sex_b = 'male', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_lung_2im_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_i_Lung_male_2', sex_b = 'male', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_2im_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_i_Stomach_male_2', sex_b = 'male', adjust = TRUE, string_search = vset3, ncores = NTHREADS)

sex_bb = 'female'
smerc_lung_3if_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Lung_female_3', sex_b = 'female', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_3if_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Stomach_female_3', sex_b = 'female', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_lung_2if_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_i_Lung_female_2', sex_b = 'female', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_2if_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_i_Stomach_female_2', sex_b = 'female', adjust = TRUE, string_search = vset3, ncores = NTHREADS)

# par(mfcol = c(1,2))
# plot(smerc_lung_3dt)
# plot(smerc_lung_3dtu)


# Mortality (d)
sex_bb = 'total'
smerc_lung_3dt_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Lung_total_3', sex_b = 'total', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_lung_2dt_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_d_Lung_total_2', sex_b = 'total', adjust = TRUE, string_search = vset3, ncores = NTHREADS)


smerc_stom_3dt_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Stomach_total_3', sex_b = 'total', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_2dt_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_total_2', yvar = 'n_d_Stomach_total_2', sex_b = 'total', adjust = TRUE, string_search = vset3, ncores = NTHREADS)

sex_bb = 'male'
smerc_lung_3dm_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_lung_2dm_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_d_Lung_male_2', sex_b = 'male', adjust = TRUE, string_search = vset3, ncores = NTHREADS)

smerc_stom_3dm_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Stomach_male_3', sex_b = 'male', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_2dm_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_male_2', yvar = 'n_d_Stomach_male_2', sex_b = 'male', adjust = TRUE, string_search = vset3, ncores = NTHREADS)

sex_bb = 'female'
smerc_lung_3df_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Lung_female_3', sex_b = 'female', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_lung_2df_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_d_Lung_female_2', sex_b = 'female', adjust = TRUE, string_search = vset3, ncores = NTHREADS)

smerc_stom_3df_v3 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Stomach_female_3', sex_b = 'female', adjust = TRUE, string_search = vset3, ncores = NTHREADS)
smerc_stom_2df_v3 = run_smerc_cancertype(data = covar_origin_05_fc, population = 'n_p_female_2', yvar = 'n_d_Stomach_female_2', sex_b = 'female', adjust = TRUE, string_search = vset3, ncores = NTHREADS)





# period 3 only ####


# Available sociodemographic and environmental variables for all periods ####
# Incidence (i)
sex_bb = 'total'
vset4 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^r_', '^n_pw', '^ap_', '^NDVI_', sep = '|')
vset4_st = str_c(str_c('^p_*.*_', sex_bb, '$'), '^r_', '^p_candiag', '^n_pw', '^ap_', '^NDVI_', sep = '|')
smerc_lung_3it_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Lung_total_3', sex_b = 'total', adjust = TRUE, string_search = vset4, ncores = NTHREADS)
smerc_stom_3it_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Stomach_total_3', sex_b = 'total', adjust = TRUE, string_search = vset4_st, ncores = NTHREADS)

sex_bb = 'male'
smerc_lung_3im_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Lung_male_3', sex_b = 'male', adjust = TRUE, string_search = vset4, ncores = NTHREADS)
smerc_stom_3im_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Stomach_male_3', sex_b = 'male', adjust = TRUE, string_search = vset4_st, ncores = NTHREADS)

sex_bb = 'female'
smerc_lung_3if_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Lung_female_3', sex_b = 'female', adjust = TRUE, string_search = vset4, ncores = NTHREADS)
smerc_stom_3if_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Stomach_female_3', sex_b = 'female', adjust = TRUE, string_search = vset4_st, ncores = NTHREADS)


# Mortality (d)
sex_bb = 'total'
smerc_lung_3dt_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Lung_total_3', sex_b = 'total', adjust = TRUE, string_search = vset4, ncores = NTHREADS)
smerc_stom_3dt_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Stomach_total_3', sex_b = 'total', adjust = TRUE, string_search = vset4_st, ncores = NTHREADS)

sex_bb = 'male'
smerc_lung_3dm_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male', adjust = TRUE, string_search = vset4, ncores = NTHREADS)
smerc_stom_3dm_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Stomach_male_3', sex_b = 'male', adjust = TRUE, string_search = vset4_st, ncores = NTHREADS)

sex_bb = 'female'
smerc_lung_3df_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Lung_female_3', sex_b = 'female', adjust = TRUE, string_search = vset4, ncores = NTHREADS)
smerc_stom_3df_v4 = run_smerc_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Stomach_female_3', sex_b = 'female', adjust = TRUE, string_search = vset4_st, ncores = NTHREADS)


save(list = ls()[grep('^smerc_', ls())],
     file = str_c(dbdir, 'Manuscript/Scan_SMERC_periods_1_3_vsets_1_4_Results_35p_p001.RData'))
