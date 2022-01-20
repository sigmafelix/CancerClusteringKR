### Run Clustering Analysis
### 01/19/22

source('./base_functions.R')
source('./Cleaning_Population_011922.r')

# knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, error = FALSE, fig.height = 7.5)

covar_origin_10 = get_basecovar(target_year = 2010)
covar_origin_10_fc = clean_consolidated(cleaned_df = covar_origin_10) %>%
    mutate_at(.vars = vars(-geom),
              .funs = list(~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    left_join(morinc_clw, by = c('sgg_cd_c'))

# Naming: smerc_[cancertype]_[period][incidence/mortality][sex][uncontrolled]
# Incidence (i)
smerc_lung_3it = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_i_Lung_total_3', sex_b = 'total', ncores = 16)
smerc_stom_3it = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_i_Stomach_total_3', sex_b = 'total', ncores = 16)
smerc_lung_3im = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_i_Lung_male_3', sex_b = 'male', ncores = 16)
smerc_stom_3im = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_i_Stomach_male_3', sex_b = 'male', ncores = 16)
smerc_lung_3if = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_i_Lung_female_3', sex_b = 'female', ncores = 16)
smerc_stom_3if = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_i_Stomach_female_3', sex_b = 'female', ncores = 16)


smerc_lung_3itu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_total_3', yvar = 'n_i_Lung_total_3', sex_b = 'total', ncores = 16)
smerc_stom_3itu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_total_3', yvar = 'n_i_Stomach_total_3', sex_b = 'total', ncores = 16)
smerc_lung_3imu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_male_3', yvar = 'n_i_Lung_male_3', sex_b = 'male', ncores = 16)
smerc_stom_3imu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_male_3', yvar = 'n_i_Stomach_male_3', sex_b = 'male', ncores = 16)
smerc_lung_3ifu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_female_3', yvar = 'n_i_Lung_female_3', sex_b = 'female', ncores = 16)
smerc_stom_3ifu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_female_3', yvar = 'n_i_Stomach_female_3', sex_b = 'female', ncores = 16)


# Mortality (d)
smerc_lung_3dt = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_d_Lung_total_3', sex_b = 'total', ncores = 16)
smerc_stom_3dt = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_d_Stomach_total_3', sex_b = 'total', ncores = 16)
smerc_lung_3dm = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_d_Lung_male_3', sex_b = 'male', ncores = 16)
smerc_stom_3dm = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_d_Stomach_male_3', sex_b = 'male', ncores = 16)
smerc_lung_3df = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_d_Lung_female_3', sex_b = 'female', ncores = 16)
smerc_stom_3df = run_smerc_cancertype(data = covar_origin_10_fc, yvar = 'n_d_Stomach_female_3', sex_b = 'female', ncores = 16)


smerc_lung_3dtu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_total_3', yvar = 'n_d_Lung_total_3', sex_b = 'total', ncores = 16)
smerc_stom_3dtu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_total_3', yvar = 'n_d_Stomach_total_3', sex_b = 'total', ncores = 16)
smerc_lung_3dmu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male', ncores = 16)
smerc_stom_3dmu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_male_3', yvar = 'n_d_Stomach_male_3', sex_b = 'male', ncores = 16)
smerc_lung_3dfu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_female_3', yvar = 'n_d_Lung_female_3', sex_b = 'female', ncores = 16)
smerc_stom_3dfu = run_smerc_cancertype_pl(data = covar_origin_10_fc, pop = 'n_p_female_3', yvar = 'n_d_Stomach_female_3', sex_b = 'female', ncores = 16)


covar_origin_10_fc_df = covar_origin_10_fc %>%
    bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
    st_drop_geometry

run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Lung_total_3', sex_b = 'total', run_glm = TRUE)
run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Stomach_total_3', sex_b = 'total', run_glm = TRUE)
run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male', run_glm = TRUE)

# covariate_adjusted
dclust_lung_3it = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Lung_total_3', sex_b = 'total')
dclust_lung_3dt = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Lung_total_3', sex_b = 'total')
dclust_lung_3im = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Lung_male_3', sex_b = 'male')
dclust_lung_3dm = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male')
dclust_lung_3if = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Lung_female_3', sex_b = 'female')
dclust_lung_3df = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Lung_female_3', sex_b = 'female')

dclust_stom_3it = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Stomach_total_3', sex_b = 'total')
dclust_stom_3dt = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Stomach_total_3', sex_b = 'total')
dclust_stom_3im = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Stomach_male_3', sex_b = 'male')
dclust_stom_3dm = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Stomach_male_3', sex_b = 'male')
dclust_stom_3if = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Stomach_female_3', sex_b = 'female')
dclust_stom_3df = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Stomach_female_3', sex_b = 'female')

# covariate_unadjusted
dclust_lung_3itu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Lung_total_3', sex_b = 'total', adjust = FALSE)
dclust_lung_3dtu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Lung_total_3', sex_b = 'total', adjust = FALSE)
dclust_lung_3imu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Lung_male_3', sex_b = 'male', adjust = FALSE)
dclust_lung_3dmu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Lung_male_3', sex_b = 'male', adjust = FALSE)
dclust_lung_3ifu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Lung_female_3', sex_b = 'female', adjust = FALSE)
dclust_lung_3dfu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Lung_female_3', sex_b = 'female', adjust = FALSE)

dclust_stom_3itu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_i_Stomach_total_3', sex_b = 'total', adjust = FALSE)
dclust_stom_3dtu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_total_3', yvar = 'n_d_Stomach_total_3', sex_b = 'total', adjust = FALSE)
dclust_stom_3imu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_i_Stomach_male_3', sex_b = 'male', adjust = FALSE)
dclust_stom_3dmu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_male_3', yvar = 'n_d_Stomach_male_3', sex_b = 'male', adjust = FALSE)
dclust_stom_3ifu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_i_Stomach_female_3', sex_b = 'female', adjust = FALSE)
dclust_stom_3dfu = run_dclust_cancertype(data = covar_origin_10_fc, population = 'n_p_female_3', yvar = 'n_d_Stomach_female_3', sex_b = 'female', adjust = FALSE)


save(list = ls()[grep('^smerc_*.*_[3]*.*', ls())],
     file = 'Scan_SMERC_2009_2013_Results.RData')
save(list = ls()[grep('^dclust_*.*_[3]*.*', ls())],
     file = 'Scan_DCLUST_2009_2013_Results.RData')
save(list = ls()[grep('^(smerc|dclust|covar_origin_10_)', ls())],
     file = 'Scan_2009_2013_Results.RData')
