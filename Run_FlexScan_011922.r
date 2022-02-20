## FlexScan Run
### 02/19/22
username = 'sigma'
if (!require(pacman)) {install.packages('pacman')}
p_load(tidyverse, sf, spdep, DCluster, rflexscan, tmap, dtplyr, fuzzyjoin, readxl, here)

source('./base_functions.R')
drive = sprintf('/mnt/c/Users/%s/OneDrive/NCC_Project/CancerClustering/', username)
#source('./Cleaning_Population_011922.r')
load(str_c(drive, "Manuscript/Clustering_Base_sf_021722.RData"))



## Tan-go's spatial scan statistic based on flexible windows
flexrun = function(case_mat, coord_mat, adj_mat, name_vec, pop_mat) {
    expected = (sum(as.vector(case_mat)/sum(pop_mat))) 
    rflexscan(x = coord_mat[,1],
              y = coord_mat[,2],
              nb = adj_mat,
              name = name_vec,
              observed = as.vector(case_mat),
              expected = expected * pop_mat,
              clustersize = 100,
              stattype = 'RESTRICTED',
              scanmethod = 'FLEXIBLE',
              rantype = 'POISSON')
}
data = covar_origin_10_fc
case_col = 'n_d_Stomach_total_3'
            pop_col = 'n_p_total_3'
            name_col = 'sgg_cd_c'
            

flexrun_sf = function(data, case_col, pop_col, name_col, covar_list = NULL, covar_control = FALSE, cl_size = 100) {
    
    coord_mat = st_coordinates(st_centroid(data))
    data_df = as.data.frame(st_drop_geometry(data))
    adj_mat = nb2mat(poly2nb(data, queen = TRUE), style = 'B', zero.policy = TRUE)
    case_mat = data_df[, case_col]
    pop_mat = data_df[, pop_col]
    name_vec = as.character(data_df[, name_col])

    if (!covar_control) {
        expected = (sum(as.vector(case_mat))/sum(as.vector(pop_mat))) 
        comment_message = str_c('No covariate control.')
    } else {
        covar_concat = paste(covar_list, collapse = '+')
        model_form = as.formula(paste(case_col, '~', covar_concat, sep = ""))
        model_n = glm(formula = model_form, data = data_df, family = poisson)
        pop_mat = model_n$fitted.values
        expected = (sum(as.vector(case_mat))/sum(pop_mat)) 
        comment_message = str_c('Outcome: ', case_col, "; covariate control: ", covar_control)
    }

    rflex = 
    rflexscan(x = coord_mat[,1],
              y = coord_mat[,2],
              nb = adj_mat,
              name = name_vec,
              observed = as.vector(case_mat),
              expected = expected * pop_mat,
              stattype = 'RESTRICTED',
              scanmethod = 'FLEXIBLE',
              rantype = 'POISSON',
              clustersize = cl_size,
              secondary = 30,
              verbose = TRUE,
              comments = comment_message)

    return(rflex)
}

##
sex_bb = "total"
vset4_st = str_c(str_c('^p_*.*_', sex_bb, '$'), '^r_', '^p_candiag', '^n_pw', '^ap_', '^NDVI_', sep = '|')
cn_vset4 = colnames(covar_origin_10_fc)[grep(vset4_st, colnames(covar_origin_10_fc), perl = TRUE)]
co10fc = covar_origin_10_fc %>%
    mutate(n_d_Stomach_total_3 = n_d_Stomach_total_3,
           n_p_total_3 = n_p_total_3 / 5)

rflex_stomach_dt_3_v4 = 
    flexrun_sf(data = covar_origin_10_fc,
            case_col = 'n_d_Stomach_total_3',
            pop_col = 'n_p_total_3',
            name_col = 'sgg_cd_c',
            covar_list = cn_vset4,
            covar_control = TRUE)
rflex_stomach_dt_3_v4 %>% plot
summary(rflex_stomach_dt_3_v4)
choropleth(st_geometry(covar_origin_10_fc), rflex_stomach_dt_3_v4, pval = 0.05)

rflex_stomach_dt_3_u = 
    flexrun_sf(data = covar_origin_10_fc,
            case_col = 'n_d_Stomach_total_3',
            pop_col = 'n_p_total_3',
            name_col = 'sgg_cd_c',
            covar_control = FALSE,
            cl_size = 50)
rflex_stomach_dt_3_u %>% plot
choropleth(st_geometry(covar_origin_10_fc), rflex_stomach_dt_3_u, pval = 0.05)

## Base data for flexrun
cmat = covar_origin_10_fc %>% st_drop_geometry
pmat_t = covar_origin_10_fc %>% .$T_2010
pmat_m = mor_dc2010_sf %>% .$M_2010
pmat_f = mor_dc2010_sf %>% .$F_2010
comat = as.matrix(st_coordinates(st_centroid(covar_origin_10_fc)))
adjmat = nb2mat(poly2nb(covar_origin_10_fc), style = 'B', zero.policy = TRUE)
nvec = mor_to2010_sf$sgg_nm


# Flexrun: All sexes
mor_to2010_stomach_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Stomach']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_t)

mor_to2010_colorectal_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Colorectal']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_t)

mor_to2010_liver_flex = 
    flexrun(case_mat = as.matrix(mor_to2010_sf %>% st_drop_geometry %>% .[,'Liver']),
            coord_mat = as.matrix(st_coordinates(st_centroid(mor_to2010_sf))),
            adj_mat = nb2mat(poly2nb(mor_to2010_sf), style = 'B', zero.policy = TRUE),
            name_vec = mor_to2010_sf$sgg_nm,
            pop_mat = pmat_t)

mor_to2010_lung_flex = 
    flexrun(case_mat = as.matrix(mor_to2010_sf %>% st_drop_geometry %>% .[,'Lung']),
            coord_mat = as.matrix(st_coordinates(st_centroid(mor_to2010_sf))),
            adj_mat = nb2mat(poly2nb(mor_to2010_sf), style = 'B', zero.policy = TRUE),
            name_vec = mor_to2010_sf$sgg_nm,
            pop_mat = pmat_t)

# Sex-specific 
mor_to2010_breast_flex = 
    flexrun(case_mat = as.matrix(mor_to2010_sf %>% st_drop_geometry %>% .[,'Breast']),
            coord_mat = as.matrix(st_coordinates(st_centroid(mor_to2010_sf))),
            adj_mat = nb2mat(poly2nb(mor_to2010_sf), style = 'B', zero.policy = TRUE),
            name_vec = mor_to2010_sf$sgg_nm,
            pop_mat = pmat_t)


## Male
mor_me2010_stomach_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Stomach_1']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_m)
mor_me2010_colorectal_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Colorectal_1']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_m)
mor_me2010_liver_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Liver_1']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_m)
mor_me2010_lung_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Lung_1']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_m)
mor_me2010_breast_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Breast_1']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_m)
mor_me2010_prostate_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Prostate_1']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_m)

## Male
mor_fe2010_stomach_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Stomach_2']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_f)
mor_fe2010_colorectal_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Colorectal_2']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_f)
mor_fe2010_liver_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Liver_2']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_f)
mor_fe2010_lung_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Lung_2']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_f)
mor_fe2010_breast_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Breast_2']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_f)
mor_fe2010_cervical_flex = 
    flexrun(case_mat = as.matrix(cmat %>%.[,'Cervical_2']),
            coord_mat = comat,
            adj_mat = adjmat,
            name_vec = nvec,
            pop_mat = pmat_f)
