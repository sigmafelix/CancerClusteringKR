## FlexScan Run


## Tan-go's spatial scan statistic based on flexible windows
flexrun = function(case_mat, coord_mat, adj_mat, name_vec, pop_mat) {
    expected = (sum(as.vector(case_mat)/sum(pop_mat))) 
    rflexscan(x = coord_mat[,1],
              y = coord_mat[,2],
              nb = adj_mat,
              name = name_vec,
              observed = as.vector(case_mat),
              expected = expected * pop_mat,
              stattype = 'RESTRICTED',
              scanmethod = 'FLEXIBLE',
              rantype = 'POISSON')
}

## Base data for flexrun
cmat = mor_to2010_sf %>% st_drop_geometry
cmat = mor_dc2010_sf %>% st_drop_geometry
pmat_t = mor_dc2010_sf %>% .$T_2010
pmat_m = mor_dc2010_sf %>% .$M_2010
pmat_f = mor_dc2010_sf %>% .$F_2010
comat = as.matrix(st_coordinates(st_centroid(mor_to2010_sf)))
adjmat = nb2mat(poly2nb(mor_to2010_sf), style = 'B', zero.policy = TRUE)
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
