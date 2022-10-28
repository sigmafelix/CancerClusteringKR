### SaTScan standalone batch
### June 22, 2022
### Last revised Oct 27, 2022
### Insang Song (sigmafelix@hotmail.com)


library(pacman)
p_load(stringr, readr, sf, dplyr)
source("./Code/Base/base_functions.R")

cldir = "/mnt/c/Users/sigma/"
"%s%" = function(x, y) str_c(x, y)

load(cldir %s% "OneDrive/NCC_Project/CancerClustering/Manuscript/Clustering_Base_sf_091522.RData")

covar_origin_00_fcd = covar_origin_00_fc %>% bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
    filter(!grepl('^(23320|37430|39)', sgg_cd_c)) %>%
    st_drop_geometry %>%
    mutate(case_normal = 1, dummyyear = 2000)
covar_origin_05_fcd = covar_origin_05_fc %>% bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
    filter(!grepl('^(23320|37430|39)', sgg_cd_c)) %>%
    st_drop_geometry %>%
    mutate(case_normal = 1, dummyyear = 2005)
covar_origin_10_fcd = covar_origin_10_fc %>% bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
    filter(!grepl('^(23320|37430|39)', sgg_cd_c)) %>%
    st_drop_geometry %>%
    mutate(case_normal = 1, dummyyear = 2010)

# write_csv(covar_origin_00_fcd, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period1.csv")
# write_csv(covar_origin_05_fcd, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period2.csv")
# write_csv(covar_origin_10_fcd, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period3.csv")

dat1 = read_csv(cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_base_data_period1.csv")
dat2 = read_csv(cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_base_data_period2.csv")
dat3 = read_csv(cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_base_data_period3.csv")

# file.copy("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period1.csv", "/home/felix/Documents/satscan_base_data_period1.csv", overwrite = TRUE)
# file.copy("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period2.csv", "/home/felix/Documents/satscan_base_data_period2.csv", overwrite = TRUE)
# file.copy("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period3.csv", "/home/felix/Documents/satscan_base_data_period3.csv", overwrite = TRUE)


dbase = cldir %s% "OneDrive/NCC_Project/CancerClustering/"
dtarg = cldir %s% "OneDrive/NCC_Project/CancerClustering/SaTScan_Results/"
geoid = "sgg_cd_c"

all_models = expand.grid(
    mid_title = c("_i_", "_d_"),
    cancertype = c("Lung", "Stomach"),
    sex = c("_male_", "_female_"), # excluded _total_ June 2022
    valuetype = c("ragest", "n"),
    period = 1:3
) %>%
    mutate(basename = str_c(valuetype, mid_title, cancertype, sex, period),
           modeltitle = str_c(cancertype, "_p", period, mid_title, valuetype, sex))


indx_p1 = 1:8# 16
indx_p2 = 17:24# 32
indx_p3 = 33:40# 48

# p1_all = mapply(function(x, y) {
#     generate_satscan_prm(data = dat1,
#                     title.analysis = y,
#                     dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
#                     filename.input = "satscan_base_data_period1.csv",
#                     filename.output = str_c(y, ".txt"),
#                     col.var = x,
#                     col.case = str_c('n', x),
#                     prm.path = "/home/felix/Documents/test.prm"
#                     )
# }, all_models$basename[10], all_models$modeltitle[10], SIMPLIFY = FALSE)


# run for the period 1 (normal) ####
p1_all = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = y,
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(y, ".txt"),
                    col.var = x,
                    weighted = T,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm"
                    )
}, all_models$basename[indx_p1], all_models$modeltitle[indx_p1], SIMPLIFY = FALSE)

p1_all_df = Reduce(rbind, p1_all)

# run for the period 2 (normal) ####
p2_all = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = y,
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(y, ".txt"),
                    col.var = x,
                    weighted = T,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm"
                    )
}, all_models$basename[indx_p2], all_models$modeltitle[indx_p2], SIMPLIFY = FALSE)

p2_all_df = Reduce(rbind, p2_all)

# run for the period 3 (normal) ####
p3_all = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = y,
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(y, ".txt"),
                    col.var = x,
                    weighted = T,
                    col.case = x, # str_c('n', x)
                    prm.path = "/home/felix/Documents/test.prm"
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_df = Reduce(rbind, p3_all)

p13_all_df = bind_rows(p1_all_df, p2_all_df) %>%
    bind_rows(p3_all_df)
write_rds(p13_all_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_uncontrolled_resid.rds")




# vset1 = str_c('^(p_65p|p_hbac)*.*_', sex_bb, '$')
# vset2 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^ap_', '^NDVI_', sep = '|')
# vset3 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^r_(?!physmid)', '^ap_', '^NDVI_', sep = '|')
# vset4 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^r_', '^n_pw', '^ap_', '^NDVI_', sep = '|')

# Covariate control (age only) ####
# run for the period 1 (normal) ####
p1_all_va = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = str_c(y, 'va'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(str_c(y, 'va'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = "age")
}, all_models$basename[indx_p1], all_models$modeltitle[indx_p1], SIMPLIFY = FALSE)

p1_all_va_df = Reduce(rbind, p1_all_va)

# run for the period 2 (normal) ####
p2_all_va = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'va'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'va'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'age'
                    )
}, all_models$basename[indx_p2], all_models$modeltitle[indx_p2], SIMPLIFY = FALSE)

p2_all_va_df = Reduce(rbind, p2_all_va)

# run for the period 3 (normal) ####
p3_all_va = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'va'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'va'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'age'
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_va_df = Reduce(rbind, p3_all_va)

p13_all_va_df = bind_rows(p1_all_va_df, p2_all_va_df) %>%
    bind_rows(p3_all_va_df)
write_rds(p13_all_va_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_age_resid.rds")


# Covariate control (excl. age in M1: educ) ####
# run for the period 1 (normal) ####
p1_all_ve = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = str_c(y, 've'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(str_c(y, 've'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = "educ")
}, all_models$basename[indx_p1], all_models$modeltitle[indx_p1], SIMPLIFY = FALSE)

p1_all_ve_df = Reduce(rbind, p1_all_ve)

# run for the period 2 (normal) ####
p2_all_ve = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 've'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 've'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'educ'
                    )
}, all_models$basename[indx_p2], all_models$modeltitle[indx_p2], SIMPLIFY = FALSE)

p2_all_ve_df = Reduce(rbind, p2_all_ve)

# run for the period 3 (normal) ####
p3_all_ve = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 've'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 've'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'educ'
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_ve_df = Reduce(rbind, p3_all_ve)

p13_all_ve_df = bind_rows(p1_all_ve_df, p2_all_ve_df) %>%
    bind_rows(p3_all_ve_df)
write_rds(p13_all_ve_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_educ_resid.rds")


# Covariate control (environmental variables only) ####
# run for the period 1 (normal) ####
p1_all_ven = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = str_c(y, 'ven'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(str_c(y, 'ven'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = "environ")
}, all_models$basename[indx_p1], all_models$modeltitle[indx_p1], SIMPLIFY = FALSE)

p1_all_ven_df = Reduce(rbind, p1_all_ven)

# run for the period 2 (normal) ####
p2_all_ven = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'ven'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'ven'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'environ'
                    )
}, all_models$basename[indx_p2], all_models$modeltitle[indx_p2], SIMPLIFY = FALSE)

p2_all_ven_df = Reduce(rbind, p2_all_ven)

# run for the period 3 (normal) ####
p3_all_ven = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'ven'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'ven'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'environ'
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_ven_df = Reduce(rbind, p3_all_ven)

p13_all_ven_df = bind_rows(p1_all_ven_df, p2_all_ven_df) %>%
    bind_rows(p3_all_ven_df)
write_rds(p13_all_ven_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_environ_resid.rds")



# Covariate control (behavioral variables only) ####
# run for the period 2 (normal) ####
p2_all_veb = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'veb'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'veb'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'behav'
                    )
}, all_models$basename[indx_p2], all_models$modeltitle[indx_p2], SIMPLIFY = FALSE)

p2_all_veb_df = Reduce(rbind, p2_all_veb)

# run for the period 3 (normal) ####
p3_all_veb = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'veb'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'veb'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'behav'
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_veb_df = Reduce(rbind, p3_all_veb)

p13_all_veb_df = bind_rows(p2_all_veb_df, p3_all_veb_df)
write_rds(p13_all_veb_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_behavior_resid.rds")


# Covariate control (Model 2: variable set 1) ####
# run for the period 1 (normal) ####
p1_all_v1 = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = str_c(y, 'v1'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(str_c(y, 'v1'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = "set1")
}, all_models$basename[indx_p1], all_models$modeltitle[indx_p1], SIMPLIFY = FALSE)

p1_all_v1_df = Reduce(rbind, p1_all_v1)

# run for the period 2 (normal) ####
p2_all_v1 = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'v1'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'v1'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set1'
                    )
}, all_models$basename[indx_p2], all_models$modeltitle[indx_p2], SIMPLIFY = FALSE)

p2_all_v1_df = Reduce(rbind, p2_all_v1)

# run for the period 3 (normal) ####
p3_all_v1 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v1'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v1'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set1'
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_v1_df = Reduce(rbind, p3_all_v1)

p13_all_v1_df = bind_rows(p1_all_v1_df, p2_all_v1_df) %>%
    bind_rows(p3_all_v1_df)
write_rds(p13_all_v1_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_vset1_resid.rds")


## Covariate control (Model 3: variable set 2) ####
# run for the period 1 (normal) ####
p1_all_v2 = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = str_c(y, 'v2'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(str_c(y, 'v2'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = "set2")
}, all_models$basename[indx_p1], all_models$modeltitle[indx_p1], SIMPLIFY = FALSE)

p1_all_v2_df = Reduce(rbind, p1_all_v2)

# run for the period 2 (normal) ####
p2_all_v2 = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'v2'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'v2'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set2'
                    )
}, all_models$basename[indx_p2], all_models$modeltitle[indx_p2], SIMPLIFY = FALSE)

p2_all_v2_df = Reduce(rbind, p2_all_v2)

# run for the period 3 (normal) ####
p3_all_v2 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v2'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v2'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set2'
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_v2_df = Reduce(rbind, p3_all_v2)

p13_all_v2_df = bind_rows(p1_all_v2_df, p2_all_v2_df) %>%
    bind_rows(p3_all_v2_df)
write_rds(p13_all_v2_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_vset2_resid.rds")


## Covariate control (Model 4: variable set 3) ####
# run for the period 2 (normal) ####
p2_all_v3 = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'v3'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'v3'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set3'
                    )
}, all_models$basename[indx_p2], all_models$modeltitle[indx_p2], SIMPLIFY = FALSE)

p2_all_v3_df = Reduce(rbind, p2_all_v3)

# run for the period 3 (normal) ####
p3_all_v3 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v3'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v3'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set3'
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_v3_df = Reduce(rbind, p3_all_v3)

p13_all_v3_df = bind_rows(p2_all_v3_df, p3_all_v3_df)
write_rds(p13_all_v3_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_vset3_resid.rds")

## Covariate control (Model 5: variable set 4) ####
# run for the period 3 (normal) ####
p3_all_v4 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v4'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v4'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set4'
                    )
}, all_models$basename[indx_p3], all_models$modeltitle[indx_p3], SIMPLIFY = FALSE)

p3_all_v4_df = Reduce(rbind, p3_all_v4)
write_rds(p3_all_v4_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_vset4_resid.rds")


## Covariate control (Model 5: variable set 3 + Stomach diagnosis rate) ####
# run for the period 3 (normal) ####
p3_all_v5 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v5'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v5'), ".txt"),
                    col.var = x,
                    col.case = x,
                    prm.path = "/home/felix/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set3',
                    add_var = "p_candiag_sto"
                    )
}, all_models$basename[c(35:36, 39:40)], all_models$modeltitle[c(35:36, 39:40)], SIMPLIFY = FALSE)

p3_all_v5_df = Reduce(rbind, p3_all_v5)
write_rds(p3_all_v5_df, cldir %s% "OneDrive/NCC_Project/CancerClustering/satscan_ASMR_stomach_diagnosis_resid.rds")
