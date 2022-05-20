### SaTScan standalone batch
### May 18, 2022
### Insang Song (sigmafelix@hotmail.com)


library(pacman)
p_load(stringr, readr, sf, dplyr)
source("base_functions.R")
load("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/Manuscript/Clustering_Base_sf_042022.RData")

covar_origin_00_fcd = covar_origin_00_fc %>% bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
    filter(!grepl('^(23320|37430|39)', sgg_cd_c)) %>%
    st_drop_geometry
covar_origin_05_fcd = covar_origin_05_fc %>% bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
    filter(!grepl('^(23320|37430|39)', sgg_cd_c)) %>%
    st_drop_geometry
covar_origin_10_fcd = covar_origin_10_fc %>% bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
    filter(!grepl('^(23320|37430|39)', sgg_cd_c)) %>%
    st_drop_geometry

# write_csv(covar_origin_00_fcd, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period1.csv")
# write_csv(covar_origin_05_fcd, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period2.csv")
# write_csv(covar_origin_10_fcd, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period3.csv")

dat1 = read_csv("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period1.csv")
dat2 = read_csv("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period2.csv")
dat3 = read_csv("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period3.csv")

# file.copy("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period1.csv", "/home/felix/Documents/satscan_base_data_period1.csv")
# file.copy("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period2.csv", "/home/felix/Documents/satscan_base_data_period2.csv")
# file.copy("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_base_data_period3.csv", "/home/felix/Documents/satscan_base_data_period3.csv")


dbase = "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/"
dtarg = "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/SaTScan_Results/"
geoid = "sgg_cd_c"

all_models = expand.grid(
    mid_title = c("_i_", "_d_"),
    cancertype = c("Lung", "Stomach"),
    sex = c("_total_", "_male_", "_female_"),
    valuetype = c("asmr", "count"),
    period = 1:3
) %>%
    mutate(basename = str_c(mid_title, cancertype, sex, period),
           modeltitle = str_c(cancertype, "_p", period, mid_title, valuetype, sex))

# run for the period 1 (normal) ####
p1_all = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = y,
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(y, ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/home/felix/Documents/test.prm"
                    )
}, all_models$basename[1:24], all_models$modeltitle[1:24], SIMPLIFY = FALSE)

p1_all_df = Reduce(rbind, p1_all)

# run for the period 2 (normal) ####
p2_all = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = y,
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(y, ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/home/felix/Documents/test.prm"
                    )
}, all_models$basename[25:48], all_models$modeltitle[25:48], SIMPLIFY = FALSE)

p2_all_df = Reduce(rbind, p2_all)

# run for the period 3 (normal) ####
p3_all = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = y,
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(y, ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm"
                    )
}, all_models$basename[49:72], all_models$modeltitle[49:72], SIMPLIFY = FALSE)

p3_all_df = Reduce(rbind, p3_all)

p13_all_df = bind_rows(p1_all_df, p2_all_df) %>%
    bind_rows(p3_all_df)
write_rds(p13_all_df, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_ASMR_uncontrolled_resid.rds")




vset1 = str_c('^(p_65p|p_hbac)*.*_', sex_bb, '$')
vset2 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^ap_', '^NDVI_', sep = '|')
vset3 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^r_(?!physmid)', '^ap_', '^NDVI_', sep = '|')
vset4 = str_c(str_c('^p_*.*_', sex_bb, '$'), '^r_', '^n_pw', '^ap_', '^NDVI_', sep = '|')

# Covariate control (Set 1) ####
# run for the period 1 (normal) ####
p1_all_v1 = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = str_c(y, 'v1'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(str_c(y, 'v1'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = "set1")
}, all_models$basename[1:24], all_models$modeltitle[1:24], SIMPLIFY = FALSE)

p1_all_v1_df = Reduce(rbind, p1_all_v1)

# run for the period 2 (normal) ####
p2_all_v1 = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'v1'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'v1'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set1'
                    )
}, all_models$basename[25:48], all_models$modeltitle[25:48], SIMPLIFY = FALSE)

p2_all_v1_df = Reduce(rbind, p2_all_v1)

# run for the period 3 (normal) ####
p3_all_v1 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v1'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v1'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set1'
                    )
}, all_models$basename[49:72], all_models$modeltitle[49:72], SIMPLIFY = FALSE)

p3_all_v1_df = Reduce(rbind, p3_all_v1)

p13_all_v1_df = bind_rows(p1_all_v1_df, p2_all_v1_df) %>%
    bind_rows(p3_all_v1_df)
write_rds(p13_all_v1_df, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_ASMR_vset1_resid.rds")


## VSet 2 ####
# run for the period 1 (normal) ####
p1_all_v2 = mapply(function(x, y) {
    generate_satscan_prm(data = dat1,
                    title.analysis = str_c(y, 'v2'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period1.csv",
                    filename.output = str_c(str_c(y, 'v2'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = "set2")
}, all_models$basename[1:24], all_models$modeltitle[1:24], SIMPLIFY = FALSE)

p1_all_v2_df = Reduce(rbind, p1_all_v2)

# run for the period 2 (normal) ####
p2_all_v2 = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'v2'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'v2'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set2'
                    )
}, all_models$basename[25:48], all_models$modeltitle[25:48], SIMPLIFY = FALSE)

p2_all_v2_df = Reduce(rbind, p2_all_v2)

# run for the period 3 (normal) ####
p3_all_v2 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v2'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v2'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set2'
                    )
}, all_models$basename[49:72], all_models$modeltitle[49:72], SIMPLIFY = FALSE)

p3_all_v2_df = Reduce(rbind, p3_all_v2)

p13_all_v2_df = bind_rows(p1_all_v2_df, p2_all_v2_df) %>%
    bind_rows(p3_all_v2_df)
write_rds(p13_all_v2_df, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_ASMR_vset2_resid.rds")


## VSet 3 ####
# run for the period 2 (normal) ####
p2_all_v3 = mapply(function(x, y) {
    generate_satscan_prm(data = dat2,
                    title.analysis = str_c(y, 'v3'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period2.csv",
                    filename.output = str_c(str_c(y, 'v3'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set3'
                    )
}, all_models$basename[25:48], all_models$modeltitle[25:48], SIMPLIFY = FALSE)

p2_all_v3_df = Reduce(rbind, p2_all_v3)

# run for the period 3 (normal) ####
p3_all_v3 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v3'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v3'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set3'
                    )
}, all_models$basename[49:72], all_models$modeltitle[49:72], SIMPLIFY = FALSE)

p3_all_v3_df = Reduce(rbind, p3_all_v3)

p13_all_v3_df = bind_rows(p2_all_v3_df, p3_all_v3_df)
write_rds(p13_all_v3_df, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_ASMR_vset3_resid.rds")

## VSet 4 ####
# run for the period 3 (normal) ####
p3_all_v4 = mapply(function(x, y) {
    generate_satscan_prm(data = dat3,
                    title.analysis = str_c(y, 'v4'),
                    dir.base = dbase, dir.target = "/home/felix/workcache/", name.idcol = geoid,
                    filename.input = "satscan_base_data_period3.csv",
                    filename.output = str_c(str_c(y, 'v4'), ".txt"),
                    col.var = str_c('ragest', x),
                    col.case = str_c('n', x),
                    prm.path = "/mnt/c/Users/sigma/Documents/test.prm",
                    adjust = TRUE,
                    vset = 'set4'
                    )
}, all_models$basename[49:72], all_models$modeltitle[49:72], SIMPLIFY = FALSE)

p3_all_v4_df = Reduce(rbind, p3_all_v4)
write_rds(p3_all_v4_df, "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/satscan_ASMR_vset4_resid.rds")

