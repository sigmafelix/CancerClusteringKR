### Run Clustering Analysis
### 01/04/22

source('./base_functions.R')

knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, error = FALSE, fig.height = 7.5)

options(repos = "https://cran.seoul.go.kr/")
if (!require(pacman)) {
    install.packages("pacman")
}
p_load(tidyverse, sf, spdep, DCluster, tmap, smerc, knitr, readxl, kableExtra, DClusterm, patchwork)

#homedir <- "/home/felix/"
drive <- str_c(basedir, "OneDrive/NCC_Project/CancerClustering/")
geopath <- str_c(basedir, "OneDrive/Data/Korea/")




# incidence
inc <- read.csv(paste(drive, "cancerInc_sgg.csv", sep = ""), fileEncoding = "EUC-KR")
mor_to <- read.csv(paste(drive, "cancerMor_sgg_total.csv", sep = ""), fileEncoding = "EUC-KR")
mor_me <- read.csv(paste(drive, "cancerMor_sgg_male.csv", sep = ""), fileEncoding = "EUC-KR")
mor_fe <- read.csv(paste(drive, "cancerMor_sgg_female.csv", sep = ""), fileEncoding = "EUC-KR")

conv_table = read.csv(paste(geopath, 'SGG_1995_2018_Conversion_Table_201108.csv', sep = ''), fileEncoding = 'EUC-KR')
# total population for the expected values
# sex code: 0(all), 1(male), 2(female)
pop <- read.csv(paste(drive, "Midyear_Population_1998_2019.csv", sep = ""), fileEncoding = "UTF-8") %>%
    mutate(
        sex0 = plyr::mapvalues(sex0, c(0, 1, 2), c("T", "M", "F")),
        sex_e = plyr::mapvalues(sex0, c("T", "M", "F"), c("Total", "Male", "Female")),
        population = as.numeric(population)
    )
pop_wide = pop %>%
    pivot_wider(names_from = c("sex0", "year"), values_from = population)

colnames(inc) <- c("sex0", "sex", "cancer_type0", "cancer_type", "sgg_cd", "sgg_nm", "year0", "year", "inc0", "type_inc", "unit", "n", "x")
colnames(mor_to) <- c("cause0", "cause", "sgg_cd", "sgg_nm", "sex0", "sex", "type0", "type", "unit", paste("Y", 1998:2019, sep = ""), "x")
colnames(mor_me) <- c("cause0", "cause", "sgg_cd", "sgg_nm", "sex0", "sex", "type0", "type", "unit", paste("Y", 1998:2019, sep = ""), "x")
colnames(mor_fe) <- c("cause0", "cause", "sgg_cd", "sgg_nm", "sex0", "sex", "type0", "type", "unit", paste("Y", 1998:2019, sep = ""), "x")


## Data preprocessing ####
sgg <- st_read(paste(geopath, "SGG_Merge_2000_2016.shp", sep = ""))
sgg2010 <- sgg %>%
    filter(BASE_YEAR == 2010)
conv_table_e = conv_table %>%
    filter(from_year >= 1999 & from_year <= 2013) %>%
    dplyr::select(fromcode, tocode)

pop_cl = pop %>%
    mutate(sgg_cd = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode),
           year_agg = cut(year, breaks = c(1998,2003,2008,2013,2018,2019), labels = c('1999-2003', '2004-2008', '2009-2013', '2014-2018', '2019'), right = TRUE)) %>%
    group_by(year_agg, sgg_cd, sex0, sex_e) %>%
    summarize(population = sum(population, na.rm = TRUE)) %>%
    ungroup


inc_cl <- inc %>%
    mutate(
        n_n = as.numeric(n),
        n_n = ifelse(is.na(n_n), 0, n_n),
        cancer_type_e = plyr::mapvalues(cancer_type, unique(cancer_type), c("Stomach", "Colorectal", "Liver", "Lung", "Breast", "Cervical", "Prostate", "Thyroid")),
        sex_e = plyr::mapvalues(sex, unique(sex), c("Total", "Male", "Female")),
        type_inc_e = plyr::mapvalues(type_inc, unique(type_inc), c("N", "r_crude", "r_agest")),
        sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)
    ) %>%
    # not sum
    pivot_wider(id_cols = c(year, sgg_cd_c, cancer_type_e, sex_e), names_from = type_inc_e, values_from = n_n, values_fn = sum) %>%
    # note that sex-specific cancers have incorrect crude rates
    left_join(pop_cl, by = c('year' = 'year_agg', 'sgg_cd_c' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    group_by(year, sgg_cd_c, cancer_type_e, sex_e) %>%
    summarize(Ntotal = sum(N, na.rm = T),
              r_crude = sum(r_crude * (population / sum(population, na.rm = T)), na.rm = T)) %>%
    ungroup %>%
    pivot_longer(cols = 5:6, values_to = 'n_n', names_to = 'type_inc_e')

# summary table
inc_cl_summary <- inc_cl %>%
    group_by(cancer_type_e, year, sex_e, type_inc_e) %>%
    summarize(
        n_nmin = min(n_n),
        n_nmedian = median(n_n),
        n_nmean = mean(n_n),
        n_nmax = max(n_n)
    ) %>%
    ungroup() %>%
    pivot_longer(cols = 5:8) %>%
    pivot_wider(names_from = c("sex_e", "name"))

mor_cl = bind_rows(mor_to, mor_me) %>%
    bind_rows(mor_fe) %>%
    pivot_longer(cols = Y1998:Y2019) %>%
    mutate(year = as.integer(str_sub(name, 2, 5)),
           year_agg = cut(year, breaks = c(1998,2003,2008,2013,2018,2019), labels = c('1999-2003', '2004-2008', '2009-2013', '2014-2018', '2019'), right = TRUE),
           cancer_type_e = plyr::mapvalues(cause, unique(cause), c("Stomach", "Colorectal", "Liver", "Lung/Bronchus", "Breast", "Cervical/Uterine", "Prostate")),
           sex_e = plyr::mapvalues(sex, unique(sex), c('Total', 'Male', 'Female')),
           type_mor_e = plyr::mapvalues(type0, c('T1', 'T4', 'T7'), c('N', 'r_crude', 'r_agest')),
           sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)
        ) %>%
    # not sum
    pivot_wider(id_cols = c(year_agg, sgg_cd_c, cancer_type_e, sex_e), names_from = type_mor_e, values_from = value, values_fn = sum) %>%
    left_join(pop_cl, by = c('year_agg' = 'year_agg', 'sgg_cd_c' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    group_by(year_agg, sgg_cd_c, cancer_type_e, sex_e) %>%
    summarize(Ntotal = sum(N, na.rm = T),
              r_crude = 1e5 * sum(N, na.rm = T)/sum(population, na.rm = T)) %>%
    ungroup %>%
    filter(year_agg %in% c('1999-2003', '2004-2008', '2009-2013'))



inc_clgg_n = inc_cl %>%
    filter(cancer_type_e %in% c('Lung', 'Stomach') & type_inc_e == 'Ntotal') %>%
    ggplot(data = .,
           mapping = aes(x = year, y = n_n, group = interaction(year, cancer_type_e, sex_e))) +
        facet_grid(sex_e ~ cancer_type_e) +
        geom_boxplot() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        #scale_y_continuous(limits = c(0, 225)) +
        labs(subtitle = 'Crude cancer incidence counts by five-year periods and types')

inc_clgg_r = inc_cl %>%
    filter(cancer_type_e %in% c('Lung', 'Stomach') & type_inc_e == 'r_crude') %>%
    ggplot(data = .,
           mapping = aes(x = year, y = n_n, group = interaction(year, cancer_type_e, sex_e))) +
        facet_grid(sex_e ~ cancer_type_e) +
        geom_boxplot() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        #scale_y_continuous(limits = c(0, 225)) +
        labs(subtitle = 'Crude cancer incidence rates per 100,000 persons by five-year periods and types')


mor_clgg_n = mor_cl %>%
    filter(cancer_type_e %in% c('Lung/Bronchus', 'Stomach')) %>%
    #left_join(pop, by = c('sgg_cd' = 'sgg_cd', 'year' = 'year', 'sex_e' = 'sex_e')) %>%
    #group_by(sgg_cd, cause, cancer_type_e, sex_e, year_agg) %>%
    #summarize(rate_100k_5yr = 1e5 * (sum(value, na.rm = T)/sum(population, na.rm = T))) %>%
    #ungroup %>%
    ggplot(data = .,
           mapping = aes(x= year_agg, y = Ntotal, group = interaction(year_agg, cancer_type_e, sex_e))) +
        facet_grid(sex_e ~ cancer_type_e) +
        geom_boxplot() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        #scale_y_continuous(limits = c(0, 225)) +
        labs(subtitle = 'Cancer mortality counts by years and types')


mor_clgg_r = mor_cl %>%
    filter(cancer_type_e %in% c('Lung/Bronchus', 'Stomach')) %>%
    #left_join(pop, by = c('sgg_cd' = 'sgg_cd', 'year' = 'year', 'sex_e' = 'sex_e')) %>%
    #group_by(sgg_cd, cause, cancer_type_e, sex_e, year_agg) %>%
    #summarize(rate_100k_5yr = 1e5 * (sum(value, na.rm = T)/sum(population, na.rm = T))) %>%
    #ungroup %>%
    ggplot(data = .,
           mapping = aes(x= year_agg, y = r_crude, group = interaction(year_agg, cancer_type_e, sex_e))) +
        facet_grid(sex_e ~ cancer_type_e) +
        geom_boxplot() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        #scale_y_continuous(limits = c(0, 225)) +
        labs(subtitle = 'Crude cancer mortality rates per 100,000 persons by years and types')



inc_cl %>% 
    dplyr::select(year, year, sgg_cd, sgg_nm) %>%
    pivot_wider(id_cols = c(sgg_cd, sgg_nm), names_from = year, values_from = year0, values_fn = length) %>%
    data.frame

mor_cl %>% 
    dplyr::select(year_agg, value, sgg_cd, sgg_nm) %>%
    pivot_wider(id_cols = c(sgg_cd, sgg_nm), names_from = year_agg, values_from = value, values_fn = function(x) sum(!is.na(x))) %>%
    data.frame



covar_origin_10 = get_basecovar(target_year = 2010)
#covar_origin_o = covar_origin_10
#colnames(covar_origin_10)[grep('^ap_sum_em_', colnames(covar_origin_10))] =
#    str_c('ap_sum_em_', c('co', 'nox', 'sox', 'tsp', 'pm10', 'voc', 'nh3'))
covar_origin_10_consol = clean_consolidated(cleaned_df = covar_origin_10) %>%
    mutate_at(.vars = vars(-geom),
              .funs = list(~ifelse(is.na(.), median(., na.rm = TRUE), .)))



run_smerc_cancertype = function(data = sgg2015, yvar = "Lung_total", sex_b = 'total', ncores = 8) {
  library(parallel)
  cns = colnames(data)[grep(str_c(str_c('(^p_*.*_', sex_b, '$'), '^(r_|ap_)', '^NDVI_)', sep = '|'), colnames(data))]
  print(cns)
  form_pois = as.formula(str_c(yvar, '~', str_c(cns, collapse = '+')))
  reg_pois = glm(formula = form_pois, data= data, family = poisson(link = 'log'))
  cls = parallel::makeCluster(spec = ncores, type = 'PSOCK')
  data_df = st_drop_geometry(data)
  eltest = smerc::elliptic.test(st_coordinates(st_centroid(data)), 
            cases = unlist(data_df[, yvar]), 
            pop = reg_pois$fitted.values,
            shape = c(1, 1.5, 2, 2.5, 3, 4, 5, 6),
            nangle = c(1, 4, 6, 12, 12, 12, 15, 18),
            cl = cls)
  parallel::stopCluster(cls)
  return(eltest)
}

run_smerc_cancertype_pl = function(data = sgg2015, population = 'n_pop_total', yvar = "Lung_total", sex_b = 'total', ncores = 8) {
  library(parallel)
  #cns = colnames(data)[grep(str_c(str_c('(^p_*.*_', sex_b, '$'), '^(r_|ap_)', '^NDVI_)', sep = '|'), colnames(data))]
  data_df = st_drop_geometry(data)
  cls = parallel::makeCluster(spec = ncores, type = 'PSOCK')
  eltest = smerc::elliptic.test(st_coordinates(st_centroid(data)), 
            cases = unlist(data_df[, yvar]), 
            pop = unlist(data_df[, population]),
            shape = c(1, 1.5, 2, 2.5, 3, 4, 5, 6),
            nangle = c(1, 4, 6, 12, 12, 12, 15, 18),
            cl = cls)
  parallel::stopCluster(cls)
  return(eltest)
}
doParallel::stopImplicitCluster()


# smerc cluster to general maps with a tmap object
tmap_smerc = function(basemap, smc, threshold = 2) {
    library(tmap)
    basemap$cluster = NA
    smc_ncl = sapply(smc$clusters, function(x) length(x$locids))
    smc_ncl_addr = (smc_ncl >= threshold)
    # p-value extraction
    smc_pval = sapply(smc$clusters, function(x) x$pvalue)
    smc_pval_addr = (smc_pval <= 0.05)
    if (!is.null(threshold)) {
        smc_pval_addr = smc_pval_addr * smc_ncl_addr
    }
    smc_pval_addr = grep(1, as.integer(smc_pval_addr))

    # loop through
    for (i in seq_len(length(smc_pval_addr))) {
        cl_ilocs = smc$clusters[[smc_pval_addr[i]]]$locids
        basemap$cluster[cl_ilocs] = i
    }
    basemap$cluster = as.factor(basemap$cluster)
    # tmap
    tm_cluster = tm_shape(basemap) +
        tm_fill('cluster', pal = 'Set3', colorNA = 'transparent', showNA = FALSE) +
        tm_borders(col = 'dark grey', lwd = 0.3)
    return(tm_cluster)
}


### Run main code
smerc_lung_t = run_smerc_cancertype(data = covar_origin_10_consol, yvar = 'n_Lung_total', sex_b = 'total', ncores = 16)
smerc_stom_t = run_smerc_cancertype(data = covar_origin_10_consol, yvar = 'n_Stomach_total', sex_b = 'total', ncores = 16)
smerc_lung_m = run_smerc_cancertype(data = covar_origin_10_consol, yvar = 'n_Lung_male', sex_b = 'male', ncores = 16)
smerc_stom_m = run_smerc_cancertype(data = covar_origin_10_consol, yvar = 'n_Stomach_male', sex_b = 'male', ncores = 16)
smerc_lung_f = run_smerc_cancertype(data = covar_origin_10_consol, yvar = 'n_Lung_female', sex_b = 'female', ncores = 16)
smerc_stom_f = run_smerc_cancertype(data = covar_origin_10_consol, yvar = 'n_Stomach_female', sex_b = 'female', ncores = 16)


smerc_lung_tu = run_smerc_cancertype_pl(data = covar_origin_10_consol, pop = 'n_pop_total', yvar = 'n_Lung_total', sex_b = 'total', ncores = 16)
smerc_stom_tu = run_smerc_cancertype_pl(data = covar_origin_10_consol, pop = 'n_pop_total', yvar = 'n_Stomach_total', sex_b = 'total', ncores = 16)
smerc_lung_mu = run_smerc_cancertype_pl(data = covar_origin_10_consol, pop = 'n_pop_male', yvar = 'n_Lung_male', sex_b = 'male', ncores = 16)
smerc_stom_mu = run_smerc_cancertype_pl(data = covar_origin_10_consol, pop = 'n_pop_male', yvar = 'n_Stomach_male', sex_b = 'male', ncores = 16)
smerc_lung_fu = run_smerc_cancertype_pl(data = covar_origin_10_consol, pop = 'n_pop_female', yvar = 'n_Lung_female', sex_b = 'female', ncores = 16)
smerc_stom_fu = run_smerc_cancertype_pl(data = covar_origin_10_consol, pop = 'n_pop_female', yvar = 'n_Stomach_female', sex_b = 'female', ncores = 16)


par(mfcol = c(1,2))
plot(smerc_stom_f)
plot(smerc_stom_fu)

par(mfcol = c(1,2))
plot(smerc_lung_t)
plot(smerc_lung_tu)




### Incidence-Mortality consolidation ####
## Mortality
mor_clw = mor_cl %>%
    filter(cancer_type_e %in% c('Stomach', 'Lung/Bronchus')) %>%
    mutate(sex_e = tolower(sex_e),
           cancer_type_e = ifelse(cancer_type_e == 'Lung/Bronchus', 'Lung', cancer_type_e)) %>%
    rename(n = Ntotal) %>%
    mutate(period = plyr::mapvalues(year_agg, unique(year_agg), seq_len(length(unique(year_agg))))) %>%
    dplyr::select(-year_agg, -r_crude) %>%
    filter(period %in% 1:3) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'n_d_',
                values_from = n)    

## Incidence
inc_clw = inc_cl %>%
    filter(cancer_type_e %in% c('Stomach', 'Lung') & type_inc_e == "Ntotal") %>%
    mutate(sex_e = tolower(sex_e),
           type_inc_e = 'n') %>%
    mutate(period = plyr::mapvalues(year, unique(year), seq_len(length(unique(year))))) %>%
    dplyr::select(-year, -type_inc_e) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'n_i_',
                values_from = n_n)    

## Population (5-year aggregation)
pop_clw = pop_cl %>%
    mutate(sex_e = tolower(sex_e),
           period = plyr::mapvalues(year_agg, unique(year_agg), seq_len(length(unique(year_agg))))) %>%
    dplyr::select(-sex0, -year_agg) %>%
    filter(period %in% 1:3) %>%
    pivot_wider(names_from = c(sex_e, period),
                names_sep = "_",
                names_prefix = 'n_p_',
                values_from = population) %>%
    dplyr::select(-ends_with('_NA'))


## Joined
morinc_clw = mor_clw %>%
    full_join(inc_clw) %>%
    full_join(pop_clw, by = c('sgg_cd_c' = 'sgg_cd')) %>%
    mutate(sgg_cd_c = plyr::mapvalues(sgg_cd_c, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    group_by(sgg_cd_c) %>%
    summarize_all(list(~sum(., na.rm = TRUE))) %>%
    ungroup


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

par(mfcol = c(1,2))
plot(smerc_lung_3dt)
plot(smerc_lung_3dtu)


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


tmap_smerc(covar_origin_10_fc, smerc_lung_3it)
tmap_smerc(covar_origin_10_fc, smerc_lung_3itu)

tmap::tmap_arrange(
tmap_smerc(covar_origin_10_fc, smerc_lung_3dt),
tmap_smerc(covar_origin_10_fc, smerc_lung_3dtu)
)

tmap::tmap_arrange(
    tmap_smerc(covar_origin_10_fc, smerc_lung_3im),
    tmap_smerc(covar_origin_10_fc, smerc_lung_3imu),
    ncol = 2
)

tmap::tmap_arrange(
    tmap_smerc(covar_origin_10_fc, smerc_lung_3if),
    tmap_smerc(covar_origin_10_fc, smerc_lung_3ifu),
    ncol = 2
)
tmap::tmap_arrange(
    tmap_smerc(covar_origin_10_fc, smerc_stom_3if),
    tmap_smerc(covar_origin_10_fc, smerc_stom_3ifu),
    ncol = 2
)

tmap::tmap_arrange(
    tmap_smerc(covar_origin_10_fc, smerc_stom_3dt),
    tmap_smerc(covar_origin_10_fc, smerc_stom_3dtu),
    ncol = 2
)


run_dclust_cancertype = function(
        data, 
        population = 'n_p_total_3', 
        yvar = "n_d_Lung_total_3", 
        sex_b = 'total',
        ncores = 20,
        adjust = TRUE,
        run_glm = FALSE) {
  cns = colnames(data)[grep(str_c(str_c('(^p_*.*_', sex_b, '$'), '^(r_|ap_)', '^NDVI_)', sep = '|'), colnames(data))]
  print(cns)
  options(mc.cores = ncores)
  form_pois = as.formula(str_c(yvar, '~ offset(log(Expected)) +', str_c(cns, collapse = '+')))
  if (!adjust) {
      form_pois = as.formula(str_c(yvar, '~ offset(log(Expected)) + 1'))
  }
  data_df = data %>%
    bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
    st_drop_geometry %>%
    mutate(Expected = !!sym(population) * sum(!!sym(yvar)) / sum(!!sym(population)))
  reg_pois = glm(formula = form_pois, data= data_df, family = poisson(link = 'log'))
  if (run_glm) { return(reg_pois)}
  
  eltest <- DetectClustersModel(data %>% as('Spatial'), 
    #radius = 1000000,
    thegrid = data_df %>% dplyr::select(X, Y), 
    fractpop = 0.5,
    alpha = 0.05, typeCluster = "S", 
    R = NULL, 
    model0 = reg_pois,
    ClusterSizeContribution = population)
  return(eltest)
}


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


# dclust cluster to general maps with a tmap object
tmap_dclust = function(basemap, dclust, threshold = 2, upto_n = NULL, alpha = 0.5) {
    library(tmap)
    dclust_f = dclust %>%
        arrange(-risk) %>%
        dplyr::filter(size >= threshold) %>%
        group_by(size, statistic) %>%
        dplyr::filter(!duplicated(risk)) %>%
        ungroup %>%
        as.data.frame
    basemap_df = basemap %>%
        bind_cols(as.data.frame(st_coordinates(st_centroid(.)))) %>%
        st_drop_geometry %>%
        rename(x = X, y = Y)
    bdclust = get.knclusters(basemap_df, dclust_f)
    basemap$cluster = NA
    
    if (is.null(upto_n)) {
        # p-value extraction
        bdclust_n = bdclust %>%
            unlist %>%
            table %>%
            as.data.frame
        colnames(bdclust_n) = c('id', 'freq')
        basemap$cluster[bdclust_n$id] = bdclust_n$freq
        # tmap
        tm_cluster = tm_shape(basemap) +
            tm_fill('cluster', pal = 'viridis', style = 'cont', colorNA = 'transparent', showNA = FALSE, alpha = alpha) +
            tm_borders(col = 'dark grey', lwd = 0.3)
    }
    
    if (!is.null(upto_n)) {
        bdclust_n = bdclust[seq_len(upto_n)]
        # loop through
        for (i in seq_len(upto_n)) {
            cl_ilocs = bdclust_n[[i]]
            basemap$cluster[cl_ilocs] = i
        }
        basemap$cluster = as.factor(basemap$cluster)
            # tmap
        tm_cluster = tm_shape(basemap) +
            tm_fill('cluster', pal = 'Set3', colorNA = 'transparent', showNA = FALSE, alpha = 0.3) +
            tm_borders(col = 'dark grey', lwd = 0.3)

    }
    return(tm_cluster)
}

tmap_arrange(
    tmap_dclust(covar_origin_10_fc, dclust_lung_3it, upto_n = NULL, alpha = 0.7),
    tmap_dclust(covar_origin_10_fc, dclust_lung_3itu, upto_n = NULL, alpha = 0.7), ncol = 2)
tmap_arrange(
    tmap_dclust(covar_origin_10_fc, dclust_lung_3dt, upto_n = NULL, alpha = 0.7),
    tmap_dclust(covar_origin_10_fc, dclust_lung_3dtu, upto_n = NULL, alpha = 0.7), ncol = 2)


tmap_dclust(covar_origin_10_fc, dclust_lung_3it, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_lung_3itu, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_lung_3dt, upto_n = NULL, alpha = 0.5)
tmap_dclust(covar_origin_10_fc, dclust_lung_3im, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_lung_3dm, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_lung_3if, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_lung_3df, upto_n = NULL, alpha = 0.7)

tmap_dclust(covar_origin_10_fc, dclust_stom_3it, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_stom_3dt, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_stom_3im, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_stom_3dm, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_stom_3if, upto_n = NULL, alpha = 0.7)
tmap_dclust(covar_origin_10_fc, dclust_stom_3df, upto_n = NULL, alpha = 0.7)



dclust_lung_3it %>%
        arrange(-risk) %>%
        filter(size != 1) %>%
        group_by(size, statistic) %>%
        filter(!duplicated(risk)) %>%
        ungroup

par(mfcol = c(2,2))
plot(smerc_lung_3it)
plot(smerc_lung_3itu)
plot(smerc_lung_3dt)
plot(smerc_lung_3dtu)

save(list = ls()[grep('^smerc_*.*_[3]*.*', ls())],
     file = 'Scan_SMERC_2009_2013_Results.RData')
save(list = ls()[grep('^dclust_*.*_[3]*.*', ls())],
     file = 'Scan_DCLUST_2009_2013_Results.RData')
save(list = ls()[grep('^(smerc|dclust|covar_origin_10_)', ls())],
     file = 'Scan_2009_2013_Results.RData')
