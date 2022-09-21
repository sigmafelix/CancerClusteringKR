## Population covariates cleaning
## Revised: 09/16/2022
options(repos = "https://cran.seoul.go.kr/")
if (!require(pacman)) {
    install.packages("pacman")
}
p_load(tidyverse, sf, DCluster, tmap, smerc, knitr, readxl, kableExtra, DClusterm, patchwork)
p_load(stars, raster, starsExtra, here, stargazer)

username = 'sigma'
basedir = sprintf('/mnt/c/Users/%s/', username)
rdatafiles = list.files(path = str_c(basedir, 'Documents/GP/'), pattern = '*.RData', full.names = TRUE)
geopath = str_c(basedir, "OneDrive/Data/Korea/")
#homedir <- "/home/felix/"
drivebase = str_c(basedir, "OneDrive/NCC_Project/CancerClustering/")
drive = str_c(drivebase, "Data/Cancer/")
drivepop = str_c(drivebase, "Data/Population/")
geopath = str_c(basedir, "OneDrive/Data/Korea/")
dbdir = drive  # here::here()

# incidence
inc <- read.csv(paste(drive, "cancerInc_sgg.csv", sep = ""), fileEncoding = "EUC-KR")
# mortality
mor_to <- read.csv(paste(drive, "cancerMor_sgg_total.csv", sep = ""), fileEncoding = "EUC-KR")
mor_me <- read.csv(paste(drive, "cancerMor_sgg_male.csv", sep = ""), fileEncoding = "EUC-KR")
mor_fe <- read.csv(paste(drive, "cancerMor_sgg_female.csv", sep = ""), fileEncoding = "EUC-KR")

conv_table = read.csv(paste(drivebase, 'SGG_before_2022_Conversion_Table_220916_CConly.csv', sep = ''), fileEncoding = 'EUC-KR')
# total population for the expected values
# sex code: 0(all), 1(male), 2(female)
pop <- read.csv(paste(drivepop, "Midyear_Population_1998_2019.csv", sep = ""), fileEncoding = "UTF-8") %>%
    mutate(
        sex0 = plyr::mapvalues(sex0, c(0, 1, 2), c("T", "M", "F")),
        sex_e = plyr::mapvalues(sex0, c("T", "M", "F"), c("total", "male", "female")),
        population = as.numeric(population)
    ) %>%
    filter()
pop_wide = pop %>%
    pivot_wider(names_from = c("sex0", "year"), values_from = population)

colnames(inc) <- c("sex0", "sex", "cancer_type0", "cancer_type", "sgg_cd", "sgg_nm", "year0", "year", "inc0", "type_inc", "unit", "n", "x")
colnames(mor_to) <- c("cause0", "cause", "sgg_cd", "sgg_nm", "sex0", "sex", "type0", "type", "unit", paste("Y", 1998:2019, sep = ""), "x")
colnames(mor_me) <- c("cause0", "cause", "sgg_cd", "sgg_nm", "sex0", "sex", "type0", "type", "unit", paste("Y", 1998:2019, sep = ""), "x")
colnames(mor_fe) <- c("cause0", "cause", "sgg_cd", "sgg_nm", "sex0", "sex", "type0", "type", "unit", paste("Y", 1998:2019, sep = ""), "x")


## Data preprocessing ####
#sgg <- st_read(paste(geopath, "SGG_Merge_2000_2016.shp", sep = ""))
#sgg2010 <- sgg %>%
#    filter(BASE_YEAR == 2010)

# sgg mort, pop, and inc data filtering using a conversion table and the spatial dataset
#' Consolidate sgg data per year
#' 
#' @param table_target long data.frame object to be cleaned that should include year and sgg_cd
#' @param sp_data sf object to be referred to select the relevant (or prioritized) sgg code
#' @param table_conv data.frame object including from_year, to_year, fromcode, and tocode.
#' @param year the target year. integer.
consolidate_sgg = function(table_target, sp_data = NULL, table_conv, year = NULL) {
    # flow
    # using table_target,
    # sp_data: 
    table_target_l = table_target %>%
        split(., .$year) %>%
        lapply(function(d) {
            dx = d
            year_u = unique(dx$year)
            if (is.numeric(year_u)) {
                tconv = table_conv %>%
                    filter(from_year <= 2013)
            } else {
                year_us = str_extract_all(year_u, '\\d{4,4}')
                year_u1 = as.integer(year_us[[1]][1])
                year_u2 = as.integer(year_us[[1]][2])
                tconv = table_conv %>%
                    filter(from_year <= 2013)
            }
            # print(tconv)
            for (k in seq_len(nrow(tconv))) {
                tconv_c = tconv[k,]
                tconv_fcode = tconv_c$fromcode
                tconv_tcode = tconv_c$tocode

                dsgg_codes = unique(dx$sgg_cd)
                if (sum(dsgg_codes %in% c(tconv_fcode, tconv_tcode)) >= 2) {
                    dx = dx %>%
                        filter(sgg_cd != tconv_fcode)
                } else {
                    dx = dx
                }
            }
            return(dx)
        })
    table_target_ldf = table_target_l %>%
        do.call(rbind, .)
    return(table_target_ldf)
}

# ttl = consolidate_sgg(inc, table_conv = conv_table_e)
# ttm = consolidate_sgg(mor_to, table_conv = conv_table_e)
# ttp = consolidate_sgg(pop, table_conv = conv_table_e)

# pop: no sgg conversion
# inc : join with pop, then convert sgg code and group_by(sgg_cd_c) and take weighted mean for rates
# mort: join with pop, split N and rates; group_by(sgg_cd, year_agg) then N(sum), rates(mean) pop(sum); 
#       convert sgg code and group_by(sgg_cd_c) and take sum(N) and weighted mean for rates
conv_table_e = conv_table %>%
    filter(to_year >= 1999 & from_year <= 2013) %>%
    dplyr::select(from_year, to_year, fromcode, tocode)

# pre-aggregated for incidence
pop_cla = pop %>%
    filter(year >= 1999 & year <= 2013) %>%
    consolidate_sgg(., table_conv = conv_table_e) %>%
    mutate(sgg_cd = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode),
           year_agg = cut(year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE)) %>%
    group_by(year_agg, sgg_cd, sex0, sex_e) %>%
    summarize(population = sum(population, na.rm = TRUE)) %>%
    ungroup %>%
    dplyr::select(-sex0)

# pre-cleaned for mortality
pop_clm = pop %>%
    filter(year >= 1999 & year <= 2013) %>%
    consolidate_sgg(., table_conv = conv_table_e) %>%
    mutate(sgg_cd = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode),
           year_agg = cut(year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE)) %>%
    group_by(year, sgg_cd, sex0, sex_e) %>%
    summarize(population = sum(population, na.rm = TRUE)) %>%
    ungroup %>%
    dplyr::select(-sex0)


inc_clo <- inc %>%
    filter(n != '-') %>%
    mutate(
        n_n = as.numeric(n),
        n_n = ifelse(is.na(n_n), 0, n_n),
        cancer_type_e = plyr::mapvalues(cancer_type, unique(cancer_type), c("Stomach", "Colorectal", "Liver", "Lung", "Breast", "Cervical", "Prostate", "Thyroid")),
        sex_e = plyr::mapvalues(sex, unique(sex), c("total", "male", "female")),
        type_inc_e = plyr::mapvalues(type_inc, unique(type_inc), c("n", "r_crude", "r_agest"))
    ) %>%
    filter(cancer_type_e %in% c('Stomach', 'Lung')) %>%
    consolidate_sgg(., table_conv = conv_table_e) %>%
    mutate(sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    left_join(pop_cla, by = c('year' = 'year_agg', 'sgg_cd_c' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    dplyr::select(-sgg_cd)

inc_cln = inc_clo %>%
    filter(type_inc_e == 'n') %>%
    # not sum
    # pivot_wider(id_cols = c(year, sgg_cd_c, cancer_type_e, sex_e), names_from = type_inc_e, values_from = n_n, values_fn = function(x) sum(x, na.rm = TRUE)) %>%
    # note that sex-specific cancers have incorrect crude rates
    # left_join(pop_cl, by = c('year' = 'year_agg', 'sgg_cd' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    group_by(year, sgg_cd_c, cancer_type_e, sex_e) %>%
    summarize(n_n = sum(n_n, na.rm = T)) %>%
    ungroup %>%
    mutate(type_inc_e = 'n')

inc_clrate = inc_clo %>%
    filter(type_inc_e != 'n') %>%
    # not sum
    # pivot_wider(id_cols = c(year, sgg_cd, cancer_type_e, sex_e), names_from = type_inc_e, values_from = n_n, values_fn = function(x) mean(x, na.rm = TRUE)) %>%
    # mutate(sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    # note that sex-specific cancers have incorrect crude rates
    # left_join(pop_cla, by = c('year' = 'year_agg', 'sgg_cd_c' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    # pivot_wider(id_cols = c(year, sgg_cd_c, cancer_type_e, sex_e), names_from = type_inc_e, values_from = n_n,
    #             value_fn = ~weighted.mean(.x, (population/sum(population, na.rm = T)))) %>%
    group_by(year, sgg_cd_c, cancer_type_e, sex_e, type_inc_e) %>%
    summarize(n_n = weighted.mean(n_n, population)) %>%
    ungroup #%>%
    #pivot_longer(cols = 6, values_to = 'n_n', names_to = 'type_inc_e')

inc_cl = bind_rows(inc_cln, inc_clrate)

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


## Mortality data cleaning
# precleaning
replace_na = function(v, targ) { v[is.na(v)] = targ ; return(v) }

mor_to = 
    mor_to %>%
        as_tibble %>%
        arrange(sgg_cd, cause0, type0) %>%
        dplyr::select(-x) %>%
        pivot_longer(cols= Y1998:Y2019, names_to = 'year') %>%
        arrange(sgg_cd, year, cause0, type0) %>%
        group_by(sgg_cd, year, cause0) %>%
        # summarize(na2 = sum(is.na(value))) %>%
        mutate(value = replace_na(value, 0)) %>% #ifelse(sum(is.na(value)) == 2, 0, value)) %>%
        ungroup %>%
        pivot_wider(names_from = year, values_from = value)
mor_fe = 
    mor_fe %>%
        as_tibble %>%
        arrange(sgg_cd, cause0, type0) %>%
        dplyr::select(-x) %>%
        pivot_longer(cols= Y1998:Y2019, names_to = 'year') %>%
        arrange(sgg_cd, year, cause0, type0) %>%
        group_by(sgg_cd, year, cause0) %>%
        # summarize(na2 = sum(is.na(value))) %>%
        mutate(value = replace_na(value, 0)) %>% #ifelse(sum(is.na(value)) == 2, 0, value)) %>%
        ungroup %>%
        pivot_wider(names_from = year, values_from = value)
        # filter(sgg_cd == 23320 & cause0 == 26)
mor_me = 
    mor_me %>%
        as_tibble %>%
        arrange(sgg_cd, cause0, type0) %>%
        dplyr::select(-x) %>%
        pivot_longer(cols= Y1998:Y2019, names_to = 'year') %>%
        arrange(sgg_cd, year, cause0, type0) %>%
        group_by(sgg_cd, year, cause0) %>%
        # summarize(na2 = sum(is.na(value))) %>%
        mutate(value = replace_na(value, 0)) %>% #ifelse(sum(is.na(value)) == 2, 0, value)) %>%
        ungroup %>%
        pivot_wider(names_from = year, values_from = value)
        # filter(sgg_cd == 23320 & cause0 == 26)



## Main
mor_clo = bind_rows(mor_to, mor_me) %>%
    bind_rows(mor_fe) %>%
    pivot_longer(cols = Y1998:Y2019) %>%
    mutate(year = as.integer(str_sub(name, 2, 5))) %>%
    filter(year >= 1999 & year <= 2013) %>%
    mutate(year_agg = cut(year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE),
           cancer_type_e = plyr::mapvalues(cause, unique(cause), c("Stomach", "Colorectal", "Liver", "Lung", "Breast", "Cervical/Uterine", "Prostate")),
           sex_e = plyr::mapvalues(sex, unique(sex), c('total', 'male', 'female')),
           type_mor_e = plyr::mapvalues(type0, c('T1', 'T4', 'T7'), c('n', 'r_crude', 'r_agest'))
        )  %>%
    filter(cancer_type_e %in% c('Stomach', 'Lung')) %>%
    group_by(sgg_cd, sex_e, cancer_type_e, type_mor_e, year, name) %>%
    filter(sum(value) != 0) %>%
    ungroup %>%
    consolidate_sgg(., table_conv = conv_table_e) %>%
    mutate(sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    left_join(pop_clm, by = c('year' = 'year', 'sgg_cd_c' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    dplyr::select(-sgg_cd)


# mor_clo = bind_rows(mor_to, mor_me) %>%
#     bind_rows(mor_fe) %>%
#     pivot_longer(cols = Y1998:Y2019) %>%
#     mutate(year = as.integer(str_sub(name, 2, 5))) %>%
#     filter(year >= 1999 & year <= 2013) %>%
#     mutate(year_agg = cut(year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE),
#            cancer_type_e = plyr::mapvalues(cause, unique(cause), c("Stomach", "Colorectal", "Liver", "Lung", "Breast", "Cervical/Uterine", "Prostate")),
#            sex_e = plyr::mapvalues(sex, unique(sex), c('total', 'male', 'female')),
#            type_mor_e = plyr::mapvalues(type0, c('T1', 'T4', 'T7'), c('n', 'r_crude', 'r_agest'))
#         )  %>%
#     filter(cancer_type_e %in% c('Stomach', 'Lung')) %>%
#     consolidate_sgg(., table_conv = conv_table_e)

mor_cln = mor_clo %>%
    filter(type_mor_e == "n") %>%
    # not sum
    pivot_wider(id_cols = c(year_agg, sgg_cd_c, cancer_type_e, sex_e), names_from = type_mor_e, values_from = value, values_fn = function(x) sum(x, na.rm = TRUE)) %>%
    # mutate(sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    # left_join(pop_cl, by = c('year_agg' = 'year_agg', 'sgg_cd' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    group_by(year_agg, sgg_cd_c, cancer_type_e, sex_e) %>%
    summarize(n_n = sum(n, na.rm = T)) %>%
    ungroup %>%
    filter(year_agg %in% c('1999-2003', '2004-2008', '2009-2013')) %>%
    # pivot_longer(cols = 5, values_to = 'n_n', names_to = 'type_mor_e')
    mutate(type_mor_e = 'n')

mor_clrate = mor_clo %>%
    filter(type_mor_e != "n") %>%
    # not sum
    group_by(year_agg, sgg_cd_c, cancer_type_e, sex_e, type_mor_e) %>%
    summarize(n_n = weighted.mean(value, population)) %>%
    ungroup #%>%

    # pivot_wider(id_cols = c(year_agg, sgg_cd, cancer_type_e, sex_e), names_from = type_mor_e, values_from = value, values_fn = function(x) mean(x, na.rm = TRUE)) %>%
    # # mutate(sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    # left_join(pop_clm, by = c('year_agg' = 'year_agg', 'sgg_cd_c' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    # group_by(year_agg, sgg_cd_c, cancer_type_e, sex_e) %>%
    # summarize(r_crude = sum(r_crude * (population/sum(population, na.rm = T)), na.rm = TRUE),
    #           r_agest = sum(r_agest * (population/sum(population, na.rm = T)), na.rm = TRUE)) %>%
    # ungroup %>%
    # filter(year_agg %in% c('1999-2003', '2004-2008', '2009-2013')) %>%
    # pivot_longer(cols = 5:6, values_to = 'n_n', names_to = 'type_mor_e')

mor_cl = bind_rows(mor_cln, mor_clrate)

# inc_clgg_n = inc_cl %>%
#     filter(cancer_type_e %in% c('Lung', 'Stomach') & type_inc_e == 'Ntotal') %>%
#     ggplot(data = .,
#            mapping = aes(x = year, y = n_n, group = interaction(year, cancer_type_e, sex_e))) +
#         facet_grid(sex_e ~ cancer_type_e) +
#         geom_boxplot() +
#         theme_bw() +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#         #scale_y_continuous(limits = c(0, 225)) +
#         labs(subtitle = 'Crude cancer incidence counts by five-year periods and types')

# inc_clgg_r = inc_cl %>%
#     filter(cancer_type_e %in% c('Lung', 'Stomach') & type_inc_e == 'r_crude') %>%
#     ggplot(data = .,
#            mapping = aes(x = year, y = n_n, group = interaction(year, cancer_type_e, sex_e))) +
#         facet_grid(sex_e ~ cancer_type_e) +
#         geom_boxplot() +
#         theme_bw() +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#         #scale_y_continuous(limits = c(0, 225)) +
#         labs(subtitle = 'Crude cancer incidence rates per 100,000 persons by five-year periods and types')


# mor_clgg_n = mor_cl %>%
#     filter(cancer_type_e %in% c('Lung/Bronchus', 'Stomach')) %>%
#     #left_join(pop, by = c('sgg_cd' = 'sgg_cd', 'year' = 'year', 'sex_e' = 'sex_e')) %>%
#     #group_by(sgg_cd, cause, cancer_type_e, sex_e, year_agg) %>%
#     #summarize(rate_100k_5yr = 1e5 * (sum(value, na.rm = T)/sum(population, na.rm = T))) %>%
#     #ungroup %>%
#     ggplot(data = .,
#            mapping = aes(x= year_agg, y = Ntotal, group = interaction(year_agg, cancer_type_e, sex_e))) +
#         facet_grid(sex_e ~ cancer_type_e) +
#         geom_boxplot() +
#         theme_bw() +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#         #scale_y_continuous(limits = c(0, 225)) +
#         labs(subtitle = 'Cancer mortality counts by years and types')


# mor_clgg_r = mor_cl %>%
#     filter(cancer_type_e %in% c('Lung/Bronchus', 'Stomach')) %>%
#     #left_join(pop, by = c('sgg_cd' = 'sgg_cd', 'year' = 'year', 'sex_e' = 'sex_e')) %>%
#     #group_by(sgg_cd, cause, cancer_type_e, sex_e, year_agg) %>%
#     #summarize(rate_100k_5yr = 1e5 * (sum(value, na.rm = T)/sum(population, na.rm = T))) %>%
#     #ungroup %>%
#     ggplot(data = .,
#            mapping = aes(x= year_agg, y = r_crude, group = interaction(year_agg, cancer_type_e, sex_e))) +
#         facet_grid(sex_e ~ cancer_type_e) +
#         geom_boxplot() +
#         theme_bw() +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#         #scale_y_continuous(limits = c(0, 225)) +
#         labs(subtitle = 'Crude cancer mortality rates per 100,000 persons by years and types')



# inc_cl %>% 
#     dplyr::select(year, year, sgg_cd, sgg_nm) %>%
#     pivot_wider(id_cols = c(sgg_cd, sgg_nm), names_from = year, values_from = year0, values_fn = length) %>%
#     data.frame

# mor_cl %>% 
#     dplyr::select(year_agg, value, sgg_cd, sgg_nm) %>%
#     pivot_wider(id_cols = c(sgg_cd, sgg_nm), names_from = year_agg, values_from = value, values_fn = function(x) sum(!is.na(x))) %>%
#     data.frame


### Incidence-Mortality consolidation ####
## Mortality (number)
mor_clw = mor_cl %>%
    filter(type_mor_e == "n") %>%
    # mutate(sex_e = tolower(sex_e),
    #        cancer_type_e = ifelse(cancer_type_e == 'Lung/Bronchus', 'Lung', cancer_type_e)) %>%
    # rename(n = Ntotal) %>%
    mutate(period = plyr::mapvalues(year_agg, unique(year_agg), seq_len(length(unique(year_agg))))) %>%
    dplyr::select(-year_agg, -type_mor_e) %>%
    filter(period %in% 1:3) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'n_d_',
                values_from = n_n)    
## Mortality (age standardized rate)
mor_clwasr = mor_cl %>%
    filter(type_mor_e == "r_agest") %>%
    # mutate(sex_e = tolower(sex_e),
    #        cancer_type_e = ifelse(cancer_type_e == 'Lung/Bronchus', 'Lung', cancer_type_e)) %>%
    # rename(n = Ntotal) %>%
    mutate(period = plyr::mapvalues(year_agg, unique(year_agg), seq_len(length(unique(year_agg))))) %>%
    dplyr::select(-year_agg, -type_mor_e) %>%
    filter(period %in% 1:3) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'ragest_d_',
                values_from = n_n)    
## Mortality (crude rate)
mor_clwcrr = mor_cl %>%
    filter(type_mor_e == "r_crude") %>%
    # mutate(sex_e = tolower(sex_e),
    #        cancer_type_e = ifelse(cancer_type_e == 'Lung/Bronchus', 'Lung', cancer_type_e)) %>%
    # rename(n = Ntotal) %>%
    mutate(period = plyr::mapvalues(year_agg, unique(year_agg), seq_len(length(unique(year_agg))))) %>%
    dplyr::select(-year_agg, -type_mor_e) %>%
    filter(period %in% 1:3) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'rcrude_d_',
                values_from = n_n)    



## Incidence (Number)
inc_clw = inc_cl %>%
    filter(type_inc_e == "n") %>%
    # mutate(sex_e = tolower(sex_e),
    #        type_inc_e = 'n') %>%
    mutate(period = plyr::mapvalues(year, unique(year), seq_len(length(unique(year))))) %>%
    dplyr::select(-year, -type_inc_e) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'n_i_',
                values_from = n_n)    

inc_clwasr = inc_cl %>%
    filter(type_inc_e == "r_agest") %>%
    # mutate(sex_e = tolower(sex_e),
    #        type_inc_e = 'r_agest') %>%
    mutate(period = plyr::mapvalues(year, unique(year), seq_len(length(unique(year))))) %>%
    dplyr::select(-year, -type_inc_e) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'ragest_i_',
                values_from = n_n)    

inc_clwcrr = inc_cl %>%
    filter(type_inc_e == "r_crude") %>%
    # mutate(sex_e = tolower(sex_e),
    #        type_inc_e = 'r_crude') %>%
    mutate(period = plyr::mapvalues(year, unique(year), seq_len(length(unique(year))))) %>%
    dplyr::select(-year, -type_inc_e) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'rcrude_i_',
                values_from = n_n)    



## Population (5-year aggregation)
pop_clw = pop_cla %>%
    mutate(sex_e = tolower(sex_e),
           period = plyr::mapvalues(year_agg, unique(year_agg), seq_len(length(unique(year_agg))))) %>%
    dplyr::select(-year_agg) %>%
    #filter(period %in% 1:3) %>%
    pivot_wider(names_from = c(sex_e, period),
                names_sep = "_",
                names_prefix = 'n_p_',
                values_from = population) %>%
    dplyr::select(-ends_with('_NA')) %>% ## this point was the last point
    mutate(sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    dplyr::select(-sgg_cd) %>%
    # group_by(sgg_cd_c) %>%
    # summarize_all(list(~sum(., na.rm = TRUE))) %>%
    ungroup 
    

## Joined
morinc_clw = mor_clw %>%
    full_join(inc_clw) %>%
    full_join(mor_clwcrr) %>%
    full_join(inc_clwcrr) %>%
    full_join(mor_clwasr) %>%
    full_join(inc_clwasr) %>%
    full_join(pop_clw) %>%
    # Ulleung-gun period 2 incidence fix (012222)
    mutate(n_i_Lung_female_2 = ifelse(sgg_cd_c == 37430, 3, n_i_Lung_female_2)) 


load(str_c("/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering", '/Manuscript/Clustering_Base_sf_062422.RData'))

fix_covar_fc = function(cfc, period, indf) {
    indf = indf %>%
        dplyr::select(sgg_cd_c, ends_with(str_c("_", period)))
    cfc_sub = cfc %>%
        dplyr::select(-matches(str_c("^(n_|ragest_)*.*(total|male|female)*.*(_", period, ")$"))) %>%
        left_join(indf)
    return(cfc_sub)
}

covar_origin_00_fco = covar_origin_00_fc
covar_origin_05_fco = covar_origin_05_fc
covar_origin_10_fco = covar_origin_10_fc

covar_origin_00_fc = fix_covar_fc(covar_origin_00_fc, 1, morinc_clw)
covar_origin_05_fc = fix_covar_fc(covar_origin_05_fc, 2, morinc_clw)
covar_origin_10_fc = fix_covar_fc(covar_origin_10_fc, 3, morinc_clw)

# ## Fix NDVI
# p_load(MODIStsp, terra, exactextractr)
# # MODIStsp()
# mtpath10 = "/mnt/d/MODIStsp/VI_16Days_250m_v6/Time_Series/RData/Terra/NDVI/MOD13Q1_NDVI_1_2010_353_2010_RData.RData"
# mtpath05 = "/mnt/d/MODIStsp2005/VI_16Days_250m_v6/Time_Series/RData/Terra/NDVI/MOD13Q1_NDVI_1_2005_353_2005_RData.RData" 
# mtpath00 = "/mnt/d/MODIStsp2000/VI_16Days_250m_v6/Time_Series/RData/Terra/NDVI/MOD13Q1_NDVI_49_2000_353_2000_RData.RData" 

# process_ndvi = function(rdpath, targetpath) {
#     load(rdpath)
#     spr = rast(raster_ts)
#     cat("Computing median cells...\n")
#     spr_med = terra::tapp(spr, rep(1, nlyr(spr)), 'median', cores = 8)
#     terra::writeRaster(spr_med, targetpath, overwrite = TRUE)
#     cat("Finished processing.\n")
#     rm(raster_ts, spr, spr_med)
#     cat("Finished postprocessing sanitization.\n")
# }

# get_ndvi_avg_median = function(fc, ndvipath) {
#     ndvi_pnt = rast(ndvipath)
#     ndvi_avg = exactextractr::exact_extract(ndvi_pnt, fc, 'mean')
#     cat("Average value is computed.\n")
#     return(ndvi_avg)
# }

# ndvi_path_10 = str_c(drivebase, "Data/NDVI_2010_Y1_median.tif")
# ndvi_path_05 = str_c(drivebase, "Data/NDVI_2005_Y1_median.tif")
# ndvi_path_00 = str_c(drivebase, "Data/NDVI_2000_Y1_median.tif")

# # Process NDVI time series
# process_ndvi(mtpath10, ndvi_path_10)
# process_ndvi(mtpath05, ndvi_path_05)
# process_ndvi(mtpath00, ndvi_path_00)

# # Average
# ndvi_avg_10 = get_ndvi_avg_median(covar_origin_10_fc, ndvi_path_10)
# ndvi_avg_05 = get_ndvi_avg_median(covar_origin_05_fc, ndvi_path_05)
# ndvi_avg_00 = get_ndvi_avg_median(covar_origin_00_fc, ndvi_path_00)

# covar_origin_10_fc = covar_origin_10_fc %>%
#     mutate(NDVI_mean = ndvi_avg_10)
# covar_origin_05_fc = covar_origin_05_fc %>%
#     mutate(NDVI_mean = ndvi_avg_05)
# covar_origin_00_fc = covar_origin_00_fc %>%
#     mutate(NDVI_mean = ndvi_avg_00)

save(covar_origin_00_fc, covar_origin_05_fc, covar_origin_10_fc,
     file = "/mnt/c/Users/sigma/OneDrive/NCC_Project/CancerClustering/Manuscript/Clustering_Base_sf_091522.RData",
     compress = 'xz', compression_level = 9)
