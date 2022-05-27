## Population covariates cleaning
## Revised: 04/20/2022
options(repos = "https://cran.seoul.go.kr/")
if (!require(pacman)) {
    install.packages("pacman")
}
p_load(tidyverse, sf, DCluster, tmap, smerc, knitr, readxl, kableExtra, DClusterm, patchwork)
p_load(stars, raster, starsExtra, here, stargazer)

basedir = sprintf('/mnt/c/Users/%s/', username)
rdatafiles = list.files(path = str_c(basedir, 'Documents/GP/'), pattern = '*.RData', full.names = TRUE)
geopath = str_c(basedir, "OneDrive/Data/Korea/")
#homedir <- "/home/felix/"
drive = str_c(basedir, "OneDrive/NCC_Project/CancerClustering/")
geopath = str_c(basedir, "OneDrive/Data/Korea/")
dbdir = drive  # here::here()

# incidence
inc <- read.csv(paste(drive, "cancerInc_sgg.csv", sep = ""), fileEncoding = "EUC-KR")
mor_to <- read.csv(paste(drive, "cancerMor_sgg_total.csv", sep = ""), fileEncoding = "EUC-KR")
mor_me <- read.csv(paste(drive, "cancerMor_sgg_male.csv", sep = ""), fileEncoding = "EUC-KR")
mor_fe <- read.csv(paste(drive, "cancerMor_sgg_female.csv", sep = ""), fileEncoding = "EUC-KR")

conv_table = read.csv(paste(dbdir, 'SGG_before_2022_Conversion_Table_220122_CConly.csv', sep = ''), fileEncoding = 'EUC-KR')
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
#sgg <- st_read(paste(geopath, "SGG_Merge_2000_2016.shp", sep = ""))
#sgg2010 <- sgg %>%
#    filter(BASE_YEAR == 2010)
conv_table_e = conv_table %>%
    filter(to_year >= 1999 & from_year <= 2013) %>%
    dplyr::select(fromcode, tocode)

pop_cl = pop %>%
    filter(year >= 1999 & year <= 2013) %>%
    mutate(sgg_cd = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode),
           year_agg = cut(year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE)) %>%
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
    pivot_wider(id_cols = c(year, sgg_cd_c, cancer_type_e, sex_e), names_from = type_inc_e, values_from = n_n, values_fn = function(x) sum(x, na.rm = TRUE)) %>%
    # note that sex-specific cancers have incorrect crude rates
    left_join(pop_cl, by = c('year' = 'year_agg', 'sgg_cd_c' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    group_by(year, sgg_cd_c, cancer_type_e, sex_e) %>%
    summarize(Ntotal = sum(N, na.rm = T),
              r_crude = 1e5 * (sum(N, na.rm = TRUE) / sum(population, na.rm = T)),
              r_agest = sum(r_agest * (population/sum(population, na.rm = T)), na.rm = TRUE)) %>%
    ungroup %>%
    pivot_longer(cols = 5:7, values_to = 'n_n', names_to = 'type_inc_e')

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
    mutate(year = as.integer(str_sub(name, 2, 5))) %>%
    filter(year >= 1999 & year <= 2013) %>%
    mutate(year_agg = cut(year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE),
           cancer_type_e = plyr::mapvalues(cause, unique(cause), c("Stomach", "Colorectal", "Liver", "Lung/Bronchus", "Breast", "Cervical/Uterine", "Prostate")),
           sex_e = plyr::mapvalues(sex, unique(sex), c('Total', 'Male', 'Female')),
           type_mor_e = plyr::mapvalues(type0, c('T1', 'T4', 'T7'), c('N', 'r_crude', 'r_agest')),
           sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)
        ) %>%
    # not sum
    pivot_wider(id_cols = c(year_agg, sgg_cd_c, cancer_type_e, sex_e), names_from = type_mor_e, values_from = value, values_fn = function(x) sum(x, na.rm = TRUE)) %>%
    left_join(pop_cl, by = c('year_agg' = 'year_agg', 'sgg_cd_c' = 'sgg_cd', 'sex_e' = 'sex_e')) %>%
    group_by(year_agg, sgg_cd_c, cancer_type_e, sex_e) %>%
    summarize(Ntotal = sum(N, na.rm = T),
              r_crude = 1e5 * sum(N, na.rm = T)/sum(population, na.rm = T),
              r_agest = sum(r_agest * (population/sum(population, na.rm = T)), na.rm = TRUE)) %>%
    ungroup %>%
    filter(year_agg %in% c('1999-2003', '2004-2008', '2009-2013'))



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
    filter(cancer_type_e %in% c('Stomach', 'Lung/Bronchus')) %>%
    mutate(sex_e = tolower(sex_e),
           cancer_type_e = ifelse(cancer_type_e == 'Lung/Bronchus', 'Lung', cancer_type_e)) %>%
    rename(n = Ntotal) %>%
    mutate(period = plyr::mapvalues(year_agg, unique(year_agg), seq_len(length(unique(year_agg))))) %>%
    dplyr::select(-year_agg, -r_crude, -r_agest) %>%
    filter(period %in% 1:3) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'n_d_',
                values_from = n)    
## Mortality (age standardized rate)
mor_clwr = mor_cl %>%
    filter(cancer_type_e %in% c('Stomach', 'Lung/Bronchus')) %>%
    mutate(sex_e = tolower(sex_e),
           cancer_type_e = ifelse(cancer_type_e == 'Lung/Bronchus', 'Lung', cancer_type_e)) %>%
    rename(n = Ntotal) %>%
    mutate(period = plyr::mapvalues(year_agg, unique(year_agg), seq_len(length(unique(year_agg))))) %>%
    dplyr::select(-year_agg, -r_crude, -n) %>%
    filter(period %in% 1:3) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'ragest_d_',
                values_from = r_agest)    



## Incidence (Number)
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

inc_clwr = inc_cl %>%
    filter(cancer_type_e %in% c('Stomach', 'Lung') & type_inc_e == "r_agest") %>%
    mutate(sex_e = tolower(sex_e),
           type_inc_e = 'n') %>%
    mutate(period = plyr::mapvalues(year, unique(year), seq_len(length(unique(year))))) %>%
    dplyr::select(-year, -type_inc_e) %>%
    pivot_wider(names_from = c(cancer_type_e, sex_e, period),
                names_sep = "_",
                names_prefix = 'ragest_i_',
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
    dplyr::select(-ends_with('_NA')) %>% ## this point was the last point
    mutate(sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    dplyr::select(-sgg_cd) %>%
    group_by(sgg_cd_c) %>%
    summarize_all(list(~sum(., na.rm = TRUE))) %>%
    ungroup 
    

## Joined
morinc_clw = mor_clw %>%
    full_join(inc_clw) %>%
    full_join(mor_clwr) %>%
    full_join(inc_clwr) %>%
    full_join(pop_clw) %>%
    # Ulleung-gun period 2 incidence fix (012222)
    mutate(n_i_Lung_female_2 = ifelse(sgg_cd_c == 37430, 3, n_i_Lung_female_2)) 
