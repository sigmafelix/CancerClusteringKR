#' ---
#' output:
#'     html_document:
#'         toc: true
#'         toc_float: true
#'     pdf_document:
#'         toc: true
#' ---

## Population data exploration
## Objective
## - Incidence and mortality: for the entire country, you should have the total counts
## - Using 5-year age groups, identify low count group then
## - Plot district-level maps for low count groups

#+ setup, include=F, echo=F
source('../Base/base_functions.R')
options(repos = 'https://cran.seoul.go.kr')
if (!require(pacman)) { install.packages('pacman') } 

p_load(tidyverse, sf, tmap, stargazer, smerc, DClusterm, kableExtra, patchwork, rmapshaper, spdep)

username = 'sigma'
basedir = sprintf('/mnt/c/Users/%s/', username)
#rdatafiles = list.files(path = str_c(basedir, 'Documents/GP/'), pattern = '*.RData', full.names = TRUE)
geopath = str_c(basedir, "OneDrive/Data/Korea/")
drive = str_c(basedir, "OneDrive/NCC_Project/CancerClustering/")
geopath = str_c(basedir, "OneDrive/Data/Korea/")
dbdir = drive  
rdsdir = sprintf("/mnt/c/Users/%s/OneDrive/NCC_Project/CancerClustering/", username)

exceldir = str_c(drive, "/Data/Cancer/")

## Part 1: histogram for period by age and sex
age_inc = readxl::read_excel(str_c(exceldir, "Incidence_Periods_Summary.xlsx"), sheet = 1) %>%
    filter(Sex != "All") %>%
    pivot_longer(cols = 4:6)
    
age_mor = readxl::read_excel(str_c(exceldir, "Mortality_Periods_Summary.xlsx"), sheet = 1) %>%
    dplyr::select(1:3, starts_with("Period")) %>%
    filter(Sex != 'All') %>%
    pivot_longer(cols = 4:6)

# plot
age_inc_gg = 
    ggplot(data = age_inc,
           mapping = aes(x = Age, y = value, fill = Sex, group = Sex)) +
    geom_col(position = 'dodge') +
    facet_grid(name~Cancer) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16),
          axis.text.x = element_text(hjust = 1, angle = 90),
          strip.text = element_text(size = 15),
          legend.position = 'top') +
    coord_flip() +
    labs(title = "Incidence by age and sex",
         caption = "Periods 1-3 are 1999-2003, 2004-2008, and 2009-2013, respectively")
age_mor_gg = 
    ggplot(data = age_mor,
           mapping = aes(x = Age, y = value, fill = Sex, group = Sex)) +
    geom_col(position = 'dodge') +
    facet_grid(name~Cancer) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16),
          axis.text.x = element_text(hjust = 1, angle = 90),
          strip.text = element_text(size = 15),
          legend.position = 'top') +
    coord_flip() +
    labs(title = "Mortality by age and sex",
         caption = "Periods 1-3 are 1999-2003, 2004-2008, and 2009-2013, respectively")

#+ setup, include=F, echo=F, fig.width=7, fig.height=9
age_inc_gg
age_mor_gg




#+ sgg cleaning, echo=F, fig.width=7, fig.height=9
age_pop = readxl::read_excel(str_c(drive, "Data/Population/MidPopulation_1999_2013_age.xlsx"))
colnames(age_pop) = c('sggcd', 'sggnm', 'sexcd', 'sex', 'agecd', 'age', 'year', 'population')
age_pop_clean = age_pop %>%
    mutate(year = as.numeric(substr(year, 1, 4)),
           period = cut(year, c(1998, 2003, 2008, 2013), labels = 1:3, right = TRUE))

# conversion table
conv_table = read.csv(paste(geopath, 'SGG_1995_2018_Conversion_Table_201108.csv', sep = ''), fileEncoding = 'EUC-KR')
conv_table_e = conv_table %>%
    filter(from_year >= 1999 & from_year <= 2013) %>%
    dplyr::select(fromcode, tocode)

# Excluding integrated cities
sgg_excl = c(38110, 37010, 35010, 31100, 31050, 31040, 31020, 31010)

age_pop_clean_agg = age_pop_clean %>%
    mutate(
        sgg_cd_c = plyr::mapvalues(sggcd, conv_table_e$fromcode, conv_table_e$tocode)
        ) %>%
    group_by(sgg_cd_c, period, sex, sexcd, age, agecd) %>%
    summarize(population = sum(population, na.rm = T)) %>%
    ungroup %>%
    mutate(sex = plyr::mapvalues(sex, unique(sex), c('Male', 'Female')))


load(str_c(drive, '/Manuscript/Clustering_Base_sf_042022.RData'))

sgg_poly = covar_origin_10_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430) & !grepl('^(39)', sgg_cd_c)) %>%
    ms_simplify(keep = 0.125, keep_shapes = TRUE)
age_pop_clean_agg %>% dplyr::select(age, agecd) %>% unique %>% arrange(agecd) %>% data.frame
age_pop_clean_w = age_pop_clean_agg %>%
    mutate(agesex = str_c(sex, "_", agecd, "_", period)) %>%
    filter(agecd %in% c('160','180','190','210', '370')) %>%
    dplyr::select(sgg_cd_c, agesex, population) %>%
    pivot_wider(names_from = agesex, values_from = population)
