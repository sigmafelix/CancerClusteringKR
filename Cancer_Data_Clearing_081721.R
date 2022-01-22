## Tidyverse-version of covariate cleaning
## 
options(repos = 'https://cran.seoul.go.kr/')
if (!require(pacman)) {install.packages('pacman')}
p_load(tidyverse, sf, spdep, DCluster, rflexscan, tmap, dtplyr, fuzzyjoin, readxl, here)

userid = 'isong'
drive = sprintf('/mnt/c/Users/%s/OneDrive/NCC_Project/CancerClustering/', userid)
geopath = sprintf('/mnt/c/Users/%s/OneDrive/Data/Korea/', userid)

# # incidence
# inc <- read.csv(paste(drive, 'cancerInc_sgg.csv', sep = ''), fileEncoding = 'EUC-KR')
# mor_to <- read.csv(paste(drive, "cancerMor_sgg_total.csv", sep = ''), fileEncoding = 'EUC-KR')
# mor_me <- read.csv(paste(drive, "cancerMor_sgg_male.csv", sep = ''), fileEncoding = 'EUC-KR')
# mor_fe <- read.csv(paste(drive, "cancerMor_sgg_female.csv", sep = ''), fileEncoding = 'EUC-KR')

# # total population for the expected values
# # sex code: 0(all), 1(male), 2(female)
# pop = read.csv(paste(drive, 'Midyear_population_1998_2019.csv', sep = ''), fileEncoding = 'UTF-8') %>%
#     mutate(sex0 = plyr::mapvalues(sex0, c(0, 1, 2), c('T', 'M', 'F')),
#            population = as.numeric(population)) %>%
#     pivot_wider(names_from = c('sex0', 'year'), values_from = population)

# colnames(inc) <- c("sex0","sex","cancer_type0","cancer_type","sgg_cd","sgg_nm","year0","year","inc0","type_inc","unit","n","x")
# colnames(mor_to) = c('cause0', 'cause', 'sgg_cd', 'sgg_nm', 'sex0', 'sex', 'type0', 'type', 'unit', paste('Y', 1998:2019, sep=''), 'x')
# colnames(mor_me) = c('cause0', 'cause', 'sgg_cd', 'sgg_nm', 'sex0', 'sex', 'type0', 'type', 'unit', paste('Y', 1998:2019, sep=''), 'x')
# colnames(mor_fe) = c('cause0', 'cause', 'sgg_cd', 'sgg_nm', 'sex0', 'sex', 'type0', 'type', 'unit', paste('Y', 1998:2019, sep=''), 'x')


# inc_cl = inc %>%
#     mutate(n_n = as.numeric(n),
#            n_n = ifelse(is.na(n_n), 0, n_n),
#            cancer_type_e = plyr::mapvalues(cancer_type, unique(cancer_type), c('Stomach', 'Colorectal', 'Liver', 'Lung', 'Breast', 'Cervical', 'Prostate', 'Thyroid')),
#            sex_e = plyr::mapvalues(sex, unique(sex), c('Total', 'Male' ,'Female')),
#            type_inc_e = plyr::mapvalues(type_inc, unique(type_inc), c('N', 'r_crude', 'r_agest')))

# # summary table
# inc_cl_summary = inc_cl %>%
#     group_by(cancer_type_e, year, sex_e, type_inc_e) %>%
#     summarize(n_nmin = min(n_n),
#             n_nmedian = median(n_n),
#             n_nmean = mean(n_n),
#             n_nmax = max(n_n)) %>%
#     ungroup() %>%
#     pivot_longer(cols = 5:8) %>%
#     pivot_wider(names_from = c('sex_e', 'name'))

# write.csv(inc_cl_summary, paste(drive, 'Cancer_Incidence_summary_081821.csv', sep = ''), fileEncoding = 'EUC-KR')

# # code structure check: sgg
# unique(inc$sgg_cd)
# unique(mor_to$sgg_cd)
# unique(mor_me$sgg_cd)
# unique(mor_fe$sgg_cd)

# unique(inc[,c('sgg_cd' ,'sgg_nm', 'year')])
# unique(mor_to[,c('sgg_cd' ,'sgg_nm')])
# unique(inc[,c('sgg_cd' ,'sgg_nm', 'year')])
# unique(inc[,c('sgg_cd' ,'sgg_nm', 'year')])

# # summary by sgg_cd 
# mor_tor = mor_to %>%
#     pivot_longer(cols = Y1998:Y2019) %>%
#     mutate(year = as.integer(gsub('Y', '', name))) %>%
#     filter(!is.na(value)) %>%
#     group_by(sgg_cd, sgg_nm) %>%
#     summarize(N = sum(!is.na(value)),
#               year_min = min(year),
#               year_max = max(year)) %>%
#     ungroup
# vroom::vroom_write(mor_tor, paste(drive, 'Cancer_Mortality_summary_All_081821.csv', sep = ''), delim = ',')

# # Will we analyze the consolidated spatial data for the analysis or spatial data by year?

# ## Data preprocessing ####
# sgg = st_read(paste(geopath, 'SGG_Merge_2000_2016.shp', sep = ''))
# sgg2010 = sgg %>%
#     filter(BASE_YEAR == 2010)


# ## Analysis: pilot in 2010 ####
# # All
# mor_to2010 = mor_to %>%
#     mutate(cancer_type_e = plyr::mapvalues(cause, unique(cause), c('Stomach', 'Colorectal', 'Liver', 'Lung', 'Breast', 'Cervical', 'Prostate'))) %>%
#     pivot_longer(cols = Y1998:Y2019) %>%
#     mutate(year = as.integer(gsub('Y', '', name))) %>%
#     filter(year == 2010)
# # Male
# mor_me2010 = mor_me %>%
#     mutate(cancer_type_e = plyr::mapvalues(cause, unique(cause), c('Stomach', 'Colorectal', 'Liver', 'Lung', 'Breast', 'Cervical', 'Prostate'))) %>%
#     pivot_longer(cols = Y1998:Y2019) %>%
#     mutate(year = as.integer(gsub('Y', '', name))) %>%
#     filter(year == 2010)
# # Female
# mor_fe2010 = mor_fe %>%
#     mutate(cancer_type_e = plyr::mapvalues(cause, unique(cause), c('Stomach', 'Colorectal', 'Liver', 'Lung', 'Breast', 'Cervical', 'Prostate'))) %>%
#     pivot_longer(cols = Y1998:Y2019) %>%
#     mutate(year = as.integer(gsub('Y', '', name))) %>%
#     filter(year == 2010)


# # per 100000 (age-standardized; all sexes)
# mor_to2010_ast0 = mor_to2010 %>%
#     filter(type0 == 'T7' & sex0 == 0) %>%
#     dplyr::select(sgg_cd, sgg_nm, cancer_type_e, value) %>%
#     pivot_wider(names_from = c('cancer_type_e'))
# # Male
# mor_me2010_ast = mor_me2010 %>%
#     filter(type0 == 'T7') %>%
#     dplyr::select(sgg_cd, sgg_nm, cancer_type_e, sex0, value) %>%
#     pivot_wider(names_from = c('cancer_type_e', sex0))
# # Female
# mor_fe2010_ast = mor_fe2010 %>%
#     filter(type0 == 'T7') %>%
#     dplyr::select(sgg_cd, sgg_nm, cancer_type_e, sex0, value) %>%
#     pivot_wider(names_from = c('cancer_type_e', sex0))

# # Counts
# mor_to2010_dc0 = mor_to2010 %>%
#     filter(type0 == 'T1' & sex0 == 0) %>%
#     dplyr::select(sgg_cd, sgg_nm, cancer_type_e, value) %>%
#     pivot_wider(names_from = c('cancer_type_e'))
# # Male
# mor_me2010_dc = mor_me2010 %>%
#     filter(type0 == 'T1') %>%
#     dplyr::select(sgg_cd, sgg_nm, cancer_type_e, sex0, value) %>%
#     pivot_wider(names_from = c('cancer_type_e', sex0))
# # Female
# mor_fe2010_dc = mor_fe2010 %>%
#     filter(type0 == 'T1') %>%
#     dplyr::select(sgg_cd, sgg_nm, cancer_type_e, sex0, value) %>%
#     pivot_wider(names_from = c('cancer_type_e', sex0))




# mor_to2010_sf = sgg2010 %>%
#     left_join(mor_to2010_ast0, by = c('SIGUNGU_CD' = 'sgg_cd')) %>%
#     left_join(mor_me2010_ast, by = c('SIGUNGU_CD' = 'sgg_cd')) %>%
#     left_join(mor_fe2010_ast, by = c('SIGUNGU_CD' = 'sgg_cd')) %>%
#     left_join(pop %>% dplyr::select(1, ends_with('2010')), by = c('SIGUNGU_CD' = 'sgg_cd'))

# mor_dc2010_sf = sgg2010 %>%
#     left_join(mor_to2010_dc0, by = c('SIGUNGU_CD' = 'sgg_cd')) %>%
#     left_join(mor_me2010_dc, by = c('SIGUNGU_CD' = 'sgg_cd')) %>%
#     left_join(mor_fe2010_dc, by = c('SIGUNGU_CD' = 'sgg_cd')) %>%
#     left_join(pop %>% dplyr::select(1, ends_with('2010')), by = c('SIGUNGU_CD' = 'sgg_cd'))



# # missing rate
# # rate
# mor_to2010_sf %>%
#     st_drop_geometry %>%
#     sapply(function(x) 100 * sum(is.na(x))/length(x)) %>%
#     round(., 2)
# # counts
# mor_dc2010_sf %>%
#     st_drop_geometry %>%
#     sapply(function(x) 100 * sum(is.na(x))/length(x)) %>%
#     round(., 2)




# Insurance Premium ####
dbdir = drive
## KCDC code to others
code_conv = read_xlsx(str_c(dbdir, '/행정자치부_심평원_통계청_코드_연계표_211119.xlsx'), sheet = 5)
code_conv_d = code_conv %>% 
  dplyr::select(1, 4, 6, SIGUNGU_PSEUDO, ends_with('2015')) %>%
  rename(SGIS = SGIS_2015,
         SIGUNGU_KCDC = SIGUNGU_PSEUDO) %>%
  mutate(SGIS = str_sub(SGIS, 1, 5)) %>%
  filter(!is.na(SGIS)) %>%
  arrange(SGIS)

pr_filters = c("^(31|32|11|23)", "^(33|34|25|29)", "^(35|36|39|24)", "^(37|38|21|22|26)")
pr_filters_df = pr_filters %>%
    split(.,.) %>%
    lapply(function(x) code_conv_d %>% filter(grepl(x, SGIS)))

# A0211: Gyeonggi ^(31|32|11|23)
# A05011: Chungcheong ^(33|34|25|29) 33*|34*|25*|29*
# A0790: Jeolla-Jeju ^(35|36|39|24) 35*|36*|39*|24*
# A1081: Gyeongsang ^(37|38|21|22|26) 37*|38*|21*|22*|26*
pr_dir = str_c(drive, "Covariates/Premium/")
pr_csvs = list.files(path = pr_dir, full.names = TRUE)
pr_csvs = pr_csvs %>%
    split(.,.) %>%
    lapply(function(x) read.csv(x, fileEncoding = 'CP949') %>% as_tibble %>% rename(SGGFLAG = X.A01.시군구별))


# Conversion
conv_table = read.csv(paste(sprintf("/mnt/c/Users/%s/OneDrive/Data/Korea/", userid), 'SGG_1995_2018_Conversion_Table_201108.csv', sep = ''), fileEncoding = 'EUC-KR')
conv_table_e = conv_table %>%
    filter(from_year >= 1999 & from_year <= 2013) %>%
    dplyr::select(fromcode, tocode)

ng_table1 = tribble(
    ~SGGFLAG,   ~SGIS,  ~NG,
    "A003", "23010",    1,
    "A028", "11020",    1,
    "A091", "32370",    1,
    "A098", "32410",    1,
    "A101", "32380",    1
)
ng_table2 = tribble(
    ~SGGFLAG,   ~SGIS,  ~NG,
    "A039", "33044",    1,
    "A040", "33042",    1,
    "A017", "34340",    1,
    "A026", "34050",    1,
    "A027", "34040",    1,
    "A033", "33350",    1
)
ng_table3 = tribble(
    ~SGGFLAG,   ~SGIS,  ~NG,
    "A047", "35010",    1
)
ng_table4 = tribble(
    ~SGGFLAG,   ~SGIS,  ~NG,
    "A002", "22010",    1,
    "A002", "26010",    1,
    "A003", "22030",    1,
    "A004", "22020",    1,
    "A004", "26030",    1,
    "A008", "22040",    1,
    "A008", "26020",    1,
    "A009", "22050",    1,
    "A009", "26040",    1,
    "A019", "21010",    1,
    "A019", "26010",    1,
    "A020", "21030",    1,
    "A020", "26030",    1,
    "A021", "21020",    1,
    "A022", "21070",    1,
    "A022", "26020",    1,
    "A023", "21080",    1,
    "A023", "26040",    1,
    "A028", "21010",    1,
    "A028", "22010",    1,
    "A029", "21070",    1,
    "A029", "22040",    1,
    "A030", "21030",    1,
    "A030", "22020",    1,
    "A031", "21080",    1,
    "A031", "22050",    1
)

ng_tables = list(
    ng_table1,
    ng_table2,
    ng_table3,
    ng_table4
)


pr_csvs_stringdist = mapply(
    FUN = function(x, y, z) {
        colnames(x) = c("SGGFLAG", "SGGWISE", "Decile_code", "Decile", "Affiliate_code", "Affiliate", "Item_code", "Item", "Unit", str_c("Y", 2013:2020), "XX")
        y2 = y %>% mutate(SGG = str_c(SIGUNGU_NM, ifelse(is.na(SIGUNGU_SUB_NM), "", str_c(" ", SIGUNGU_SUB_NM))))
        x %>%
            stringdist_join(., y2, by = c("SGGWISE" = "SGG"), mode = 'left', method = "jaccard", useBytes = T, max_dist = 0.2) %>%
            filter(SGGWISE == SGG) %>%
            left_join(z) %>%
            filter(is.na(NG)) %>%
            mutate(Y2013 = as.integer(Y2013)) %>%
            dplyr::select(1:10, 19:24) -> xy
        return(xy)
    }, pr_csvs, pr_filters_df, ng_tables, SIMPLIFY = FALSE
)

pr_csvs_con = pr_csvs_stringdist %>%
    do.call(bind_rows, .)

pr_csvs_con_base = pr_csvs_con %>%
    filter(Affiliate_code == "C001") %>%
    #filter(Affiliate_code == "C001" & Item_code == "T01") %>%
    dplyr::select(2, SGIS, SGG, ends_with("_code"), Item, Y2013) %>%
    mutate(Decile_code = plyr::mapvalues(Decile_code, unique(Decile_code), sprintf("D%02d", 0:10)))# %>%
    #pivot_wider(names_from = c(Affiliate_code, Item_code, Decile_code),
    #            names_sep = "_",
    #            values_from = Y2013)

# population weights
pr_csvs_con_pwprem = pr_csvs_con_base %>%
    dplyr::select(2, 4, 6, 8) %>%
    pivot_wider(names_from = Item_code,
                values_from = Y2013) %>%
    group_by(SGIS) %>%
    summarize(n_pwpremium = sum((T04/T04[1])[-1] * T03[-1])) %>%
    ungroup
pr_csvs_con_totalpop = pr_csvs_con_base %>%
    dplyr::select(2, 4, 6, 8) %>%
    pivot_wider(names_from = Item_code,
                values_from = Y2013) %>%
    filter(Decile_code == "D00") %>%
    dplyr::select(1, 5) %>%
    rename(n_insured = T04)

# join population-weighted mean premium and total number of insured
pr_csvs_con_prempop = pr_csvs_con_pwprem %>%
    left_join(pr_csvs_con_totalpop) %>%
    mutate(SGIS = plyr::mapvalues(SGIS, conv_table_e$fromcode, conv_table_e$tocode)) %>%
    group_by(SGIS) %>%
    summarize(n_pwpremium = sum(n_pwpremium * (n_insured / sum(n_insured)))) %>%
    ungroup



# pr_csvs_con_base %>% 
#     filter_at(.vars = vars(4:ncol(.)),
#               any_vars(.>1))
# strategy: weighted mean of weighted means
# - Weighted mean by decile population ratio and average premium
# - Weighted mean of weighted means by total population weights

