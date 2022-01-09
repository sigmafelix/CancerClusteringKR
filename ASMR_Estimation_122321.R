## Age-Standardized Mortality Ratio Estimation
## VERSION: 01092022
## DESCRIPTION: get district-level incidence and death counts by age stratum
## METHOD: integer programming
## TODO: validate the approach; ultimately to make a national ("one-shot") model
library(pacman)
p_load(tidyverse, lpSolve, lpsymphony, sf, spdep, DCluster, tmap, smerc, knitr, readxl, kableExtra, DClusterm, patchwork)

homedir <- "/home/felix/"
#homedir = '/mnt/c/Users/sigma/'

drive <- str_c(homedir, "OneDrive/NCC_Project/CancerClustering/")
geopath <- str_c(homedir, "OneDrive/Data/Korea/")

# incidence
inc <- read.csv(paste(drive, "cancerInc_sgg.csv", sep = ""), fileEncoding = "EUC-KR")
mor_to <- read.csv(paste(drive, "cancerMor_sgg_total.csv", sep = ""), fileEncoding = "EUC-KR")
mor_me <- read.csv(paste(drive, "cancerMor_sgg_male.csv", sep = ""), fileEncoding = "EUC-KR")
mor_fe <- read.csv(paste(drive, "cancerMor_sgg_female.csv", sep = ""), fileEncoding = "EUC-KR")
colnames(mor_to) = colnames(mor_me) = colnames(mor_fe) = 
    c("cause0", "cause", "sgg_cd", "sgg_nm", "sex0", "sexnm", "type0", "type", "unit", paste("Y", 1998:2019, sep = ""), "x")

spop <- read.csv(paste(drive, 'StandardPopulation_1998_2020.csv', sep = ""), fileEncoding = 'EUC-KR')
colnames(spop) = c('sggnm', 'sexnm', 'age', paste('Y', 2020:1998, sep = ''))
spop = spop %>%
    pivot_longer(cols = Y2020:Y1998, names_to = 'year') %>%
    mutate(year = as.integer(gsub('Y', '', year)),
           sex_e = plyr::mapvalues(sexnm, c('계', '남자', '여자'), c('Total', 'Male', 'Female')))

# sgg_pop
sgg_pop_to <- read.csv(paste(drive, "MidPopulation_Total_1998_2013.csv", sep = ""), fileEncoding = "EUC-KR")
sgg_pop_me <- read.csv(paste(drive, "MidPopulation_Male_1998_2013.csv", sep = ""), fileEncoding = "EUC-KR")
sgg_pop_fe <- read.csv(paste(drive, "MidPopulation_Female_1998_2013.csv", sep = ""), fileEncoding = "EUC-KR")
colnames(sgg_pop_to) = colnames(sgg_pop_me) = colnames(sgg_pop_fe) =
    c('sgg_cd', 'sggnm', 'sexcd', 'sexnm', 'agecd', 'agenm', 'itemcd', 'itemnm', 'unit', paste('Y', 1998:2013, sep = ""), "x")

sgg_pop_to = sgg_pop_to %>%
    dplyr::select(-c(ncol(.))) %>%
    pivot_longer(cols = Y1998:Y2013, names_to = 'year') %>%
    mutate(year = as.integer(gsub('Y', '', year)))
sgg_pop_me = sgg_pop_me %>%
    dplyr::select(-c(ncol(.))) %>%
    pivot_longer(cols = Y1998:Y2013, names_to = 'year') %>%
    mutate(year = as.integer(gsub('Y', '', year)))
sgg_pop_fe = sgg_pop_fe %>%
    dplyr::select(-c(ncol(.))) %>%
    pivot_longer(cols = Y1998:Y2013, names_to = 'year') %>%
    mutate(year = as.integer(gsub('Y', '', year)))
sgg_pop_cl = bind_rows(sgg_pop_to, sgg_pop_me) %>%
    bind_rows(sgg_pop_fe) %>%
    mutate(sex_e = plyr::mapvalues(sexcd, 0:2, c('Total', 'Male', 'Female')))



#conversion table
conv_table = read.csv(paste(geopath, 'SGG_1995_2018_Conversion_Table_201108.csv', sep = ''), fileEncoding = 'EUC-KR')
conv_table_e = conv_table %>%
    filter(from_year >= 1999 & from_year <= 2013) %>%
    dplyr::select(fromcode, tocode)

# Excluding integrated cities
sgg_excl = c(38110, 37010, 35010, 31100, 31050, 31040, 31020, 31010)

mor_cl = bind_rows(mor_to, mor_me) %>%
    bind_rows(mor_fe) %>%
    pivot_longer(cols = Y1998:Y2019) %>%
    mutate(year = as.integer(str_sub(name, 2, 5)),
           year_agg = cut(year, breaks = c(1998,2003,2008,2013,2018,2019), labels = c('1999-2003', '2004-2008', '2009-2013', '2014-2018', '2019'), right = TRUE),
           cancer_type_e = plyr::mapvalues(cause, unique(cause), c("Stomach", "Colorectal", "Liver", "Lung/Bronchus", "Breast", "Cervical/Uterine", "Prostate")),
           sex_e = plyr::mapvalues(sexnm, unique(sexnm), c('Total', 'Male', 'Female')),
           type_mor_e = plyr::mapvalues(type0, c('T1', 'T4', 'T7'), c('N', 'r_crude', 'r_agest')),
           sgg_cd_c = plyr::mapvalues(sgg_cd, conv_table_e$fromcode, conv_table_e$tocode)
        ) %>%
    filter(!sgg_cd %in% sgg_excl)

find_mort_n_unit = function(mort, sex, sggcd, pop_ref, pop_sgg, year_f, cancertype, solve_engine = 'lpsolve') {
    # input: total population (reference) by year, age-stratified population by age range
    # output: found integers by age ranges by district
    mort_year = mort %>%
        filter(cancer_type_e == cancertype & sex_e == sex & year == year_f) %>%
        filter(sgg_cd == sggcd) %>%
        split(., .$type0)
    #print(mort_year)
    pop_ref_a = pop_ref %>% 
        filter(year == year_f & sex_e == sex)
    pop_sgg = pop_sgg %>%
        filter(sex_e == sex & year == year_f & sgg_cd == sggcd)
    print(mort_year)
    print(pop_ref_a)
    print(pop_sgg)
}    

#find_mort_n_unit(mor_cl, 'Total', 11010, spop, sgg_pop_cl, 2000, 'Stomach', 'lpsymphony')


# Main function ####
find_mort_n = function(mort, sex, sggcd, pop_ref, pop_sgg, year_f, cancertype, solve_engine = 'lpsolve', varwise = FALSE, setmax = TRUE) {
    # input: total population (reference) by year, age-stratified population by age range
    # output: found integers by age ranges by district
    mort_year = mort %>%
        filter(cancer_type_e == cancertype & sex_e == sex & year == year_f) %>%
        filter(sgg_cd == sggcd) %>%
        split(., .$type0) %>% #T1: N; T4: r_crude; T7: r_agestd
        lapply(function(x) x$value)
    #print(mort_year)
    # 2005 instead of year_f , sex_e 'Total' instead of sex
    pop_ref_a = pop_ref %>% 
        filter(year == 2005 & sex_e == 'Total') %>%
        .$value
    # 2005 instead of year_f 
    pop_sgg = pop_sgg %>%
        filter(year == year_f & sex_e == sex & sgg_cd == sggcd) %>%
        .$value
    
    # recheck
    pop_w = pop_ref_a / sum(pop_ref_a)
    pop_wc = (1e5) / sum(pop_sgg)
    pop_sgg_coef = (1e5) / pop_sgg
    pop_sgg_coef_agest = pop_sgg_coef * pop_w
    print(pop_wc)
    #print(pop_sgg_coef)
    print(pop_sgg_coef_agest)
    
    # maximize
    cmat = rbind(rep(1, 18),
                #pop_wc, 
                #pop_wc, 
                pop_sgg_coef_agest,
                pop_sgg_coef_agest)
    cmat_ext = rbind(cmat,
                diag(nrow = 18))
    # constraints of all-time average death counts by 5-year ranges
    if (cancertype == 'Stomach') {
        const_agebound = c(1,1,0,1,1,2,4,5,6,7,15,20,20,15,15,12,12,10) 
        if (sex != 'Total') {
            const_agebound = ceiling(const_agebound * 0.66)
        }
    } else if (cancertype == 'Lung/Bronchus') {
        const_agebound = c(1,1,1,1,1,1,1,2,3,4,5,12,20,30,30,36,30,20)
        if (sex != 'Total') {
            const_agebound = ceiling(const_agebound * 0.66)
        }
    }

    directions = c('==', 
                   #'>=', '<=', 
                   '>=', '<=')
    directions_ext = c(directions,
                    rep('<=', 18))
    rhss = c(mort_year[[1]], 
            #mort_year[[2]] - 0.04, mort_year[[2]] + 0.04, 
            mort_year[[3]] - 0.04, mort_year[[3]] + 0.04)
    rhss_ext = c(rhss, const_agebound)
    boundlist = list(upper = list(ind = 1:18, val = c(rep(10,5), rep(36,13))), 
                     lower = list(ind = 1:18, val = rep(0, 18)))
    #print(cmat)
    if (solve_engine == 'lpsolve') {
        lpref = lpSolve::lp(direction = if (setmax) "max" else "min",
                            all.int = TRUE,
                            objective.in = rep(1, 18),
                            const.mat = if (varwise) cmat_ext else cmat,
                            const.dir = if (varwise) directions_ext else directions,
                            const.rhs = if (varwise) rhss_ext else rhss)
        nvec = lpref$solution
    } else if (solve_engine == 'lpsymphony') {
        lpref = lpsymphony::lpsymphony_solve_LP(
                            max = setmax,
                            types = rep("I", 18),
                            #bounds = boundlist,
                            obj = rep(1, 18),
                            mat = if (varwise) cmat_ext else cmat,
                            dir = if (varwise) directions_ext else directions,
                            rhs = if (varwise) rhss_ext else rhss)
        #return(lpref)
        nvec = lpref$solution

    } else if (solve_engine == 'glpk') {
        lpref = Rglpk::Rglpk_solve_LP(
                            max = setmax,
                            types = rep("I", 18),
                            #bounds = boundlist,
                            obj = rep(1, 18),
                            mat = if (varwise) cmat_ext else cmat,
                            dir = if (varwise) directions_ext else directions,
                            rhs = if (varwise) rhss_ext else rhss,
                            control = list(presolve = TRUE))
        #return(lpref)
        nvec = lpref$solution
    }
    res = list(n_original = mort_year[[1]],
                r_agest_original = mort_year[[3]],
                r_crude_original = mort_year[[2]],
                r_agest = sum(nvec * pop_sgg_coef_agest),
                r_crude = sum(nvec * pop_wc),
                sol = nvec,
                #cmat = cmat,
                fitobj = lpref)

    return(res)
}

mor_cll = mor_cl %>% filter(!is.na(value))

total_2000 = 
unique(mor_cl$sgg_cd) %>%
    split(.,.) %>%
    lapply(function(x) { 
        tryCatch({
            find_mort_n(mor_cll, 'Total', x, spop, sgg_pop_cl, 2000, 'Stomach', 'glpk', varwise = T, setmax =T)$sol},
            error = function(e) return(rep(0,18)))
    })
Reduce(`+`, total_2000) 
Reduce(rbind, total_2000)

total_2010 = 
unique(mor_cl$sgg_cd) %>%
    split(.,.) %>%
    lapply(function(x) { 
        tryCatch({
            find_mort_n(mor_cll, 'Total', x, spop, sgg_pop_cl, 2010, 'Stomach', 'glpk', varwise = T, setmax =T)$sol},
            error = function(e) return(rep(0,18)))
    })
Reduce(`+`, total_2010) 



jj = find_mort_n(mor_cl, 'Male', 21030, spop, sgg_pop_cl, 1999, 'Stomach', 'glpk', varwise = T, setmax =T)
kk = find_mort_n(mor_cl, 'Total', 35020, spop, sgg_pop_cl, 2011, 'Stomach', 'glpk', varwise = T, setmax =T)
kk = find_mort_n(mor_cl, 'Total', 23050, spop, sgg_pop_cl, 2003, 'Lung/Bronchus', 'glpk', varwise = T, setmax = T)
qw = find_mort_n(mor_cl, 'Male', 38030, spop, sgg_pop_cl, 2013, 'Lung/Bronchus', 'glpk', varwise = T, setmax = T)
qw = find_mort_n(mor_cl, 'Female', 38020, spop, sgg_pop_cl, 2005, 'Lung/Bronchus', 'glpk', varwise = T, setmax = T)
qw = find_mort_n(mor_cl, 'Male', 39010, spop, sgg_pop_cl, 2009, 'Stomach', 'glpk', varwise = T, setmax = F)

fmn_11210 = 
    1999:2013 %>% 
        split(.,.) %>%
        lapply(function(x) find_mort_n(mor_cll, 'Total', 11210, spop, sgg_pop_cl, x, 'Stomach', 'glpk', varwise = T)$sol)
Reduce(`+`, fmn_11210)

fmn_34040 = 
    1999:2013 %>% 
        split(.,.) %>%
        lapply(function(x) find_mort_n(mor_cl, 'Total', 34040, spop, sgg_pop_cl, x, 'Stomach', 'glpk', varwise = T)$sol)
Reduce(`+`, fmn_34040)

fmn_38030 = 
    1999:2013 %>% 
        split(.,.) %>%
        lapply(function(x) find_mort_n(mor_cl, 'Total', 38030, spop, sgg_pop_cl, x, 'Lung/Bronchus', 'glpk', varwise = T)$sol)
Reduce(`+`, fmn_38030)

fmn_38040 = 
    1999:2013 %>% 
        split(.,.) %>%
        lapply(function(x) find_mort_n(mor_cl, 'Total', 38040, spop, sgg_pop_cl, x, 'Lung/Bronchus', 'glpk', varwise = T)$sol)
Reduce(`+`, fmn_38040)


spop_t = spop %>% filter(성별 == "계")
p13 = spop_t$X2013
p13_r = p13 / sum(p13)
p13_coef = p13_r * (1000 / p13)


