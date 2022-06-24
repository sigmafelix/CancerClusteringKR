### Covariate processing (generalized)
### 01/22/2022
PCTempDir <- Sys.getenv("TEMP")

#detect and delete folders with pattern "Rtmp"
folders <- dir(PCTempDir, pattern = "Rtmp", full.names = TRUE)
unlink(folders, recursive = TRUE, force = TRUE, expand = TRUE)


# Clean covariates: soon to be obsolete
clean_covar = function(db_dir = dbdir,
                       geo_path = geopath,
                       target_year = 2010,
                       gp_base_dir = rdatafiles
                        ) {

    load(gp_base_dir[grep(target_year, gp_base_dir)])

    # Cancer mortality (number only)
    cmorts = list.files(pattern = '^cancerMor',
                        path = db_dir,
                        full.names = TRUE)
    cmorts = lapply(cmorts, function(x) read.csv(x, fileEncoding = 'CP949'))
    cmorts = lapply(cmorts,
                    function(x) {
                    colnames(x) = c('cause0', 'cause', 'sgg_cd', 'sgg_nm', 'sex0', 'sex', 'type0', 'type', 'unit', paste('Y', 1998:2019, sep=''), 'x')
                    x = x %>%
                    mutate(cancer_type_e = plyr::mapvalues(cause, unique(cause), c('Stomach', 'Colorectal', 'Liver', 'Lung', 'Breast', 'Cervical', 'Prostate'))) 
                    return(x)
                    }
    )
    cmorts_df = do.call(rbind, cmorts) %>%
        dplyr::select(sgg_cd, cancer_type_e, sex0, type0, str_c('Y', target_year)) %>%
        filter(type0 == 'T1') %>%
        mutate(sex = plyr::mapvalues(sex0, 0:2, c('total', 'male', 'female'))) %>%
        dplyr::select(-sex0, -type0) %>%
        pivot_wider(names_from = c(cancer_type_e, sex), values_from = str_c('Y', target_year))

    # Total Population by sex
    midpop = read.csv(str_c(db_dir, '/MidPopulation_2000_2020.csv'), fileEncoding = 'CP949')
    midpop = midpop %>%
        dplyr::select(1, 3, 6, 10:14)
    colnames(midpop) = c('sggcd', 'sex0', 'agegroup', 'Y2000', 'Y2005', 'Y2010', 'Y2015', 'Y2020')
    midpop = midpop %>%
        dplyr::select(1:3, 7) %>%
        filter(agegroup == '계') %>%
        pivot_wider(names_from = sex0, values_from = sym(str_c('Y', target_year))) %>%
        dplyr::select(-2) %>%
        rename(pop_total = `0`,
            pop_male = `1`,
            pop_female = `2`)


    # NDVI
    st_crs(carreg) = 5179
    ndvi_stars = st_as_stars(ndvi_mean)
    #plot(carreg)
    if (target_year > 2010) {
        sgg = carreg#
    } else {
        if (target_year == 2010) {
            carreg = carreg %>%
                group_by(SIGUNGU_CD) %>%
                summarize_at(.vars = vars(X01:Car_Mean),
                             .funs = list(~unique(.))) %>%
                ungroup
        }
        if (sum(grepl('^(sigungu|sgg|SGG)', colnames(carreg))) != 0) {
            colnames(carreg)[grep('^(sigungu|sgg|SGG)', colnames(carreg))] = 'SGGCD'
        }
        sgg = st_transform(carreg, 5174)
        st_crs(ndvi_stars) = 5174
    }
    sgg_ndvi = starsExtra::extract2(ndvi_stars, sgg, function(x) mean(x, na.rm = T))

    # Air pollution exposure
    airpol = read.csv(str_c(db_dir, '/Covariates/Airpollution/', target_year, '_new.csv'), fileEncoding = 'CP949')
    airpol = airpol %>%
        transmute(sggcd = SGG_cd,
                ap_PM10_pred = PM10_pred,
                ap_NO2_pred = NO2_pred)

    # Emission
    emission = st_as_stars(get(str_c('em.', target_year)), crs = 5174)
    st_crs(emission) = 5174
    sgg_5174 = st_transform(sgg, 5174)
    sgg_emission = starsExtra::extract2(emission, sgg_5174, function(x) sum(x, na.rm = T))
    sgg_emission = sgg_emission %>%
        as.data.frame
    colnames(sgg_emission) = str_c('sum_', attr(emission, 'dimensions')$band$values)

    # Education attainment
    educ = read.csv(str_c(db_dir, '/Covariates/Education_', target_year, '.csv'), fileEncoding = 'CP949') %>%
        .[,c(1, 4:7)]
    colnames(educ) = c('sggcd', 'total_6yo', 'bachelor', 'masters', 'doctorate')
    educ = educ %>%
        transmute(sggcd = as.integer(str_split(sggcd, ' ', simplify = TRUE)[,1]),
                p_eduhigh = 100 * (bachelor + masters + doctorate) / total_6yo)


    # KCDC Community Health Data: separate tables: 2010, 2015
    ## KCDC code to others
    code_conv = read_xlsx(str_c(db_dir, '/행정자치부_심평원_통계청_코드_연계표_211119.xlsx'), sheet = 4)
    code_conv = code_conv %>% 
        dplyr::select(SIGUNGU_PSEUDO, ends_with(target_year)) %>%
        rename(SGIS = sym('SGIS_', target_year),
                SIGUNGU_KCDC = SIGUNGU_PSEUDO) %>%
        mutate(SGIS = str_sub(SGIS, 1, 5)) %>%
        filter(!is.na(SGIS))

    ## Main
    kcdc_csvs = 
    list.files(pattern = '시군구별_*.*.csv$', 
            path = str_c(db_dir, "/Covariates/"),
            full.names = TRUE)
    kcdc_list = lapply(kcdc_csvs, 
        function(x) read_csv(x, skip = 1, locale = locale(encoding = 'CP949')) %>%
            mutate_at(.vars = vars(-1:-3), .funs = list(~as.numeric(.))))

    kcdc_cns = 
    c('sido', 'sigungu', 'sub_sigungu',
                str_c(rep(c('N_', 'CR_', 'CR_SE_', 'SR_', 'SR_SE_'), 13), rep(2008:2020, each = 5)))
    kcdc_cn4 = c('sido', 'sigungu', 'sub_sigungu',
                str_c(rep(c('N_', 'CR_', 'CR_SE_', 'SR_', 'SR_SE_'), 12), rep(2009:2020, each = 5)))
    colnames(kcdc_list[[1]]) = kcdc_cns
    colnames(kcdc_list[[2]]) = kcdc_cns[-59:-63]
    colnames(kcdc_list[[3]]) = kcdc_cns
    colnames(kcdc_list[[4]]) = kcdc_cn4[-grep('2018$', kcdc_cn4)]
    colnames(kcdc_list[[5]]) = kcdc_cns


    kcdc_covars = kcdc_list %>%
    lapply(function(x) x %>% 
        mutate_at(.vars = vars(1:3), .funs = list(~gsub('[[:blank:]]|[가-힣]', '', .))) %>%
        transmute(SIGUNGU_KCDC = ifelse(sub_sigungu == '', sigungu, sub_sigungu),
                    SR = sym(str_c('SR_', target_year)))) %>%
        mapply(function(x, y) {colnames(x)[2] = y; return(x)},
                ., c('r_walking', 'r_obesity', 'r_alcoholmonth', 'r_physmid', 'r_smoking') %>% split(., 1:5), SIMPLIFY = FALSE) %>%
        plyr::join_all(.)
    kcdc_covars = kcdc_covars %>%
        left_join(code_conv, .) %>%
        mutate(SGIS = as.integer(SGIS))


    # KCDC Community Health Data
    kcdc = readxl::read_xlsx(str_c(db_dir, '/Covariates/CommunityHealth_Covariates_KCDC.xlsx'), 
                            sheet = str_c('Data_', target_year))
    kcdc_sub = kcdc %>% 
        dplyr::select(2:5, #contains('비만율'), contains('음주율'), contains('신체활동'), contains(' 흡연율'), contains('현재흡연율'), 
        contains('암검진'), #contains('걷기'),
        contains('건강검진수검율'), contains('암검진율'), contains('실업률'), contains('지가변동률'), contains('산림면적비율'),
        contains('암검진_대상인원'), contains('암검진_수검인원')) %>%
        transmute(sggcd = 코드,
            sido = 시도,
            sigungu = 시군구,
            sdsgg = 지역,
            #p_obese_std = 비만율_자가보고_표준화율,
            #p_obese_cru = 비만율_자가보고_조율,
            #p_alcl_std = 평생음주율_표준화율,
            #p_alcl_cru = 평생음주율_조율,
            #p_alcy_std = 연간음주율_표준화율,
            #p_alcy_cru = 연간음주율_조율,
            #p_smk_std = 현재흡연율_표준화율,
            #p_smk_cru = 현재흡연율_조율,
            #p_active_std = `중등도신체활동실천율_표준화율`,
            #p_active_cru = `중등도신체활동실천율_조율`,
            #p_walking_std = `걷기실천율_표준화율`,
            #p_walking_cru = `걷기실천율_조율`,
            p_candiag_std = 암검진율_표준화율,
            p_candiag_cru = 암검진율_조율,
            p_candiag_sto = 위암검진_수검인원/위암검진_대상인원 * 100,
            p_candiag_col = 대장암검진_수검인원/대장암검진_대상인원 * 100,
            p_candiag_liv = 간암검진_수검인원/간암검진_대상인원 * 100,
            p_candiag_bre = 유방암검진_수검인원/유방암검진_대상인원 * 100,
            p_candiag_cer = 자궁경부암검진_수검인원/자궁경부암검진_대상인원 * 100,
            r_unemp = 실업률,
            r_landprice = 지가변동률,
            r_forest = 산림면적비율) %>% # 산림면적비율 = 산림면적/도시면적 * 100 (may exceed 100)
        # for data-specific problem: different data granularity of availability
        mutate(sggcd_pseudo = ifelse(grepl('^29$', sggcd), 29010, ifelse(nchar(sggcd) == 5, str_sub(sggcd, 1, 4), sggcd)),
            sggcd_unemp = ifelse(grepl('^[1-2][0-9]|^39', sggcd), str_sub(sggcd, 1, 2), ifelse(grepl('[3][0-8]..[0-9]', sggcd), str_sub(sggcd, 1, 4), sggcd))) %>%
        mutate(sigungu_b = if_else(grepl('[[:blank:]]', sigungu), 
                                str_split(sigungu, ' ', simplify = TRUE)[,1], 
                                str_c(sido, sigungu))
            )
    kcdc_sub_fctrs = sapply(kcdc_sub, function(x) !is.character(x) & !is.factor(x)) %>%
                    as.logical %>%
                    which %>%
                    .[-1]
    kcdc_sub = kcdc_sub %>%
        group_by(sggcd_pseudo) %>%
        mutate_at(vars(-group_cols(), -sggcd, c(kcdc_sub_fctrs), -r_unemp), list(~ifelse(any(!is.na(.)), .[which(!is.na(.))], .))) %>%
        ungroup %>%
        group_by(sggcd_unemp) %>%
        mutate_at(vars(-group_cols(), -sggcd, c(kcdc_sub_fctrs)), list(~ifelse(any(!is.na(.)), .[which(!is.na(.))], .))) %>%
        #ungroup %>%
        #group_by(sggcd) %>%
        # take the first row (sigungu-representative) from each sigungu code
        #summarize_at(vars(-group_cols()), list(~.[which(!is.na(.))])) %>%
        ungroup %>%
        mutate(sggcd = ifelse(sggcd == 29, 29010, sggcd))

    sgg_covars_cleaned = sgg %>%
        mutate(NDVI_mean = sgg_ndvi) %>%
        bind_cols(sgg_emission) %>%
        left_join(cmorts_df, by = c('SGGCD' = 'sgg_cd')) %>%
        left_join(midpop, by = c('SGGCD' = 'sggcd')) %>%
        left_join(kcdc_sub %>% filter(!duplicated(sggcd)), by = c('SGGCD' = 'sggcd')) %>%
        left_join(kcdc_covars, by = c('SGGCD' = 'SGIS')) %>%
        left_join(airpol, by = c('SGGCD' = 'sggcd')) %>%
        left_join(educ, by = c('SGGCD' = 'sggcd')) 
   

    # Conversion
    conv_table = read.csv(paste(db_dir, 'SGG_before_2022_Conversion_Table_220122_CConly.csv', sep = ''), fileEncoding = 'EUC-KR')

    conv_table_e = conv_table %>%
        filter(to_year >= 1999 & from_year <= 2013) %>%
        dplyr::select(fromcode, tocode)

    sgg_covars_cleaned = sgg_covars_cleaned %>%
        mutate(sgg_cd_c = plyr::mapvalues(SGGCD, conv_table_e$from_code, conv_table_e$to_code)) %>%
        group_by(sgg_cd_c) %>%


    return(sgg_covars_cleaned)
    }



## To aggregate by 5-year periods
get_basecovar = function(db_dir = dbdir,
                       geo_path = geopath,
                       target_year = 2010,
                       target_year_kcdc = 2010,
                       gp_base_dir = rdatafiles
                        ) {
    ##
    if (target_year == 2000) {
        target_year_r = 2001
    } else {
        target_year_r = target_year
    }

    load(gp_base_dir[grep(target_year_r, gp_base_dir)])


    ## remove after the GP data fix
    assign(str_c('em.', target_year_r),
            raster::mosaic(get(str_c('seed.raster.a.', target_year_r)),
                           get(str_c('seed.raster.l.', target_year_r)),
                           get(str_c('seed.raster.p.', target_year_r)),
                           fun = sum))

    # Cancer mortality (number only)
    cmorts = list.files(pattern = '^cancerMor',
                        path = db_dir,
                        full.names = TRUE)
    cmorts = lapply(cmorts, function(x) read.csv(x, fileEncoding = 'CP949'))
    cmorts = lapply(cmorts,
                    function(x) {
                    colnames(x) = c('cause0', 'cause', 'sgg_cd', 'sgg_nm', 'sex0', 'sex', 'type0', 'type', 'unit', paste('Y', 1998:2019, sep=''), 'x')
                    x = x %>%
                        mutate(cancer_type_e = plyr::mapvalues(cause, unique(cause), c('Stomach', 'Colorectal', 'Liver', 'Lung', 'Breast', 'Cervical', 'Prostate'))) 
                    return(x)
                    }
    )
    cmorts_df = do.call(rbind, cmorts) %>%
        dplyr::select(sgg_cd, cancer_type_e, sex0, type0, str_c('Y', target_year)) %>%
        filter(type0 == 'T1') %>%
        mutate(sex = plyr::mapvalues(sex0, 0:2, c('total', 'male', 'female'))) %>%
        dplyr::select(-sex0, -type0) %>%
        pivot_wider(names_from = c(cancer_type_e, sex), 
                    names_prefix = 'n_', 
                    values_from = str_c('Y', target_year))

    # Total Population by sex
    midpop = read.csv(str_c(db_dir, '/MidPopulation_2000_2020.csv'), fileEncoding = 'CP949')
    midpop = midpop %>%
        dplyr::select(1, 3, 6, 5, 10:14)
    colnames(midpop) = c('sggcd', 'sex0', 'agegroup', 'agegroup_f', 'Y2000', 'Y2005', 'Y2010', 'Y2015', 'Y2020')
    midpop_elderly = midpop %>%
        filter(agegroup_f %in% c(0, 280, 310, 330, 340)) %>%
        group_by(sggcd, sex0) %>%
        summarize(n_calc_pop = (!!sym(str_c("Y", target_year)))[1],
                  n_calc_pop_65p = sum((!!sym(str_c("Y", target_year)))[-1])) %>%
        ungroup %>%
        mutate(sex0 = plyr::mapvalues(sex0, 0:2, c('total', 'male', 'female'))) %>%
        pivot_wider(names_from = sex0,
                    values_from = c(n_calc_pop, n_calc_pop_65p))

    midpop_sex = midpop %>%
        dplyr::select(1:3, ends_with(as.character(target_year))) %>%
        filter(agegroup == '계') %>%
        pivot_wider(names_from = sex0, values_from = !!sym(str_c('Y', target_year))) %>%
        dplyr::select(-2) %>%
        rename(n_pop_total = `0`,
            n_pop_male = `1`,
            n_pop_female = `2`)

    # NDVI
    st_crs(carreg) = 5179
    ndvi_stars = st_as_stars(ndvi_mean)
    #plot(carreg)
    if (target_year > 2010) {
        sgg = carreg#
    } else {
        if (target_year == 2010) {
            carreg = carreg %>%
                group_by(SIGUNGU_CD) %>%
                summarize_at(.vars = vars(X01:Car_Mean),
                             .funs = list(~unique(.))) %>%
                ungroup %>%
                rename(SGGCD = SIGUNGU_CD)
        }
        if (sum(grepl('^(sigungu|sgg|SGG)', colnames(carreg))) != 0) {
            colnames(carreg)[grep('^(sigungu|sgg|SGG)', colnames(carreg))] = 'SGGCD'
        }
        sgg = st_transform(carreg, 5174)
        st_crs(ndvi_stars) = 5174
    }
    sgg_ndvi = starsExtra::extract2(ndvi_stars, sgg, function(x) mean(x, na.rm = T))

    # Air pollution exposure
    airpol = read.csv(str_c(db_dir, '/Covariates/Airpollution/', target_year_r, '_new.csv'), fileEncoding = 'CP949')
    airpol = airpol %>%
        transmute(sggcd = SGG_cd,
                ap_PM10_pred = PM10_pred,
                ap_NO2_pred = NO2_pred)

    # Emission
    emission = st_as_stars(get(str_c('em.', target_year_r)), crs = 5174)
    st_crs(emission) = 5174
    sgg_5174 = st_transform(sgg, 5174)
    sgg_emission = starsExtra::extract2(emission, sgg_5174, function(x) sum(x, na.rm = T))
    sgg_emission = sgg_emission %>%
        as.data.frame
    sgg_emission_dims = attr(emission, 'dimensions')$band$values
    if (length(sgg_emission_dims) == 7) {
        colnames(sgg_emission) = str_c('ap_sum_em_', c('co', 'nox', 'sox', 'tsp', 'pm10', 'voc', 'nh3'))
    } else {
        colnames(sgg_emission) = str_c('ap_sum_em_', attr(emission, 'dimensions')$band$values)
    }

    # Education attainment
    educ = read.csv(str_c(db_dir, '/Covariates/Education_', target_year, '.csv'), fileEncoding = 'CP949') %>%
        .[,c(1:2, 4:7)]
    colnames(educ) = c('sggcd', 'sex_cd', 'total_6yo', 'bachelor', 'masters', 'doctorate')
    educ = educ %>%
        transmute(sggcd = as.integer(str_split(sggcd, ' ', simplify = TRUE)[,1]),
                sex_cd = plyr::mapvalues(sex_cd, unique(sex_cd), c('total', 'male', 'female')),
                bachelor = bachelor,
                masters = masters, 
                doctorate = doctorate,
                total_6yo = total_6yo) %>%
        pivot_longer(cols = 3:6) %>%
        pivot_wider(names_from = c(name, sex_cd), names_prefix = 'n_', names_sep = '_')


    sgg_covars_cleaned = sgg %>%
        mutate(SGGCD = as.integer(as.character(SGGCD))) %>%
        mutate(NDVI_mean = sgg_ndvi) %>%
        bind_cols(sgg_emission) %>%
        #left_join(cmorts_df, by = c('SGGCD' = 'sgg_cd')) %>%
        left_join(airpol, by = c('SGGCD' = 'sggcd')) %>%
        left_join(educ, by = c('SGGCD' = 'sggcd')) %>%
        left_join(midpop_sex, by = c('SGGCD' = 'sggcd')) %>%
        left_join(midpop_elderly, by = c('SGGCD' = 'sggcd'))



    if (target_year_kcdc >= 2008) {
        # KCDC Community Health Data: separate tables: 2010, 2015
        ## KCDC code to others
        code_conv = read_xlsx(str_c(db_dir, '/행정자치부_심평원_통계청_코드_연계표_211119.xlsx'), sheet = 4)
        code_conv = code_conv %>% 
            dplyr::select(SIGUNGU_PSEUDO, ends_with(as.character(target_year_kcdc))) %>%
            rename(SGIS = !!sym(str_c('SGIS_', target_year_kcdc)),
                    SIGUNGU_KCDC = SIGUNGU_PSEUDO) %>%
            mutate(SGIS = str_sub(SGIS, 1, 5)) %>%
            filter(!is.na(SGIS))
        ## Main
        kcdc_csvs = 
        list.files(pattern = '시군구별_*.*.csv$', 
                path = str_c(db_dir, "/Covariates/"),
                full.names = TRUE)
        kcdc_list = lapply(kcdc_csvs, 
            function(x) read_csv(x, skip = 1, locale = locale(encoding = 'CP949')) %>%
                mutate_at(.vars = vars(-1:-3), .funs = list(~as.numeric(.))))

        kcdc_cns = 
            c('sido', 'sigungu', 'sub_sigungu',
                    str_c(rep(c('N_', 'CR_', 'CR_SE_', 'SR_', 'SR_SE_'), 13), rep(2008:2020, each = 5)))
        kcdc_cn4 = c('sido', 'sigungu', 'sub_sigungu',
                    str_c(rep(c('N_', 'CR_', 'CR_SE_', 'SR_', 'SR_SE_'), 12), rep(2009:2020, each = 5)))
        colnames(kcdc_list[[1]]) = kcdc_cns
        colnames(kcdc_list[[2]]) = kcdc_cns[-59:-63]
        colnames(kcdc_list[[3]]) = kcdc_cns
        colnames(kcdc_list[[4]]) = kcdc_cn4[-grep('2018$', kcdc_cn4)]
        colnames(kcdc_list[[5]]) = kcdc_cns


        kcdc_covars = kcdc_list %>%
            lapply(function(x) x %>% 
                mutate_at(.vars = vars(1:3), .funs = list(~gsub('[[:blank:]]|[가-힣]', '', .))) %>%
                mutate(SR_2008 = if (sum(grepl('SR_2008', colnames(x))) ==0)  NA else SR_2008) %>%
                transmute(SIGUNGU_KCDC = ifelse(sub_sigungu == '', sigungu, sub_sigungu),
                            SR = !!sym(str_c('SR_', target_year_kcdc)))) %>%
            mapply(function(x, y) {colnames(x)[2] = y; return(x)},
                    ., c('r_walking', 'r_obesity', 'r_alcoholmonth', 'r_physmid', 'r_smoking') %>% split(., 1:5), SIMPLIFY = FALSE) %>%
            plyr::join_all(.)
        kcdc_covars = kcdc_covars %>%
            left_join(code_conv, .) %>%
            mutate(SGIS = as.integer(SGIS))


        # KCDC Community Health Data
        kcdc = readxl::read_xlsx(str_c(db_dir, '/Covariates/CommunityHealth_Covariates_KCDC.xlsx'), 
                                sheet = str_c('Data_', target_year_kcdc))
        kcdc_sub = kcdc %>% 
            dplyr::select(2:5, #contains('비만율'), contains('음주율'), contains('신체활동'), contains(' 흡연율'), contains('현재흡연율'), 
            contains('암검진'), #contains('걷기'),
            contains('건강검진수검율'), contains('암검진율'), contains('실업률'), contains('지가변동률'), contains('산림면적비율'),
            contains('암검진_대상인원'), contains('암검진_수검인원')) %>%
            transmute(sggcd = 코드,
                sido = 시도,
                sigungu = 시군구,
                sdsgg = 지역,
                n_candiag_std = 암검진율_표준화율,
                n_candiag_cru = 암검진율_조율,
                n_candiag_sto_denom = ifelse(target_year_kcdc >= 2010, 위암검진_수검인원, NA),
                n_candiag_sto_nom = ifelse(target_year_kcdc >= 2010, 위암검진_대상인원, NA),
                n_candiag_col_denom = ifelse(target_year_kcdc >= 2010, 대장암검진_수검인원, NA),
                n_candiag_col_nom = ifelse(target_year_kcdc >= 2010, 대장암검진_대상인원, NA),
                n_candiag_liv_denom = ifelse(target_year_kcdc >= 2010, 간암검진_수검인원, NA),
                n_candiag_liv_nom = ifelse(target_year_kcdc >= 2010, 간암검진_대상인원, NA),
                n_candiag_bre_denom = ifelse(target_year_kcdc >= 2010, 유방암검진_수검인원, NA),
                n_candiag_bre_nom = ifelse(target_year_kcdc >= 2010, 유방암검진_대상인원, NA),
                n_candiag_cer_denom = ifelse(target_year_kcdc >= 2010, 자궁경부암검진_수검인원, NA),
                n_candiag_cer_nom = ifelse(target_year_kcdc >= 2010, 자궁경부암검진_대상인원, NA),
                r_unemp = 실업률,
                r_landprice = 지가변동률) %>% # 산림면적비율 = 산림면적/도시면적 * 100 (may exceed 100)
            # for data-specific problem: different data granularity of availability
            mutate(sggcd_pseudo = ifelse(grepl('^29$', sggcd), 29010, ifelse(nchar(sggcd) == 5, str_sub(sggcd, 1, 4), sggcd)),
                sggcd_unemp = ifelse(grepl('^[1-2][0-9]|^39', sggcd), str_sub(sggcd, 1, 2), ifelse(grepl('[3][0-8]..[0-9]', sggcd), str_sub(sggcd, 1, 4), sggcd))) %>%
            mutate(sigungu_b = if_else(grepl('[[:blank:]]', sigungu), 
                                    str_split(sigungu, ' ', simplify = TRUE)[,1], 
                                    str_c(sido, sigungu))
                )
        kcdc_sub_fctrs = sapply(kcdc_sub, function(x) !is.character(x) & !is.factor(x)) %>%
                        as.logical %>%
                        which %>%
                        .[-1]
        kcdc_sub = kcdc_sub %>%
            group_by(sggcd_pseudo) %>%
            mutate_at(vars(-group_cols(), -sggcd, all_of(kcdc_sub_fctrs), -r_unemp), list(~ifelse(any(!is.na(.)), .[which(!is.na(.))], .))) %>%
            ungroup %>%
            group_by(sggcd_unemp) %>%
            mutate_at(vars(-group_cols(), -sggcd, all_of(kcdc_sub_fctrs)), list(~ifelse(any(!is.na(.)), .[which(!is.na(.))], .))) %>%
            #ungroup %>%
            #group_by(sggcd) %>%
            # take the first row (sigungu-representative) from each sigungu code
            #summarize_at(vars(-group_cols()), list(~.[which(!is.na(.))])) %>%
            ungroup %>%
            mutate(sggcd = ifelse(sggcd == 29, 29010, sggcd))

    sgg_covars_cleaned = sgg_covars_cleaned %>%
        left_join(kcdc_sub %>% filter(!duplicated(sggcd)), by = c('SGGCD' = 'sggcd')) %>%
        left_join(kcdc_covars, by = c('SGGCD' = 'SGIS'))

    }
   
    return(sgg_covars_cleaned)
    }


# to make unified district data
clean_consolidated = function(db_dir = dbdir, cleaned_df, target_year = 2010) {

    # Conversion
    conv_table = read.csv(paste(db_dir, 'SGG_before_2022_Conversion_Table_220122_CConly.csv', sep = ''), fileEncoding = 'EUC-KR')

    conv_table_e = conv_table %>%
        filter(to_year >= 1999 & from_year <= 2013) %>%
        dplyr::select(fromcode, tocode)

    cleaned_df = cleaned_df %>%
        mutate(sgg_cd_c = plyr::mapvalues(SGGCD, conv_table_e$fromcode, conv_table_e$tocode))
    # population weighted rates
    cleaned_df_rates = cleaned_df %>%
        st_drop_geometry %>%
        group_by(sgg_cd_c) %>%
        summarize_at(.vars = vars(starts_with('r_'), ends_with('_pred'), NDVI_mean),
                     .funs = list(~weighted.mean(., w = n_pop_total/sum(n_pop_total, na.rm = T), na.rm = T))) %>%
        ungroup
    # counts and compound counts
    cleaned_df_counts = cleaned_df %>%
        st_drop_geometry %>%
        group_by(sgg_cd_c) %>%
        summarize_at(.vars = vars(starts_with('n_')),
                     .funs = list(~sum(., na.rm = TRUE))) %>%
        ungroup %>%
        mutate(p_65p_total = 100 * n_calc_pop_65p_total / n_calc_pop_total,
               p_65p_male = 100 * n_calc_pop_65p_male / n_calc_pop_male,
               p_65p_female = 100 * n_calc_pop_65p_female / n_calc_pop_female,
               p_hbac_total = 100 * (n_bachelor_total + n_masters_total + n_doctorate_total) / n_total_6yo_total,
               p_hbac_male = 100 * (n_bachelor_male + n_masters_male + n_doctorate_male) / n_total_6yo_male,
               p_hbac_female = 100 * (n_bachelor_female + n_masters_female + n_doctorate_female) / n_total_6yo_female,
               p_candiag_sto = if (target_year >= 2008) 100 * n_candiag_sto_denom / n_candiag_sto_nom else NA) %>%
        dplyr::select(sgg_cd_c, starts_with('n_Stomach'), starts_with('n_Lung'), starts_with('p_hbac'), p_candiag_sto, starts_with('n_pop'), starts_with('p_65p'))
    # air pollution
    cleaned_df_apsum = cleaned_df %>%
        st_drop_geometry %>%
        group_by(sgg_cd_c) %>%
        summarize_at(.vars = vars(starts_with('ap_sum')),
                     .funs = list(~sum(., na.rm = TRUE))) %>%
        ungroup

    # Consolidation
    cleaned_df_consol =
    cleaned_df %>%
        group_by(sgg_cd_c) %>%
        summarize(Car_Mean = sum(Car_Mean, na.rm = T)) %>%
        ungroup %>%
        left_join(cleaned_df_counts) %>%
        left_join(cleaned_df_rates) %>%
        left_join(cleaned_df_apsum)

    if (target_year == 2010) {   
        source("Cancer_Data_Clearing_081721.R")
        pr_csvs_con_prempop = pr_csvs_con_prempop %>%
            mutate(SGIS = as.integer(SGIS))

        cleaned_df_consol = cleaned_df_consol %>%
            left_join(pr_csvs_con_prempop, by = c('sgg_cd_c' = 'SGIS'))
    }
    return(cleaned_df_consol)
}

#covar_origin_10 = get_basecovar(target_year = 2010)
#covar_origin_o = covar_origin_10
#colnames(covar_origin_10)[grep('^ap_sum_em_', colnames(covar_origin_10))] =
#    str_c('ap_sum_em_', c('co', 'nox', 'sox', 'tsp', 'pm10', 'voc', 'nh3'))
#covar_origin_10_consol = clean_consolidated(cleaned_df = covar_origin_10)

# Run estimate model (aka covariate adjustment)
regress_counts = function(data, 
                          population,
                          yvar,
                          sex_b,
                          string_search = str_c(str_c('(^p_*.*_', sex_b, '$'), '^(r_|ap_)', '^NDVI_)', sep = '|'),
                          add_var = NULL) {
    if (!is.null(add_var)) {
        string_search = str_replace(string_search, '\\^NDVI_\\)', str_c('^NDVI_|', add_var, ')'))
    }
    #print(string_search)
    cns = colnames(data)[grep(string_search, colnames(data), perl = TRUE)]
    #print(cns)
    form_pois = as.formula(str_c(yvar, '~offset(log(', population, '))+', str_c(cns, collapse = '+')))
    reg_pois = glm(formula = form_pois, data= data, family = poisson(link = 'log'))
    return(reg_pois)
}

regress_rates = function(data, 
                          population = NULL,
                          yvar,
                          sex_b,
                          string_search = str_c(str_c('(^p_*.*_', sex_b, '$'), '^(r_|ap_)', '^NDVI_)', sep = '|'),
                          add_var = NULL) {
    if (!is.null(add_var)) {
        string_search = str_replace(string_search, '\\^NDVI_\\)', str_c('^NDVI_|', add_var, ')'))
    }
    #print(string_search)
    cns = colnames(data)[grep(string_search, colnames(data), perl = TRUE)]
    #print(cns)
    form_norm = as.formula(str_c(yvar, '~', str_c(cns, collapse = '+')))
    reg_norm = lm(formula = form_norm, data= data)
    return(reg_norm)
}

generate_satscan_prm = function(data,
                    title.analysis,
                    dir.base,
                    dir.target,
                    name.idcol,
                    filename.input, filename.output,
                    col.var, col.case, #col.weight = '',
                    coord.x = "X", coord.y = "Y",
                    prm.path,
                    #model = "normal",
                    weighted = TRUE,
                    adjust = FALSE,
                    string_search = str_c(str_c('(^p_*.*_', sex_b, '$'), '^(r_|ap_)', '^NDVI_)', sep = '|'),
                    add_var = NULL,
                    vset = NULL
                    ) {
    if (!file.exists(prm.path)) {
        file.create(prm.path)
    }
    
    dcols = colnames(data)
    iscount = grepl('_n_', title.analysis)
    for_residuals = !iscount
    #fullpath.input0 = str_c(dir.base, filename.input)
    fullpath.input = str_c(dir.base, filename.input)
    fullpath.output = str_c(dir.target, filename.output)
    indx.idcol = grep(name.idcol, dcols)
    indx.case = ifelse(iscount, grep(col.var, dcols), 
        ifelse(weighted, grep(str_replace(col.var, "ragest_", "n_"), dcols),
                grep('case_normal', dcols)))
    # Soon to be deprecated: indx.weight and col.weight
    indx.weight = ifelse(iscount | weighted, '', 
                        str_c(",", grep(str_replace(col.var, "ragest_", "n_"), dcols)))
    indx.var = ifelse(iscount, '', str_c(",,", grep(col.var, dcols)))
    indx.xcoord = grep(str_c("^", coord.x, "$"), dcols)
    indx.ycoord = grep(str_c("^", coord.y, "$"), dcols)
    indx.year = grep("dummyyear", dcols)
    filename.temp = str_replace(filename.input, "\\.csv", str_c("_", col.var, ".csv"))
    fullpath.temp = str_c(dir.target, filename.temp)
    file.copy(fullpath.input, fullpath.temp)

    # weighted normal: 0 or 1
    #indx.casenormal = grep('case_normal', dcols)
    #file.copy(fullpath.input, str_c("/home/felix/Documents/", filename.input), overwrite = TRUE)
    #file.copy(fullpath.input, )
    if (iscount) {
        indx.pop = grep(str_c("n_p", str_extract(title.analysis, "_(female|male|total)")), dcols)
        line.popdef = str_c("PopulationFile=", fullpath.temp)
        line.pop = "PopulationFile-SourceType=0"
        line.popmap = str_c("PopulationFile-SourceFieldMap=", indx.idcol, ",", indx.year, ",", indx.pop)
    } else {
        line.popdef = "PopulationFile="
        line.pop = "PopulationFile-SourceType="
        line.popmap = "PopulationFile-SourceFieldMap="
    }

    # If adjust, we replace the population with the "covariate-controlled" estimates
    if (adjust) {
        sex_t = str_extract(title.analysis, "(male|female|total)")
        string_search_n =
            switch(vset,
                    set1 = str_c('^(p_65p|p_hbac)*.*_', sex_t, '$'),
                    set2 = str_c(str_c('^p_*.*_', sex_t, '$'), '^ap_', '^NDVI_', sep = '|'),
                    set3 = str_c(str_c('^p_*.*_', sex_t, '$'), '^r_(?!physmid)', '^ap_', '^NDVI_', sep = '|'),
                    set4 = str_c(str_c('^p_*.*_', sex_t, '$'), '^r_', '^n_pw', '^ap_', '^NDVI_', sep = '|'))

        if (iscount) {
            lms = regress_counts(data = data,
                        yvar = col.var,
                        population = dcols[indx.pop],
                        sex_b = sex_t,
                        string_search = string_search_n)
            col.est = str_c(dcols[indx.pop], "_est")
            lmsfit = lms$fitted.values
            
        } else {
            # for rates, covariate control is always assumed to return residuals
            lms = regress_rates(data = data,
                        yvar = col.var,
                        sex_b = sex_t,
                        string_search = string_search_n)
            lmsfit = residuals(lms)
            col.est = col.var
        }
        data = data %>%
            mutate({{col.est}} := lmsfit)
        write_csv(data, fullpath.temp)

        newdata = read_csv(fullpath.temp)
        dcols.new = colnames(newdata)
        if (iscount) {
            indx.pop = grep(col.est, dcols.new)
            line.popmap = str_c("PopulationFile-SourceFieldMap=", indx.idcol, ",", indx.year, ",", indx.pop)
        }
    }
    
    # print(fullpath.input)
    # print(fullpath.output)
    # print(indx.idcol)
    # print(indx.case)
    # print(indx.var)
    # print(indx.xcoord)
    # print(indx.ycoord)
    modeltype = ifelse(iscount, 0, 5)

    xml_analysis_block = str_glue(
        '[Input]
;case data filename
CaseFile={fullpath.temp}
;source type (CSV=0, DBASE=1, SHAPE=2)
CaseFile-SourceType=0
;source field map (comma separated list of integers, oneCount, generatedId, shapeX, shapeY) (id, case, date, attr, weight)
CaseFile-SourceFieldMap={indx.idcol},{indx.case}{indx.var}{indx.weight}
;csv source delimiter (leave empty for space or tab delimiter)
CaseFile-SourceDelimiter=,
;csv source group character
CaseFile-SourceGrouper="
;csv source skip initial lines (i.e. meta data)
CaseFile-SourceSkip=0
;csv source first row column header
CaseFile-SourceFirstRowHeader=y
;control data filename
ControlFile=
;time precision (0=None, 1=Year, 2=Month, 3=Day, 4=Generic)
PrecisionCaseTimes=0
;study period start date (YYYY/MM/DD)
StartDate=2000/1/1
;study period end date (YYYY/MM/DD)
EndDate=2000/12/31
;population data filename
{line.popdef}
;population source type
{line.pop}
;csv source first row column header
PopulationFile-SourceFirstRowHeader=y
;population source fieldmap (comma separated list)
{line.popmap}
;coordinate data filename
CoordinatesFile={fullpath.temp}
;source type (CSV=0, DBASE=1, SHAPE=2)
CoordinatesFile-SourceType=0
;source field map (comma separated list of integers, oneCount, generatedId, shapeX, shapeY)
CoordinatesFile-SourceFieldMap={indx.idcol},{indx.xcoord},{indx.ycoord}
;csv source delimiter (leave empty for space or tab delimiter)
CoordinatesFile-SourceDelimiter=,
;csv source group character
CoordinatesFile-SourceGrouper="
;csv source skip initial lines (i.e. meta data)
CoordinatesFile-SourceSkip=0
;csv source first row column header
CoordinatesFile-SourceFirstRowHeader=y
;use grid file? (y/n)
UseGridFile=n
;grid data filename
GridFile=
;coordinate type (0=Cartesian, 1=latitude/longitude)
CoordinatesType=0

[Analysis]
;analysis type (1=Purely Spatial, 2=Purely Temporal, 3=Retrospective Space-Time, 4=Prospective Space-Time, 5=Spatial Variation in Temporal Trends, 6=Prospective Purely Temporal, 7=Seasonal Temporal)
AnalysisType=1
;model type (0=Discrete Poisson, 1=Bernoulli, 2=Space-Time Permutation, 3=Ordinal, 4=Exponential, 5=Normal, 6=Continuous Poisson, 7=Multinomial, 8=Rank, 9=UniformTime)
ModelType={modeltype}
;scan areas (1=High Rates(Poison,Bernoulli,STP); High Values(Ordinal,Normal); Short Survival(Exponential); Higher Trend(Poisson-SVTT), 2=Low Rates(Poison,Bernoulli,STP); Low Values(Ordinal,Normal); Long Survival(Exponential); Lower Trend(Poisson-SVTT), 3=Both Areas)
ScanAreas=3
;time aggregation units (0=None, 1=Year, 2=Month, 3=Day, 4=Generic)
TimeAggregationUnits=1
;time aggregation length (Positive Integer)
TimeAggregationLength=1

[Output]
;analysis main results output filename
ResultsFile={fullpath.output}
;output Google Earth KML file (y/n)
OutputGoogleEarthKML=n
;output shapefiles (y/n)
OutputShapefiles=n
;output cartesian graph file (y/n)
OutputCartesianGraph=y
;output cluster information in ASCII format? (y/n)
MostLikelyClusterEachCentroidASCII=y
;output cluster information in dBase format? (y/n)
MostLikelyClusterEachCentroidDBase=n
;output cluster case information in ASCII format? (y/n)
MostLikelyClusterCaseInfoEachCentroidASCII=y
;output cluster case information in dBase format? (y/n)
MostLikelyClusterCaseInfoEachCentroidDBase=n
;output location information in ASCII format? (y/n)
CensusAreasReportedClustersASCII=n
;output location information in dBase format? (y/n)
CensusAreasReportedClustersDBase=n
;output risk estimates in ASCII format? (y/n)
IncludeRelativeRisksCensusAreasASCII=n
;output risk estimates in dBase format? (y/n)
IncludeRelativeRisksCensusAreasDBase=n
;output simulated log likelihoods ratios in ASCII format? (y/n)
SaveSimLLRsASCII=n
;output simulated log likelihoods ratios in dBase format? (y/n)
SaveSimLLRsDBase=n
;generate Google Maps output (y/n)
OutputGoogleMaps=n

[Multiple Data Sets]
; multiple data sets purpose type (0=Multivariate, 1=Adjustment)
MultipleDataSetsPurposeType=0

[Data Checking]
;study period data check (0=Strict Bounds, 1=Relaxed Bounds)
StudyPeriodCheckType=0
;geographical coordinates data check (0=Strict Coordinates, 1=Relaxed Coordinates)
GeographicalCoordinatesCheckType=0

[Locations Network]
;locations network filename
LocationsNetworkFilename=
;use locations network file
UseLocationsNetworkFile=n
;purpose of locations network file (0=Coordinates File Override, 1=Network Definition)
PurposeLocationsNetworkFile=1

[Spatial Neighbors]
;use neighbors file (y/n)
UseNeighborsFile=n
;neighbors file
NeighborsFilename=
;use meta locations file (y/n)
UseMetaLocationsFile=n
;meta locations file
MetaLocationsFilename=
;multiple coordinates type (0=OnePerLocation, 1=AtLeastOneLocation, 2=AllLocations)
MultipleCoordinatesType=0

[Spatial Window]
;maximum spatial size in population at risk (<=50%)
MaxSpatialSizeInPopulationAtRisk=50
;restrict maximum spatial size - max circle file? (y/n)
UseMaxCirclePopulationFileOption=n
;maximum spatial size in max circle population file (<=50%)
MaxSpatialSizeInMaxCirclePopulationFile=50
;maximum circle size filename
MaxCirclePopulationFile=
;restrict maximum spatial size - distance? (y/n)
UseDistanceFromCenterOption=n
;maximum spatial size in distance from center (positive integer)
MaxSpatialSizeInDistanceFromCenter=1
;include purely temporal clusters? (y/n)
IncludePurelyTemporal=n
;window shape (0=Circular, 1=Elliptic)
SpatialWindowShapeType=1
;elliptic non-compactness penalty (0=NoPenalty, 1=MediumPenalty, 2=StrongPenalty)
NonCompactnessPenalty=1
;isotonic scan (0=Standard, 1=Monotone)
IsotonicScan=0

[Temporal Window]
;minimum temporal cluster size (in time aggregation units)
MinimumTemporalClusterSize=1
;how max temporal size should be interpretted (0=Percentage, 1=Time)
MaxTemporalSizeInterpretation=0
;maximum temporal cluster size (<=90%)
MaxTemporalSize=50
;include purely spatial clusters? (y/n)
IncludePurelySpatial=n
;temporal clusters evaluated (0=All, 1=Alive, 2=Flexible Window)
IncludeClusters=0
;flexible temporal window start range (YYYY/MM/DD,YYYY/MM/DD)
IntervalStartRange=2000/1/1,2000/12/31
;flexible temporal window end range (YYYY/MM/DD,YYYY/MM/DD)
IntervalEndRange=2000/1/1,2000/12/31

[Cluster Restrictions]
;risk limit high clusters (y/n)
RiskLimitHighClusters=n
;risk threshold high clusters (1.0 or greater)
RiskThresholdHighClusters=1
;risk limit low clusters (y/n)
RiskLimitLowClusters=n
;risk threshold low clusters (0.000 - 1.000)
RiskThresholdLowClusters=1
;minimum cases in low rate clusters (positive integer)
MinimumCasesInLowRateClusters=0
;minimum cases in high clusters (positive integer)
MinimumCasesInHighRateClusters=2

[Space and Time Adjustments]
;time trend adjustment type (0=None, 2=LogLinearPercentage, 3=CalculatedLogLinearPercentage, 4=TimeStratifiedRandomization, 5=CalculatedQuadratic, 1=TemporalNonparametric)
TimeTrendAdjustmentType=0
;time trend adjustment percentage (>-100)
TimeTrendPercentage=0
;time trend type - SVTT only (Linear=0, Quadratic=1)
TimeTrendType=0
;adjust for weekly trends, nonparametric
AdjustForWeeklyTrends=n
;spatial adjustments type (0=None, 1=SpatiallyStratifiedRandomization, 2=SpatialNonparametric)
SpatialAdjustmentType=0
;use adjustments by known relative risks file? (y/n)
UseAdjustmentsByRRFile=n
;adjustments by known relative risks file name (with HA Randomization=1)
AdjustmentsByKnownRelativeRisksFilename=

[Inference]
;p-value reporting type (Default p-value=0, Standard Monte Carlo=1, Early Termination=2, Gumbel p-value=3) 
PValueReportType=3
;early termination threshold
EarlyTerminationThreshold=50
;report Gumbel p-values (y/n)
ReportGumbel=n
;Monte Carlo replications (0, 9, 999, n999)
MonteCarloReps=999
;adjust for earlier analyses(prospective analyses only)? (y/n)
AdjustForEarlierAnalyses=n
;prospective surveillance start date (YYYY/MM/DD)
ProspectiveStartDate=2000/12/31
;perform iterative scans? (y/n)
IterativeScan=n
;maximum iterations for iterative scan (0-32000)
IterativeScanMaxIterations=100
;max p-value for iterative scan before cutoff (0.000-1.000)
IterativeScanMaxPValue=0.05

[Cluster Drilldown]
;perform detected cluster standard drilldown (y/n)
PerformStandardDrilldown=n
;perform detected cluster Bernoulli drilldown (y/n)
PerformBernoulliDrilldown=n
;minimum number of locations in detected cluster to perform drilldown (positive integer)
DrilldownMinimumClusterLocations=2
;minimum number of cases in detected cluster to perform drilldown (positive integer)
DrilldownMinimumClusterCases=10
;p-value cutoff of detected cluster to perform drilldown (0.000-1.000)
DrilldownClusterPvalueCutoff=0.05
;adjust for weekly trends, purely spatial Bernoulli drilldown
DrilldownAdjustForWeeklyTrends=n

[Miscellaneous Analysis]
CalculateOliveira=n
;number of bootstrap replications for Oliveira calculation (minimum=100, multiple of 100)
NumBootstrapReplications=1000
;p-value cutoff for clusters in Oliveira calculation (0.000-1.000)
OliveiraPvalueCutoff=0.05
;frequency of prospective analyses type (0=Same Time Aggregation, 1=Daily, 2=Weekly, 3=Monthy, 4=Quarterly, 5=Yearly)
ProspectiveFrequencyType=0
;frequency of prospective analyses  (positive integer)
ProspectiveFrequency=1

[Power Evaluation]
;perform power evaluation - Poisson only (y/n)
PerformPowerEvaluation=n
;power evaluation method (0=Analysis And Power Evaluation Together, 1=Only Power Evaluation With Case File, 2=Only Power Evaluation With Defined Total Cases)
PowerEvaluationsMethod=0
;total cases in power evaluation
PowerEvaluationTotalCases=600
;critical value type (0=Monte Carlo, 1=Gumbel, 2=User Specified Values)
CriticalValueType=0
;power evaluation critical value .05 (> 0)
CriticalValue05=0
;power evaluation critical value .001 (> 0)
CriticalValue01=0
;power evaluation critical value .001 (> 0)
CriticalValue001=0
;power estimation type (0=Monte Carlo, 1=Gumbel)
PowerEstimationType=0
;number of replications in power step
NumberPowerReplications=1000
;power evaluation alternative hypothesis filename
AlternativeHypothesisFilename=
;power evaluation simulation method for power step (0=Null Randomization, 1=N/A, 2=File Import)
PowerEvaluationsSimulationMethod=0
;power evaluation simulation data source filename
PowerEvaluationsSimulationSourceFilename=
;report power evaluation randomization data from power step (y/n)
ReportPowerEvaluationSimulationData=n
;power evaluation simulation data output filename
PowerEvaluationsSimulationOutputFilename=

[Spatial Output]
;automatically launch map viewer - gui only (y/n)
LaunchMapViewer=n
;create compressed KMZ file instead of KML file (y/n)
CompressKMLtoKMZ=n
;whether to include cluster locations kml output (y/n)
IncludeClusterLocationsKML=n
;threshold for generating separate kml files for cluster locations (positive integer)
ThresholdLocationsSeparateKML=1000
;report hierarchical clusters (y/n)
ReportHierarchicalClusters=y
;criteria for reporting secondary clusters(0=NoGeoOverlap, 1=NoCentersInOther, 2=NoCentersInMostLikely,  3=NoCentersInLessLikely, 4=NoPairsCentersEachOther, 5=NoRestrictions)
CriteriaForReportingSecondaryClusters=2
;report gini clusters (y/n)
ReportGiniClusters=n
;gini index cluster reporting type (0=optimal index only, 1=all values)
GiniIndexClusterReportingType=0
;spatial window maxima stops (comma separated decimal values[<=50%] )
SpatialMaxima=5,10,12,15,20,25,30,33,35,40,45,50
;max p-value for clusters used in calculation of index based coefficients (0.000-1.000)
GiniIndexClustersPValueCutOff=0.05
;report gini index coefficents to results file (y/n)
ReportGiniIndexCoefficents=y
;restrict reported clusters to maximum geographical cluster size? (y/n)
UseReportOnlySmallerClusters=n
;maximum reported spatial size in population at risk (<=50%)
MaxSpatialSizeInPopulationAtRisk_Reported=35
;restrict maximum reported spatial size - max circle file? (y/n)
UseMaxCirclePopulationFileOption_Reported=n
;maximum reported spatial size in max circle population file (<=50%)
MaxSizeInMaxCirclePopulationFile_Reported=35
;restrict maximum reported spatial size - distance? (y/n)
UseDistanceFromCenterOption_Reported=n
;maximum reported spatial size in distance from center (positive integer)
MaxSpatialSizeInDistanceFromCenter_Reported=1

[Temporal Output]
;output temporal graph HTML file (y/n)
OutputTemporalGraphHTML=n
;temporal graph cluster reporting type (0=Only most likely cluster, 1=X most likely clusters, 2=Only significant clusters)
TemporalGraphReportType=0
;number of most likely clusters to report in temporal graph (positive integer)
TemporalGraphMostMLC=1
;significant clusters p-value cutoff to report in temporal graph (0.000-1.000)
TemporalGraphSignificanceCutoff=0.05

[Other Output]
;report critical values for .01 and .05? (y/n)
CriticalValue=n
;report cluster rank (y/n)
ReportClusterRank=n
;print ascii headers in output files (y/n)
PrintAsciiColumnHeaders=y
;user-defined title for results file
ResultsTitle={analysis.title2}

[Elliptic Scan]
;elliptic shapes - one value for each ellipse (comma separated decimal values)
EllipseShapes=1.5,2,2.5,3,4,5
;elliptic angles - one value for each ellipse (comma separated integer values)
EllipseAngles=4,6,8,9,12,15

[Power Simulations]
;simulation methods (0=Null Randomization, 1=N/A, 2=File Import)
SimulatedDataMethodType=0
;simulation data input file name (with File Import=2)
SimulatedDataInputFilename=
;print simulation data to file? (y/n)
PrintSimulatedDataToFile=n
;simulation data output filename
SimulatedDataOutputFilename=

[Run Options]
;number of parallel processes to execute (0=All Processors, x=At Most X Processors)
NumberParallelProcesses=12
;suppressing warnings? (y/n)
SuppressWarnings=n
;log analysis run to history file? (y/n)
LogRunToHistoryFile=n
;analysis execution method  (0=Automatic, 1=Successively, 2=Centrically)
ExecutionType=0

[System]
;system setting - do not modify
Version=10.0.2\\n',
    fullpath.temp = fullpath.temp,
    indx.idcol = indx.idcol,
    indx.case = indx.case,
    indx.var = indx.var,
    #indx.year = indx.year,
    indx.casenormal = indx.casenormal,
    modeltype = modeltype,
    line.popdef = line.popdef,
    line.pop = line.pop,
    line.popmap = line.popmap,
    fullpath.output = fullpath.output,
    indx.xcoord = indx.xcoord,
    indx.ycoord = indx.ycoord,
    analysis.title2 = str_c(title.analysis, "out.txt")
    )
    print(xml_analysis_block)
    sink(file = prm.path)
    cat(xml_analysis_block)
    sink()

    system(str_c("/home/felix/SaTScan/SaTScanBatch64 ", prm.path))
    cat(str_c("Remove the artifacts of the current analysis...\n"))
    #file.copy(prm.path, str_c(dir.base, title.analysis, ".prm"))
    file.remove(prm.path)
    #file.remove(fullpath.temp)

    ## Part 2: get tsv file
    if (modeltype == 0) {
        dist.info = "Poisson"
        cn = c('cluster', 'locid', 'x', 'y', 'minor', 'major', 'angle', 'shape',
               'start_date', 'end_date', 'number_locs', 'LLR', 'pvalue', 'obs', 'ex', 'obsoverex', 'RR', 'pop')
    } else if (modeltype == 5) {
        dist.info = "Weighted normal"
        cn = c('cluster', 'locid', 'x', 'y', 'minor', 'major', 'angle', 'shape',
               'start_date', 'end_date', 'number_locs', 'LLR', 'stat_test', 'pvalue', 'observed', 'weight_in', 'mean_in', 'mean_out', 'variance', 'std', 'w_mean_in', 'w_mean_out', 'w_variance', 'w_std', 'GINI')
    }
    tss = read_fwf(str_c(dir.target, title.analysis, ".col.txt"), skip = 1)
    colnames(tss) = cn
    tss = tss %>% 
        mutate(analysis_title = title.analysis) %>%
        dplyr::select(analysis_title, 1:(ncol(.)-1))
    cat(str_glue("Modeltype={modeltype} ({dist.info})\nTarget variable={col.var}\nFile written:{fullpath.temp}\n",
                modeltype = modeltype, dist.info = dist.info, col.var = col.var, fullpath.temp = fullpath.temp))

    system(str_glue("rm {dir.target}/*", dir.target = dir.target))
    return(tss)
    }


# Run smerc elliptic.test
run_smerc_cancertype = function(data = sgg2015, 
                                population = 'n_pop_total', 
                                yvar = "Lung_total", 
                                sex_b = 'total',
                                type = 'poisson',
                                alpha = 0.01,
                                ubpop = 0.35,
                                string_search = str_c(str_c('(^p_*.*_', sex_b, '$'), '^(r_|ap_)', '^NDVI_)', sep = '|'),
                                add_var = NULL,
                                adjust = FALSE, 
                                ncores = 8) {
    library(parallel)
    data_df = st_drop_geometry(data)
 
    if (adjust) {
        reg_pois = regress_counts(data_df, population = population, yvar = yvar, 
                                    sex_b = sex_b, add_var = add_var, string_search = string_search)
        pop_in = reg_pois$fitted.values
    } 
    if (!adjust) {
        pop_in = unlist(data_df[, population])       
    }
    css = unlist(data_df[, yvar])
    cls = parallel::makeCluster(spec = ncores, type = 'PSOCK')
    eltest = smerc::elliptic.test(st_coordinates(st_centroid(data)), 
            cases = css, 
            pop = pop_in,
            shape = c(1, 1.5, 2, 2.5, 3, 4, 6),
            nangle = c(1, 4, 6, 8, 9, 12, 18),
            ubpop = ubpop,
            type = type,
            #min.cases = min(css) + 1,
            nsim = 499,
            alpha = alpha,
            cl = cls)
    parallel::stopCluster(cls)
    return(eltest)
}


# smerc cluster to general maps with a tmap object
tmap_smerc = function(basemap, smc, threshold = 2, significance = 0.01, alpha = 0.4, return_ellipses = FALSE) {
    rotate = function(a) {a = a*pi/180; matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)}
    
    library(tmap)
    basemap$cluster = NA
    smc_ncl = sapply(smc$clusters, function(x) length(x$locids))
    smc_ncl_addr = (smc_ncl >= threshold)
    # p-value extraction
    smc_pval = sapply(smc$clusters, function(x) x$pvalue)
    smc_pval_addr = (smc_pval <= significance)
    if (!is.null(threshold)) {
        smc_pval_addr = smc_pval_addr * smc_ncl_addr
    }
    smc_pval_addr = grep(1, as.integer(smc_pval_addr))

    if (length(smc_pval_addr) == 0) {
        tm_cluster = tm_shape(basemap) +
            tm_borders(col = 'dark grey', lwd = 0.8)
        return(tm_cluster)
    }

    smc_shps = vector('list', length = length(smc_pval_addr))
    # loop through
    if (!return_ellipses) {
        for (i in seq_len(length(smc_pval_addr))) {
            cl_ilocs = smc$clusters[[smc_pval_addr[i]]]$locids
            basemap$cluster[cl_ilocs] = i
        }
        basemap$cluster = as.factor(basemap$cluster)
        tm_cluster = tm_shape(basemap) +
            #tm_fill('cluster', pal = 'Set3', colorNA = 'transparent', showNA = FALSE, alpha = alpha) +
            tm_borders(col = 'dark grey', lwd = 0.3) +
            tm_polygons('cluster', pal = 'Set3', colorNA = 'transparent', lwd = 1.5, showNA = FALSE, alpha = alpha)
    }
    if (return_ellipses) {
        smc_shps = 
        smc$clusters %>%
            .[smc_pval_addr] %>%
            lapply(function(x) {
                cntr = st_geometry(st_point(x$centroid))
                ell = (nngeo::st_ellipse(pnt = cntr,
                                  ex = x$semiminor_axis,
                                  ey = x$semimajor_axis,
                                  res = 80) - cntr) * rotate(x$angle) + cntr
                ell = ell %>%
                    st_sf %>%
                    mutate(
                           RR = x$rr,
                           statistic = x$test_statistic,
                           p_value = x$pvalue)
                return(ell)})
        smc_shps = do.call(rbind, smc_shps) %>%
            mutate(cluster = seq_len(length(smc_pval_addr)),
                   class_cl = c('#C41B40', rep('#D980AF', nrow(.)-1)))
        st_crs(smc_shps) = st_crs(basemap)
        #smc_shps = st_transform(smc_shps, st_crs(basemap))
        # mf_init(x = basemap)
        # mf_map(basemap, add = TRUE)
        # mf_map(smc_shps, type = 'typo', var = 'class_cl', leg_title = 'Primary', leg_pos = NA, alpha = 0.5,add = TRUE)
        
        tm_cluster = tm_shape(basemap) +
            tm_borders(col = 'dark grey', lwd = 0.8) +
            tm_shape(smc_shps) +
            tm_fill('class_cl', colorNA = 'transparent', showNA = FALSE, alpha = alpha) +
            tm_borders('#FF0000', lwd = 1.5)
        #tm_cluster = smc_shps
    }

    # tmap
    return(tm_cluster)
}

## Plot SaTScan results directly
tmap_satscan = function(basemap, sats, threshold = 2, significance = 0.01, alpha = 0.4, return_ellipses = TRUE, gini = FALSE) {
    rotate = function(a) {a = a*pi/180; matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)}
    
    library(tmap)
    basemap$cluster = NA

    if (nrow(sats) == 0 | (gini & sum(sats$GINI) == 0)) {
        tm_cluster = tm_shape(basemap) +
            tm_borders(col = 'dark grey', lwd = 0.8)
        return(tm_cluster)
    }
    if (gini) {
        sats = sats %>% filter(GINI)
    }
    # test if the angle was weirdly reflected
    sats = sats %>% mutate(angle = if_else(angle < 0, angle + 360, angle))
    nclust = nrow(sats)
    smc_shps = vector('list', length = length(sats$pvalue))
    # loop through
    if (return_ellipses) {
        smc_shps = 
        sats %>%
            split(., .$cluster) %>%
            lapply(function(x) {
                cntr = st_geometry(st_point(c(x$x, x$y)))
                ell = (nngeo::st_ellipse(pnt = cntr,
                                  ex = x$minor,
                                  ey = x$major,
                                  res = 80) - cntr) * rotate(x$angle) + cntr
                ell = ell %>%
                    st_sf %>%
                    mutate(
                           RR = x$LLR,
                           statistic = x$stat_test,
                           p_value = x$pvalue)
                return(ell)})
        smc_shps = do.call(rbind, smc_shps) %>%
            mutate(cluster = seq_len(nclust),
                   class_cl = c('#C41B4050', rep(NA, nrow(.)-1))) #'#D980AF'
        st_crs(smc_shps) = st_crs(basemap)
        #smc_shps = st_transform(smc_shps, st_crs(basemap))
        # mf_init(x = basemap)
        # mf_map(basemap, add = TRUE)
        # mf_map(smc_shps, type = 'typo', var = 'class_cl', leg_title = 'Primary', leg_pos = NA, alpha = 0.5,add = TRUE)
        
        if (nrow(smc_shps) == 1) {
            tm_cluster = tm_shape(basemap) +
                    tm_borders(col = 'dark grey', lwd = 0.8) +
                    tm_shape(smc_shps[1,]) +
                    tm_borders('#FF0000', lwd = 1.2) +
                    tm_fill('class_cl', colorNA = 'transparent', showNA = FALSE)
        } else {
            tm_cluster = tm_shape(basemap) +
                tm_borders(col = 'dark grey', lwd = 0.8) +
                tm_shape(smc_shps[-1,]) +
                tm_borders("#F200AA") +
                tm_shape(smc_shps[1,]) +
                tm_borders('#FF0000', lwd = 1.2) +
                tm_fill('class_cl', colorNA = 'transparent', showNA = FALSE)
        }
    }

    # tmap
    return(tm_cluster)
}



# Run dclustm analysis
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



# dclust cluster to general maps with a tmap object
# Plot dclustm results with a basemap
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



## LISA mapping
map_lisa = function(map, lisa.var, vset = NULL, degree.s = 1, signif = 0.05, title = "", special_label = TRUE, adjust = FALSE, for_residuals = TRUE) {

    if (adjust) {
        sex_t = str_extract(lisa.var, "(male|female|total)")
        string_search_n =
            switch(vset,
                    set1 = str_c('^(p_65p|p_hbac)*.*_', sex_t, '$'),
                    set2 = str_c(str_c('^p_*.*_', sex_t, '$'), '^ap_', '^NDVI_', sep = '|'),
                    set3 = str_c(str_c('^p_*.*_', sex_t, '$'), '^r_(?!physmid)', '^ap_', '^NDVI_', sep = '|'),
                    set4 = str_c(str_c('^p_*.*_', sex_t, '$'), '^r_', '^n_pw', '^ap_', '^NDVI_', sep = '|'))

        if (model == "normal") {
            lms = regress_rates(data = data,
                        yvar = col.var,
                        sex_b = sex_t,
                        string_search = string_search_n)
        } else if (model == "poisson") {
            lms = regress_counts(data = data,
                        yvar = col.var,
                        population = col.case,
                        sex_b = sex_t,
                        string_search = string_search_n)
        }
        lmsfit = lms$fitted.values
        if (for_residuals) {
            lmsfit = residuals(lms)
        }
        map_df = map_df %>%
            mutate({{lisa.var}} := lmsfit)
    }

    map_df = st_drop_geometry(map)
    lw = nb2listw(poly2nb(map), zero.policy = TRUE)
    if (degree.s > 1) {
        lw = nblag_cumul(nblag(poly2nb(map), degree.s))
        lw = nb2listw(map[,id.col])
    }
    v.o = as.vector(scale(unlist(map_df[,lisa.var])))
    v.lag = lag.listw(lw, v.o, zero.policy = TRUE)
    lisa_res = 
        localmoran(v.o,
                listw = lw,
                zero.policy = TRUE)
    lisa_resdf = as.data.frame(lisa_res)
    colnames(lisa_resdf) = c("I", "EI", 'VarI', 'ZI', 'pvalue')

    quadrant = rep(NA, length(v.o))
    quadrant[v.o >0 & v.lag>0] <- "HH"      
    quadrant[v.o <0 & v.lag<0] <- "LL"     
    quadrant[v.o <0 & v.lag>0] <- "LH"
    quadrant[v.o >0 & v.lag<0] <- "HL"
    quadrant[!is.na(quadrant) & (lisa_resdf[,5]>signif)] <- NA

    colpal = c("blue", "deepskyblue1", "pink", "red")
    namepal = c("HH", "HL", "LH", "LL")
    vals_real = sort(unique(quadrant))
    pal_fin = colpal[grep(str_c("(", str_c(vals_real, collapse = "|"), ")"), namepal)]

    map_lisa = bind_cols(map, lisa_resdf)
    map_lisa$LISA = factor(quadrant, levels = vals_real)
    
    maptitle = str_c(title)
    map_result = tm_shape(map_lisa) +
        tm_fill("LISA", style = 'cat', palette = pal_fin, colorNA = 'transparent', textNA = "Insignificant") +
        tm_borders('grey', lwd = 0.66) +
        tm_layout(frame = FALSE, title = maptitle)
    map_result
 
}
