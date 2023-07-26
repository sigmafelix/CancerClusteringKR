source('./Code/Base/base_functions.R')
#options(repos = 'https://cran.seoul.go.kr')
if (!require(pacman)) { install.packages('pacman') ; library(pacman)} 

p_load(tidyverse, sf, stars, raster, starsExtra, readxl, here, tmap, stargazer, smerc, DClusterm, kableExtra, patchwork, rmapshaper, spdep)

dirpattern = "/home/%s/"
username = 'felix'
basedir = sprintf(dirpattern, username)
rdatafiles = list.files(path = str_c(basedir, 'Documents/GP/'), pattern = '*.RData', full.names = TRUE)
geopath = str_c(basedir, "OneDrive/Data/Korea/")
drive = str_c(basedir, "OneDrive/NCC_Project/CancerClustering/")
geopath = str_c(basedir, "OneDrive/Data/Korea/")
dbdir = drive  
rdsdir = drive #sprintf("/mnt/c/Users/%s/OneDrive/NCC_Project/CancerClustering/", username)

sf_use_s2(F)

load(str_c(drive, '/Manuscript/Clustering_Base_sf_091522.RData'))

covar_origin_10_fc = covar_origin_10_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430, 39010, 39020)) %>%
    ms_simplify(keep = 0.1, keep_shapes = TRUE)
covar_origin_05_fc = covar_origin_05_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430, 39010, 39020)) %>%
    ms_simplify(keep = 0.1, keep_shapes = TRUE)
covar_origin_00_fc = covar_origin_00_fc %>%
    filter(!sgg_cd_c %in% c(23320, 37430, 39010, 39020)) %>%
    ms_simplify(keep = 0.1, keep_shapes = TRUE)

clear_input = function(filepath, alpha = 0.01, locmin = 2) {
    # alpha: significance level
    # locmin: how many points should be included at least?
    res = read_rds(filepath)
    labels = as.data.frame(str_split_fixed(res$analysis_title, "_", 6))
    colnames(labels) = c("cancertype", "period", "measure", "target", "sex", "vset")
    res = res %>%
          bind_cols(labels, .) %>%
          filter(pvalue <= alpha & number_locs >= locmin)
    resl = res %>% split(., .$analysis_title)
    return(resl)
}


### FigS1
target_sdname = c("Seoul (Capital)", "Busan", "Daegu", "Incheon", "Gwangju", "Daejeon", "Ulsan",
        "Sejong", "Gyeonggi", "Gangwon", "North Chungcheong", "South Chungcheong",
        "North Jeolla", "South Jeolla", "North Gyeongsang", "South Gyeongsang")

covar_origin_10_fcc = covar_origin_10_fc %>%
    dplyr::filter(!sgg_cd_c %in% c(23320, 37430, 39010, 39020)) %>%
    sf::st_cast("MULTIPOLYGON") %>%
    sf::st_buffer(0)
covar_origin_10_fcc = # %>%
    rmapshaper::ms_simplify(input = covar_origin_10_fcc, keep = 0.5, keep_shapes = TRUE)
covar_origin_10_fcc = covar_origin_10_fcc %>%
    mutate(sdcd = str_sub(sgg_cd_c, 1, 2)) %>%
    mutate(sdname = plyr::mapvalues(sdcd, unique(sdcd),
        target_sdname)) %>%
    mutate(sdname = factor(sdname, levels = target_sdname))

mf_map(covar_origin_10_fcc, var = "sdname", leg_title = NA, type = "typo")


korea_metros = covar_origin_10_fcc %>%
    group_by(sdname, sdcd) %>%
    summarize(N = n(), geom = st_union(geom)) %>%
    ungroup %>%
    filter(as.numeric(sdcd) < 29)

covar_origin_10_fccc = covar_origin_10_fcc %>%
    st_centroid %>%
    group_by(sdname) %>%
    summarize(N = n()) %>%
    ungroup %>%
    st_centroid
covar_origin_10_fccc$geom[[9]] = covar_origin_10_fccc$geom[[9]] + c(3e3, -2e4)

geodata::geodata_path("./Data/")
prk = geodata::gadm("PRK", 0)
prksf = prk %>%
    st_as_sf %>%
    st_transform("EPSG:5179")
kor = geodata::gadm("KOR", 1)
ggisf = kor %>%
    st_as_sf %>%
    st_transform("EPSG:5179") %>%
    filter(NAME_1 == "Gyeonggi-do")

sdb = read_sf("/home/felix/OneDrive/NCC_Project/CancerClustering/Data/SD_boundaries.gpkg")
# prkgsf = st_read("/mnt/c/Users/sigma/Documents/PRKG_SF.gpkg")
prkgsf = prksf

tm_sido = 
    tm_shape(covar_origin_10_fcc) +
    tm_borders("transparent") +
    tm_shape(prkgsf) +
    tm_fill('grey', alpha = 0.5) +
    tm_shape(covar_origin_10_fcc) +
    tm_fill("sdname", pal = "Set3", legend.show = FALSE) +
    tm_borders(lwd = 0.3, col = "dark grey") +
    tm_scale_bar(breaks = c(0, 25, 50), text.size = 0.8) +
    tm_layout(title = NA,#"Stomach incidence rate\n(female)",
                  frame = FALSE,
                  legend.title.color = 'transparent',
                  legend.title.size = 0.01,
                  #legend.title.size = 1,
                  legend.text.size = 1,
                  panel.label.size = 1.4,
                  fontfamily = "Pretendard",
                  outer.margins = c(0.005, 0.005, 0.005, 0.005)) +
    tm_shape(sdb) +
    tm_lines(col = "black", lwd = 1.5, lty = "solid", alpha = 0.5) +
    tm_shape(korea_metros) +
    tm_borders("orange", lwd = 2, alpha = 0.5) +
    tm_shape(covar_origin_10_fccc) +
    tm_text("sdname", size = 1, col = 'black', fontface = 'bold', fontfamily = "Pretendard", shadow = TRUE) 

tmap_save(tm_sido, filename = str_c(drive, "/Manuscript/Provincial_Map_07232023.png"), 
        width = 10, height = 12, units = 'in', dpi = 300)

cofc3= covar_origin_10_fc %>% mutate(period = "2009-2013") %>% st_transform(5179) %>%
    dplyr::select(sgg_cd_c, period, matches("ragest_(i|d)_(Lung|Stomach)_*"), geom) %>%
    dplyr::select(-contains("_total_")) %>%
    pivot_longer(cols = seq(3, ncol(.)-1)) %>%
    mutate(valuetype = ifelse(grepl("_i_", name), "Incidence", "Mortality"),
           sex = ifelse(grepl("_(female)_", name), "Female", "Male"),
           cancertype = ifelse(grepl("_Lung_", name), "Lung", "Stomach"),
           cancerlabel = paste(cancertype, " ", valuetype, sep = ""))
cofc2= covar_origin_05_fc %>% mutate(period = "2004-2008") %>% st_transform(5179) %>%
    dplyr::select(sgg_cd_c, period, matches("ragest_(i|d)_(Lung|Stomach)_*"), geom) %>%
    dplyr::select(-contains("_total_")) %>%
    pivot_longer(cols = seq(3, ncol(.)-1)) %>%
    mutate(valuetype = ifelse(grepl("_i_", name), "Incidence", "Mortality"),
           sex = ifelse(grepl("_(female)_", name), "Female", "Male"),
           cancertype = ifelse(grepl("_Lung_", name), "Lung", "Stomach"),
           cancerlabel = paste(cancertype, " ", valuetype, sep = ""))
cofc1= covar_origin_00_fc %>% mutate(period = "1999-2003") %>% st_transform(5179) %>%
    dplyr::select(sgg_cd_c, period, matches("ragest_(i|d)_(Lung|Stomach)_*"), geom) %>%
    dplyr::select(-contains("_total_")) %>%
    pivot_longer(cols = seq(3, ncol(.)-1)) %>%
    mutate(valuetype = ifelse(grepl("_i_", name), "Incidence", "Mortality"),
           sex = ifelse(grepl("_(female)_", name), "Female", "Male"),
           cancertype = ifelse(grepl("_Lung_", name), "Lung", "Stomach"),
           cancerlabel = paste(cancertype, " ", valuetype, sep = ""))



cofcs = bind_rows(
    cofc1, cofc2, cofc3
) %>%
    st_as_sf(sf_column_name = "geom") %>%
    mutate(cancerlabel = factor(cancerlabel,
        levels = c("Lung Incidence", "Lung Mortality", "Stomach Incidence", "Stomach Mortality")))




concat_tmaps = function(sfd, cn = 'value', title_custom = "Age standardized rate"){
    
    tmap_style("white")
    # tm_shape(sfd) +
    #         tm_fill(cn, n = 5, style = 'quantile', pal = "PuRd",
    #                 title = title_custom) +
    #         tm_borders(lwd = 0.05, col = 'grey') +
    #         tm_facets(c("period", 'cancerlabel')) +
    #         tm_layout(
    #             legend.position = F,
    #             title.position = F,
    #         legend.outside = TRUE,
    #         legend.outside.position = c('bottom'), 
    #         outer.margins = c(-0.1, 0.01, 0.01, 0.01),
    #         panel.label.bg.color = NA,
    #         frame.lwd = NA,
    #         panel.label.size = 1.4, panel.label.height = 1.5,
    #         legend.title.size = 1.2,
    #         legend.text.size = 1)

    mapsep = lapply(split(sfd, sfd$cancerlabel),
    function(k) {
        # Old: PuRd, New: Reds (colorblind safe); 07/23/2023
        tm_shape(k) +
            tm_fill(cn, n = 5, style = 'quantile', pal = "Reds",
                    title = title_custom) +
            tm_borders(lwd = 0.05, col = 'grey') +
            tm_facets(c("period", "cancerlabel"), nrow = 3) +
            tm_layout(#legend.position = "bottom",
            legend.outside = TRUE,
            legend.outside.position = "bottom", 
            outer.margins = c(0.1, 0.01, 0.01, 0.01),
                panel.show = T,
                frame.lwd = NA,
                panel.label.bg.color = NA,
                frame = F,
            panel.label.size = 1.2, panel.label.height = 1.5,
            legend.title.size = 1.5,
            legend.text.size = 1)
    }    )
    tmap_arrange(mapsep[[1]], mapsep[[2]], mapsep[[3]], mapsep[[4]],
        ncol = 4, outer.margins = c(-0.125, 0.02,0.0,0.02))
}


cofcsmap1 = concat_tmaps(cofcs %>% filter(sex == "Male"), 'value')
cofcsmap2 = concat_tmaps(cofcs %>% filter(sex == "Female"), 'value',
                title_custom = "Age standardized rate")

tmap_save(cofcsmap1, filename = str_c(drive, "/Manuscript/ASR_Map_Male_3periods.png"),
        width = 12, height = 15, units = 'in', dpi = 300, pointsize = 16)
tmap_save(cofcsmap2, filename = str_c(drive, "/Manuscript/ASR_Map_Female_3periods.png"),
        width = 12, height = 15, units = 'in', dpi = 300, pointsize = 16)


## Separate legends
cofcsmap1e = cofcs %>% filter(sex == "Male") %>%
    tm_shape(.) +
        tm_fill('value', n = 5, style = 'quantile', pal = "Reds",
                title = NA) +
        tm_borders(lwd = 0.05, col = 'grey') +
        tm_facets(c("period", "cancerlabel"), ncol = 4, nrow = 3, free.scales.fill = TRUE) +
        tm_layout(#legend.position = "bottom",
            #legend.outside = TRUE,
            #legend.outside.position = "bottom", 
            legend.position = c(-0.02, 0.7),
            outer.margins = c(0.01, 0.01, 0.01, 0.01),
            panel.show = T,
            frame.lwd = NA,
            panel.label.bg.color = NA,
            frame = F,
            panel.label.size = 1.2, 
            panel.label.height = 1.5,
            legend.format = list(text.separator = "–"),
            legend.title.color = 'transparent',
            legend.title.size = 0.01,
            legend.text.size = 0.7)
cofcsmap2e = cofcs %>% filter(sex == "Female") %>%
    tm_shape(.) +
        tm_fill('value', n = 5, style = 'quantile', pal = "Reds",
                title = NA) +
        tm_borders(lwd = 0.05, col = 'grey') +
        tm_facets(c("period", "cancerlabel"), ncol = 4, nrow = 3, free.scales.fill = TRUE) +
        tm_layout(#legend.position = "bottom",
            #legend.outside = TRUE,
            #legend.outside.position = "bottom", 
            legend.position = c(-0.02, 0.7),
            outer.margins = c(0.01, 0.01, 0.01, 0.01),
            panel.show = T,
            frame.lwd = NA,
            panel.label.bg.color = NA,
            frame = F,
            panel.label.size = 1.2, 
            panel.label.height = 1.5,
            legend.format = list(text.separator = "–"),
            legend.title.color = 'transparent',
            legend.title.size = 0.01,
            legend.text.size = 0.7)

tmap_save(cofcsmap1e, filename = str_c(drive, "/Manuscript/ASR_Map_Male_seplegend_072323.png"),
        width = 15, height = 12, units = 'in', dpi = 300, pointsize = 24)
tmap_save(cofcsmap2e, filename = str_c(drive, "/Manuscript/ASR_Map_Female_seplegend_072323.png"),
        width = 15, height = 12, units = 'in', dpi = 300, pointsize = 24)



## boxplots


c0010_ggm = cofcs %>% filter(sex == "Male") %>%
    ggplot(data = ., 
           mapping = aes(x = period, y = value)) +
        theme_minimal() +
        geom_boxplot() +
        facet_wrap(sex ~ cancerlabel, scales = "free_y", nrow = 1) +
        theme(text = element_text(family = 'Pretendard'),
              axis.title = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank())
c0010_ggf = cofcs %>% filter(sex == "Female") %>%
    ggplot(data = ., 
           mapping = aes(x = period, y = value)) +
        theme_minimal() +
        geom_boxplot() +
        xlab("Period") +
        facet_wrap(sex ~ cancerlabel, scales = "free_y", nrow = 1) +
        theme(text = element_text(family = 'Pretendard'),
              axis.title = element_blank(),
              axis.text.x = element_text(size = rel(1), hjust = 1, angle = 60))

c0010_gg = 
patchwork::wrap_elements(wrap_plots(c0010_ggm, c0010_ggf,
            nrow = 2, byrow = TRUE, heights=c(0.48, 0.52))) +
    labs(tag = "Age-standardized rate (per 100,000)") +
    theme(plot.tag = element_text(size = rel(1), angle = 90, family = "Pretendard"),
          plot.tag.position = 'left')
            #grid::textGrob('Period', x = 0.52, y = 0.02)
# grid::grid.draw(grid::textGrob('Age-standardized rate (per 100,000)', x = 0.03, rot = 90))
# grid::grid.draw(grid::textGrob('Period', x = 0.52, y = 0.02))

        # ylab('Age-standardized rate (per 100,000)') +
        # xlab("Period")
ggsave(c0010_gg, 
       filename = str_c(drive, "/Manuscript/Fig1_period_rate_sex.png"),
       width = 10, height = 6, units = 'in', dpi = 508, scale = 0.75)




###
fulllist = c(
    "p_65p_male", "p_65p_female",
    "p_hbac_male", "p_hbac_female",
    "r_smoking", "r_alcoholmonth", "r_obesity", "r_walking",
    "NDVI_mean", 
    "ap_PM10_pred", "ap_NO2_pred",
    "ap_sum_em_co", "ap_sum_em_nox", "ap_sum_em_sox", "ap_sum_em_tsp", "ap_sum_em_pm10", "ap_sum_em_voc", "ap_sum_em_nh3"
)

fulllist_name = 
    c("65 years old or older\n(male, %)",
    "65 years old or older\n(female, %)",
    "Education attainment\nbachelor's degree and higher\n(male, %)",
    "Education attainment\nbachelor's degree and higher\n(female, %)",
    "People who have smoked 100 cigarettes\nfor the entire life and currently smoke\n(%)",
    "People who consumed alcohol\nat least once a month in one recent year\n(%)",
    "People who deemed themselves\nas a little or very obese\n(%)",
    "People who have 30 minutes or more of\nwalk for at least five days a week\n(%)",
    "Average annual pixel-level\nmedian Normalized Difference Vegetation Index",
    "Annual average concentrations of\nparticles 10 micrometer or less\nin diameter (PM10)\n(μg/m3)",
    "Annual average concentrations of\nnitrogen oxide (NO2)\n(ppm)",
    "Total amount of\ncarbon monoxide emissions estimated from\npoint, line, and area sources\n(kg)",
    "Total amount of\nnitrogen oxide emissions\n(kg)",
    "Total amount of\nsulfur oxide emissions\n(kg)",
    "Total amount of\ntotal suspended particles emissions\n(kg)",
    "Total amount of\nPM10 emissions\n(kg)",
    "Total amount of\nvolatile organic compounds emissions\n(kg)",
    "Total amount of\nammonia emissions\n(kg)")

covar_origin_00_fcdf = read_csv(str_c(drive, "/satscan_base_data_period1.csv"))
covar_origin_05_fcdf = read_csv(str_c(drive, "/satscan_base_data_period2.csv"))
covar_origin_10_fcdf = read_csv(str_c(drive, "/satscan_base_data_period3.csv"))


covarlist00 = colnames(covar_origin_00_fcdf)[c(4:5,11:12,13:22)]
covarlist05 = colnames(covar_origin_05_fcdf)[c(4:5,11:12,15:17, 19:29)]
covarlist10 = colnames(covar_origin_10_fcdf)[c(4:5,11:12,15,16,17,19,20,21,22,23,24,25,26,27,28,29)]

par(mfrow = c(5,4))

export_mm = function(insf, covarlist, covarnamelist, filename) {
    #oldpar = par()
    png(filename, 15, 15, "in", 16, res = 508)
    par(mfrow = c(5,4))
    for (ce in seq_len(length(covarlist))) {
        mf_init(insf)
        mf_map(insf[,covarlist[ce]], covarlist[ce], 
                leg_title = covarnamelist[ce],
                type = "choro", nbreaks = 5, border = 'grey', lwd = 0.05, add = TRUE)
        # mf_title(ce)
    }
    dev.off()
    par(new=TRUE)
}




export_mm(covar_origin_00_fc, fulllist[-5:-8], fulllist_name[-5:-8], str_c(drive,"/Manuscript/Covar_dist_00_Re.png"))
export_mm(covar_origin_05_fc, fulllist, fulllist_name, str_c(drive,"/Manuscript/Covar_dist_05_Re.png"))
export_mm(covar_origin_10_fc, fulllist, fulllist_name, str_c(drive,"/Manuscript/Covar_dist_10_Re.png"))





####

p_load(corrplot)

# y
yvec_i = c("ragest_i_Lung_male_", "ragest_i_Lung_female_", "ragest_i_Stomach_male_", "ragest_i_Stomach_female_")
yvec_d = c("ragest_d_Lung_male_", "ragest_d_Lung_female_", "ragest_d_Stomach_male_", "ragest_d_Stomach_female_")
yvec_p1 = c(str_c(yvec_i, 1), str_c(yvec_d, 1))
yvec_p2 = c(str_c(yvec_i, 2), str_c(yvec_d, 2))
yvec_p3 = c(str_c(yvec_i, 3), str_c(yvec_d, 3))

plabs = str_c("Period ", 1:3)

## ypairs
yvec_ilm = str_c(yvec_i[1], 1:3)
yvec_ilf = str_c(yvec_i[2], 1:3)
yvec_ism = str_c(yvec_i[3], 1:3)
yvec_isf = str_c(yvec_i[4], 1:3)
yvec_mlm = str_c(yvec_d[1], 1:3)
yvec_mlf = str_c(yvec_d[2], 1:3)
yvec_msm = str_c(yvec_d[3], 1:3)
yvec_msf = str_c(yvec_d[4], 1:3)


# x
xvec_m1 = c(
    "p_65p_male", "p_65p_female",
    "p_hbac_male", "p_hbac_female",
    #"r_smoking", "r_alcoholmonth", "r_obesity", "r_walking",
    "NDVI_mean", 
    "ap_PM10_pred", "ap_NO2_pred",
    "ap_sum_em_co", "ap_sum_em_nox", "ap_sum_em_sox", "ap_sum_em_tsp", "ap_sum_em_pm10", "ap_sum_em_voc", "ap_sum_em_nh3"
)
xvec_f1 = c(
    "p_65p_female",
    "p_hbac_female",
    #"r_smoking", "r_alcoholmonth", "r_obesity", "r_walking",
    "NDVI_mean", 
    "ap_PM10_pred", "ap_NO2_pred",
    "ap_sum_em_co", "ap_sum_em_nox", "ap_sum_em_sox", "ap_sum_em_tsp", "ap_sum_em_pm10", "ap_sum_em_voc", "ap_sum_em_nh3"
)
xvec_m2 = c(
    "p_65p_male", "p_65p_female",
    "p_hbac_male", "p_hbac_female",
    "r_smoking", "r_alcoholmonth", "r_obesity", "r_walking",
    "NDVI_mean", 
    "ap_PM10_pred", "ap_NO2_pred",
    "ap_sum_em_co", "ap_sum_em_nox", "ap_sum_em_sox", "ap_sum_em_tsp", "ap_sum_em_pm10", "ap_sum_em_voc", "ap_sum_em_nh3"
)
xvec_f2 = c(
    "p_65p_female",
    "p_hbac_female",
    "r_smoking", "r_alcoholmonth", "r_obesity", "r_walking",
    "NDVI_mean", 
    "ap_PM10_pred", "ap_NO2_pred",
    "ap_sum_em_co", "ap_sum_em_nox", "ap_sum_em_sox", "ap_sum_em_tsp", "ap_sum_em_pm10", "ap_sum_em_voc", "ap_sum_em_nh3"
)


fulllist_name2 = 
    c("65 years old or older\n(male, %)",
    "65 years old or older\n(female, %)",
    "Education attainment\nbachelor's degree and higher (male, %)",
    "Education attainment\nbachelor's degree and higher (female, %)",
    "People who have smoked 100 cigarettes\nfor the entire life and currently smoke (%)",
    "People who consumed alcohol\nat least once a month in one recent year (%)",
    "People who deemed themselves\nas a little or very obese (%)",
    "People who have 30 minutes or more of\nwalk for at least five days a week (%)",
    "Average annual pixel-level\nmedian Normalized Difference Vegetation Index",
    "Annual average concentrations of\nparticles 10 micrometer or less in diameter (PM10) (μg/m3)",
    "Annual average concentrations of\nnitrogen oxide (NO2) (ppm)",
    "Total amount of\ncarbon monoxide emissions estimated from point, line, and area sources (kg)",
    "Total amount of\nnitrogen oxide emissions (kg)",
    "Total amount of\nsulfur oxide emissions (kg)",
    "Total amount of\ntotal suspended particles emissions (kg)",
    "Total amount of\nPM10 emissions (kg)",
    "Total amount of\nvolatile organic compounds emissions (kg)",
    "Total amount of\nammonia emissions (kg)")



# xlab
xvec_lab1 = c('Senior (male)', 'Senior (female)', '>Bachelor (male)', '>Bachelor (female)', 
            'PM10_ex', 'NO2_ex',
            'Average NDVI', 'CO', 'NOx', 'SOx', 'TSP', 'PM10_em',
            'VOC', 'NH3')
xvec_lab2 = c('Senior (male)', 'Senior (female)', '>Bachelor (male)', '>Bachelor (female)', 
            'Smoking (%)', 'Alcohol consumption (%)', 'Self-reported obesity (%)', 'Walking (%)', 
            'Average NDVI', 'PM10_ex', 'NO2_ex',
            'CO', 'NOx', 'SOx', 'TSP', 'PM10_em',
            'VOC', 'NH3')
# xvec_lab1 = fulllist_name2[-5:-8]
# xvec_lab2 = fulllist_name2

# main
corr_single = function(dat, y, x, ylab, xlab, title) {
    subdat = st_drop_geometry(dat) %>% 
        dplyr::select(all_of(c(y, x))) %>%
        data.frame
    colnames(subdat) = c(ylab, xlab)
    subcor = subdat %>%
        cor.mtest(., conf.level = 0.95)
    
    corrplot(cor(subdat), #p.mat = subcor$p, 
            method = 'number', diag = F, #insig = 'pch', 
            is.corr = TRUE,
            insig = 'pch',
            pch.col = 'black',
            pch.cex = 0.5,
            addCoef.col = 'black',
            title = title, type = 'upper', tl.cex = 0.8, mar = c(0, 0, 1.2, 0), number.cex = 0.9, tl.col = 'dark blue')#, tl.pos = 'td', cl.pos = 'r', tl.srt = 60)
}


ylab_all = 
c("ASIR (Lung, Male)", "ASIR (Lung, Female)", "ASIR (Stomach, Male)", "ASIR (Stomach, Female)", 
"ASMR (Lung, Male)", "ASMR (Lung, Female)", "ASMR (Stomach, Male)", "ASMR (Stomach, Female)")



png(filename = str_c(rdsdir, "Manuscript/Total_corrplot_061823.png"), 
        width = 15, height = 15, units = 'in', pointsize = 10, res = 508,
        type = "cairo-png")
par(mfrow = c(2, 2), mar = c(0.1, 0.1, 0.1, 0.1), mai = rep(0.01, 4))
corr_single(covar_origin_00_fc, yvec_p1, xvec_m1, ylab_all, xvec_lab1, 'Period 1')
corr_single(covar_origin_05_fc, yvec_p2, xvec_m2, ylab_all, xvec_lab2, 'Period 2')
corr_single(covar_origin_10_fc, yvec_p3, xvec_m2, ylab_all, xvec_lab2, 'Period 3')
dev.off()

png(filename = str_c(rdsdir, "Manuscript/Total_corrplot_panel1_052223.png"), 
        width = 10, height = 10, units = 'in', pointsize = 15, res = 508,
        type = "cairo-png")
par(mfcol = c(1, 1), mar = c(0.1, 0.1, 0.1, 0.1))
corr_single(covar_origin_00_fc, yvec_p1, xvec_m1, ylab_all, xvec_lab1, 'Period 1')
dev.off()

png(filename = str_c(rdsdir, "Manuscript/Total_corrplot_panel2_052223.png"), 
        width = 10, height = 10, units = 'in', pointsize = 15, res = 508,
        type = "cairo-png")
par(mfcol = c(1, 1), mar = c(0.1, 0.1, 0.1, 0.1))
corr_single(covar_origin_05_fc, yvec_p2, xvec_m2, ylab_all, xvec_lab2, 'Period 2')
dev.off()

png(filename = str_c(rdsdir, "Manuscript/Total_corrplot_panel3_052223.png"), 
        width = 10, height = 10, units = 'in', pointsize = 15, res = 508,
        type = "cairo-png")
par(mfcol = c(1, 1), mar = c(0.1, 0.1, 0.1, 0.1))
corr_single(covar_origin_10_fc, yvec_p3, xvec_m2, ylab_all, xvec_lab2, 'Period 3')
dev.off()

