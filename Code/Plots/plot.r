source('../Base/base_functions.R')
options(repos = 'https://cran.seoul.go.kr')
if (!require(pacman)) { install.packages('pacman') ; library(pacman)} 

p_load(tidyverse, sf, stars, raster, starsExtra, readxl, here, tmap, stargazer, smerc, DClusterm, kableExtra, patchwork, rmapshaper, spdep)

dirpattern = "/mnt/c/Users/%s/"
username = 'sigma'
basedir = sprintf(dirpattern, username)
rdatafiles = list.files(path = str_c(basedir, 'Documents/GP/'), pattern = '*.RData', full.names = TRUE)
geopath = str_c(basedir, "OneDrive/Data/Korea/")
drive = str_c(basedir, "OneDrive/NCC_Project/CancerClustering/")
geopath = str_c(basedir, "OneDrive/Data/Korea/")
dbdir = drive  
rdsdir = drive #sprintf("/mnt/c/Users/%s/OneDrive/NCC_Project/CancerClustering/", username)



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


cofc3= covar_origin_10_fc %>% mutate(period = "2009-2013") %>% st_transform(5179) %>%
    dplyr::select(sgg_cd_c, period, matches("ragest_(i|d)_(Lung|Stomach)_*"), geom) %>%
    pivot_longer(cols = seq(3, ncol(.)-1)) %>%
    mutate(valuetype = ifelse(grepl("_i_", name), "Incidence", "Mortality"),
           sex = ifelse(grepl("_(female)_", name), "female", "male"),
           cancertype = ifelse(grepl("_Lung_", name), "Lung", "Stomach"),
           cancerlabel = paste(cancertype, " ", valuetype, sep = ""))
cofc2= covar_origin_05_fc %>% mutate(period = "2004-2008") %>% st_transform(5179) %>%
    dplyr::select(sgg_cd_c, period, matches("ragest_(i|d)_(Lung|Stomach)_*"), geom) %>%
    pivot_longer(cols = seq(3, ncol(.)-1)) %>%
    mutate(valuetype = ifelse(grepl("_i_", name), "Incidence", "Mortality"),
           sex = ifelse(grepl("_(female)_", name), "female", "male"),
           cancertype = ifelse(grepl("_Lung_", name), "Lung", "Stomach"),
           cancerlabel = paste(cancertype, " ", valuetype, sep = ""))
cofc1= covar_origin_00_fc %>% mutate(period = "1999-2003") %>% st_transform(5179) %>%
    dplyr::select(sgg_cd_c, period, matches("ragest_(i|d)_(Lung|Stomach)_*"), geom) %>%
    pivot_longer(cols = seq(3, ncol(.)-1)) %>%
    mutate(valuetype = ifelse(grepl("_i_", name), "Incidence", "Mortality"),
           sex = ifelse(grepl("_(female)_", name), "female", "male"),
           cancertype = ifelse(grepl("_Lung_", name), "Lung", "Stomach"),
           cancerlabel = paste(cancertype, " ", valuetype, sep = ""))



cofcs = bind_rows(
    cofc1, cofc2, cofc3
) %>%
    st_as_sf(sf_column_name = "geom") %>%
    mutate(cancerlabel = factor(cancerlabel,
        levels = c("Lung Incidence", "Lung Mortality", "Stomach Incidence", "Stomach Mortality")))




concat_tmaps = function(sfd, cn, title_custom = "Age standardized rate (Male)"){
    
    tmap_style("white")
        
    tm_shape(sfd) +
        tm_fill(cn, n = 7, style = 'jenks', pal = "PuRd",
                title = title_custom) +
        tm_borders(lwd = 0.05, col = 'grey') +
        tm_facets(c("period", "cancerlabel")) +
        tm_layout(#legend.position = "bottom",
        legend.outside = TRUE,
        legend.outside.position = "bottom", 
        outer.margins = c(-0.1, 0.01, 0.01, 0.01),
        panel.label.bg.color = 'white',
        panel.label.size = 1.5, panel.label.height = 2,
        legend.title.size = 1.6,
        legend.text.size = 1)
    
}


cofcsmap1 = concat_tmaps(cofcs %>% filter(sex == "male"), 'value')
cofcsmap2 = concat_tmaps(cofcs %>% filter(sex == "female"), 'value',
                title_custom = "Age standardized rate (Female)")

tmap_save(cofcsmap1, filename = str_c(drive, "/Manuscript/ASR_Map_Male.png"),
        width = 12, height = 15, units = 'in', dpi = 300, pointsize = 16)
tmap_save(cofcsmap2, filename = str_c(drive, "/Manuscript/ASR_Map_Female.png"),
        width = 12, height = 15, units = 'in', dpi = 300, pointsize = 16)
