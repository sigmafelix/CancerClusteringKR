source('./Code/Base/base_functions.R')
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
           sex = ifelse(grepl("_(female)_", name), "Female", "Male"),
           cancertype = ifelse(grepl("_Lung_", name), "Lung", "Stomach"),
           cancerlabel = paste(cancertype, " ", valuetype, sep = ""))
cofc2= covar_origin_05_fc %>% mutate(period = "2004-2008") %>% st_transform(5179) %>%
    dplyr::select(sgg_cd_c, period, matches("ragest_(i|d)_(Lung|Stomach)_*"), geom) %>%
    pivot_longer(cols = seq(3, ncol(.)-1)) %>%
    mutate(valuetype = ifelse(grepl("_i_", name), "Incidence", "Mortality"),
           sex = ifelse(grepl("_(female)_", name), "Female", "Male"),
           cancertype = ifelse(grepl("_Lung_", name), "Lung", "Stomach"),
           cancerlabel = paste(cancertype, " ", valuetype, sep = ""))
cofc1= covar_origin_00_fc %>% mutate(period = "1999-2003") %>% st_transform(5179) %>%
    dplyr::select(sgg_cd_c, period, matches("ragest_(i|d)_(Lung|Stomach)_*"), geom) %>%
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

        tm_shape(k) +
            tm_fill(cn, n = 5, style = 'quantile', pal = "PuRd",
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

tmap_save(cofcsmap1, filename = str_c(drive, "/Manuscript/ASR_Map_Male2.png"),
        width = 12, height = 15, units = 'in', dpi = 300, pointsize = 16)
tmap_save(cofcsmap2, filename = str_c(drive, "/Manuscript/ASR_Map_Female2.png"),
        width = 12, height = 15, units = 'in', dpi = 300, pointsize = 16)

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
covarlist00 = colnames(covar_origin_00_fc)[c(4:5,11:12,13:22)]
covarlist05 = colnames(covar_origin_05_fc)[c(4:5,11:12,15:17, 19:29)]
covarlist10 = colnames(covar_origin_10_fc)[c(4:5,11:12,15,16,17,19,20,21,22,23,24,25,26,27,28,29)]
par(mfrow = c(5,4))

export_mm = function(insf, covarlist, filename) {
    #oldpar = par()
    png(filename, 15, 15, "in", 16, res = 508)
    par(mfrow = c(5,4))
    for (ce in seq_len(length(covarlist))) {
        mf_init(insf)
        mf_map(insf[,covarlist[ce]], covarlist[ce], type = "choro", nbreaks = 5, border = 'grey', lwd = 0.05, add = TRUE)
        # mf_title(ce)
    }
    dev.off()
    par(new=TRUE)
}

export_mm(covar_origin_00_fc, covarlist00, str_c(drive,"/Manuscript/Covar_dist_00.png"))
export_mm(covar_origin_05_fc, covarlist05, str_c(drive,"/Manuscript/Covar_dist_05.png"))
export_mm(covar_origin_10_fc, covarlist10, str_c(drive,"/Manuscript/Covar_dist_10.png"))



    png(str_c(drive,"/Manuscript/Covar_dist_00.png"), width=15, height=15, units="in", pointsize=16, res = 508)
    par(mfrow = c(5,4))
    for (ce in seq_len(length(covarlist00))) {
        mf_init(covar_origin_00_fc)
        mf_map(covar_origin_00_fc[,covarlist00[ce]], covarlist00[ce], type = "choro", nbreaks = 5, border = 'grey', lwd = 0.05, add = TRUE)
        # mf_title(ce)
    }
    dev.off()
    par(new=TRUE)
