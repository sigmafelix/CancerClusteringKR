## Timeline: light review

##
library(pacman)
p_load(tidyverse, readxl, patchwork, gridExtra)
rev_tab = read_excel('/mnt/c/Users/isong/OneDrive/NCC_Project/CancerClustering/Cancer_Clustering_Results_PubMed_Re_032222.xlsx', sheet = 2)

rev_tab_l = rev_tab %>%
	dplyr::select(PMID, `Create Date`, SaTScan_C:Others) %>%
	pivot_longer(cols = SaTScan_C:Others) %>%
	rename(Cluster_Method = name,
		   occurrence = value) %>%
	filter(!is.na(occurrence)) %>%
    mutate(Cluster_Method = factor(Cluster_Method, 
                                levels = rev(c('Knox','GlobalMoran','LocalMoran','SaTScan_C','SaTScan_E','FleXScan','Others'))))


rev_tab_s = rev_tab %>%
	dplyr::select(PMID, `Create Date`, SaTScan_C:Others) %>%
	pivot_longer(cols = SaTScan_C:Others) %>%
	rename(Cluster_Method = name,
		   occurrence = value) %>%
    mutate(Cluster_Method = factor(Cluster_Method, 
                                levels = c('Knox','GlobalMoran','LocalMoran','SaTScan_C','SaTScan_E','FleXScan','Others'))) %>%
	filter(!is.na(occurrence)) %>%
    group_by(Cluster_Method) %>%
    summarize(N = n()) %>%
    ungroup
 


rev_plot = rev_tab_l %>%
	ggplot(data = .,
		   mapping = aes(x = `Create Date`, y = Cluster_Method, color = Cluster_Method)) +
		geom_point(size = 3, pch = 19, alpha = 0.5) +
    ylab('Clustering method') +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 13, face = 'bold'),
          legend.position = 'none') 

rev_plot + tableGrob(rev_tab_s) + plot_layout(ncol = 2, widths = c(4, 1))
