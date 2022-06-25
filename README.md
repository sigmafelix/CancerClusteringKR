# CancerClusteringKR
Spatial clustering analysis of cancer incidence and mortality in South Korea with spatial scan statistic

## Major application
- Using the standalone command line software SaTScanBatch64, I batched the purely spatial cluster analysis by automatically generating SaTScan setting files (*.prm)
- Incidence and mortality of Lung (male only) and Stomach (male and female) cancers in South Korea 
- Applying weighted and plain normal spatial scan statistic

## Change log
- The initial version relied on the R package [`smerc`](https://cran.r-project.org/web/packages/smerc/index.html)
- Due to the difference in the processing speed, I shifted from R package to standalone [SaTScan software](https://www.satscan.org)
- The directory structure was revised
- Data will be shared around the time when the manuscript will be submitted
