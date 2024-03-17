

library(MODIStsp)

MODIStsp(
  gui             = FALSE, 
  out_folder      = "/home/felix/Documents/MODIS", 
  selprod         = "Vegetation Indexes_16Days_250m (M*D13Q1)",
  bandsel         = c("NDVI"), 
  prod_version = "061",
  quality_bandsel = NULL, 
  indexes_bandsel = NULL, 
  user            = "sigmafelix" ,
  password        = PASS,
  sensor = "Terra",
  start_x = 27,
  end_x = 28,
  start_y = 5,
  end_y = 5,
  start_date      = "2014.01.01", 
  end_date        = "2014.12.31", 
  verbose         = FALSE,
  parallel        = FALSE
)


MODIStsp(gui = TRUE)
