source("0_functions_plot.R")
res.data <- read.csv("complete_data.csv",check.names=FALSE)
res.data <- filter(res.data,
                   Study_midyear>=1990,
                   Meas_method %in% c('IRGA','Gas Chromatography') 
                   &Ecosystem_type != 'Agriculture'
                   & Outlier == 0)

dv = c('Rs_annual')
idv = c('Study_midyear','Latitude','Biome',
        'MAT','MAP', 'delta.Mean.Temp.','MAP','delta.Mean.Precip.', 'SOC_stock',
        'Meas_method','Ecosystem_type2','Partition_method','Stage','Elevation')

work.data = na.omit(res.data[,c(idv,dv)])
work.data <- droplevels(work.data)

work.data$Year = as.character(lapply(work.data$Study_midyear,labelYear1)) %>% as.factor()
work.data$Year2 = as.character(lapply(work.data$Study_midyear,labelYear2)) %>% as.factor()
boreal = filter(work.data, Biome == "boreal"|Biome == "Arctic")
Temperate = filter(work.data, Biome == "Temperate"|Biome == "Mediterranean" )
Tropical = filter(work.data, Biome == "Tropical" | Biome == "Subtropical")
Bef2008 = filter(work.data, Year2 == "1990-2010" )
Post2008 = filter(work.data, Year2 == "2010-2022" )
