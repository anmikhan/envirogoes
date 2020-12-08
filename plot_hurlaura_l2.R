library (rgdal)
library (raster)
library (sp)
library (viridis)
library (ncdf4)
library (RColorBrewer)

#Read in shapefiles
extent = readOGR('/Users/badgrs/data/hurr_laura', 'small_ext')
track = readOGR('/Users/badgrs/data/hurr_laura', 'advisory_forecast_track')
states = readOGR('/Users/badgrs/data/hurr_laura', 'states')

#Find rasters
rast_dir = '/Users/badgrs/data/goes16/hur_laura'
rast_files = list.files (rast_dir)

#Hurricane: doy 239 237 242 244 245 246
pre_file = file.path (rast_dir, grep ('*s20202351601.*', rast_files, value = T))
post_file = file.path (rast_dir, grep ('*s202024616.*', rast_files, value = T)[1])

#-----Load veg indices pre and post-------
pre_ndvi = raster (pre_file, varname = 'NDVI')
post_ndvi = raster (post_file, varname = 'NDVI')

pre_nirv = raster (pre_file, varname = 'NIRV')
post_nirv = raster (post_file, varname = 'NIRV')

pre_ndii = raster (pre_file, varname = 'NDII')
post_ndii = raster (post_file, varname = 'NDII')

#--------Make EVI raster--------------
pre_blue = raster (pre_file, varname = 'CMI_C01')
post_blue = raster (post_file, varname = 'CMI_C01')
pre_red = raster (pre_file, varname = 'CMI_C02')
post_red = raster (post_file, varname = 'CMI_C02')
pre_nir = raster (pre_file, varname = 'CMI_C03')
post_nir = raster (post_file, varname = 'CMI_C03')

#Calculate pre and post EVI
pre_evi = 2.5 * ((pre_nir - pre_red) / 
                   (pre_nir + (6 * pre_red) - (7.5 * pre_blue) + 1)) 
pre_evi [pre_evi < 0] = NA
pre_evi [pre_evi > 1] = NA

post_evi = 2.5 * ((post_nir - post_red) / 
                    (post_nir + (6 * post_red) - (7.5 * post_blue) + 1))  
post_evi [post_evi < 0] = NA
post_evi [post_evi > 1] = NA

#Set extent and crs for all rasters
rast_ext = c(-1941892.833,
             -370743.263,
             2525062.054,
             3535086.510
             )
goes_crs = '+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x'

extent(pre_ndvi) = rast_ext #ndvi
crs(pre_ndvi) = goes_crs
pre_ndvi = crop (pre_ndvi, extent(extent))

extent(post_ndvi) = rast_ext #ndvi
crs(post_ndvi) = goes_crs
post_ndvi = crop (post_ndvi, extent(extent))

extent(pre_nirv) = rast_ext #nirv
crs(pre_nirv) = goes_crs
pre_nirv = crop (pre_nirv, extent(extent))

extent(post_nirv) = rast_ext #nirv
crs(post_nirv) = goes_crs
post_nirv = crop (post_nirv, extent(extent))

extent(pre_evi) = rast_ext #ndii
crs(pre_evi) = goes_crs
pre_evi = crop (pre_evi, extent(extent))

extent(post_evi) = rast_ext #ndii
crs(post_evi) = goes_crs
post_evi = crop (post_evi, extent(extent))

extent(pre_ndii) = rast_ext #ndii
crs(pre_ndii) = goes_crs
pre_ndii = crop (pre_ndii, extent(extent))

extent(post_ndii) = rast_ext #ndii
crs(post_ndii) = goes_crs
post_ndii = crop (post_ndii, extent(extent))

#Mark NIRv rasters as Na where NDVI is less than 0 to get NIR reflectance of veg
pre_nirv [pre_ndvi < 0] = NA
post_nirv [post_ndvi < 0] = NA

#Mark evi rasters as NA where NDVI is NA to cloud mask
pre_evi [is.na(pre_ndvi)] = NA
post_evi [is.na(post_ndvi)] = NA

#Get differences for each veg index
diff_ndvi = pre_ndvi - post_ndvi
diff_nirv = pre_nirv - post_nirv
diff_evi = pre_evi - post_evi
diff_ndii = pre_ndii - post_ndii

#-----Make NA rasters to plot under actual rasters to dsiplay NA as gray----
pre_ndvi_na <- pre_ndvi
pre_ndvi_na [is.na(pre_ndvi_na)] <- 999

post_ndvi_na <- post_ndvi
post_ndvi_na [is.na(post_ndvi_na)] <- 999

diff_ndvi_na <- diff_ndvi
diff_ndvi_na [is.na(diff_ndvi_na)] <- 999
#nirv
pre_nirv_na <- pre_nirv
pre_nirv_na [is.na(pre_nirv_na)] <- 999

post_nirv_na <- post_nirv
post_nirv_na [is.na(post_nirv_na)] <- 999

diff_nirv_na <- diff_nirv
diff_nirv_na [is.na(diff_nirv_na)] <- 999
#evi
pre_evi_na <- pre_evi
pre_evi_na [is.na(pre_evi_na)] <- 999

post_evi_na <- post_evi
post_evi_na [is.na(post_evi_na)] <- 999

diff_evi_na <- diff_evi
diff_evi_na [is.na(diff_evi_na)] <- 999
#ndii
pre_ndii_na <- pre_ndii
pre_ndii_na [is.na(pre_ndii_na)] <- 999

post_ndii_na <- post_ndii
post_ndii_na [is.na(post_ndii_na)] <- 999

diff_ndii_na <- diff_ndii
diff_ndii_na [is.na(diff_ndii_na)] <- 999



#-----Plot----
par(mfrow = c(4, 3),
    mar = c(2, 3, 4, 10), #c(bottom, left, top, right)
    cex.axis = 1.3, 
    cex.main = 1.2)

pal <- colorRampPalette(c("darkgreen", "yellow", "purple", "black"))

#pre NDVI
plot(pre_ndvi_na, 
     col = 'gray', 
     legend = F)
plot(pre_ndvi, 
     col = viridis(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T)
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

title ('DOY 235')
title ('TOA NDVI', adj = 0)

#post NDVI
plot(post_ndvi_na, 
     col = 'gray', 
     legend = F)
plot(post_ndvi, 
     col = viridis(256), 
     legend = T,
     #breaks = ndvi_brks,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T)
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

title ('DOY 246')

#create breaks for difference
diff_ndvi_brks = seq(-min(max(values(diff_ndvi), na.rm = T), abs(min(values(diff_ndvi), na.rm = T))),
                     min(max(values(diff_ndvi), na.rm = T), abs(min(values(diff_ndvi), na.rm = T))),
                     length.out = 255)
diff_ndvi_at = c(seq (range(diff_ndvi_brks)[1], 
                      0, 
                      length.out = 4),
                 seq (0,
                      range(diff_ndvi_brks)[2],
                      length.out = 4))
#Take out extra zeros
diff_ndvi_at = unique(diff_ndvi_at)

diff_ndvi_labels = c(formatC(diff_ndvi_at[1:3], digits = 2), 
                     '0',
                     formatC(diff_ndvi_at[5:6], digits = 2),
                     paste('>=', 
                           formatC(diff_ndvi_at[7], digits = 2),
                           sep = ' ')
                     )
#Make new difference raster where all values above or equal to range(diff_ndvi_brks)[2] and equal to range(diff_ndvi_brks)[2]
diff_ndvi_plot = diff_ndvi
diff_ndvi_plot [diff_ndvi_plot > range(diff_ndvi_brks)[2]] = range(diff_ndvi_brks)[2]

# plot difference NDVI
plot(diff_ndvi_na, 
     col = 'gray', 
     legend = F)
plot(diff_ndvi_plot, 
     col = pal(256), 
     legend = T,
     breaks = diff_ndvi_brks,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T,
     axis.args = list (at = diff_ndvi_at,
                       labels = diff_ndvi_labels,
                       cex.axis = 1.2))
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)
title ('DOY 235 - DOY 246')

#pre NIRv
plot(pre_nirv_na, 
     col = 'gray', 
     legend = F)
plot(pre_nirv, 
     col = viridis(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T,
     xlab = 'TOA NIRv')
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

title ('TOA NIRv', adj = 0)

#post NIRv
plot(post_nirv_na, 
     col = 'gray', 
     legend = F)
plot(post_nirv, 
     col = viridis(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T)
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

#create breaks
diff_nirv_brks = seq(-min(max(values(diff_nirv), na.rm = T), 
                          abs(min(values(diff_nirv), na.rm = T))),
                     min(max(values(diff_nirv), na.rm = T), 
                         abs(min(values(diff_nirv), na.rm = T))),
                     length.out = 255)
diff_nirv_at = c(seq (range(diff_nirv_brks)[1], 
                      0, 
                      length.out = 4),
                 seq (0,
                      range(diff_nirv_brks)[2],
                      length.out = 4))
#Take out extra zeros
diff_nirv_at = unique(diff_nirv_at)

diff_nirv_labels = c(formatC(diff_nirv_at[1:3], digits = 1), 
                     '0',
                     formatC(diff_nirv_at[5:6], digits = 1),
                     paste('>=', 
                           formatC(diff_nirv_at[7], digits = 1), 
                           sep = ' '))
#Make new difference raster where all values above or equal to range(diff_ndvi_brks)[2] are equal to range(diff_ndvi_brks)[2]
diff_nirv_plot = diff_nirv
diff_nirv_plot [diff_nirv_plot > range(diff_nirv_brks)[2]] = range(diff_nirv_brks)[2]

#difference NIRv
plot(diff_nirv_na, 
     col = 'gray', 
     legend = F)
plot(diff_nirv_plot, 
     col = pal(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 0.2,
     add = T,
     axis.args = list (at = diff_nirv_at,
                       labels = diff_nirv_labels,
                       cex.axis = 1.2))
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

#-----Plot EVI and NDII---------------------
pal <- colorRampPalette(c("darkgreen", "yellow", "purple", "black"))
#pre EVI
plot(pre_evi_na, 
     col = 'gray', 
     legend = F)
plot(pre_evi, 
     col = inferno(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T,
     xlab = 'TOA EVI')
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)


title ('TOA EVI', adj = 0)

#post evi
plot(post_evi_na, 
     col = 'gray', 
     legend = F)
plot(post_evi, 
     col = inferno(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T)
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

#create breaks
diff_evi_brks = seq(-min(max(values(diff_evi), na.rm = T), 
                          abs(min(values(diff_evi), na.rm = T))),
                     min(max(values(diff_evi), na.rm = T), 
                         abs(min(values(diff_evi), na.rm = T))),
                     length.out = 255)
diff_evi_at = c(seq (range(diff_evi_brks)[1], 
                      0, 
                      length.out = 4),
                 seq (0,
                      range(diff_evi_brks)[2],
                      length.out = 4))
#Take out extra zeros
diff_evi_at = unique(diff_evi_at)

diff_evi_labels = c(formatC(diff_evi_at[1:3], digits = 2), 
                     '0',
                     formatC(diff_evi_at[5:6], digits = 2),
                     paste('>=', 
                           formatC(diff_evi_at[7], digits = 2), 
                           sep = ' '))
#Make new difference raster where all values above or equal to range(diff_evi_brks)[2] are equal to range(diff_evi_brks)[2]
diff_evi_plot = diff_evi
diff_evi_plot [diff_evi_plot > range(diff_evi_brks)[2]] = range(diff_evi_brks)[2]

#difference evi
plot(diff_evi_na, 
     col = 'gray', 
     legend = F)
plot(diff_evi_plot, 
     col = pal(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T,
     axis.args = list (at = diff_evi_at,
                       labels = diff_evi_labels,
                       cex.axis = 1.2))
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

#pre NDII
plot(pre_ndii_na, 
     col = 'gray', 
     legend = F)
plot(pre_ndii, 
     col = viridis(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T,
     xlab = 'TOA NDII',
     ext = extent(extent))
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

title ('TOA NDII', adj = 0)

#post NDII
plot(post_ndii_na, 
     col = 'gray', 
     legend = F)
plot(post_ndii, 
     col = viridis(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T)
plot(states, add = T)
plot(track, add = T, col = 'red', lwd = 2)

#create breaks
diff_ndii_brks = seq(-min(max(values(diff_ndii), na.rm = T), 
                         abs(min(values(diff_ndii), na.rm = T))),
                    min(max(values(diff_ndii), na.rm = T), 
                        abs(min(values(diff_ndii), na.rm = T))),
                    length.out = 255)
diff_ndii_at = c(seq (range(diff_ndii_brks)[1], 
                     0, 
                     length.out = 4),
                seq (0,
                     range(diff_ndii_brks)[2],
                     length.out = 4))

#Take out extra zeros
diff_ndii_at = unique(diff_ndii_at)

diff_ndii_labels = c(formatC(diff_ndii_at[1:3], digits = 2), 
                    '0',
                    formatC(diff_ndii_at[5:6], digits = 2),
                    paste('>=', 
                          formatC(diff_ndii_at[7], digits = 2), 
                          sep = ' '))
#Make new difference raster where all values above or equal to range(diff_ndii_brks)[2] and equal to range(diff_ndii_brks)[2]
diff_ndii_plot = diff_ndii
diff_ndii_plot [diff_ndii_plot > range(diff_ndii_brks)[2]] = range(diff_ndii_brks)[2]

#difference NDII
plot(diff_ndii_na, 
     col = 'gray', 
     legend = F)
plot(diff_ndii_plot, 
     col = pal(256), 
     legend = T,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T,
     axis.args = list (at = diff_ndii_at,
                       labels = diff_ndii_labels,
                       cex.axis = 1.2))
plot(states, add = T, border = 'white')
plot(track, add = T, col = 'red', lwd = 2)