library (rgdal)
library (raster)
library (sp)
library (viridis)

#Read in NDVI rasters
ndvi_295 = raster ('/Users/badgrs/data/califire_maps/goes_imagery/OR_ABI-L2-MCMIPC-M6_G16_s20192952101323_e20192952104097_c20192952104259.nc',
                   varname = 'NDVI')
ndvi_312 = raster ('/Users/badgrs/data/califire_maps/goes_imagery/OR_ABI-L2-MCMIPC-M6_G16_s20193122141211_e20193122143583_c20193122144141.nc',
                   varname = 'NDVI')
#Set extent and crs
extent(ndvi_295) <- c(-3398813.3659705431200564, 
                      -3358732.9255137387663126, 
                      3621259.3909196606837213, 
                      3663343.7944607832469046)

crs(ndvi_295) <- '+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x'

extent(ndvi_312) <- c(-3398813.3659705431200564, 
                      -3358732.9255137387663126, 
                      3621259.3909196606837213,
                      3663343.7944607832469046)
crs(ndvi_312) <- '+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x'

#NDVI difference
ndvi_diff = ndvi_295 - ndvi_312

#Read in Kincade fire perimiter
kincade = readOGR('/Users/badgrs/data/fire', 'kincade')
kin_geos = spTransform(kincade, '+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x')

#Read in Cali borders and raster extent
states = readOGR('/Users/badgrs/data/tiger/tl_2019_us_state', 
                 'tl_2019_us_state')
raster_ext = readOGR('/Users/badgrs/data/califire_maps',
                     'raster_extent')

cali = states[states$STUSPS == 'CA',]

#Read in points of time series extract
kincade_in = readOGR('/Users/badgrs/data/califire_maps', 'kincade_point_infire')
kincade_in = spTransform(kincade_in, '+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x')

kincade_out = readOGR('/Users/badgrs/data/califire_maps', 'kincade_point_outfire')
kincade_out = spTransform(kincade_out, '+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x')


#Make NA rasters to plot under actual rasters to dsiplay NA as gray
ndvi_295_na <- ndvi_295
ndvi_295_na [is.na(ndvi_295_na)] <- 999

ndvi_312_na <- ndvi_312
ndvi_312_na [is.na(ndvi_312_na)] <- 999

ndvi_diff_na <- ndvi_diff
ndvi_diff_na [is.na(ndvi_diff_na)] <- 999

###Make rasters with same range to make same legends####
#Set breaks and labels according to the maximum minimum value between two days
#and the minimum maximum values between the two days
ndvi_brks = seq(max(min(values(ndvi_295), na.rm = T), min(values(ndvi_312), na.rm = T)),
                min(max(values(ndvi_295), na.rm = T), max(values(ndvi_312), na.rm = T)),
                length.out = 255)
ndvi_at = c(seq(range(ndvi_brks)[1],
                range(ndvi_brks)[2],
                length.out = 4))
ndvi_295_labels = c(formatC(ndvi_at[1:3], digits = 2),
                    paste('>=', formatC(ndvi_at[4], digits = 2))
)

ndvi_312_labels = c(paste('<=', formatC(ndvi_at[1], digits = 2)),
                    formatC(ndvi_at[2:3], digits = 2),
                    paste('>=', formatC(ndvi_at[4], digits = 2)))

#Adjust raster ranges
ndvi_295_plot = ndvi_295
ndvi_295_plot[ndvi_295_plot >= range(ndvi_brks)[2]] = range(ndvi_brks)[2]

ndvi_312_plot = ndvi_312
ndvi_312_plot[ndvi_312_plot >= range(ndvi_brks)[2]] = range(ndvi_brks)[2]
ndvi_312_plot[ndvi_312_plot <= range(ndvi_brks)[1]] = range(ndvi_brks)[1]




#--------------PLOT---------------------------
#ndvi_brks = seq(0, 1, by = 0.1)

par(mfrow = c(1, 4),
    mar = c(5, 5, 4, 4), #c(bottom, left, top, right)
    cex.axis = 1.3, 
    cex.main = 1.2)
plot(ndvi_295_na, 
     col = 'gray', 
     legend = F)
plot(ndvi_295_plot, 
     col = viridis(256), 
     legend = T,
     breaks = ndvi_brks,
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T,
     axis.args = list (at = ndvi_at,
                       labels = ndvi_295_labels,
                       cex.axis = 1.2))
plot(kin_geos, border = 'gray', lwd = 2, add = T)
plot(kincade_in[1,], col = 'red', pch = 20, cex = 1.5, add = T)
plot(kincade_out, col = 'red', pch = 8, cex = 1.5, add = T)
title ('A', adj = 0)

plot(ndvi_312_na, 
     col = 'gray', 
     legend = F)
plot(ndvi_312_plot, 
     col = viridis(256), 
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T,
     breaks = ndvi_brks,
     axis.args = list (at = ndvi_at,
                       labels = ndvi_312_labels,
                       cex.axis = 1.2))
plot(kin_geos, border = 'gray', lwd = 2, add = T)
plot(kincade_in[1,], col = 'red', pch = 20, cex = 1.5, add = T)
plot(kincade_out, col = 'red', pch = 8, cex = 1.5, add = T)
title ('B', adj = 0)

plot(ndvi_diff_na, 
     col = 'gray', 
     legend = F)
plot(ndvi_diff, 
     col = inferno(256), 
     legend.width = 2, 
     legend.shrink = 1,
     cex = 1.2,
     add = T)
plot(kin_geos, border = 'gray', lwd = 2, add = T)
plot(kincade_in[1,], col = 'red', pch = 20, add = T, cex = 1.5)
plot(kincade_out, col = 'red', pch = 8, add = T, cex = 1.5)
title ('C', adj = 0)

plot(cali, col = 'darkgray', border = 'black', axes = T)
plot(raster_ext, border = 'red', lwd = 2, add = T)
title ('D', adj = 0)

legend ('topright',
        legend = c('Fire perimeter',
                   'Raster extent',
                   'Inside perimeter',
                   'Outside perimeter'),
        col = c(NA, NA, 'red', 'red'),
        pch = c(NA, NA, 20, 8),
        fill = c(NA, NA, NA, NA),
        border = c('gray', 'red', NA, NA),
        cex = 1.4,
        bty = 'n'
        )


