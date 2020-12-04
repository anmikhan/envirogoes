library (data.table)
library (lubridate)

#Function for calculating solar zenith
#--------------------------------------------------------->
# This routine used to estimate the zenith angle (Degree)
# Modified version of code from Jehn-Yih Juang, Jul 2004
# ref: Campbell and Norman,1998 eq. 11.1 - 11.4
# input
# Hr: hour in decimal form, ex: 14:30-->14.5
# LAT: latitude, in degree
# LON: longitude, in degree
# DOY: day of year
# output
# ZA: zenith angel in degree
# cosZA: cosine of zenith angle
#--------------------------------------------------------->
  
zenith_calc <- function  (LAT,LON,DOY,Hr){
  A = (356.6 + 0.9856 * DOY) * pi/180
  A = (278.97 + 0.9856 * DOY + 1.9165 * sin(A)) * pi/180
  sinSD = 0.39785 * sin(A)
  SD = asin (sinSD)

  STM = seq(-165, 180, by = 15) #Standard meridians
  LCa = (LON-STM) / 15
  LC = min(abs(LCa))
  LC = LCa[which(abs(LCa) == LC)] #Longitudinal correction
  
  f = 279.575 + 0.9856 * DOY
  A = f * pi/180
  ET = (-104.7 * sin(A) 
        + 596.2 * sin(2 * A) 
        + 4.3 * sin(3 * A) 
        - 12.7 * sin(4 * A) 
        - 429.3 * cos(A) 
        - 2.0 * cos(2 * A) 
        + 19.3 * cos (3 * A)) / 3600
  t0 = 12 - LC - ET
  
  A = LAT * pi/180.;
  cosZA = (sin(A) * sin(SD) 
           + cos(A) * cos(SD) * cos (15 * (Hr - t0) * pi/180))
  
  ZA = acos(cosZA) * 180/pi
  return (c(ZA, cosZA))
  
}

#--------------------------------------------------------------------------------------------- 
#Read in data from goes
fires <- fread('/Users/badgrs/data/fire/califires2019_goes_cmi2km_mar2019_march2020_cmask.csv')

#convert stimestamp to utc POSIXct datatype
fires$stimestamp_utc = as.POSIXct(fires$stimestamp, 
                                  tz = 'UTC',
                                  format="%Y-%m-%dT%H:%M:%OSZ")

#convert to pacific
fires$stimestamp_ldt = format(fires$stimestamp_utc, tz = 'US/Pacific')
fires$stimestamp_ldt = as.POSIXct(fires$stimestamp_ldt, tz = 'US/Pacific')

#Get date as just day to calculate stats by day
fires$date_utc = date(fires$stimestamp_utc)
fires$date_ldt = date(fires$stimestamp_ldt)

#Calculate decimal hour
fires[, dec_hour_utc := (hour(stimestamp_utc) 
                         + (minute(stimestamp_utc) / 60) 
                         + (second(stimestamp_utc) / 3600)
                         )]
fires[, dec_hour_ldt := (hour(stimestamp_ldt) 
                         + (minute(stimestamp_ldt) / 60) 
                         + (second(stimestamp_ldt) / 3600)
                         )]

#Calculate solar zenith angle and cos(sza) for each unique fire timestamp
fires[, sza := zenith_calc(lat, lon, doy, dec_hour_ldt)[1], 
      by = list(fire, stimestamp, inside_perim)]
fires[, cos_sza := zenith_calc(lat, lon, doy, dec_hour_ldt)[2], 
      by = list(fire, stimestamp, inside_perim)]

#Clear Kincade fire observations
kin = fires[fire == 'Kincade' & cloud == 0]

#Calculate NDVI and NDII
kin[, NDVI := (CMI_C03 - CMI_C02) / (CMI_C03 + CMI_C02)]
kin[, NDII_5 := (CMI_C03 - CMI_C05) / (CMI_C03 + CMI_C05)]
kin[, NDII_6 := (CMI_C03 - CMI_C06) / (CMI_C03 + CMI_C06)]

kin_day = kin[sza < 70]

#Calculate mean daytime NDVI, NDII_5, NDII_6
in_mean_ndvi = kin_day[inside_perim == T, 
                      list(mean_NDVI = mean(NDVI),
                           mean_NDII5 = mean(NDII_5),
                           mean_NDII6 = mean(NDII_6),
                           stimestamp_ldt, 
                           stimestamp_utc), 
                      by = date_ldt]
out_mean_ndvi = kin_day[inside_perim == F, 
                        list(mean_NDVI = mean(NDVI), 
                             mean_NDII5 = mean(NDII_5),
                             mean_NDII6 = mean(NDII_6),
                             stimestamp_ldt, 
                             stimestamp_utc), 
                       by = date_ldt]


#Calculate median NDVI,  NDII_5, NDII_6 surrounding noon hours
in_mednoon_ndvi = kin_day[inside_perim == T 
                          & dec_hour_ldt > 10 
                          & dec_hour_ldt < 14, 
                      list(median_NDVI = median(NDVI), 
                           median_NDII5 = median(NDII_5),
                           median_NDII6 = median(NDII_6),
                           stimestamp_ldt, 
                           stimestamp_utc), 
                      by = date_ldt]
out_mednoon_ndvi = kin_day[inside_perim == F 
                           & dec_hour_ldt > 10 
                           & dec_hour_ldt < 14, 
                           list(median_NDVI = median(NDVI), 
                                median_NDII5 = median(NDII_5),
                                median_NDII6 = median(NDII_6),
                                stimestamp_ldt, 
                                stimestamp_utc),  
                       by = date_ldt]

#Read in modis data at fire point
fires.modis <- fread('/Users/badgrs/data/fire/califires2019_modis_l1b1km_mar2019_march2020.csv')

#Calculate NDVI
fires.modis [, NDVI := (band2 - band1) / (band2 + band1)]

#convert stimestamp to utc POSIXct datatype
fires.modis$timestamp_utc = as.POSIXct(fires.modis$timestamp, 
                                       tz = 'UTC',
                                       format="%Y-%m-%dT%H:%M:%OSZ")

#convert to pacific
fires.modis$timestamp_ldt = format(fires.modis$timestamp_utc, 
                                   tz = 'US/Pacific')
fires.modis$timestamp_ldt = as.POSIXct(fires.modis$timestamp_ldt, 
                                       tz = 'US/Pacific')

#Get date as just day 
fires.modis$date_utc = date(fires.modis$timestamp_utc)
fires.modis$date_ldt = date(fires.modis$timestamp_ldt)

#Calculate decimal hour
fires.modis[, dec_hour_utc := (hour(timestamp_utc) 
                               + (minute(timestamp_utc) / 60) 
                               + (second(timestamp_utc) / 3600))]
fires.modis[, dec_hour_ldt := (hour(timestamp_ldt) 
                               + (minute(timestamp_ldt) / 60) 
                               + (second(timestamp_ldt) / 3600))]

#Use cloud mask from goes data to cloud mask modis data
fires.modis$cloud = 1 #Add column to be populated

for (row in 1:nrow(fires.modis)){
  mod_hour = fires.modis[row,]$dec_hour_ldt #Get modis local hour
  mod_date_ldt = fires.modis[row,]$date_ldt #Get modis local date
  gfire = fires[fire == 'Kincade' & inside_perim == T & date_ldt == mod_date_ldt] #get goes fire rows for modis date and hour
  if (nrow(gfire) > 0){
    gfire [,goes_mod_hr_diff := dec_hour_ldt - mod_hour] #Calculate difference between modis hour and all goes hours for the day
    diff_index = which(abs(gfire$goes_mod_hr_diff) == min(abs(gfire$goes_mod_hr_diff)))
    if (is.na(gfire [diff_index,]$cloud)){
      fires.modis[row,]$cloud = NA
    }else{
      fires.modis[row,]$cloud = gfire [diff_index,]$cloud
    }
  }
}

#----------------------PLOT NDVI inside fire--------
par(mfrow = c(2, 1))

ylims = c(min(kin_day[inside_perim == T]$NDVI, 
              fires.modis[Name == 'kincade_infire']$NDVI, 
              na.rm = T), 
          max(kin_day[inside_perim == T]$NDVI, 
              fires.modis[Name == 'kincade_infire']$NDVI, 
              na.rm = T)
)

#Set up plot
plot(NDVI ~ stimestamp_ldt, 
     kin_day[inside_perim == T], 
     xlab = '', 
     xaxt = 'n',
     ylab = 'TOA NDVI',
     ylim = ylims,
     pch = '',
     cex.axis = 1)

#Shade fire dates
x = c(min(kin_day[inside_perim == T & doy == 296]$stimestamp_ldt), 
      min(kin_day[inside_perim == T & doy == 310]$stimestamp_ldt), 
      min(kin_day[inside_perim == T & doy == 310]$stimestamp_ldt), 
      min(kin_day[inside_perim == T & doy == 296]$stimestamp_ldt))
y = c(-2, -2, 2, 2)
polygon(x, y, col = 'gray88', border = 'gray88')

#Add NDVI
points(NDVI ~stimestamp_ldt,
       kin_day[inside_perim == T],
       pch = 20, cex = 0.5, col = 'gray50')

#Add modis NDVI
points(NDVI ~timestamp_ldt,
       fires.modis[Name == 'kincade_infire' & cloud == 0],
       pch = 20, cex = 0.7, col = 'red')

#Add GOES noon time median NDVI
points(median_NDVI ~ stimestamp_ldt,
       in_mednoon_ndvi,
       pch = 20, cex = 0.6, col = 'goldenrod2')

#Add GOES mean NDVI
points(mean_NDVI ~ stimestamp_ldt,
       in_mean_ndvi,
       #pch = 5, 
       cex = 0.5, col = 'purple')

#Add x-axis
date_seq = seq(
  min(kin_day[inside_perim == T]$stimestamp_ldt), 
  max(kin_day[inside_perim == T]$stimestamp_ldt), 
  by = 'week')

axis.POSIXct(1, 
             kin_day[inside_perim == T]$stimestamp_ldt, 
             at = date_seq, 
             format = '%m-%d', 
             las = 3,
             cex = 0.8)
title ("Inside fire perimeter", col = 'gray')


######-------Plot NDVI outside fire perimeter--------##########
ylims_out = c(min(kin_day[inside_perim == F]$NDVI, 
              fires.modis[Name == 'kincade_outfire']$NDVI, 
              na.rm = T), 
          max(kin_day[inside_perim == F]$NDVI, 
              fires.modis[Name == 'kincade_outfire']$NDVI, 
              na.rm = T)
)

plot(NDVI ~ stimestamp_ldt, 
     kin_day[inside_perim == F & NDVI > -1 & NDVI < 1], 
     xlab = '', 
     ylab = 'TOA NDVI',
     ylim = ylims_out,
     pch = '',
     xaxt = 'n')

#Shade fire dates
x = c(min(kin_day[inside_perim == F & doy == 296]$stimestamp_ldt), 
      min(kin_day[inside_perim == F & doy == 310]$stimestamp_ldt), 
      min(kin_day[inside_perim == F & doy == 310]$stimestamp_ldt), 
      min(kin_day[inside_perim == F & doy == 296]$stimestamp_ldt))
y = c(-2, -2, 2, 2)
polygon(x, y, col = 'gray88', border = 'gray88')

#Add NDVI
points(NDVI ~stimestamp_ldt,
       kin_day[inside_perim == F],
       pch = 20, cex = 0.5, col = 'gray50')

#Add modis NDVI
points(NDVI ~timestamp_ldt,
       fires.modis[Name == 'kincade_outfire' & cloud == 0],
       pch = 20, cex = 0.7, col = 'red')

#Add median noon NDVI
points(median_NDVI ~ stimestamp_ldt,
       out_mednoon_ndvi,
       pch = 20, cex = 0.6, col = 'goldenrod2')

#Add daily mean NDVI
points(mean_NDVI ~ stimestamp_ldt,
       out_mean_ndvi,
       cex = 0.5, col = 'purple')

#Add x-axis
date_seq = seq(
  min(kin_day[inside_perim == F]$stimestamp_ldt), 
  max(kin_day[inside_perim == F]$stimestamp_ldt), 
  by = 'week')

axis.POSIXct(1, 
             kin_day[inside_perim == F]$stimestamp_ldt, 
             at = date_seq, 
             format = '%m-%d', 
             las = 3)

title ('Outside fire perimeter', col = 'gray')

#Add legend
legend('bottomright',
       legend=c("MODIS TOA NDVI",
                "GOES TOA NDVI",
                "GOES median (local hours 10 - 14)", 
                "GOES daily daytime mean",
                "Fire duration"),
       pch = c(20, 20, 20, 1, NA),
       col=c("red", "gray50", "goldenrod2", "purple", NA),
       border = c(NA, NA, NA, NA, NA),
       fill = c(NA, NA, NA, NA, 'gray'),
       ncol=1, 
       bty="n")


#######-------Plot NDII inside fire perimeter--------##########
par(mfrow = c(2, 1))
ylims_out = c(min(kin_day$NDII_5, 
                  na.rm = T), 
              max(kin_day$NDII_5, 
                  na.rm = T)
)

plot(NDII_5 ~ stimestamp_ldt, 
     kin_day[inside_perim == T & NDII_5 > -1 & NDII_5 < 1], 
     xlab = '', 
     ylab = 'TOA NDII',
     ylim = ylims_out,
     pch = '',
     xaxt = 'n')

#Add fire dates
x = c(min(kin_day[inside_perim == T & doy == 296]$stimestamp_ldt), 
      min(kin_day[inside_perim == T & doy == 310]$stimestamp_ldt), 
      min(kin_day[inside_perim == T & doy == 310]$stimestamp_ldt), 
      min(kin_day[inside_perim == T & doy == 296]$stimestamp_ldt))
y = c(-2, -2, 2, 2)
polygon(x, y, col = 'gray88', border = 'gray88')

#Add NDII
points(NDII_5 ~ stimestamp_ldt,
       kin_day[inside_perim == T],
       pch = 20, cex = 0.5, col = 'gray50')

#Add median noon NDII
points(median_NDII5 ~ stimestamp_ldt,
       in_mednoon_ndvi,
       pch = 20, cex = 0.6, col = 'goldenrod2')

#Add daily mean NDII
points(mean_NDII5 ~ stimestamp_ldt,
       in_mean_ndvi,
       cex = 0.5, col = 'purple')

#Add x-axis
date_seq = seq(
  min(kin_day[inside_perim == T]$stimestamp_ldt), 
  max(kin_day[inside_perim == T]$stimestamp_ldt), 
  by = 'week')

axis.POSIXct(1, kin_day[inside_perim == T]$stimestamp_ldt, at = date_seq, format = '%m-%d', las = 3)
title ('Inside fire perimeter', col = 'gray')

#Add legend
legend('topright',
       legend=c("GOES TOA NDII",
                "GOES median (local hours 10 - 14)", 
                "GOES daily daytime mean",
                "Fire duration"),
       pch = c(20, 20, 1, NA),
       col=c("gray50", "goldenrod2", "purple", NA),
       border = c(NA, NA, NA, NA),
       fill = c(NA, NA, NA, 'gray'),
       ncol=1, 
       #xpd=NA, 
       bty="n")

#######-------Plot NDII outside of fire perimeter--------##########
ylims_out = c(min(kin_day$NDII_5, 
                  na.rm = T), 
              max(kin_day$NDII_5, 
                  na.rm = T)
)

plot(NDII_5 ~ stimestamp_ldt, 
     kin_day[inside_perim == F & NDII_5 > -1 & NDII_5 < 1], 
     xlab = '', 
     ylab = 'TOA NDII',
     ylim = ylims_out,
     pch = '',
     xaxt = 'n')

#Shade fire dates
x = c(min(kin_day[inside_perim == F & doy == 296]$stimestamp_ldt), 
      min(kin_day[inside_perim == F & doy == 310]$stimestamp_ldt), 
      min(kin_day[inside_perim == F & doy == 310]$stimestamp_ldt), 
      min(kin_day[inside_perim == F & doy == 296]$stimestamp_ldt))
y = c(-2, -2, 2, 2)
polygon(x, y, col = 'gray88', border = 'gray88')

#Add NDII
points(NDII_5 ~ stimestamp_ldt,
       kin_day[inside_perim == F],
       pch = 20, cex = 0.5, col = 'gray50')

#Add median noon NDII
points(median_NDII5 ~ stimestamp_ldt,
       out_mednoon_ndvi,
       pch = 20, cex = 0.6, col = 'goldenrod2')

#Add daily mean NDII
points(mean_NDII5 ~ stimestamp_ldt,
       out_mean_ndvi,
       #pch = 4, 
       cex = 0.5, col = 'purple')

#Add x-axis
date_seq = seq(
  min(kin_day[inside_perim == F]$stimestamp_ldt), 
  max(kin_day[inside_perim == F]$stimestamp_ldt), 
  by = 'week')

axis.POSIXct(1, 
             kin_day[inside_perim == F]$stimestamp_ldt, 
             at = date_seq, 
             format = '%m-%d', 
             las = 3)

title ('Outside fire perimeter', col = 'gray')