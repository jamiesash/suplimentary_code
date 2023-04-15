source("code//functions.R")
source("code//libraries.R")

# ------------------------------------------------------------------------------
# loads the summary table created in stat_table.R file and makes some fiugures 
# from the information. 

sumtbl = read.csv("data/summary_statistics/full_sum_20230218.csv")

# making a pretty table --------------------------------------------------------

tbl = sumtbl
rownames(tbl) = NULL
#idx = !year(tbl$sdate) == 2008
#tbl = tbl[idx,]
kbl(tbl, 
    booktabs = T, 
    latex_options = c("striped"),
    caption = "Summary statistics of each bloom around St. ALOHA and 30N between 2002 to 2019. Magnitude taken as the maximum CHL $mg/m^3$ value the blooms reached, and area is the maximum area the bloom reached.") %>%
  kable_styling(full_width = F, latex_options = "HOLD_position") %>%
  column_spec(1:2, bold = F, color = "black") %>%
  column_spec(1:7, width = "10em")
# row_spec(c(4, 5, 8, 11, 12, 17), bold = F, color = "black", background = "#add8e6")
# add_header_above(c("Table of (sub)mesoscale contributions: St. ALOHA" = 5), bold=T)

# ------------------------------------------------------------------------------
# Big Bloom Magnitudes

png(filename = paste("figures\\bloom_magnitudes_", Sys.Date(), ".png",sep = ""),
    width = 3.5,
    height = 6,
    units = "in",
    res = 300,
    pointsize = 10)

tbl = sumtbl
tbl[y == 2010, ]$mag = 0.142
x = tbl$mag
y = year(tbl$sdate)

main = ""
ylab = TeX(r'(CHL $[mg/mg^3]$)')
xlab = "Bloom year"
# this changes the linewidth outsidet he plot
opar = par(lwd = 2)

par(mar = c(5,5,3,3))
barplot(x, 
        names.arg = y, 
        main = "",
        ylab = "", 
        xlab = "",
        ylim = c(0, max(x) + max(x)*0.05),
        axes = FALSE,
        cex.names = 1.25,
        col = "grey80",
        #las=2,
        cex.axis=0.2,
        bourder = colvect("grey22", alpha = 0.9)
        #density = 1
)
box(which = "plot", lty = "solid", lwd = 3, col = colvect("grey22", alpha = 0.9))
grid(nx = 6, # X-axis divided in two sections
     ny = 6, # Y-axis divided in three sections
     lty = 2, 
     col = colvect(c("gray69"), alpha = 0.5), lwd = 1)
axis(side = 2,
     las  = 2,
     lwd  = 2,
     mgp  = c(1, 0.75, 0),
     cex.axis = 1.25,
     col = colvect("grey22", alpha = 0.9))
title(main = main,
      cex.lab = 2,
      line= 0.75,
      adj = 0)
title(ylab = ylab, cex.lab = 1.5, line = 2.75)
title(xlab = xlab, cex.lab = 1.5, line = 2.5)
dev.off()
# ------------------------------------------------------------------------------
# Big Bloom Areas

png(filename = paste("figures\\bloom_areas_", Sys.Date(), ".png",sep = ""),
    width = 3.5,
    height = 6,
    units = "in",
    res = 300,
    pointsize = 10)

tbl = sumtbl
#idx = tbl$area > 550000
#tbl = tbl[idx,]
y = year(tbl$sdate)
tbl[y == 2022, ]$area = 682983.1
x = tbl$area/100000
y = year(tbl$sdate)
main = ""
ylab = TeX(r'(Area $[km \times 10^5]$)')
xlab = "Bloom year"

# this changes the linewidth outsidet he plot
par(mar = c(5,5,3,3))
opar = par(lwd = 2)
barplot(x, 
        axes = FALSE,
        names.arg = y, 
        main = "",
        ylab = "", 
        xlab = "",
        ylim = c(0, max(x) + max(x)*0.05))
box(which = "plot", lty = "solid", lwd = 3, col = colvect("grey22", alpha = 0.9))
grid(nx = 6, # X-axis divided in two sections
     ny = 6, # Y-axis divided in three sections
     lty = 2, 
     col = colvect(c("gray69"), alpha = 0.5), lwd = 1)
axis(side = 2,
     las  = 2, 
     lwd  = 2, 
     mgp  = c(1, 0.75, 0), 
     cex.axis = 1.25,
     col = colvect("grey22", alpha = 0.9))
title(main = main,
      cex.lab = 2,
      line= 0.75,
      adj = 0)
title(ylab = ylab, cex.lab = 1.5, line = 2.75)
title(xlab = xlab, cex.lab = 1.5, line = 2.5)
dev.off()

# ------------------------------------------------------------------------------
# Bloom Duration
tbl = sumtbl
#idx = tbl$area > 550000
#tbl = tbl[idx,]

png(filename = paste("figures\\bloom_durations_", Sys.Date(), ".png",sep = ""),
    width = 3.5,
    height = 6,
    units = "in",
    res = 300,
    pointsize = 10)

start = yday(tbl$sdate)
end   = yday(tbl$edate)
# end[end < 100] = 365
# end[length(end)] = end[length(end)] + 60

durations = data.frame(id = as.character(year(tbl$sdate)), start, end)
durations$middle = apply(durations[, c("start", "end")], 1, FUN = mean, na.rm = TRUE)
durations$year = year(tbl$sdate)

# seperating years with two blooms to plot correctly
idx_dup   = duplicated(durations$year)
twobloom  = durations[idx_dup,]
durations = durations[!idx_dup,]

idx_dup    = duplicated(twobloom$year)
threebloom = twobloom[idx_dup,]
twobloom   = twobloom[!idx_dup,]

ggplot(durations) +
  geom_boxplot(
    fill = "grey100",
    colour = colvect("grey22", alpha = 0.9),
    size = 0.75,
    stat = "identity",
    aes(x = id, 
        lower  = start, 
        middle = middle, 
        upper  = end, 
        ymin   = start, 
        ymax   = end
    )) +
  geom_boxplot(data = twobloom,
               fill = "grey100",
               colour = colvect("grey22", alpha = 0.9),
               size = 0.75,
               stat = "identity",
               aes(x = id, 
                   lower  = start, 
                   middle = middle, 
                   upper  = end, 
                   ymin   = start, 
                   ymax   = end
               )) +
  geom_boxplot(data = threebloom,
               fill = "grey100",
               colour = colvect("grey22", alpha = 0.9),
               size = 0.75,
               stat = "identity",
               aes(x = id, 
                   lower  = start, 
                   middle = middle, 
                   upper  = end, 
                   ymin   = start, 
                   ymax   = end
               )) +
  xlab("Bloom year") +
  ylab("Day of year") +
  labs(title = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(colour = colvect("grey22", alpha = 0.9), 
                                fill = NA, 
                                size = 2),
    axis.title.x = element_text(size = rel(1.25)),
    axis.title.y = element_text(size = rel(1.25)),
    axis.text.x = element_text(size = rel(1.25)),
    axis.text.y = element_text(size = rel(1.25)),
    plot.title = element_text(size  = rel(1.25))
  ) +
  coord_flip()
dev.off()
# ------------------------------------------------------------------------------
# Bloom Centers

png(filename = paste("figures\\bloom_centers_", Sys.Date(), ".png",sep = ""),
    width = 3.5,
    height = 6,
    units = "in",
    res = 300,
    pointsize = 10)

tbl = sumtbl
e = extent(-173, -125, 11, 46)
wdmap <- getMap(resolution = "high")
temp = tbl
idx = tbl$lon > -165
tbl = tbl[idx,]
#idx = !year(tbl$sdate) == 2015
#tbl = tbl[idx,]

# radius in km
a = tbl$area
r = sqrt(a/pi) 
r = r/110.574 
tbl$rad = r 
tbl$y = as.character(year(tbl$sdate))

idx = order(tbl$area, decreasing = TRUE)
txt = tbl[idx,]
txt = txt[1:3,]

jamie_theme(e[1:2], e[3:4],
            main = "", 
            ylab = "Latitude", 
            xlab = "Longitude",
            dt = FALSE,
            xaxes = TRUE,
            yaxes = TRUE,
            mar = c(5,5,2,3)
            #asp = 1.08
)
points(tbl$lon, tbl$lat, pch = 19)
for(i in 1:nrow(tbl)) draw.circle(tbl$lon[i], 
                                  tbl$lat[i], 
                                  tbl$rad[i],
                                  border = colvect("grey22", alpha = 0.6),
                                  lwd = 2,
                                  nv  = 500,
                                  lty = 2)
plot(wdmap, 
     xlim = e[1:2], 
     ylim = e[3:4], 
     asp = 1, 
     bg = "black", 
     border = "black", 
     col = "black", 
     #wrap = c(-180, 180), 
     add = TRUE)
text(txt$lon, txt$lat, label = txt$y, adj = 1.2)
box(which = "plot", lty = "solid", lwd = 3, col = colvect("grey22", alpha = 0.9))
#legend(-140, 22, legend = "centers", pch = 20)
dev.off()

















