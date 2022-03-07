library(magrittr)
library(data.table)
library(ggpubr)
library(ggrepel)
library(vegan)
library(googlesheets4)
library(adespatial)
library(iNEXT)
library(future.apply)
library(plyr)
library(pheatmap)
library(geosphere)
library(BiodiversityR)
library(gmodels)
library(ncdf4)
library(ape)
library(factoextra)

Sys.getlocale()
# "English_United States.1252"
Sys.setlocale(locale = "cht")

pn <- 19999 # caution: may be very time consuming
d0 <- 
d <- fread("data.csv")
d[, site.ordered := ordered(site, levels = d[, .(area, site, lat, lon)] %>% unique %>% .[order(lat, decreasing = T)] %>% .$site)]
d[, substrate.ordered := ordered(substrate, levels = c("living snail", "hermit crab", "cobble", "reef"))]

d.x <- d[, !grepl("^mOTU[0-9]+", names(d)), with = F]
d.y <- d[, grepl("^mOTU[0-9]+", names(d)), with = F] %>% .[, colSums(.) > 0, with = F]

#### SCBD (22 samples * 26 mOTUs)
plan(multisession(workers = 10))
SCBD.0 <- beta.div(d.y, "hellinger", nperm = 1)$SCBD
this.d <-
  future_replicate(pn,
                   d.y %>% apply(., 1, sample) %>% t %>% set_colnames(colnames(d.y)),
                   simplify = F)
SCBD.pval <-
  future_sapply(this.d, function(x) {
    beta.div(x, method = "hellinger", nperm = 1)$SCBD
  }) %>%
  cbind(SCBD.0, .) %>%
  apply(., 1, function(x) {
    sum(x >= x[1]) / (pn + 1)
  })
SCBD <- rbind(SCBD = SCBD.0, pval = SCBD.pval)[, order(-SCBD.0)]
plan(sequential)
SCBD.pick <- which(cumsum(SCBD["SCBD", ]) < 0.95) %>% names %T>% print
#  [1] "mOTU027" "mOTU023" "mOTU139" "mOTU075" "mOTU097" "mOTU101" "mOTU116" "mOTU025" "mOTU098" "mOTU104" "mOTU103" "mOTU114"
# [13] "mOTU117" "mOTU058" "mOTU041"
SCBD.dt <-
  SCBD %>% {
    foo <- .
    data.table(mOTU = foo %>% colnames,
               SCBD = foo[1,],
               pval = foo[2,])
  } %>%
  .[order(SCBD, decreasing = T)] %>%
  .[, SCBD.cumsum := SCBD %>% cumsum] %>%
  .[, x := 1:.N] %>%
  print
#        mOTU        SCBD    pval SCBD.cumsum  x
#  1: mOTU027 0.258987217 0.00005   0.2589872  1
#  2: mOTU023 0.182777611 0.00005   0.4417648  2
#  3: mOTU139 0.095904537 0.01985   0.5376694  3
#  4: mOTU075 0.079387412 0.06380   0.6170568  4
#  5: mOTU097 0.069536258 0.11685   0.6865930  5
#  6: mOTU101 0.048636024 0.30725   0.7352291  6
#  7: mOTU116 0.042124862 0.40245   0.7773539  7
#  8: mOTU025 0.039185949 0.43435   0.8165399  8
#  9: mOTU098 0.035215857 0.49540   0.8517557  9
# 10: mOTU104 0.027084466 0.61665   0.8788402 10
# 11: mOTU103 0.015762360 0.81110   0.8946026 11
# 12: mOTU114 0.015079417 0.82120   0.9096820 12
# 13: mOTU117 0.013752658 0.84530   0.9234346 13
# 14: mOTU058 0.012324150 0.87125   0.9357588 14
# 15: mOTU041 0.012296258 0.87190   0.9480550 15
# 16: mOTU065 0.008449264 0.93130   0.9565043 16
# 17: mOTU100 0.008178611 0.93355   0.9646829 17
# 18: mOTU040 0.007102765 0.94955   0.9717857 18
# 19: mOTU002 0.006288709 0.96025   0.9780744 19
# 20: mOTU110 0.005197295 0.96880   0.9832717 20
# 21: mOTU013 0.003681572 0.98210   0.9869533 21
# 22: mOTU112 0.003510336 0.98320   0.9904636 22
# 23: mOTU147 0.003510336 0.98115   0.9939739 23
# 24: mOTU039 0.002515741 0.99130   0.9964897 24
# 25: mOTU106 0.001755168 0.99265   0.9982448 25
# 26: mOTU107 0.001755168 0.99240   1.0000000 26
#        mOTU        SCBD    pval SCBD.cumsum  x



#### dbRDA
permute.free <- how(nperm = pn)
m00 <- 
  capscale(
    d.y %>% decostand("hellinger") ~ 
      substrate.ordered + site.ordered + SST.year.month + SSS.year.month, 
    d.x, 
    distance = "bray")
m1.f <-
  formula(
    d.y %>% decostand("hellinger") ~ 
      (substrate.ordered + site.ordered + SST.year.month + SSS.year.month) ^ 2 + poly(SST.year.month, 2) + poly(SSS.year.month, 2)
  )
ordistep(m00, m1.f, permutations = permute.free) 
#                                    Df    AIC      F Pr(>F)  
# + poly(SST.year.month, 2)           1 17.373 1.8715 0.0898 .
# + SST.year.month:SSS.year.month     1 17.576 1.7535 0.1089  
# + poly(SSS.year.month, 2)           1 18.419 1.2739 0.2648  
# + site.ordered:SST.year.month       2 19.176 0.9032 0.5449  
# + substrate.ordered:SST.year.month  2 19.714 0.7606 0.6914  
# + site.ordered:SSS.year.month       2 19.798 0.7387 0.6974  
# + substrate.ordered:SSS.year.month  2 20.249 0.6223 0.8247  
# + substrate.ordered:site.ordered    0 18.829                

mbest <- update(m00, . ~ substrate.ordered + site.ordered + SST.year.month + SSS.year.month)
anova(mbest, by = "margin", permutations = permute.free)
#                   Df SumOfSqs      F  Pr(>F)    
# substrate.ordered  2  1.09904 6.6385 0.00010 ***
# site.ordered       5  1.25844 3.0405 0.00015 ***
# SST.year.month     1  0.07334 0.8859 0.41505    
# SSS.year.month     1  0.08163 0.9861 0.35850    
# Residual          12  0.99334                   

summary(mbest)
# Importance of components:
#                         CAP1    CAP2    CAP3    CAP4    CAP5    CAP6    CAP7     CAP8      CAP9    MDS1    MDS2    MDS3    MDS4
# Eigenvalue            3.2370 0.41312 0.39859 0.24433 0.20316 0.10139 0.06288 0.056375 0.0049982 0.33677 0.20019 0.18796 0.08891
# Proportion Explained  0.5664 0.07228 0.06974 0.04275 0.03555 0.01774 0.01100 0.009864 0.0008746 0.05893 0.03503 0.03289 0.01556
# Cumulative Proportion 0.5664 0.63867 0.70841 0.75116 0.78671 0.80445 0.81545 0.825318 0.8261929 0.88512 0.92015 0.95303 0.96859
#                         MDS5     MDS6     MDS7     MDS8     MDS9    MDS10     MDS11     MDS12
# Eigenvalue            0.0663 0.035659 0.032312 0.019805 0.012300 0.009376 0.0035296 2.218e-04
# Proportion Explained  0.0116 0.006239 0.005654 0.003465 0.002152 0.001641 0.0006176 3.881e-05
# Cumulative Proportion 0.9802 0.986432 0.992086 0.995551 0.997703 0.999344 0.9999612 1.000e+00


#### dbRDA plot (3 panels)
mbest1 <- mbest
d.spp <-
  scores(mbest)$species %>% as.data.table(keep.rownames = "mOTU") %>%
  .[, length := CAP1 ^ 2 + CAP2 ^ 2] %>%
  .[order(length, decreasing = T)]

color8 <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')[c(1:6, 9, 11)] %>% rev
color8.1 <- c('#00429d', '#4771b2', '#aaaaaa', '#a5d5d8', '#ffbcaf', '#f4777f', '#cf3759', '#93003a')
color8.2 <- "white" %>% rep(8)

f.site <-
  cbind(d, scores(mbest1)$site) %>%
  ggplot(aes(CAP1, CAP2)) +
  geom_point(
    aes(
      color = site.ordered,
      shape = site.ordered,
    ),
    stroke = 1,
    size = 1,
    position = position_jitter(width = 0, height = 0)
  ) +
  coord_fixed() +
  theme_pubr(8, border = T, legend = "right") + 
  theme(legend.text = element_text(size = 8)) +
  scale_color_manual("site.ordered", values = color8.1) +
  scale_shape_manual("site.ordered", values = c(0,1,2,3,4,5)) +
  geom_segment(
    data = d.spp,
    aes(
      x = 0,
      y = 0,
      xend = CAP1,
      yend = CAP2
    ),
    arrow = arrow(length = unit(0, "mm"))
  ) +
  geom_text(data = d.spp[mOTU %in% SCBD.pick[1:9]],
            aes(x = CAP1*1, y = CAP2*1, label = mOTU),
            size = 0.352777778 * 8) +
  scale_x_continuous("CAP1 (56.6%)", expand = expansion(c(0.1, 0.1))) + 
  scale_y_continuous("CAP2 (7.2%)", expand = expansion(c(0.1, 0.1)))
f.substrate <-
  cbind(d, scores(mbest1)$site) %>%
  ggplot(aes(CAP1, CAP2)) +
  geom_point(
    aes(
      color = substrate.ordered,
      shape = substrate.ordered
    ),
    stroke = 1,
    position = position_jitter(width = 0, height = 0),
    size = 1,
  ) +
  coord_fixed() +
  theme_pubr(8, border = T, legend = "right") + 
  theme(legend.text = element_text(size = 8)) +
  scale_color_manual("substrate.ordered", values = c('#00429d', '#73a2c6', '#f4777f', '#93003a')) +
  geom_segment(
    data = d.spp,
    aes(
      x = 0,
      y = 0,
      xend = CAP1,
      yend = CAP2
    ),
    arrow = arrow(length = unit(0, "mm"))
  ) +
  geom_text(data = d.spp[mOTU %in% SCBD.pick[1:9]],
            aes(x = CAP1*1, y = CAP2*1, label = mOTU),
            size = 0.352777778 * 8) +
  scale_x_continuous("CAP1 (56.6%)", expand = expansion(c(0.1, 0.1))) + 
  scale_y_continuous("CAP2 (7.2%)", expand = expansion(c(0.1, 0.1)))
f.season <-
  cbind(d, scores(mbest1)$site) %>%
  ggplot(aes(CAP1, CAP2)) +
  geom_point(
    aes(
      color = SST.year.month,
      size = SSS.year.month
    ),
    stroke = 1,
    position = position_jitter(width = 0, height = 0)
  ) +
  coord_fixed() +
  theme_pubr(8, border = T, legend = "right") + 
  theme(legend.text = element_text(size = 8)) +
  scale_colour_steps("SST (°C)", n.breaks = 6) +
  scale_size_binned("SSS (PSU)", range = c(0.2, 3)) +
  geom_segment(
    data = d.spp,
    aes(
      x = 0,
      y = 0,
      xend = CAP1,
      yend = CAP2
    ),
    arrow = arrow(length = unit(0, "mm"))
  ) +
  geom_text(data = d.spp[mOTU %in% SCBD.pick[1:9]],
            aes(x = CAP1*1, y = CAP2*1, label = mOTU),
            size = 0.352777778 * 8) +
  geom_segment(
    data = mbest1$CCA$biplot %>% as.data.table(keep.rownames = "IV") %>% .[grep("^SS", IV)],
    aes(
      x = 0,
      y = 0,
      xend = CAP1*4,
      yend = CAP2*4
    ),
    arrow = arrow(length = unit(3, "mm")),
    color = "red"
  ) +
  geom_text(data = mbest1$CCA$biplot %>% as.data.table(keep.rownames = "IV") %>% .[grep("^SS", IV)],
            aes(x = CAP1*4, y = CAP2*4, label = IV),
            size = 0.352777778 * 8,
            color = "red") +
  scale_x_continuous("CAP1 (56.6%)", expand = expansion(c(0.1, 0.1))) + 
  scale_y_continuous("CAP2 (7.2%)", expand = expansion(c(0.1, 0.1)))

ggarrange(
  f.site + theme_pubr(8, border = T, legend = "top") + theme(legend.text = element_text(size = 8), legend.direction = "vertical") + geom_hline(yintercept = 0, color = gray(0.8), linetype = 2) + geom_vline(xintercept = 0, color = gray(0.8), linetype = 2),
  f.substrate + theme_pubr(8, border = T, legend = "top") + theme(legend.text = element_text(size = 8), legend.direction = "vertical") + geom_hline(yintercept = 0, color = gray(0.8), linetype = 2) + geom_vline(xintercept = 0, color = gray(0.8), linetype = 2),
  f.season + theme_pubr(8, border = T, legend = "top") + theme(legend.text = element_text(size = 8), legend.direction = "vertical") + geom_hline(yintercept = 0, color = gray(0.8), linetype = 2) + geom_vline(xintercept = 0, color = gray(0.8), linetype = 2),
  align = "hv",
  nrow = 1
)

#### Variance partition
RsquareAdj(mbest)
# $r.squared
# [1] 0.884005
# 
# $adj.r.squared
# [1] 0.7970087

d.x.mat <-
  model.matrix(
    ~ substrate + site + scale(SST.year.month) + scale(SSS.year.month),
    data = d.x,
    contrasts.arg = list(
      substrate = contr.helmert,
      site = contr.helmert
    )
  )[,-1]
colnames(d.x.mat)
varpart(d.y %>% decostand("hellinger") %>% vegdist("bray"),
        d.x.mat[, 3:7], # site
        d.x.mat[, 1:2], # substrate
        d.x.mat[, 8:9]  # SST & SSS
        ) %T>% 
  print %>% 
  plot(cutoff = -1, digits = 2)
# Partition table:
#                       Df R.square Adj.R.square Testable
# [a+d+f+g] = X1         5  0.62821      0.51203     TRUE
# [b+d+e+g] = X2         2  0.60309      0.56131     TRUE
# [c+e+f+g] = X3         2  0.03156     -0.07038     TRUE
# [a+b+d+e+f+g] = X1+X2  7  0.82742      0.74112     TRUE
# [a+c+d+e+f+g] = X1+X3  7  0.64960      0.47441     TRUE
# [b+c+d+e+f+g] = X2+X3  4  0.63297      0.54661     TRUE
# [a+b+c+d+e+f+g] = All  9  0.85199      0.74097     TRUE
# Individual fractions                                   
# [a] = X1 | X2+X3       5               0.19436     TRUE
# [b] = X2 | X1+X3       2               0.26657     TRUE
# [c] = X3 | X1+X2       2              -0.00015     TRUE
# [d]                    0               0.35043    FALSE
# [e]                    0              -0.03747    FALSE
# [f]                    0              -0.01455    FALSE
# [g]                    0              -0.01821    FALSE
# [h] = Residuals                        0.25903    FALSE
# Controlling 1 table X                                  
# [a+d] = X1 | X3        5               0.54479     TRUE
# [a+f] = X1 | X2        5               0.17981     TRUE
# [b+d] = X2 | X3        2               0.61700     TRUE
# [b+e] = X2 | X1        2               0.22910     TRUE
# [c+e] = X3 | X1        2              -0.03762     TRUE
# [c+f] = X3 | X2        2              -0.01470     TRUE

#### heatmap
map.mat <- 
  cbind(
    d.y %>% decostand("total", margin = 1),
    d.x[, .(site.ordered, substrate.ordered, plot, SST.year.month, SSS.year.month)]
  ) %>%
  as.data.table %>%
  melt(
    .,
    id.vars = c("site.ordered", "substrate.ordered", "plot", "SST.year.month", "SSS.year.month"),
    variable.name = "mOTU",
    value.name = "ratio"
  )
map.mat.3factor <- 
  map.mat %>%
  dcast(site.ordered + substrate.ordered + plot + SST.year.month + SSS.year.month ~ mOTU, value.var = "ratio") %>% {
    foo <- .
    foo[, -1:-5] %>%
      as.matrix %>%
      set_rownames(paste(foo$site.ordered, foo$substrate.ordered, round(foo$SST.year.month,1), round(foo$SSS.year.month,1), sep = "-"))
  } %>%
  .[, SCBD.pick]
drows <- map.mat.3factor %>% cbind(Others = 1 - rowSums(.)) %>% decostand("hellinger") %>% vegdist
map.mat.3factor %>%
  cbind(Others = 1 - rowSums(.)) %>% {
    foo <- .
    foo[abs(foo) < 0.000001] <- NA
    foo
  } %>%
  pheatmap(
    .,
    na_col = 8,
    breaks = seq(0, 1, 0.1),
    color = c('#94003a', '#a33047', '#b14c55', '#bf6664', '#cd8074', '#d99984', '#e5b296', '#f0cbaa', '#fae5c0', '#ffffe0') %>% rev,
    scale = "none",
    cluster_rows = T,
    cluster_cols = F,
    clustering_distance_rows = drows
  )


#### distance decay (with Baosheng)
mantel(
  distm(d.x[substrate == "reef", .(lon, lat)]),
  d.y[d.x$substrate == "reef"] %>% decostand("hellinger") %>% vegdist("bray", diag = T, upper = T),
  method = "pearson",
  permutations = 19999
)
# Mantel statistic r: 0.5512 
#       Significance: 0.00075 
ff2 <-
  data.table(
    `geographical distance (km)` = distm(d.x[substrate == "reef", .(lon, lat)]) %>% .[lower.tri(.)] %>% divide_by(1000),
    `Bray-Curtis dissimilarity` = d.y[d.x$substrate == "reef"] %>% decostand("hellinger") %>% vegdist("bray", diag = F, upper = F) %>% as.vector) %>%
  ggplot(aes(`geographical distance (km)`, 1 - `Bray-Curtis dissimilarity`)) +
  geom_point(color = 1,
             shape = 21,
             alpha = 1) +
  geom_smooth(method = MASS::rlm, se = F) +
  theme_pubr(8, border = T, legend = "right") +
  theme(legend.text = element_text(size = 8)) +
  scale_y_continuous("Community similarity", breaks = seq(0, 1, 0.2), limits = c(0, 1))

#### distance decay (without Baosheng)
d.substrate.content <- 
  fread("5sites-sedimentation.csv") %>% 
  .[, month := factor(month, levels = month.name[c(7,8,10:12)], labels = month.abb[c(7,8,10:12)], ordered = T)] %>% 
  .[, lapply(.SD, mean, na.rm = TRUE), by = .(Site, month), .SDcols = c(3:4)] %>%
  .[order(month, Site), ] %>%
  .[, Weight := Mud + Sand] %>%
  .[, logRatio := log(Sand / Mud)] %>% 
  dcast(Site ~ month, value.var = c("Weight", "logRatio")) %>%
  setnames("Site", "site") %>% {
    cbind(.[, .(site)],
          .[, lapply(.SD, scale), .SDcol = 2:11])
  } %T>%
  print
rows <- which(d.x$site %in% d.substrate.content$site & d.x$substrate == "reef")
d.xx <- merge(
  d.x[rows],
  d.substrate.content,
  by = "site",
  sort = F
)
d.yy <- d.y[rows]

mantel(
  distm(d.xx[, .(lon, lat)]),
  d.yy %>% decostand("hellinger") %>% vegdist("bray", diag = T, upper = T),
  method = "pearson",
  permutations = 19999
)
# Mantel statistic r: 0.5783 
#       Significance: 0.0014 
mantel(
  d.xx[, .(lon, lat)] %>% unique %>% distm,
  d.xx[, Weight_Jul.V1:logRatio_Dec.V1]  %>% unique %>% dist,
  # d.xx[, mud.mean_August.V1:sand.mean_October.V1]  %>% unique %>% dist,
  method = "pearson",
  permutations = 19999
)
# Mantel statistic r: 0.6494 
#       Significance: 0.025 
mantel(
  # d.xx[, mud.mean_August.V1:sand.mean_October.V1] %>% dist,
  d.xx[, Weight_Jul.V1:logRatio_Dec.V1] %>% dist,
  d.yy %>% decostand("hellinger") %>% vegdist("bray", diag = T, upper = T),
  method = "pearson",
  permutations = 19999
)
# Mantel statistic r: 0.5689 
#       Significance: 0.0033 


f.ty.1 <-
  data.table(
    `Geographical distance (km)` = distm(d.xx[, .(lon, lat)]) %>% .[lower.tri(.)] %>% divide_by(1000),
    `Bray-Curtis dissimilarity` = d.yy[] %>% decostand("hellinger") %>% vegdist("bray", diag = F, upper = F) %>% as.vector
  ) %>%
  ggplot(aes(`Geographical distance (km)`, 1 - `Bray-Curtis dissimilarity`)) +
  geom_point(color = 1,
             shape = 21,
             alpha = 1) +
  geom_smooth(method = MASS::rlm, se = F) +
  theme_pubr(8, border = T, legend = "right") +
  theme(legend.text = element_text(size = 8)) +
  scale_x_continuous(breaks = seq(0, 10, 2.5), limits = c(0, 10)) +
  scale_y_continuous("Community similarity", breaks = seq(0, 1, 0.2), limits = c(0, 1))
f.ty.2 <-
  data.table(
    `Geographical distance (km)` = d.xx[, .(lon, lat)] %>% unique %>% distm %>% .[lower.tri(.)] %>% divide_by(1000), 
    # `Sedimentation dissimilarity` = d.xx[, mud.mean_August.V1:sand.mean_October.V1] %>% unique %>%  dist %>% as.vector
    `Sedimentation dissimilarity` = d.xx[, Weight_Jul.V1:logRatio_Dec.V1] %>% unique %>%  dist %>% as.vector
  ) %>%
  ggplot(aes(`Geographical distance (km)`, `Sedimentation dissimilarity`)) +
  geom_point(color = 1,
             shape = 21,
             alpha = 1) +
  geom_smooth(method = MASS::rlm, se = F) +
  theme_pubr(8, border = T, legend = "right") +
  theme(legend.text = element_text(size = 8)) +
  scale_x_continuous(breaks = seq(0, 10, 2.5), limits = c(0, 10))
f.ty.3 <-
  data.table(
    # `Sedimentation dissimilarity` = d.xx[, mud.mean_August.V1:sand.mean_October.V1] %>%  dist %>% as.vector,
    `Sedimentation dissimilarity` = d.xx[, Weight_Jul.V1:logRatio_Dec.V1] %>%  dist %>% as.vector,
    `Bray-Curtis dissimilarity` = d.yy[] %>% decostand("hellinger") %>% vegdist("bray", diag = F, upper = F) %>% as.vector
  ) %>%
  ggplot(aes(`Sedimentation dissimilarity`, 1 - `Bray-Curtis dissimilarity`)) +
  geom_point(color = 1,
             shape = 21,
             alpha = 1) +
  geom_smooth(method = MASS::rlm, se = F) +
  theme_pubr(8, border = T, legend = "right") +
  theme(legend.text = element_text(size = 8)) +
  scale_y_continuous("Community similarity", breaks = seq(0, 1, 0.2), limits = c(0, 1))
ggarrange(
  ff2 + xlab("Geographical distance (km)"),
  f.ty.1 + xlab("Geographical distance (km)"),
  f.ty.3 + xlab("Geographical distance (km)"),
  f.ty.2,
  nrow = 2, ncol = 2, 
  align = "hv"
)

#### sedimentation plot
f.mud1 <-
  fread("5sites-sedimentation.csv") %>% 
  .[, lapply(.SD, mean, na.rm = TRUE), by = .(Site, month), .SDcols = c(3:4)] %>%
  .[order(month, Site), ] %>%
  .[, Weight := Mud + Sand] %>%
  .[, logRatio := log(Sand / Mud)] %>% 
  dcast(Site ~ month, value.var = c("Weight", "logRatio")) %>%
  setnames("Site", "site") %>%
  melt(id.vars = "site") %>%
  .[, type := tstrsplit(variable, "[\\._]", keep = 1)] %>%
  .[, month := tstrsplit(variable, "[\\._]", keep = 2)] %>% 
  .[, month := factor(month, levels = month.name[c(7,8,10:12)], labels = month.abb[c(7,8,10:12)], ordered = T)] %>% 
  .[, variable := NULL] %>%
  .[, site := factor(site, levels = c("Baiyu", "Datan G1", "Datan G2", "Yongxin", "Yongan") %>% rev, ordered = T)] %>%
  dcast(site + month ~ type) %>%
  ggplot(aes(month, site)) +
  geom_point(aes(size = logRatio, fill = logRatio), color = 8, shape = 21) +
  scale_size_binned("log sand-mud ratio", n.breaks = 5, range = c(1, 6)) +
  scale_fill_binned("log sand-mud ratio", n.breaks = 5) +
  theme_pubr(8, border = T, legend = "top") + theme(legend.text = element_text(size = 8)) +
  xlab("Month") + ylab("Site")
f.mud2 <-
  fread("5sites-sedimentation.csv") %>% 
  .[, lapply(.SD, mean, na.rm = TRUE), by = .(Site, month), .SDcols = c(3:4)] %>%
  .[order(month, Site), ] %>%
  .[, Weight := Mud + Sand] %>%
  .[, logRatio := log(Sand / Mud)] %>% 
  dcast(Site ~ month, value.var = c("Weight", "logRatio")) %>%
  setnames("Site", "site") %>%
  melt(id.vars = "site") %>%
  .[, type := tstrsplit(variable, "[\\._]", keep = 1)] %>%
  .[, month := tstrsplit(variable, "[\\._]", keep = 2)] %>% 
  .[, month := factor(month, levels = month.name[c(7,8,10:12)], labels = month.abb[c(7,8,10:12)], ordered = T)] %>% 
  .[, variable := NULL] %>%
  .[, site := factor(site, levels = c("Baiyu", "Datan G1", "Datan G2", "Yongxin", "Yongan") %>% rev, ordered = T)] %>%
  dcast(site + month ~ type) %>%
  ggplot(aes(month, site)) +
  geom_point(aes(size = log(Weight), fill = log(Weight)), color = 8, shape = 21) +
  scale_size_binned("Sedimentation rate\n(log mg cm-2 d-1)", n.breaks = 5, range = c(1, 6)) +
  scale_fill_binned("Sedimentation rate\n(log mg cm-2 d-1)", n.breaks = 5) +
  theme_pubr(8, border = T, legend = "top") + theme(legend.text = element_text(size = 8)) +
  xlab("Month") + ylab("Site")
ggarrange(f.mud2 + coord_fixed(), f.mud1 + coord_fixed(), align = "hv")


#### SST cluster
d.cluster.site <- fread("data-site.csv")
YYYY <- 2000:2020
MONTH <- 1:12
d.cluster.site.SST <-
  expand.grid(site = d.cluster.site$site, YYYY = YYYY, MONTH = MONTH, stringsAsFactors = F) %>% as.data.frame %>%
  merge(d.cluster.site %>% as.data.frame, ., by = "site") %>%
  as.data.table  # rows = 2016 = 8site*12MM*21YYYY 

d.cluster.site.SST <- fread("SST-cluster.csv")
d.cluster.site.SST[, date := as.Date(paste(YYYY, MONTH, 15, sep = "-"))]
d.cluster.site.SST[, site.ordered := ordered(site, levels = d.cluster.site[, .(area, site, lat, lon)] %>% .[order(lat, decreasing = T)] %>% .$site %>% unique)]
d.cluster.site.SST[, area.ordered := ordered(area, levels = d.cluster.site[, .(area, lat)] %>% .[order(lat, decreasing = T)] %>% .$area %>% unique)] 
setkey(d.cluster.site.SST, date, site)

f0 <-
  d.cluster.site.SST %>%
  .[, .(SSTmonthly = mean(SSTmonthly)), by = list(site, site.ordered, area, area.ordered, date, YYYY, MONTH)] %>%
  .[, .(ave = mean(SSTmonthly), max = max(SSTmonthly), min = min(SSTmonthly)), by = .(area.ordered, site.ordered, YYYY)] %>%
  melt(., id.vars = c("YYYY", "area.ordered", "site.ordered")) %>%
  ggplot(aes(YYYY, value)) +
  geom_line(aes(color = area.ordered, group = site.ordered)) +
  facet_wrap(vars(variable), scales = "free_y") +
  theme_pubr(8, border = T) +
  theme(strip.text = element_text(size = 8), legend.text = element_text(size = 8)) +
  guides(color = guide_legend(title = "Area")) +
  xlab("Year")


d.cluster.site.SST.wide <- 
  d.cluster.site.SST %>% 
  .[, .(SSTmonthly = mean(SSTmonthly)), by = list(site, site.ordered, area, area.ordered, date, YYYY, MONTH)] %>%
  dcast(site.ordered + site + area + area.ordered ~ date, value.var = "SSTmonthly")
SST.dist <- dist(d.cluster.site.SST.wide[, -1:-4] %>% as.matrix %>% set_rownames(d.cluster.site.SST.wide$site), "manhattan")
fit.pcoa <- pcoa(SST.dist, correction = "none") %T>% print
d.cluster.site.SST.wide[, PCoA1 := fit.pcoa$vector[, 1]]
d.cluster.site.SST.wide[, PCoA2 := fit.pcoa$vector[, 2]]
f1 <-
  ggplot(d.cluster.site.SST.wide, aes(PCoA1, PCoA2)) +
  geom_hline(yintercept = 0, color = 8, linetype = 2) +
  geom_vline(xintercept = 0, color = 8, linetype = 2) +
  geom_point(aes(color = area)) +
  # geom_text(aes(label = site), size = 8*25.6/72, nudge_x = 1, hjust = 0) +
  coord_fixed() +
  theme_pubr(8, border = T) +
  theme(legend.text = element_text(size = 8)) +
  xlab("PCoA1 (%)") +
  ylab("PCoA2 (%)")

f21 <- 
  fviz_nbclust(d.cluster.site.SST.wide[, -1:-4], kmeans, method = "wss", k.max = 10, diss = SST.dist) +
  theme_pubr(8, border = T)
f22 <- 
  fviz_nbclust(d.cluster.site.SST.wide[, -1:-4], kmeans, method = "silhouette", k.max = 10, diss = SST.dist) +
  theme_pubr(8, border = T)


hclust.fit <- hclust(SST.dist, method = "ward.D")
color14 <- c('#00429d', '#325da9', '#4e78b5', '#6694c1', '#80b1cc', '#9dced6', '#c0eade', '#ffdac4', '#ffb3a7', '#fb8a8c', '#eb6574', '#d5405e', '#b81b4a', '#93003a')
label.cols <- 
  data.table(site = hclust.fit$labels) %>% 
  merge(d.cluster.site.SST[, .(site.ordered, area.ordered, site, area)] %>% unique, by = "site", sort = F) %T>% print %>% {
    foo <- .
    foo[, col := area.ordered]
    levels(foo$col) <- color14 %>% rev
    foo
  } %>% 
  .$col
f3 <-
  fviz_dend(hclust.fit, k = 3, horiz = T, cex = 7/12, color_labels_by_k = F, label_cols = label.cols %>% as.character) + 
  theme_pubr(8, border = F) +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title.y = element_blank())
  

ggarrange(
  f0 + 
    theme_pubr(8, border = T, legend = "top") + 
    theme(legend.text = element_text(size = 8), strip.text = element_text(size = 8)) +
    scale_color_manual(values = color14), 
  ggarrange(
    f21,f22,
    align = "hv",
    nrow = 1,
    ncol = 2,
    widths = c(1, 1)
  ),
  f3,
  ncol = 1, nrow = 3, heights = c(0.8, 0.6, 2)
) 

#### coverage plot
d.cover <- fread("coverage.csv") %>% .[, benthos := factor(benthos, levels = c("CCA", "MA", "Turf", "ACA", "Zoantharia"))]
d.cover.ave <- d.cover[, .(cover.ave = mean(cover)), by = .(site, season, benthos)]
d.cover.wide <- 
  dcast(d.cover, site + season + replicate ~ benthos, value.var = "cover") %>% 
  .[, Others := 100 - ACA - CCA - MA - Turf - Zoantharia]
d.cover.ave.wide <-
  dcast(d.cover.ave, site + season ~ benthos, value.var = "cover.ave") %>% 
  .[, Others := 100 - ACA - CCA - MA - Turf - Zoantharia]

d.cover.ave.wide %>%
  melt(id.vars = c("site", "season"), variable.name = "benthos", value.name = "coverage") %>%
  ggplot(aes(season)) +
  geom_col(aes(y = coverage, fill = benthos), color = 1, width = 0.8, position = position_stack(reverse = T)) +
  scale_y_continuous("Coverage (%)", expand = expansion(0), limits = c(0, 120), breaks = seq(0,100,10)) +
  scale_fill_manual("Benthos", values = c('#00348c', '#535195', '#80709d', '#a792a5', '#cbb7ae', '#dee1be')) +
  coord_flip() +
  theme_pubr(legend = "right", border = T, base_size = 8) +
  theme(legend.text = element_text(size = 8)) +
  scale_x_discrete("Site-Season", expand = expansion(0)) +
  theme(legend.key.size = unit(3, "mm"))
