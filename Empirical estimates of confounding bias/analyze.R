
########################### PRELIMINARIES ########################### 

library(dplyr)
library(readxl)
library(robumeta)
library(Hmisc)
library(data.table)
library(testthat)
library(here)
library(MetaUtility)
library(msm)

prepped.data.dir = here("Empirical estimates of confounding bias/Prepped data")
code.dir = here("Empirical estimates of confounding bias")
results.dir = here("Empirical estimates of confounding bias/Results from R")

setwd(code.dir)
source("helper.R")

setwd(prepped.data.dir)
d = read.csv("prepped_merged_data.csv")


########################### CALCULATE CALIBRATED ESTIMATES WITHIN SOURCE ########################### 

# yi's are log-ratios (obs vs. RCT)
d = d %>% group_by(study) %>%
  mutate( calib = calib_ests(yi = yi,
                             sei = sqrt(vyi) ) )

# sanity check
expect_equal( d$calib[1],
              calib_ests(yi = d$yi[ d$study == d$study[1] ],
                   sei = sqrt(d$vyi[ d$study == d$study[1] ]) )[1] )

# big difference between looking at ratios with and without calibration! 
# e.g., 27.5% greater than 2-fold vs. only 6% after calibration
mean( d$yi > log(2) | d$yi < log(0.5) )
mean( d$calib > log(2) | d$calib < log(0.5) )


##### Summary Table of Calibrated Estimates #####
d %>% group_by(study) %>%
  summarise( ratioGt1.25 = mean( calib > log(1.1) | calib < log(1/1.1) ),
             ratioGt1.5 = mean( calib > log(1.5) | calib < log(1/1.5) ),
             ratioGt2 = mean( calib > log(2) | calib < log(0.5) ) )


d %>% group_by(study) %>%
  summarise( mean( calib < log(0.8) ) )


# calculate Phats for different ratios
t = d %>% group_by(study) %>%
  summarise( 
          ratioGt1.1 = quick_Phat( .q = log(1.1),  .yi = yi, .vyi = vyi )$Phat,
          ratioLt1.1.Inv = quick_Phat( .q = log(1/1.1), .yi = yi, .vyi = vyi )$Phat,
          
          ratioGt1.25 = quick_Phat( .q = log(1.25),  .yi = yi, .vyi = vyi )$Phat,
          ratioLt1.25.Inv = quick_Phat( .q = log(1/1.25), .yi = yi, .vyi = vyi )$Phat,
  
          ratioGt2 = quick_Phat( .q = log(2),  .yi = yi, .vyi = vyi  )$Phat,
          ratioLt2.Inv = quick_Phat( .q = log(1/2),  .yi = yi, .vyi = vyi  )$Phat
           )

# simple sanity check
# d %>% group_by(study) %>%
#   summarise( ratioGt1.25 = mean( calib > log(1.1) | calib < log(1/1.1) ),
#              ratioGt1.5 = mean( calib > log(1.5) | calib < log(1/1.5) ),
#              ratioGt2 = mean( calib > log(2) | calib < log(0.5) ) )


# save it (transposed)
setwd(results.dir)
write.csv( t(t), "Phat_table.csv")



########################### LOOK AT CALIBRATED ESTIMATES IN FRANKLIN ONLY ########################### 

sort( exp(d$calib[ d$study == "Franklin 2020"]) )

########################### PLOT CALIBRATED ESTIMATES WITHIN SOURCE ########################### 

# choose axis scaling
summary( exp(d$calib) )
xmin = -0.25
xmax = 1.75
tickJump = 0.5  # space between tick marks

colors = c("orange", "darkgray", "red")
levels( as.factor(d$study) )


ggplot( data = d,
        aes( x = calib,
             fill = study,
             color = study ) ) +
  
  # mean estimates from subset models
  geom_vline( xintercept = mean( d$calib[ d$study == "Franklin 2020" ] ),
              color = colors[1],
              lty = 2 ) +
  
  geom_vline( xintercept = mean( d$calib[ d$study == "Golder 2011" ] ),
              color = colors[2],
              lty = 2 ) +
  
  geom_vline( xintercept = mean( d$calib[ d$study == "Shikata 2006" ] ),
              color = colors[3],
              lty = 2 ) +
  
   geom_vline(xintercept = 0,
             color = "gray",
             lty = 1) +
  
  geom_density(alpha = 0.3) +
  
  theme_bw() +
  
  xlab("Estimated distribution of observational:RCT ratios") +
  scale_x_continuous( limits = c(xmin, xmax), breaks = seq(xmin, xmax, tickJump)) +
  
  ylab("Density") +
  
  scale_color_manual( values = rev(colors), name = "Source" ) +
  scale_fill_manual( values = rev(colors), name = "Source" ) +
  
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


my_ggsave( name = "calibrated_plot.pdf",
           width = 8,
           height = 5 )




