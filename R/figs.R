# setup ---------------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(here)
library(patchwork)
library(EBASE)
library(readr)
library(scales)

source(here('R/funcs.R'))

fwdat <- read_csv(file = url('https://raw.githubusercontent.com/fawda123/BASEmetab_script/master/data/apafwoxy2.csv'))

fwdatcmp <- fwdat %>% 
  mutate(
    DateTimeStamp = dmy_hms(datet, tz = 'America/Jamaica'),
    Date = as.Date(DateTimeStamp, tz = 'America/Jamaica'),
    DO_obs = `oxy,mmol/m3`, 
    a = `aparam,(mmolO2/m2/d)/(W/m2)`, # no conversion to volumetric, just for this plot
    R = `er,mmol/m2/d`,
    P = `gpp,mmol/m2/d`,
    D = -1 * `gasex,mmol/m2/d`,
    b = 100 * 3600 * `kw,m/s` / `wspd2,m2/s2` / (`sc,dimensionless` / 660) ^ -0.5 # (m/s)/(m2/s2) to (cm/hr) / (m2/s2)
  ) %>% 
  select(Date, DateTimeStamp, DO_obs, a, b, P, R, D)

# fwoxy for input to ebase
fwdatinp <- fwdat %>% 
  mutate(
    datet = dmy_hms(datet, tz = 'America/Jamaica')
  ) %>% 
  select(
    DateTimeStamp = datet,
    DO_obs = `oxy,mmol/m3`, 
    Temp = `temp,degC`, 
    Sal = `salt,ppt`, 
    PAR = `par,W/m2`, 
    WSpd = `wspd2,m2/s2`, 
    H = `ht,m`
  ) %>% 
  mutate(
    DO_obs = DO_obs / 1000 * 32, # to mg/L
    WSpd = sqrt(WSpd)
  )

# prior distribution plot ---------------------------------------------------------------------

p <- prior_plot(n = 3e6)
p <- p + theme(axis.text = element_text(size = 10))

png(here('figs/priorplot.png'), height = 2.5, width = 6.5, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# simulated apa data --------------------------------------------------------------------------

p1 <- ggplot(fwdatcmp, aes(x = DateTimeStamp, y = DO_obs)) + 
  geom_line(linewidth = 0.2) + 
  scale_x_datetime(breaks = "1 month", labels = date_format("%b")) +
  theme_minimal() + 
  labs(
    x = NULL, 
    y = expression(paste(O [2], ' (mmol ', m^-3, ')')),
    title = '(a) Synthetic dissolved oxygen'
  )

p2 <- ggplot(fwdatcmp, aes(x = DateTimeStamp, y = a)) + 
  geom_line(linewidth = 0.2) + 
  scale_x_datetime(breaks = "1 month", labels = date_format("%b")) +
  theme_minimal() + 
  labs(
    x = NULL, 
    y = expression(paste(italic(a), ' (mmol ', m^-2, d^-1, ') / (W ', m^{-2}, ')')), 
    title = expression(paste('(b) Synthetic ', italic(a), ' parameter'))
  )

p3 <- ebase_plot(fwdatcmp, instantaneous = F) +
  scale_x_date(breaks = "1 month", labels = date_format("%b")) +
  labs(title = '(c) Synthetic metabolic estimates')

p <- p1 + p2 + p3 + plot_layout(ncol = 1) & 
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10)
    )

png(here('figs/synapa.png'), height = 7, width = 7, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# Fwoxy apa comparison ------------------------------------------------------------------------

fl <- paste0(tempdir(), '/apasumdat.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apasumdat.RData', destfile = fl)
load(file = fl)

p1 <- priorcomp(apasumdat, met = 'nse')

png(here('figs/priorcomp.png'), height = 4, width = 6, family = 'serif', units = 'in', res = 500)
print(p1)
dev.off()

# Fwoxy apa comparison met sum --------------------------------------------

fl <- paste0(tempdir(), '/apasumdat.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apasumdat.RData', destfile = fl)
load(file = fl)

p1 <- priorsumcomp(apasumdat, met = 'nse')

png(here('figs/priorsumcomp.png'), height = 4, width = 5, family = 'serif', units = 'in', res = 500)
print(p1)
dev.off()

# best/worst fwoxy apa comp -------------------------------------------------------------------

# grid ests
fl1 <- paste0(tempdir(), '/apagrd1.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apagrd1.RData', destfile = fl1)
load(file = fl1)
fl7 <- paste0(tempdir(), '/apagrd7.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apagrd7.RData', destfile = fl7)
load(file = fl7)
fl30 <- paste0(tempdir(), '/apagrd30.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apagrd30.RData', destfile = fl30)
load(file = fl30)

apagrd <- bind_rows(apagrd1, apagrd7, apagrd30)

# summary ests
fl <- paste0(tempdir(), '/apasumdat.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apasumdat.RData', destfile = fl)
load(file = fl)

# row y axis limits
lims <- list(
  P = c(0, 600), 
  R = c(0, 600), 
  a = c(0, 11)
)

p <- optex(apagrd, fwdatcmp, apasumdat, rnkmetsum = c(1, 16), met = 'nse', lims = lims)

png(here('figs/optex.png'), height = 8.7, width = 9, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# synthetic with noise ------------------------------------------------------------------------

# load apacp data run through wtreg
fl <- paste0(tempdir(), '/apacpdtd.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apacpdtd.RData', destfile = fl)
load(file = fl)

nosdat <- apacpdtd %>% 
  mutate(
    tidnoise = DO_prd - DO_nrm, 
    obsnoise = DO_obs - DO_prd
  ) %>% 
  select(DateTimeStamp, DO_act = DO_obs, DO_prd, DO_nrm, tidnoise, obsnoise)

# add tidal noise to fwdatinp - need to figure out values less than zero
tomod <- fwdatinp %>% 
  left_join(nosdat, by = 'DateTimeStamp') %>% 
  mutate(
    DO_nos = pmax(0, DO_obs + tidnoise + obsnoise)
  )

toplo <- tomod %>% 
  filter(month(DateTimeStamp) == 6) 
toplo1 <- toplo %>% 
  select(DateTimeStamp, DO_obs, DO_nos) %>% 
  pivot_longer(-DateTimeStamp) %>% 
  mutate(
    name = case_when(
      name == 'DO_obs' ~ 'Synthetic', 
      name == 'DO_nos' ~ 'Synethic + noise'
    ), 
    value = value / 32 * 1000
  )
toplo2 <- toplo %>% 
  select(DateTimeStamp, tidnoise, obsnoise) %>%
  pivot_longer(-DateTimeStamp) %>% 
  mutate(
    name = case_when(
      name == 'obsnoise' ~ 'Instrument', 
      name == 'tidnoise' ~ 'Tidal'
    ),
    value = value / 32 * 1000
  )

p1 <- ggplot(toplo1, aes(x = DateTimeStamp, y = value, color = name)) + 
  geom_line() + 
  labs(
    color = NULL, 
    title = '(a) Time series'
  )
p2 <- ggplot(toplo2, aes(x = DateTimeStamp, y = value, color = name)) + 
  geom_line() + 
  scale_color_manual(values = c('darkgrey', 'black')) +
  labs(
    color = NULL,
    title = '(b) Noise'
  )

p <- p1 + p2 + plot_layout(ncol = 1) & scale_x_datetime(expand = c(0,0)) & theme_minimal() & labs(x = NULL, y = expression(paste(O [2], ' (mmol ', m^-3, ')'))) & theme(panel.grid.minor = element_blank(), legend.position = 'top')

png(here('figs/synapanos.png'), height = 6, width = 8, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# synthetic with noise results ----------------------------------------------------------------

# ebase results on observed and noisy observed
fl <- paste0(tempdir(), '/resobs.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/resobs.RData', destfile = fl)
load(file = fl)
fl <- paste0(tempdir(), '/resnos.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/resnos.RData', destfile = fl)
load(file = fl)

p <- syncomp_plo(resobs, resnos) 
  
png(here('figs/synapanoscmp.png'), height = 4.75, width = 6, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# actual apa comparison -----------------------------------------------------------------------

load(file = url('https://github.com/fawda123/BASEmetab_script/raw/master/data/apacmp.RData'))

toplo <- apacmp %>% 
  filter(typ != 'BASEmetab') %>% 
  mutate(
    val = case_when(
      var == 'Rt_vol' ~ -1 * val, 
      T ~ val
    ), 
    var = factor(var, 
                 levels = c('NEM', 'P', 'R', 'D')
    ), 
    seas = case_when(
      month(Date) %in% c(1:5, 10:12) ~ 'dry', 
      T ~ 'wet'
    )
  ) %>% 
  pivot_wider(names_from = 'typ', values_from = 'val')

p1a <- apacmp_plo(toplo, xlb = NULL, ylb = 'Odum', dotyp = 'observed', addtitle = T)
p2a <- apacmp_plo(toplo, xlb = 'EBASE', ylb = 'Odum', dotyp = 'detided', addtitle = T)

p <- p1a[[1]] + p1a[[2]] + p1a[[3]] + p1a[[4]] + 
  p2a[[1]] + p2a[[2]] + p2a[[3]] + p2a[[4]] + 
  plot_layout(ncol = 4, guides = 'collect') & 
  theme(legend.position = 'bottom')

png(here('figs/apacmpfig.png'), height = 4.75, width = 8, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()
