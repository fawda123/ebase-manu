library(tidyverse)
library(lubridate)
library(here)
library(patchwork)
library(EBASE)
library(readr)
library(scales)

source(here('R/funcs.R'))

# prior distribution plot ---------------------------------------------------------------------

p <- prior_plot(n = 1e5)

png(here('figs/priorplot.png'), height = 3, width = 8, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# simulated apa data --------------------------------------------------------------------------

fwdat <- read_csv(file = url('https://raw.githubusercontent.com/fawda123/BASEmetab_script/master/data/apafwoxy.csv'))

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
  ggplot2::scale_color_discrete(
    breaks = c('P', 'R', 'D')
  ) +
  scale_x_date(breaks = "1 month", labels = date_format("%b")) +
  labs(title = '(c) Synthetic metabolic estimates')

p <- p1 + p2 + p3 + plot_layout(ncol = 1) & theme(panel.grid.minor = element_blank())

png(here('figs/synapa.png'), height = 8, width = 8, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# Fwoxy apa comparison ------------------------------------------------------------------------

fl <- paste0(tempdir(), '/apasumdat.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apasumdat.RData', destfile = fl)
load(file = fl)

p1 <- priorcomp(apasumdat, 1)

png(here('figs/apasumdat.png'), height = 6, width = 6, family = 'serif', units = 'in', res = 500)
print(p1)
dev.off()

# best fwoxy apa comparison -------------------------------------------------------------------

fl <- paste0(tempdir(), '/apagrd.RData')
download.file('https://github.com/fawda123/BASEmetab_script/raw/master/data/apagrd.RData', destfile = fl)
load(file = fl)

fwdat <- read_csv(file = url('https://raw.githubusercontent.com/fawda123/BASEmetab_script/master/data/apafwoxy.csv'))

fwdatcmp <- fwdat %>% 
  mutate(
    DateTimeStamp = dmy_hms(datet, tz = 'America/Jamaica'),
    Date = as.Date(DateTimeStamp, tz = 'America/Jamaica'),
    DO_obs = `oxy,mmol/m3`, 
    a = `aparam,(mmolO2/m2/d)/(W/m2)` / `ht,m`,
    Rt_vol = `er,mmol/m2/d` / `ht,m`,
    Pg_vol = `gpp,mmol/m2/d` / `ht,m`,
    D = -1 * `gasex,mmol/m2/d` / `ht,m`,
    b = 100 * 3600 * `kw,m/s` / `wspd2,m2/s2` / (`sc,dimensionless` / 660) ^ -0.5 # (m/s)/(m2/s2) to (cm/hr) / (m2/s2)
  ) %>% 
  select(Date, DateTimeStamp, DO_obs, a, b, Pg_vol, Rt_vol, D)

subttl <- expression(paste('(a) ndays = 1, ', italic(a), ' (sd) = 1, ', italic(r), ' (sd) = 50, ', italic(b), ' (sd) = 0.001'))
p1 <- optex(apagrd, fwdatcmp, asdin = 1, rsdin = 50, bsdin = 0.001, ndaysin = 1, subttl = subttl)

subttl <- expression(paste('(b) ndays = 7, ', italic(a), ' (sd) = 0.01, ', italic(r), ' (sd) = 50, ', italic(b), ' (sd) = 0.001'))
p2 <- optex(apagrd, fwdatcmp, asdin = 0.01, rsdin = 50, bsdin = 0.001, ndaysin = 7, subttl = subttl, ylbs = F)

p <- ((p1 + plot_layout(ncol = 1)) | (p2 + plot_layout(ncol = 1)))  + plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = 'top')

png(here('figs/optex.png'), height = 8.75, width = 8, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()

# actual apa comparison -----------------------------------------------------------------------

load(file = url('https://github.com/fawda123/BASEmetab_script/raw/master/data/apacmp.RData'))

toplo <- apacmp %>% 
  mutate(
    typ = ifelse(typ == 'BASEmetab', 'BASE', typ),
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
p1b <- apacmp_plo(toplo, xlb = NULL, ylb = 'BASE', dotyp = 'observed', addtitle = F)
p2a <- apacmp_plo(toplo, xlb = NULL, ylb = 'Odum', dotyp = 'detided', addtitle = T)
p2b <- apacmp_plo(toplo, xlb = 'EBASE', ylb = 'BASE', dotyp = 'detided', addtitle = F)

p <- p1a[[1]] + p1a[[2]] + p1a[[3]] + p1a[[4]] + p1b[[1]] + p1b[[2]] + p1b[[3]] + p1b[[4]] + 
  p2a[[1]] + p2a[[2]] + p2a[[3]] + p2a[[4]] + p2b[[1]] + p2b[[2]] + p2b[[3]] + p2b[[4]] + 
  plot_layout(ncol = 4, guides = 'collect') & 
  theme(legend.position = 'bottom')

png(here('figs/apacmpfig.png'), height = 9, width = 8, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()
