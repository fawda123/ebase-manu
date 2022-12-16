library(tidyverse)
library(lubridate)
library(here)
library(patchwork)

load(file = url('https://github.com/fawda123/BASEmetab_script/raw/master/data/apacmp.RData'))

toplo <- apacmp %>% 
  mutate(
    val = case_when(
      var == 'Rt_vol' ~ -1 * val, 
      T ~ val
    ), 
    var = factor(var, 
                 levels = c('NEM', 'Pg_vol', 'Rt_vol', 'D'), 
                 labels = c('NEM', 'P', 'R', 'D')
    ), 
    seas = case_when(
      month(Date) %in% c(1:5, 10:12) ~ 'dry', 
      T ~ 'wet'
    )
  ) %>% 
  pivot_wider(names_from = 'typ', values_from = 'val')

p1a <- apacmp_plo(toplo, xlb = NULL, ylb = 'Odum', dotyp = 'observed', addtitle = T)
p1b <- apacmp_plo(toplo, xlb = NULL, ylb = 'BASEmetab', dotyp = 'observed', addtitle = F)
p2a <- apacmp_plo(toplo, xlb = NULL, ylb = 'Odum', dotyp = 'detided', addtitle = T)
p2b <- apacmp_plo(toplo, xlb = 'EBASE', ylb = 'BASEmetab', dotyp = 'detided', addtitle = F)

p <- p1a[[1]] + p1a[[2]] + p1a[[3]] + p1a[[4]] + p1b[[1]] + p1b[[2]] + p1b[[3]] + p1b[[4]] + 
  p2a[[1]] + p2a[[2]] + p2a[[3]] + p2a[[4]] + p2b[[1]] + p2b[[2]] + p2b[[3]] + p2b[[4]] + 
  plot_layout(ncol = 4, guides = 'collect') & 
  theme(legend.position = 'bottom')

png(here('figs/apacmp2.png'), height = 9, width = 8, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()