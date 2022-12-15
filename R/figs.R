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

toplo1 <- toplo %>% 
  filter(dotyp == 'observed')
toplo1a <- toplo1 %>% 
  select(-BASEmetab)
toplo1b <- toplo1 %>% 
  select(-Odum)

toplo2 <- toplo %>% 
  filter(dotyp == 'detided') 
toplo2a <- toplo2 %>% 
  select(-BASEmetab)
toplo2b <- toplo2 %>% 
  select(-Odum)

thm <- theme_minimal() + 
  theme(
    strip.background = element_blank(), 
    legend.position = 'top', 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 6)
  )

alph <- 0.8
cols <- c('darkgrey', 'black')

# obs do, odum v ebase
vars <- levels(toplo1a$var)
for(vr in vars){
  
  toplotmp <- toplo1a %>% 
    filter(var == vr)
  ind <- which(vr == vars)
  
  corv <- cor.test(toplotmp$Odum, toplotmp$EBASE) %>% 
    .$estimate %>% 
    round(., 2) %>% 
    format(., nsmall = 2)
  
  toplotmp <- toplotmp %>% 
    mutate(
      var = paste0(var, ' (', corv, ')')
    )
  
  ylb <- NULL
  ttl <- NULL
  if(ind == 1){
    ylb <- 'Odum'
    ttl <- '(a) Observed dissolved oxygen'
  }
  
  lims <- range(toplotmp[, c('Odum', 'EBASE')], na.rm = T)
  
  ptmp <- ggplot(toplotmp, aes(x = EBASE, y = Odum, colour = seas)) + 
    geom_point(alpha = alph) + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~var, ncol = 1) +
    scale_colour_manual(values = cols) + 
    scale_x_continuous(limits = lims) + 
    scale_y_continuous(limits = lims) + 
    thm + 
    labs(
      colour = 'Season', 
      y = ylb, 
      title = ttl, 
      x = NULL
    )

  assign(paste0('p1a', ind), ptmp)
  
}

# obs do, BASEmetab v ebase
vars <- levels(toplo1b$var)
for(vr in vars){
  
  toplotmp <- toplo1b %>% 
    filter(var == vr)
  ind <- which(vr == vars)
  
  corv <- cor.test(toplotmp$BASEmetab, toplotmp$EBASE) %>% 
    .$estimate %>% 
    round(., 2) %>% 
    format(., nsmall = 2)
  
  toplotmp <- toplotmp %>% 
    mutate(
      var = paste0(var, ' (', corv, ')')
    )
  
  ylb <- NULL
  if(ind == 1)
    ylb <- 'BASEmetab'
  
  lims <- range(toplotmp[, c('BASEmetab', 'EBASE')], na.rm = T)
  
  ptmp <- ggplot(toplotmp, aes(x = EBASE, y = BASEmetab, colour = seas)) + 
    geom_point(alpha = alph) + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~var, ncol = 1) +
    scale_colour_manual(values = cols) + 
    scale_x_continuous(limits = lims) + 
    scale_y_continuous(limits = lims) + 
    thm + 
    labs(
      colour = 'Season', 
      y = ylb, 
      x = NULL
    )
  
  assign(paste0('p1b', ind), ptmp)
  
}

# dtd do, odum v ebase
vars <- levels(toplo2a$var)
for(vr in vars){
  
  toplotmp <- toplo2a %>% 
    filter(var == vr)
  ind <- which(vr == vars)
  
  corv <- cor.test(toplotmp$Odum, toplotmp$EBASE) %>% 
    .$estimate %>% 
    round(., 2) %>% 
    format(., nsmall = 2)
  
  toplotmp <- toplotmp %>% 
    mutate(
      var = paste0(var, ' (', corv, ')')
    )
  
  ylb <- NULL
  ttl <- NULL
  if(ind == 1){
    ylb <- 'Odum'
    ttl <- '(a) Detided dissolved oxygen'
  }
  
  lims <- range(toplotmp[, c('Odum', 'EBASE')], na.rm = T)
  
  ptmp <- ggplot(toplotmp, aes(x = EBASE, y = Odum, colour = seas)) + 
    geom_point(alpha = alph) + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~var, ncol = 1) +
    scale_colour_manual(values = cols) + 
    scale_x_continuous(limits = lims) + 
    scale_y_continuous(limits = lims) + 
    thm + 
    labs(
      colour = 'Season', 
      y = ylb, 
      title = ttl, 
      x = NULL
    )
  
  assign(paste0('p2a', ind), ptmp)
  
}

# dtd do, BASEmetab v ebase
vars <- levels(toplo2b$var)
for(vr in vars){
  
  toplotmp <- toplo2b %>% 
    filter(var == vr)
  ind <- which(vr == vars)
  
  corv <- cor.test(toplotmp$BASEmetab, toplotmp$EBASE) %>% 
    .$estimate %>% 
    round(., 2) %>% 
    format(., nsmall = 2)
  
  toplotmp <- toplotmp %>% 
    mutate(
      var = paste0(var, ' (', corv, ')')
    )
  
  ylb <- NULL
  if(ind == 1)
    ylb <- 'BASEmetab'
  
  lims <- range(toplotmp[, c('BASEmetab', 'EBASE')], na.rm = T)
  
  ptmp <- ggplot(toplotmp, aes(x = EBASE, y = BASEmetab, colour = seas)) + 
    geom_point(alpha = alph) + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~var, ncol = 1) +
    scale_colour_manual(values = cols) + 
    scale_x_continuous(limits = lims) + 
    scale_y_continuous(limits = lims) + 
    thm + 
    labs(
      colour = 'Season', 
      y = ylb, 
      x = 'EBASE'
    )
  
  assign(paste0('p2b', ind), ptmp)
  
}

p <- p1a1 + p1a2 + p1a3 + p1a4 + p1b1 + p1b2 + p1b3 + p1b4 + 
  p2a1 + p2a2 + p2a3 + p2a4 + p2b1 + p2b2 + p2b3 + p2b4 + 
  plot_layout(ncol = 4, guides = 'collect') & 
  theme(legend.position = 'bottom')

png(here('figs/apacmp.png'), height = 9, width = 8, family = 'serif', units = 'in', res = 500)
print(p)
dev.off()