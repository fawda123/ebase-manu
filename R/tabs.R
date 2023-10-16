library(tidyverse)
library(lubridate)
library(here)

source(file = here('R/funcs.R'))

# comparison of EBASE to other methods with actual apa data -----------------------------------

load(file = url('https://github.com/fawda123/BASEmetab_script/raw/master/data/apacmp.RData'))

tosum <- apacmp %>% 
  filter(var != 'D') %>% 
  mutate(
    var = factor(var, levels = c('NEM', 'P', 'R'), 
    )
  ) %>% 
  pivot_wider(names_from = 'typ', values_from = 'val') %>% 
  na.omit()

grd <- crossing(
  dotyp = unique(tosum$dotyp), 
  var = levels(tosum$var)
  ) %>% 
  mutate(
    corv = NA, 
    int = NA,
    slo = NA, 
    slolo = NA,
    slohi = NA,
    rse = NA,
    rmsd = NA
  )

for(i in 1:nrow(grd)){
  
  cat(i, '\t')
  
  tosel <- grd[i, ]
               
  dotyp <- grd[i, ]$dotyp
  var <- grd[i, ]$var
  
  # filter data
  tocmp <- tosum %>% 
    filter(dotyp == !!dotyp) %>% 
    filter(var == !!var) %>% 
    rename(
      yval = 'Odum', 
      xval = 'EBASE'
      )
  
  # get summaries
  corv <- cor.test(tocmp$xval, tocmp$yval)
  corv <- paste0(round(corv$estimate, 2), p_ast(corv$p.value))
 
  lmmod <- lm(yval ~ xval, data = tocmp)
  lmmodsum <- summary(lmmod)
  tval <- (abs(lmmodsum$coefficients[2, 1] - 1)) / (lmmodsum$coefficients[2, 2]) 
  pval <- pt(tval, df = lmmodsum$df[2], lower.tail = F) * 2
  int <- paste0(round(lmmodsum$coefficients[1, 1], 2), p_ast(lmmodsum$coefficients[1, 4]))
  slo <- paste0(round(lmmodsum$coefficients[2, 1], 2), p_ast(pval)) # test if different from one
  slolo <- confint(lmmod)[2, 1]
  slohi <- confint(lmmod)[2, 2]
  rse <- round(sigma(lmmod), 2)
  rmsd <- as.character(round(sqrt(sum((tocmp$yval - tocmp$xval)^2) / nrow(tocmp)), 2))

  # append to output
  grd[i, 'corv'] <- corv
  grd[i, 'int'] <- int
  grd[i, 'slo'] <- slo
  grd[i, 'slolo'] <- slolo
  grd[i, 'slohi'] <- slohi
  grd[i, 'rse'] <- rse
  grd[i, 'rmsd'] <- rmsd
  
}

totab <- grd %>% 
  select(dotyp, var, corv, rmsd) %>% 
  pivot_longer(corv:rmsd) %>% 
  unite('name', name, dotyp) %>% 
  mutate(
    name = factor(name, levels = c('corv_observed', 'corv_detided', 'rmsd_observed', 'rmsd_detided')),
    var = factor(var, levels = c('P', 'R', 'NEM')), 
    value = gsub('\\*+', '', value)
  ) %>% 
  pivot_wider(names_from = 'name', values_from = 'value') %>% 
  unite('Corr.', corv_observed, corv_detided, sep = ' / ') %>% 
  unite('RMSD', rmsd_observed, rmsd_detided, sep = ' / ') %>% 
  arrange(var) %>%
  rename(
    Parameter = var
  ) 

apacmptab <- totab %>%
  knitr::kable(
    col.names = c('Parameter', '$\\rho$', 'RMSD (mmol O$^2$/m$^2$/d)')
  )

save(apacmptab, file = here('tabs/apacmptab.RData'))
