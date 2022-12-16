library(flextable)
library(tidyverse)
library(here)

source(here('R/funcs.R'))

load(file = here('data/apacmp.RData'))


tosum <- apacmp %>% 
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

grd <- crossing(
  dotyp = unique(tosum$dotyp), 
  var = levels(tosum$var), 
  seas = c('dry', 'wet', 'all'), 
  comp = c('Odum v EBASE', 'BASEmetab v EBASE')
  ) %>% 
  mutate(
    corv = NA, 
    int = NA,
    slo = NA, 
    rse = NA
  )

for(i in 1:nrow(grd)){
  
  cat(i, '\t')
  
  tosel <- grd[i, ]
               
  dotyp <- grd[i, ]$dotyp
  var <- grd[i, ]$var
  seas <- grd[i, ]$seas
  comp <- grd[i, ]$comp
  
  if(seas == 'all')
    seas <- c('dry', 'wet')
  
  comp <- strsplit(comp, ' v ')[[1]]
  
  # filter data
  tocmp <- tosum %>% 
    filter(dotyp == !!dotyp) %>% 
    filter(var == !!var) %>% 
    filter(seas %in% !!seas) %>% 
    rename(
      yval = !!comp[1], 
      xval = !!comp[2]
    )
  
  # get summaries
  corv <- cor.test(tocmp$xval, tocmp$yval)
  corv <- paste0(round(corv$estimate, 2), p_ast(corv$p.value))
  
  lmmod <- lm(yval ~ xval, data = tocmp)
  int <- paste0(round(summary(lmmod)$coefficients[1, 1], 2), p_ast(summary(lmmod)$coefficients[1, 4]))
  slo <- paste0(round(summary(lmmod)$coefficients[2, 1], 2), p_ast(summary(lmmod)$coefficients[2, 4])) # test if different from one
  rse <- round(sigma(lmmod), 2)
  
  # append to output
  grd[i, 'corv'] <- corv
  grd[i, 'int'] <- int
  grd[i, 'slo'] <- slo
  grd[i, 'rse'] <- rse
  
}

apacmptab <- grd

save(apacmptab, file = here('tabs/apacmptab.RData'))
