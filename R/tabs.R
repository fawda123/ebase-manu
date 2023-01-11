library(flextable)
library(tidyverse)
library(lubridate)
library(here)

source(file = here('R/funcs.R'))

# comparison of EBASE to other methods with actual apa data -----------------------------------

load(file = url('https://github.com/fawda123/BASEmetab_script/raw/master/data/apacmp.RData'))

tosum <- apacmp %>% 
  mutate(
    val = case_when(
      var == 'Rt_vol' ~ -1 * val, 
      T ~ val
    ), 
    var = factor(var, 
                 levels = c('NEM', 'Pg_vol', 'Rt_vol', 'D'), 
                 labels = c('NEM', 'P', 'R', 'D')
    )
  ) %>% 
  pivot_wider(names_from = 'typ', values_from = 'val')

grd <- crossing(
  dotyp = unique(tosum$dotyp), 
  var = levels(tosum$var), 
  comp = c('Odum v EBASE', 'BASEmetab v EBASE')
  ) %>% 
  mutate(
    corv = NA, 
    int = NA,
    slo = NA, 
    slolo = NA,
    slohi = NA,
    rse = NA
  )

for(i in 1:nrow(grd)){
  
  cat(i, '\t')
  
  tosel <- grd[i, ]
               
  dotyp <- grd[i, ]$dotyp
  var <- grd[i, ]$var
  comp <- grd[i, ]$comp
  
  comp <- strsplit(comp, ' v ')[[1]]
  
  # filter data
  tocmp <- tosum %>% 
    filter(dotyp == !!dotyp) %>% 
    filter(var == !!var) %>% 
    rename(
      yval = !!comp[1], 
      xval = !!comp[2]
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
  
  # append to output
  grd[i, 'corv'] <- corv
  grd[i, 'int'] <- int
  grd[i, 'slo'] <- slo
  grd[i, 'slolo'] <- slolo
  grd[i, 'slohi'] <- slohi
  grd[i, 'rse'] <- rse
  
}

totab <- grd %>% 
  select(dotyp, comp, var, corv, int, slo, rse) %>% 
  mutate(
    dotyp = factor(dotyp, levels = c('observed', 'detided'), labels = c('Observed', 'Detided')), 
    comp = factor(comp, levels = c('Odum v EBASE', 'BASEmetab v EBASE')), 
    var = factor(var, levels = c('NEM', 'P', 'R', 'D'))
  ) %>% 
  arrange(dotyp, comp, var) %>%
  group_by(dotyp) %>% 
  mutate(
    comp = ifelse(duplicated(comp), '', as.character(comp))
  ) %>% 
  rename(
    `Dissolved Oxygen` = dotyp, 
    `Comparison` = comp,
    Estimate = var, 
    `Corr.` = corv,
    Intercept = int, 
    Slope = slo
  ) %>% 
  flextable::as_grouped_data(groups = c('Dissolved Oxygen'))

apacmptab <- totab %>%
  flextable() %>% 
  fontsize(part = 'all', size = 12) %>% 
  font(part = 'all', fontname = 'Times New Roman') %>% 
  padding(padding = 1, part = 'all') %>% 
  width(width = 6.5 / ncol_keys(.)) %>% 
  flextable::compose(part = 'header', j = 'rse', value = as_paragraph(as_equation("\\sigma"))) %>% 
  flextable::compose(part = 'header', j = 'Corr.', value = as_paragraph(as_equation("\\rho"))) %>% 
  flextable::align(align = 'center', j = 4:7, part = 'all') %>% 
  flextable::add_footer_lines(value = '* p < 0.05, ** p < 0.005') %>% 
  font(part = 'footer', fontname = 'Times New Roman') %>% 
  flextable::align(align = 'right', part = 'footer')

save(apacmptab, file = here('tabs/apacmptab.RData'))
