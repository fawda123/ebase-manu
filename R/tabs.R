library(flextable)
library(tidyverse)
library(lubridate)
library(here)

source(file = here('R/funcs.R'))

# EBASE i/o -----------------------------------------------------------------------------------

totab <- tibble(
  Type = c('Input', 'Input', 'Input', 'Input', 'Input', 'Input',
           'Input-derived', 'Input-derived', 'Input-derived',
           'Output', 'Output', 'Output', 'Output', 'Output', 'Output'),
  Description = c('Dissolved oxygen', 'Water temperature', 'Salinity', 'Total photosynthetically active radiation', 'Wind speed', 'Water column depth',
                  'Wind speed at 10 meter height, squared', 'Schmidt number (from water temperature and salinity)', 'Dissolved oxygen at saturation (from water temperature and salinity)',
                  'Dissolved oxygen (modelled)', 'Production', 'Respiration', 'Gas exchange', 'Light efficiency', 'b'),
  `Model notation` =  c('C_d', '-', '-', 'PAR', '-', 'H', 
                        'U102', 'Schmidt number', 'C_sat', 
                        'C_mod', 'P', 'R',  'D', 'a', 'b'),
  Units = c('mg/L', 'C', 'psu', 'W/m2/s', 'm/s', 'm', 
            'm2', 'unitless', 'mmol/m3',
            'mmol/m3', 'mmol O2 m2/d', 'mmol O2 m2/d', 'mmol O2 m2/d', '(mmol/m3/d)/(W/m2)', '(cm/hr)/(m2/s2)')
) %>% 
  flextable::as_grouped_data(groups = 'Type')

eqsz <- fp_text_default(font.size = 9)
# col2wid <- 1
# col3wid <- 2
# col14wid <- 6.5 - (col2wid + col3wid)
ebaseiotab <- flextable(totab) %>% 
  flextable::compose(i = 2, j = 3, value = as_paragraph(as_equation('C_d', props = eqsz))) %>% 
  flextable::compose(i = 5, j = 3, value = as_paragraph(as_equation('PAR', props = eqsz))) %>% 
  flextable::compose(i = 7, j = 3, value = as_paragraph(as_equation('H', props = eqsz))) %>% 
  flextable::compose(i = 9, j = 3, value = as_paragraph(as_equation('U_{10}^2', props = eqsz))) %>% 
  flextable::compose(i = 10, j = 3, value = as_paragraph(as_equation('S_c', props = eqsz))) %>% 
  flextable::compose(i = 11, j = 3, value = as_paragraph(as_equation('C_{sat}', props = eqsz))) %>% 
  flextable::compose(i = 13, j = 3, value = as_paragraph(as_equation('C_{mod}', props = eqsz))) %>% 
  flextable::compose(i = 14, j = 3, value = as_paragraph(as_equation('P,\\,aPAR', props = eqsz))) %>% 
  flextable::compose(i = 15, j = 3, value = as_paragraph(as_equation('R,\\,r', props = eqsz))) %>% 
  flextable::compose(i = 16, j = 3, value = as_paragraph(as_equation('D,\\,\\frac{1}{H}\\left[-bU_{10}^2\\left(\\frac{s_c}{600} \\right)^{-0.5} \\left(C_{Sat} - C_d\\right )\\right]', props = eqsz))) %>% 
  flextable::compose(i = 17, j = 3, value = as_paragraph(as_equation('a', props = eqsz))) %>% 
  flextable::compose(i = 18, j = 3, value = as_paragraph(as_equation('b', props = eqsz))) %>% 
  flextable::compose(i = 2, j = 4, value = as_paragraph(as_equation('\\text{mg}\\,\\text{L}^{-1}', props = eqsz))) %>% 
  flextable::compose(i = 3, j = 4, value = as_paragraph(as_equation('\\degree\\text{C}', props = eqsz))) %>% 
  flextable::compose(i = 4, j = 4, value = as_paragraph(as_equation('\\text{psu}', props = eqsz))) %>% 
  flextable::compose(i = 5, j = 4, value = as_paragraph(as_equation('\\text{Watts}\\,\\text{m}^{-2}', props = eqsz))) %>% 
  flextable::compose(i = 6, j = 4, value = as_paragraph(as_equation('\\text{m}\\,\\text{s}^{-1}', props = eqsz))) %>% 
  flextable::compose(i = 7, j = 4, value = as_paragraph(as_equation('\\text{m}', props = eqsz))) %>% 
  flextable::compose(i = 9, j = 4, value = as_paragraph(as_equation('\\text{m}^2\\,\\text{s}^{-2}', props = eqsz))) %>% 
  flextable::compose(i = 11, j = 4, value = as_paragraph(as_equation('\\text{mmol}\\,\\text{m}^{-3}', props = eqsz))) %>% 
  flextable::compose(i = 13, j = 4, value = as_paragraph(as_equation('\\text{mmol}\\,\\text{m}^{-3}', props = eqsz))) %>% 
  flextable::compose(i = 14, j = 4, value = as_paragraph(as_equation('\\text{mmol}\\,\\text{m}^{-2}\\,\\text{d}^{-1},\\,\\text{mmol}\\,\\text{m}^{-3}\\,\\text{d}^{-1}', props = eqsz))) %>% 
  flextable::compose(i = 15, j = 4, value = as_paragraph(as_equation('\\text{mmol}\\,\\text{m}^{-2}\\,\\text{d}^{-1},\\, \\text{mmol}\\,\\text{m}^{-3}\\,\\text{d}^{-1}', props = eqsz))) %>% 
  flextable::compose(i = 16, j = 4, value = as_paragraph(as_equation('\\text{mmol}\\,\\text{m}^{-2}\\,\\text{d}^{-1},\\, \\text{mmol}\\,\\text{m}^{-3}\\,\\text{d}^{-1}', props = eqsz))) %>% 
  flextable::compose(i = 17, j = 4, value = as_paragraph(as_equation('\\left(\\text{mmol}\\,\\text{m}^{-3}\\,\\text{d}^{-1}\\right)/\\left(\\text{Watts}\\,\\text{m}^{-2}\\right)', props = eqsz))) %>% 
  flextable::compose(i = 18, j = 4, value = as_paragraph(as_equation('\\left(\\text{cm}\\,\\text{hr}^{-1}\\right)/\\left(\\text{m}^2\\,\\text{s}^{-2}\\right)', props = eqsz))) %>%
  flextable::font(part = 'all', fontname = 'Times New Roman') %>% 
  flextable::align(align = 'center', j = 3:4, part = 'all') %>% 
  width(j = 1, width = 0.5) %>% 
  width(j = 2, width = 2) %>%
  width(j = 3, width = 2) %>% 
  width(j = 4, width = 2) %>% 
  padding(padding = 1, part = 'all')  

save(ebaseiotab, file = here('tabs/ebaseiotab.RData'))

# comparison of EBASE to other methods with actual apa data -----------------------------------

load(file = url('https://github.com/fawda123/BASEmetab_script/raw/master/data/apacmp.RData'))

tosum <- apacmp %>% 
  mutate(
    typ = ifelse(typ == 'BASEmetab', 'BASE', typ),
    var = factor(var, levels = c('NEM', 'P', 'R', 'D'), 
    )
  ) %>% 
  pivot_wider(names_from = 'typ', values_from = 'val')

grd <- crossing(
  dotyp = unique(tosum$dotyp), 
  var = levels(tosum$var), 
  comp = c('Odum v EBASE', 'BASE v EBASE')
  ) %>% 
  mutate(
    corv = NA, 
    int = NA,
    slo = NA, 
    slolo = NA,
    slohi = NA,
    rse = NA,
    rmse = NA
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
  rmse <- round(modelr::rmse(lmmod, tocmp), 2)
  
  # append to output
  grd[i, 'corv'] <- corv
  grd[i, 'int'] <- int
  grd[i, 'slo'] <- slo
  grd[i, 'slolo'] <- slolo
  grd[i, 'slohi'] <- slohi
  grd[i, 'rse'] <- rse
  grd[i, 'rmse'] <- rmse
  
}

totab <- grd %>% 
  select(dotyp, comp, var, corv, int, slo, rmse) %>% 
  mutate(
    dotyp = factor(dotyp, levels = c('observed', 'detided'), labels = c('Observed', 'Detided')), 
    comp = factor(comp, levels = c('Odum v EBASE', 'BASE v EBASE')), 
    var = factor(var, levels = c('NEM', 'P', 'R', 'D'))
  ) %>% 
  arrange(dotyp, comp, var) %>%
  filter(comp != 'BASE v EBASE') %>%  # remove comparison to BASE
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
    Slope = slo,
    RMSE = rmse
  ) %>% 
  select(-Comparison) %>% 
  flextable::as_grouped_data(groups = c('Dissolved Oxygen'))

apacmptab <- totab %>%
  flextable() %>% 
  fontsize(part = 'all', size = 12) %>% 
  font(part = 'all', fontname = 'Times New Roman') %>% 
  padding(padding = 1, part = 'all') %>% 
  width(width = 6.5 / ncol_keys(.)) %>% 
  flextable::compose(part = 'header', j = 'RMSE', value = as_paragraph(as_equation("\\text{RMSE}~(\\text{mmol}~\\text{O}_2/\\text{m}^2/\\text{d})"))) %>% 
  flextable::compose(part = 'header', j = 'Intercept', value = as_paragraph(as_equation("\\text{Intercept}~(\\text{mmol}~\\text{O}_2/\\text{m}^2/\\text{d})"))) %>% 
  flextable::compose(part = 'header', j = 'Corr.', value = as_paragraph(as_equation("\\rho"))) %>% 
  flextable::align(align = 'center', j = 3:6, part = 'all') %>% 
  flextable::add_footer_lines(value = '* p < 0.05, ** p < 0.005') %>% 
  font(part = 'footer', fontname = 'Times New Roman') %>% 
  flextable::align(align = 'right', part = 'footer')

save(apacmptab, file = here('tabs/apacmptab.RData'))
