# function for formatting p-values in tables
p_ast <- function(x){
  
  sig_cats <- c('**', '*', '')
  sig_vals <- c(-Inf, 0.005, 0.05, Inf)
  
  out <- cut(x, breaks = sig_vals, labels = sig_cats, right = FALSE)
  out <- as.character(out)
  
  return(out)
  
}

# capitalization function
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# create subplots for apalachicola obs v dtd metab comparisons
apacmp_plo <- function(dat, xlb = 'EBASE', ylb = 'Odum', dotyp = 'observed', addtitle = T, alph = 0.8, cols = c('darkgrey', 'black')){

  thm <- theme_minimal() + 
    theme(
      strip.background = element_blank(), 
      legend.position = 'top', 
      panel.grid.minor = element_blank(), 
      axis.text = element_text(size = 6)
    )
  
  names(dat)[names(dat) == ylb] <- 'yval'
  
  subfigs <- list(
    observed = '(a)', 
    detided = '(b)'
  )
  
  dat <- dat %>% 
    filter(dotyp == !!dotyp) %>% 
    select(var, seas, yval, EBASE) %>% 
    na.omit()
  
  vars <- levels(dat$var)
  for(vr in vars){
    
    toplotmp <- dat %>% 
      filter(var == vr)
    ind <- which(vr == vars)
  
    if(ind != 1)
      ylb <- NULL
    
    ttl <- NULL
    if(addtitle & ind == 1){
      subfig <- subfigs[[dotyp]]
      ttl <- paste(subfig, simpleCap(dotyp), 'dissolved oxygen')
    }
    
    lims <- range(toplotmp[, c('yval', 'EBASE')], na.rm = T)
    
    ptmp <- ggplot(toplotmp, aes(x = EBASE, y = yval)) + 
      geom_point(aes(colour = seas), alpha = alph) + 
      geom_abline(intercept = 0, slope = 1) +
      geom_smooth(method = 'lm', se = F, colour = 'tomato1', size = 0.7) +
      facet_wrap(~var, ncol = 1) +
      scale_colour_manual(values = cols) + 
      scale_x_continuous(limits = lims) + 
      scale_y_continuous(limits = lims) + 
      thm + 
      labs(
        colour = 'Season', 
        y = ylb, 
        title = ttl, 
        x = xlb
      )
    
    assign(paste0('p', ind), ptmp)
    
  }
  
  out <- list(p1, p2, p3, p4)
  
  return(out)
  
}

# Fwoxy apa comparison to EBASE for different priors sd only
priorcomp <- function(dat, ind){
  
  met <- tibble(
    lbs = c('r2', 'rmse', 'aved'),
    lbspr = c('italic(R)^2', 'RMSE', 'Ave.\nDiff.'), 
    direc = c(-1, 1, 1)
  )
  
  toshw <- met$lbs[ind]
  leglb <- met$lbspr[ind]
  direc <- met$direc[ind]
  
  toplo <- dat %>% 
    select(-amean, -rmean, -bmean) %>% 
    unnest('ests') %>% 
    select(ndays, asd, rsd, bsd, var, matches(toshw)) %>%
    filter(!var %in% 'b') %>% 
    pivot_wider(names_from = 'var', values_from = !!toshw) %>% 
    mutate(
      ind = sort(rep(1: (nrow(.) / 2), times = 2)), 
      ndays = case_when(
        ndays == 1 ~ paste(ndays, 'day'), 
        T ~ paste(ndays, 'days')
      )
    )
  
  toplo1 <- toplo %>% 
    select(ind, asd, rsd, bsd) %>% 
    unique() %>% 
    mutate(
      def = case_when(
        asd == 0.1 & rsd == 5 & bsd == 0.01 ~ '*', 
        T ~ ''
      ),
      asd = factor(asd, labels = c('L', 'M', 'H')), 
      rsd = factor(rsd, labels = c('L', 'M', 'H')), 
      bsd = factor(bsd, labels = c('L', 'M', 'H')), 
    ) %>% 
    pivot_longer(-c('ind', 'def'), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('asd', 'rsd', 'bsd'), 
                   labels = c('italic(a)', 'italic(r)', 'italic(b)'))
    ) 
  
  toplo2 <- toplo %>% 
    select(-asd, -rsd, -bsd) %>%
    pivot_longer(-c(ind, ndays), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('DO_mod', 'Pg_vol', 'Rt_vol', 'D', 'a'), 
                   labels = c('DO', 'P', 'R', 'D', 'italic(a)')
      )
    )
  
  p1 <- ggplot(toplo1, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    theme(
      axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size = 12),
      axis.ticks = element_blank(), 
      legend.position = 'left', 
      legend.title = element_blank()
    ) + 
    scale_fill_brewer(palette = 'Greys') + 
    scale_x_discrete(position = 'top', expand = c(0, 0), labels = parse(text = levels(toplo1$var))) + 
    scale_y_reverse(expand = c(0, 0), breaks = toplo1$ind, labels = toplo1$def) + 
    labs(
      y = NULL, 
      x = 'Variance of prior',
      caption = '* EBASE default'
    )
  
  p2 <- ggplot(toplo2, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    theme(
      axis.text.x = element_text(face = 'italic', size = 12), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      legend.position = 'right', 
      strip.background = element_blank(), 
      strip.text = element_text(hjust = 0, size = 12, face = 'bold')
    ) + 
    facet_wrap(~ndays, ncol = 2) + 
    scale_fill_distiller(palette = 'YlOrRd', direction = direc, limits = c(0, 100)) + 
    scale_x_discrete(position = 'top', expand = c(0, 0), labels = parse(text = levels(toplo2$var))) + 
    scale_y_reverse(expand = c(0, 0)) + 
    labs(
      y = NULL, 
      fill = parse(text = leglb),
      x = 'Parameter from EBASE vs simulated,\nby optimization period'
    )
  
  out <- p1 + p2 + plot_layout(ncol = 2, widths = c(0.3, 1))
  
  return(out)
  
}

# comparison of fwoxy and ebase results for selected prior at time step of ndays
optex <- function(apagrd, fwdatcmp, asdin, rsdin, bsdin, ndaysin){
  
  res <- apagrd %>% 
    filter(
      asd == asdin & rsd == rsdin & bsd == bsdin & ndays == ndaysin
    ) %>% 
    pull(out) %>% 
    .[[1]] %>% 
    .[[1]]
  cmp <- inner_join(fwdatcmp, res, by = c('Date', 'DateTimeStamp')) %>%
    select(-converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>%
    rename(
      DO_mod.x = DO_obs.x,
      DO_mod.y = DO_mod
    ) %>%
    pivot_longer(!all_of(c('DateTimeStamp', 'Date', 'grp')), names_to = 'var', values_to = 'val') %>%
    separate(var, c('var', 'mod'), sep = '\\.') %>%
    mutate(
      mod = case_when(
        mod == 'x' ~ 'Fwoxy',
        mod == 'y' ~ 'EBASE'
      )
    ) %>%
    pivot_wider(names_from = 'mod', values_from = 'val')
  
  toplo1 <- cmp %>% 
    filter(var %in% c('Pg_vol', 'Rt_vol', 'D')) %>% 
    group_by(grp, var) %>% 
    summarise(
      Fwoxy = mean(Fwoxy, na.rm = T), 
      EBASE = mean(EBASE, na.rm = T),
      Date  = min(Date),
      .groups = 'drop'
    ) %>% 
    pivot_longer(-c(Date, grp, var), names_to = 'model', values_to = 'est') %>% 
    mutate(
      var = factor(var, 
                   levels = c('Pg_vol', 'Rt_vol', 'D'), 
                   labels = c('Pg [vol]', 'Rt [vol]', 'D')
      )
    ) %>% 
    select(-grp)
  
  ylab <- expression(paste(O [2], ' (mmol ', m^-3, ' ', d^-1, ')'))
  
  p1 <- ggplot(toplo1, aes(x = Date, y = est, group = model, color = model)) + 
    geom_line() +
    geom_point() + 
    facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y', labeller = label_parsed) + 
    theme_minimal() + 
    theme(
      strip.placement = 'outside', 
      strip.background = element_blank(), 
      legend.position = 'top', 
      legend.title = element_blank(),
      strip.text = element_text(size = rel(1))
    ) + 
    labs(
      x = NULL, 
      y = ylab
    )
  
  labs <- c('DO[mod]~(mmol~O[2]~m^{3}~d^{-1})',
            'a~(mmol~m^{-3}~d^{-1})/(W~m^{-2})', 
            'b~(cm~hr^{-1})/(m^{2}~s^{-2})'
  )
  
  toplo2 <- cmp %>% 
    filter(var %in% c('DO_mod', 'a', 'b')) %>% 
    group_by(grp, var) %>%
    summarise(
      Fwoxy = mean(Fwoxy, na.rm = T),
      EBASE = mean(EBASE, na.rm = T),
      Date = min(Date),
      .groups = 'drop'
    ) %>%
    pivot_longer(-c(Date, grp, var), names_to = 'model', values_to = 'est') %>% 
    mutate(
      var = factor(var, 
                   levels = c('DO_mod', 'a', 'b'), 
                   labels = labs
      )
    ) %>% 
    select(-grp)
  
  p2 <- ggplot(toplo2, aes(x = Date, y = est, group = model, color = model)) + 
    geom_line() +
    geom_point() +
    facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y', labeller = label_parsed) + 
    theme_minimal() + 
    theme(
      strip.placement = 'outside', 
      strip.background = element_blank(), 
      legend.position = 'top', 
      legend.title = element_blank(), 
      strip.text = element_text(size = rel(1))
    ) + 
    labs(
      x = NULL, 
      y = NULL
    )
  
  p1 + p2 + plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = 'top')
  
}
