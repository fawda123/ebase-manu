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
apacmp_plo <- function(dat, xlb = 'EBASE', ylb = 'Odum', dotyp = 'observed', addtitle = T, alph = 0.5, col = 'black'){

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
      geom_point(color = col, alpha = alph, stroke = 0, size = 2) + 
      geom_abline(intercept = 0, slope = 1) +
      geom_smooth(formula = y ~ x ,method = 'lm', se = F, colour = 'tomato1', size = 0.7) +
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
priorcomp <- function(dat, ind, topbot = 3){

  met <- tibble(
    lbs = c('r2', 'rmse', 'aved'),
    lbspr = c('italic(R)^2', 'RMSE', 'Ave.\nDiff.'), 
    direc = c(-1, 1, 1)
  )
  
  toshw <- met$lbs[ind]
  leglb <- met$lbspr[ind]
  direc <- met$direc[ind]
  
  toplo <- dat %>% 
    unnest('ests') %>% 
    select(ndays, amean, asd, rmean, rsd, bmean, bsd, var, matches(toshw)) %>%
    filter(!var %in% 'b') %>% 
    pivot_wider(names_from = 'var', values_from = !!toshw) %>% 
    mutate(
      ind = rep(1: (nrow(.) / 2), times = 2), 
      ndays = case_when(
        ndays == 1 ~ paste(ndays, 'day'), 
        T ~ paste(ndays, 'days')
      )
    )
  
  
  # ave r2 summary
  optest <- dat %>% 
    mutate(
      ind = rep(1: (nrow(.) / 2), times = 2), 
      ndays = case_when(
        ndays == 1 ~ paste(ndays, 'day'), 
        T ~ paste(ndays, 'days')
      )
    ) %>% 
    unnest('ests') %>% 
    summarise(
      r2 = mean(r2), 
      .by = c('ind', 'ndays')
    ) %>% 
    mutate(
      rnkr2 = rank(-r2), 
      .by = 'ndays'
    ) 
  
  topn <- sort(unique(optest$rnkr2))[1:topbot]
  botn <- sort(unique(optest$rnkr2), decreasing = T)[1:topbot]
  
  optest <- optest %>% 
    mutate(
      fnts = ifelse(!rnkr2 %in% c(topn, botn), 'plain', ifelse(rnkr2 %in% topn, 'bold', 'italic')), 
      cols = ifelse(rnkr2 %in% c(botn, topn), 'black', '#828282'), 
      .by = 'ndays'
    )
  
  toplo1 <- toplo %>% 
    select(ind, amean, asd, rmean, rsd, bmean, bsd) %>% 
    unique() %>% 
    mutate_at(-ind, factor, labels = c('L', 'H')) %>% 
    pivot_longer(-c('ind'), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('amean', 'asd', 'rmean', 'rsd', 'bmean', 'bsd'), 
                   labels = c('italic(a)~mu', 'italic(a)~sigma', 'italic(r)~mu', 'italic(r)~sigma', 'italic(b)~mu', 'italic(b)~sigma')), 
      fac = ''
    )
  
  toplo2 <- toplo %>% 
    select(-amean, -asd, -rmean, -rsd, -bmean, -bsd) %>%
    pivot_longer(-c(ind, ndays), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('DO_mod', 'P', 'R', 'D', 'a'), 
                   labels = c('DO', 'P', 'R', 'D', 'italic(a)')
      )
    ) %>% 
    left_join(optest, by = c('ind', 'ndays'), multiple= 'all')

  p1 <- ggplot(toplo1, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    theme(
      axis.text.x = element_text(size = 9, angle = 45, hjust = 0), 
      axis.ticks = element_blank(), 
      legend.position = 'left', 
      legend.title = element_blank(),
      axis.text.y = element_blank(), 
      strip.placement = 'outside', 
      strip.background = element_blank()
    ) + 
    facet_wrap(~fac) + 
    scale_fill_brewer(palette = 'Blues') + 
    scale_x_discrete(position = 'top', expand = c(0, 0), labels = parse(text = levels(toplo1$var))) + 
    scale_y_reverse(expand = c(0, 0)) + 
    labs(
      y = NULL, 
      x = 'Prior parameters'
    )

  p2 <- ggplot(toplo2, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    geom_text(aes(y = ind, x = 5.6, label = rnkr2), size = 2.35, hjust = 0, color = toplo2$cols, fontface = toplo2$fnts) + 
    theme(
      axis.text.x = element_text(face = 'italic', size = 12), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      legend.position = 'right', 
      strip.background = element_blank(), 
      strip.placement = 'outside', 
      strip.text = element_text(vjust = 0, size = 12, hjust = 0.5, face = "plain"), 
      panel.grid = element_blank(), 
      panel.background = element_blank()
    ) + 
    facet_wrap(~ndays, ncol = 2) + 
    scale_fill_distiller(palette = 'YlOrRd', direction = direc, limits = c(0, 100)) + 
    scale_x_discrete(position = 'top', expand = c(0, 0), labels = parse(text = levels(toplo2$var))) + 
    scale_y_reverse(expand = c(0, 0)) + 
    coord_cartesian(clip = 'off') + 
    theme(panel.spacing = unit(1.25, "lines")) +
    labs(
      y = NULL, 
      fill = parse(text = leglb),
      x = 'Parameter from EBASE vs simulated,\nby optimization period'
    )
  
  out <- p1 + p2 + plot_layout(ncol = 2, widths = c(0.45, 1))
  
  return(out)
  
}

# comparison of fwoxy and ebase results for selected prior at time step of ndays
optex <- function(apagrd, fwdatcmp, asdin, rsdin, bsdin, ndaysin, subttl, ylbs = TRUE){
  
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
                   labels = c('P~(mmol~m^{3}~d^{-1})', 'R~(mmol~m^{3}~d^{-1})', 'D~(mmol~m^{3}~d^{-1})')
      )
    ) %>% 
    select(-grp)
  
  ylab <- expression(paste(O [2], ' (mmol ', m^-3, ' ', d^-1, ')'))
  
  p1 <- ggplot(toplo1, aes(x = Date, y = est, group = model, color = model)) + 
    geom_line(alpha = 1) +
    # geom_point() + 
    facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y', labeller = label_parsed) + 
    theme_minimal() + 
    theme(
      strip.placement = 'outside', 
      strip.background = element_blank(), 
      legend.position = 'top', 
      legend.title = element_blank(),
      strip.text = element_text(size = rel(0.7)), 
      axis.text.x = element_blank()
    ) + 
    labs(
      x = NULL, 
      y = NULL,
      subtitle = subttl
    )
  
  labs <- c('DO~(mmol~m^{3}~d^{-1})',
            'italic(a)~(mmol~m^{-3}~d^{-1})/(W~m^{-2})', 
            'italic(b)~(cm~hr^{-1})/(m^{2}~s^{-2})'
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
    geom_line(alpha = 1) +
    # geom_point() +
    facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y', labeller = label_parsed) + 
    theme_minimal() + 
    theme(
      strip.placement = 'outside', 
      strip.background = element_blank(), 
      legend.position = 'top', 
      legend.title = element_blank(), 
      strip.text = element_text(size = rel(0.7))
    ) + 
    labs(
      x = NULL, 
      y = NULL
    )
  
  if(!ylbs){
    p1 <- p1 + 
      theme(strip.text = element_blank()) + 
      labs(y = NULL)
    p2 <- p2 + 
      theme(strip.text = element_blank())
    
  }
    
  p1 + p2
  
}