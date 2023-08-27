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
priorcomp <- function(dat, met, topbot = 3){

  metsel <- tibble(
      met = c('r2', 'rmse', 'mape', 'mae', 'nse'),
      lbspr = c('italic(R)^2', 'RMSE', 'MAPE', 'mae', 'NSE'), 
      direc = c(-1, 1, 1, 1, -1),
      limts = list(c(0, 100), NULL, NULL, NULL, NULL),
      trans = c('identity', 'log10', 'log10', 'log10', 'identity')
    ) %>% 
    filter(met == !!met)
 
  toshw <- metsel$met
  leglb <- metsel$lbspr
  direc <- metsel$direc
  limts <- metsel$limts[[1]]
  trans <- metsel$trans
  
  toplo <- dat %>% 
    unnest('ests') %>% 
    select(ndays, amean, asd, rmean, rsd, var, matches(toshw)) %>%
    filter(!var %in% 'b') %>% 
    pivot_wider(names_from = 'var', values_from = !!toshw) %>% 
    mutate(
      ind = rep(1: (nrow(.) / 3), times = 3), 
      ndays = case_when(
        ndays == 1 ~ paste(ndays, 'day'), 
        T ~ paste(ndays, 'days')
      )
    )
  
  # ave summary
  optest <- metsum_fun(dat, met) %>% 
    mutate(
      ndays = case_when(
        ndays == 1 ~ paste(ndays, 'day'), 
        T ~ paste(ndays, 'days')
      ), 
      ndays = factor(ndays, levels = unique(ndays))
    )
  
  topn <- sort(unique(optest$rnkmetsum))[1:topbot]
  botn <- sort(unique(optest$rnkmetsum), decreasing = T)[1:topbot]
  
  optest <- optest %>% 
    mutate(
      fnts = ifelse(!rnkmetsum %in% c(topn, botn), 'plain', ifelse(rnkmetsum %in% topn, 'bold', 'italic')), 
      cols = ifelse(rnkmetsum %in% c(botn, topn), 'black', '#828282'), 
      .by = 'ndays'
    )
  
  toplo1 <- toplo %>% 
    select(ind, amean, asd, rmean, rsd) %>% 
    unique() %>% 
    mutate_at(vars(-matches('ind')), factor, labels = c('L', 'H')) %>% 
    pivot_longer(-c('ind'), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('amean', 'asd', 'rmean', 'rsd'), 
                   labels = c('italic(a)~mu', 'italic(a)~sigma', 'italic(R)~mu', 'italic(R)~sigma')), 
      fac = ''
    )
  
  toplo2 <- toplo %>% 
    select(-amean, -asd, -rmean, -rsd) %>%
    pivot_longer(-c(ind, ndays), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('DO_mod', 'D', 'P', 'R', 'a'), 
                   labels = c('DO', 'D', 'P', 'R', 'italic(a)')
      ), 
      ndays = factor(ndays, levels = unique(ndays))
    ) %>% 
    left_join(optest, by = c('ind', 'ndays'), multiple= 'all')

  p1 <- ggplot(toplo1, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    theme(
      axis.text.x = element_text(size = 9, angle = 45, hjust = 0), 
      axis.ticks = element_blank(), 
      legend.position = 'left', 
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
      x = 'Prior parameters', 
      fill = 'Prior\nvalue'
    )

  # if(met == 'nse')
  #   toplo2 <- toplo2 %>% 
  #     mutate(
  #       val = pmax(val, -2)
  #     )
  
  p2 <- ggplot(toplo2, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    geom_text(aes(y = ind, x = 5.6, label = rnkmetsum), size = 2.35, hjust = 0, color = toplo2$cols, fontface = toplo2$fnts) + 
    theme(
      axis.text.x = element_text(face = 'italic', size = 8), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      legend.position = 'right', 
      strip.background = element_blank(), 
      strip.placement = 'outside', 
      strip.text = element_text(vjust = 0, size = 12, hjust = 0.5, face = "plain"), 
      panel.grid = element_blank(), 
      panel.background = element_blank()
    ) + 
    facet_wrap(~ndays, ncol = 3) + 
    scale_x_discrete(position = 'top', expand = c(0, 0), labels = parse(text = levels(toplo2$var))) + 
    scale_y_reverse(expand = c(0, 0)) + 
    coord_cartesian(clip = 'off') + 
    theme(panel.spacing = unit(1.25, "lines")) +
    labs(
      y = NULL, 
      fill = parse(text = leglb),
      x = 'Parameter from EBASE vs synthetic,\nby optimization period'
    )
  
  # special color palette for nse
  if(met != 'nse')
    p2 <- p2 + 
      scale_fill_distiller(palette = 'YlOrRd', direction = direc, limits = limts) 
    
  if(met == 'nse')
    p2 <- p2 + 
      scale_fill_gradient2(low = 'red', mid = 'white', high = 'green', midpoint = 1, trans = 'exp' ,
                           breaks = c(-500, -1, 0, 0.5, 0.9))
    
  out <- p1 + p2 + plot_layout(ncol = 2, widths = c(0.35, 1.2))
  
  return(out)
  
}

# Fwoxy apa comparison to EBASE for different priors, rmse summary
priorsumcomp <- function(dat, met = 'r2'){

  toplo <- dat %>% 
    unnest(ests) %>% 
    select(amean, asd, rmean, rsd, ndays, ind, var, !!met) %>% 
    rename(met = !!met) %>% 
    mutate_at(vars(matches('mean$|sd$')), factor, labels = c('L', 'H')) %>% 
    mutate(
      ndays = case_when(
        ndays == 1 ~ paste(ndays, 'day'), 
        T ~ paste(ndays, 'days')
      )
    ) %>% 
    filter(!var %in% 'DO_mod') %>% 
    pivot_longer(cols = matches('mean$|sd$'), names_to = 'prior', values_to = 'val') %>% 
    summarise(
      medv = median(met, na.rm = T), 
      iqr = IQR(met, na.rm = T),
      .by = c('ndays', 'prior', 'val')
    ) %>% 
    mutate(
      grandmed = median(medv), 
      .by = c('ndays', 'val')
    ) %>% 
    mutate(
           medv = 100 * (medv - grandmed) / ((grandmed + medv) / 2),
           param = gsub('mean$|sd$', '', prior),
           param = factor(param, levels = c('a', 'r'), labels = c('a', 'R')),
           prior = gsub('^a|^r|^b', '', prior), 
           prior = factor(prior, levels = c('mean', 'sd'), labels = c('mu', 'sigma')), 
           ndays = factor(ndays, levels = c('1 day', '7 days', '30 days'), labels = c('1~day', '7~days', '30~days'))
    )
  
  metsel <- tibble(
    met = c('r2', 'rmse', 'mape', 'mae', 'nse'),
    lbspr = c('italic(R)^2', 'RMSE', 'MAPE', 'mae', 'NSE')
  ) %>% 
    filter(met == !!met)
  ylab <- parse(text = paste0('`%`', '~change~', metsel$lbspr))
  
  wd <- 0.3
  
  p <- ggplot(toplo, aes(x = param, group = val)) + 
    geom_linerange(aes(ymin = 0, ymax = medv, x = param), position = position_dodge(width = wd), linetype = 'dashed') +
    geom_hline(yintercept = 0) + 
    geom_point(aes(y = medv, fill = val, size = iqr), pch = 21, position = position_dodge(width = wd)) +
    facet_grid(prior~ndays, labeller = label_parsed) + 
    scale_fill_brewer(palette = 'Blues') + 
    scale_size(range = c(1, 10)) + 
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    theme_bw() + 
    theme(
      axis.text.x = element_text(size = 14, face = 'italic'), 
      strip.text.y = element_text(size = 14, face = 'italic'), 
      strip.text.x = element_text(size = 12),
      strip.background = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = NULL, 
      y = ylab, 
      fill = 'Prior\nvalue', 
      size = 'IQR'
    )
  
  return(p)
  
}

# plot for comparing Fwoxy and EBASE, used in optext and defex
cmpplo <- function(res, fwdatcmp, subttl, ylbs = TRUE){
  
  cmp <- inner_join(fwdatcmp, res, by = c('Date', 'DateTimeStamp')) %>%
    mutate(
      a.x = a.x / H, # fwoxy is m-2, ebase is m-3
    ) %>% 
    select(-H, -converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>%
    rename(
      DO_mod.x = DO_obs.x,
      DO_mod.y = DO_mod
    ) %>%
    pivot_longer(!all_of(c('DateTimeStamp', 'Date', 'grp')), names_to = 'var', values_to = 'val') %>%
    separate(var, c('var', 'mod'), sep = '\\.') %>%
    mutate(
      mod = case_when(
        mod == 'x' ~ 'Synthetic',
        mod == 'y' ~ 'EBASE'
      )
    ) %>%
    pivot_wider(names_from = 'mod', values_from = 'val')
  
  toplo1 <- cmp %>% 
    filter(var %in% c('P', 'R', 'D')) %>% 
    group_by(grp, var) %>% 
    summarise(
      Synthetic = mean(Synthetic, na.rm = T), 
      EBASE = mean(EBASE, na.rm = T),
      Date  = min(Date),
      .groups = 'drop'
    ) %>% 
    pivot_longer(-c(Date, grp, var), names_to = 'model', values_to = 'est') %>% 
    mutate(
      var = factor(var, 
                   levels = c('P', 'R', 'D'), 
                   labels = c('P~(mmol~m^{2}~d^{-1})', 'R~(mmol~m^{2}~d^{-1})', 'D~(mmol~m^{2}~d^{-1})')
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
  
  labs <- c('O[2]~(mmol~m^{-3})',
            'italic(a)~(mmol~m^{-3}~d^{-1})/(W~m^{-2})', 
            'italic(b)~(cm~hr^{-1})/(m^{2}~s^{-2})'
  )
  
  toplo2 <- cmp %>% 
    filter(var %in% c('DO_mod', 'a', 'b')) %>% 
    group_by(grp, var) %>%
    summarise(
      Synthetic = mean(Synthetic, na.rm = T),
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
      axis.text.x = element_text(size = 8),
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

# comparison of fwoxy and ebase results for ranked model by prior at opt of ndays
optex <- function(apagrd, fwdatcmp, apasumdat, rnkmetsum, met, lims = NULL){

  if(!is.null(lims))
    lims <- lapply(lims, function(x) tibble(min = x[1], max = x[2])) %>% 
      tibble::enframe(name = 'var') %>% 
      unnest('value')
  
  # find the prior comparison
  rnkfnd <- metsum_fun(apasumdat, met = met, parms = T) %>% 
    filter(rnkmetsum %in% !!rnkmetsum) %>% 
    select(-ind, -metsum)
  
  res <- apagrd %>% 
    inner_join(rnkfnd, by = c("amean", "asd", "rmean", "rsd", "ndays")) %>% 
    arrange(ndays, rnkmetsum) %>% 
    mutate(
      subfig = paste0('(', letters[1:nrow(.)], ')'), 
      subfig = case_when(
        rnkmetsum == !!rnkmetsum[1] ~ paste(subfig, 'Best', ndays, 'day', sep = '~'), 
        rnkmetsum == !!rnkmetsum[2] ~ paste(subfig, 'Worst', ndays, 'day', sep = '~') 
      ), 
      subfig = factor(subfig)
    ) %>% 
    select(subfig, out) %>% 
    unnest(out) %>% 
    unnest(out)
  
  cmp <- inner_join(fwdatcmp, res, by = c('Date', 'DateTimeStamp')) %>%
    select(-H, -converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>%
    rename(
      DO_mod.x = DO_obs.x,
      DO_mod.y = DO_mod
    ) %>%
    pivot_longer(!all_of(c('subfig', 'DateTimeStamp', 'Date', 'grp')), names_to = 'var', values_to = 'val') %>%
    separate(var, c('var', 'mod'), sep = '\\.') %>%
    mutate(
      mod = case_when(
        mod == 'x' ~ 'Synthetic',
        mod == 'y' ~ 'EBASE'
      )
    ) %>%
    pivot_wider(names_from = 'mod', values_from = 'val')

  toplo <- cmp %>% 
    left_join(lims, by = 'var') %>% 
    summarise(
      Synthetic = mean(Synthetic, na.rm = T),
      EBASE = mean(EBASE, na.rm = T),
      Date = min(Date),
      .by = c('subfig', 'grp', 'var', 'min', 'max')
    ) %>% 
    pivot_longer(-c(subfig, Date, grp, var, min, max), names_to = 'model', values_to = 'est') %>% 
    mutate(
      var = factor(var, 
                   levels = c('P', 'R', 'D', 'DO_mod', 'a', 'b'), 
                   labels = c('P~(mmol~m^{2}~d^{-1})', 
                              'R~(mmol~m^{2}~d^{-1})', 
                              'D~(mmol~m^{2}~d^{-1})',
                              'O[2]~(mmol~m^{-3})',
                              'italic(a)~(mmol~m^{-3}~d^{-1})/(W~m^{-2})', 
                              'italic(b)~(cm~hr^{-1})/(m^{2}~s^{-2})'))
    ) %>% 
    select(-grp)
  
  brks <- seq.Date(min(toplo$Date), max(toplo$Date), by = '3 months')
  
  pout <- toplo %>% 
    group_nest(var) %>% 
    mutate(
      p = purrr::pmap(list(var, data), function(var, data){

        p <- ggplot(data, aes(x = Date, y = est, group = model, color = model)) + 
          geom_line(alpha = 1) +
          facet_grid(~ subfig, scales = 'free_y', labeller = label_parsed) +
          scale_x_date(date_labels = '%b', breaks = brks) + 
          theme_minimal() +
          theme(
            strip.placement = 'outside', 
            strip.background = element_blank(), 
            legend.position = 'top', 
            legend.title = element_blank(), 
            legend.text = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.title.y = element_text(size = 8.5),
            axis.text.y = element_text(size = 10),
            strip.text.x = element_text(size = rel(1.4)), 
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill = 'gray97', color = NA),
            plot.margin = unit(c(0,0,0,0), "cm")
          ) + 
          labs(
            x = NULL, 
            y = parse(text = as.character(var))
          )

        if(!grepl('P', var))
          p <- p + 
            theme(strip.text.x = element_blank())
        
        if(!grepl('b', var))
          p <- p +
            theme(axis.text.x = element_blank())
        
        minv <- unique(data$min)
        maxv <- unique(data$max)
        if(any(!is.na(c(minv, maxv))))
          p <- p + coord_cartesian(ylim = c(minv, maxv))
        
        return(p)
        
      })
    )

  p <- pout$p[[1]] + pout$p[[2]] + pout$p[[3]] + pout$p[[4]] + pout$p[[5]] + pout$p[[6]] + 
    plot_layout(ncol = 1, guides = 'collect') & 
    theme(legend.position = 'top')
    
  return(p)
  
}

# get mean summary of r2 or rmse for apasumdat, parms as T/F to return means and sd
metsum_fun <- function(apasumdat, met, parms = F){

  metsel <- tibble(
      met = c('r2', 'rmse', 'mape', 'mae', 'nse'),
      direc = c(-1, 1, 1, 1, -1)
    ) %>% 
    filter(met == !!met)

  out <- apasumdat %>% 
    mutate(
      ind = rep(1: (nrow(.) / 3), times = 3)
    ) %>% 
    tidyr::unnest('ests') %>% 
    rename(met = !!met) %>% 
    summarise(
      metsum = median(met, na.rm = T), 
      .by = c('amean', 'asd', 'rmean', 'rsd', 'ind', 'ndays')
    ) %>% 
    mutate(
      rnkdir = metsel$direc * metsum,
      rnkmetsum = rank(rnkdir), 
      .by = 'ndays'
    ) %>% 
    select(-rnkdir)
  
  if(!parms)
    out <- out %>% 
      select(-matches('mean$|sd$'))
  
  return(out)
  
}

# compare synthetic and synthetic + noise
syncomp_plo <- function(resobs, resnos){
  
  res <- list(
    'Synthetic' = resobs,
    'Synthetic + noise' = resnos
  ) %>% 
    enframe(name = 'dotyp') %>% 
    unnest('value') %>% 
    select(-matches('hi$|lo$|converge|rsq|^Date$|dDO|^H$')) %>% 
    pivot_longer(-c('dotyp', 'DateTimeStamp', 'grp')) %>% 
    summarise(
      value = mean(value, na.rm = T),
      DateTimeStamp = min(DateTimeStamp),
      .by = c('dotyp', 'grp', 'name')
    ) %>% 
    filter(!name == 'DO_obs') %>% 
    mutate(
      name = factor(name, 
                    levels = c('DO_mod', 'a', 'b', 'P', 'R', 'D'), 
                    labels = c('O[2]~(mmol~m^{-3})', 'italic(a)~(mmol~m^{-2}~d^{-1})/(W~m^{-2})', 'italic(b)~(cm~hr^{-1})/(m^{2}~s^{-2})', 'P~(mmol~m^{2}~d^{-1})', 'R~(mmol~m^{2}~d^{-1})', 'D~(mmol~m^{2}~d^{-1})')
      )
    ) %>% 
    pivot_wider(names_from = dotyp, values_from = value)
  
  vrs <- levels(res$name)
  
  for(i in seq_along(vrs)){
    
    vr <- vrs[i]

    ylab <- NULL
    xlab <- NULL
    if(i %in% c(1, 4))
      ylab <- 'Synthetic + noise'
    if(i %in% 4:6)
      xlab <- 'Synthetic'
    
    toplo <- res %>% filter(name == vr)
    
    lims <- lims <- range(toplo[, c('Synthetic', 'Synthetic + noise')], na.rm = T)
    
    thm <- theme_minimal() + 
      theme(
        strip.background = element_blank(), 
        legend.position = 'top', 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 6)
      )
    
    ptmp <- ggplot(toplo, aes(x = Synthetic, y = `Synthetic + noise`)) + 
      geom_point() + 
      facet_wrap(~name, labeller = label_parsed) + 
      geom_abline(intercept = 0, slope = 1) +
      geom_smooth(formula = y ~ x ,method = 'lm', se = F, colour = 'tomato1', linewidth = 0.7) +
      coord_cartesian(
        xlim = lims, 
        ylim = lims
      ) + 
      labs(
        x = xlab, 
        y = ylab
      ) +
      thm
    
    assign(paste0('p', i), ptmp)
    
  }
  
  p <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3)
  
  return(p)
  
}