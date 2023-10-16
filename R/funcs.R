# function for formatting p-values in tables
p_ast <- function(x){
  
  sig_cats <- c('**', '*', '')
  sig_vals <- c(-Inf, 0.005, 0.05, Inf)
  
  out <- cut(x, breaks = sig_vals, labels = sig_cats, right = FALSE)
  out <- as.character(out)
  
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

  
  metsel <- tibble(
    met = c('r2', 'rmse', 'mape', 'mae', 'nse'),
    lbspr = c('italic(R)^2', 'RMSE', 'MAPE', 'mae', 'NSE')
  ) %>% 
    filter(met == !!met)
  ylab <- parse(text = paste0('`%`', '~change~', metsel$lbspr))
  
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
      upr = medv + (iqr / 2), 
      lwr = medv - (iqr / 2),
      param = gsub('mean$|sd$', '', prior),
      param = factor(param, levels = c('a', 'r'), labels = c('a', 'R')),
      prior = gsub('^a|^r|^b', '', prior), 
      prior = factor(prior, levels = c('mean', 'sd'), labels = c('mu', 'sigma')), 
      ndays = factor(ndays, levels = c('1 day', '7 days', '30 days'), labels = c('1~day', '7~days', '30~days'))
    )
  
  ylab <- toupper(metsel$met)
  
  wd <- 0.9
  
  lvs <- levels(toplo$ndays)
  for(i in seq_along(lvs)){
    
    toploi <- toplo %>% 
      filter(ndays %in% lvs[i])
    
    ylb <- NULL
    if(i == 1)
      ylb <- ylab
    
    # yrng <- max(abs(range(toploi[, c('upr', 'lwr')])))
    
    p <- ggplot(toploi, aes(x = param, group = val, fill = val)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr, color = val), position = position_dodge(width = wd), width = 0, linewidth = 1) +
      geom_col(aes(y = medv, color = val), position = position_dodge(width = wd), color = 'grey') +
      facet_grid(prior~ndays, labeller = label_parsed) + 
      scale_fill_brewer(palette = 'Blues') +
      scale_color_brewer(palette = 'Blues') +
      # coord_cartesian(ylim = c(-1 * yrng, yrng)) +
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
        y = ylb,
        color = 'Prior\nvalue',
        fill = 'Prior\nvalue'
      ) 
    
    if(i != 3)
      p <- p + theme(strip.text.y = element_blank())
    
    assign(paste0('p', i), p)
    
  }
  
  p <- p1 + p2 + p3 + plot_layout(ncol = length(lvs), guides = 'collect')
  
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
    filter(var != 'b') %>% 
    summarise(
      Synthetic = mean(Synthetic, na.rm = T),
      EBASE = mean(EBASE, na.rm = T),
      Date = min(Date),
      .by = c('subfig', 'grp', 'var', 'min', 'max')
    ) %>% 
    pivot_longer(-c(subfig, Date, grp, var, min, max), names_to = 'model', values_to = 'est') %>% 
    mutate(
      var = factor(var, 
                   levels = c('P', 'R', 'D', 'DO_mod', 'a'), 
                   labels = c('P~(mmol~m^{2}~d^{-1})', 
                              'R~(mmol~m^{2}~d^{-1})', 
                              'D~(mmol~m^{2}~d^{-1})',
                              'O[2]~(mmol~m^{-3})',
                              'italic(a)~(mmol~d^{-1}~W^{-1})')
      )
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
            axis.title.y = element_text(size = 12),
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
        
        if(!grepl('a', var))
          p <- p +
            theme(axis.text.x = element_blank())
        
        minv <- unique(data$min)
        maxv <- unique(data$max)
        if(any(!is.na(c(minv, maxv))))
          p <- p + coord_cartesian(ylim = c(minv, maxv))
        
        return(p)
        
      })
    )

  p <- pout$p[[1]] + pout$p[[2]] + pout$p[[3]] + pout$p[[4]] + pout$p[[5]] + 
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
syncomp_plo <- function(resobs, resnos, fwdatcmp){
  
  toplo1 <- list(
    'Synthetic' = fwdatcmp,
    'EBASE recovered' = resobs,
    'EBASE recovered (+ noise)' = resnos
  ) %>% 
    enframe() %>% 
    unnest('value') %>% 
    select(matches('name|Date|DateTimeStamp|grp|^a$|^R$')) %>% 
    pivot_longer(names_to = 'var', values_to = 'val', a:R) %>% 
    group_nest(name, var) %>% 
    mutate(
      data = purrr::pmap(list(name, data), function(name, data){
        
        if(name == 'Synthetic')
          return(data)
        
        out <- data %>% 
          group_by(grp) %>% 
          reframe(
            val = unique(val, na.rm = T), 
            DateTimeStamp = min(DateTimeStamp)
          ) %>% 
          na.omit()
        
        return(out)
        
      })
    ) %>% 
    unnest('data') %>% 
    mutate(
      var = factor(var, 
                   levels = c('a', 'R'), 
                   labels = c('italic(a)~(mmol~m^{-2}~d^{-1})/(W~m^{-2})',
                              'R~(mmol~m^{-2}~d^{-1})'
                   )
      ), 
      name = factor(name, levels = c('Synthetic', 'EBASE recovered', 'EBASE recovered (+ noise)'))
    )
  
  p1 <- ggplot(toplo1, aes(x = DateTimeStamp, y = val, color = name)) + 
    geom_line(linewidth = 0.8) +
    facet_wrap(~var, scales = 'free', ncol = 1, strip.position = 'left', labeller = label_parsed) + 
    theme_minimal(base_size = 12) + 
    theme(
      strip.placement = 'outside', 
      legend.position = 'top', 
      panel.grid.minor = element_blank()
    ) + 
    labs(
      y = NULL,
      x = NULL,
      color = NULL
    )
  
  unidt <- toplo1 %>% 
    filter(name == 'EBASE recovered') %>% 
    pull(DateTimeStamp) %>% 
    unique() %>% 
    as.Date()
  
  toplo2 <- toplo1 %>% 
    select(-Date) %>% 
    mutate(DateTimeStamp = as.Date(DateTimeStamp)) %>% 
    filter(DateTimeStamp %in% unidt) %>% 
    group_by(name, var, DateTimeStamp) %>% 
    reframe(
      val = mean(val)
    ) %>% 
    pivot_wider(names_from = 'name', values_from = 'val') %>% 
    pivot_longer(names_to = 'EBASE', values_to = 'val', matches('EBASE'))
  
  toplo2a <- toplo2 %>% 
    filter(grepl('a', var))
  rnga <- range(c(toplo2a$Synthetic, toplo2a$val)) 
  p2a <- ggplot(toplo2a, aes(x = Synthetic, y = val, color = EBASE)) + 
    geom_point(show.legend = F, size = 0.5) +
    geom_smooth(method = 'lm', se = F, formula = y~x, show.legend = F) +
    scale_color_manual(values = c("#00BA38", "#619CFF")) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(
      xlim = c(rnga[1], rnga[2]), 
      ylim = c(rnga[1], rnga[2])
    ) +
    theme(legend.position = 'none') +
    facet_wrap(~var, scales = 'free', labeller = label_parsed, strip.position = 'top') + 
    labs(
      color = NULL,
      y = 'EBASE'
    )
  
  toplo2b <- toplo2 %>% 
    filter(grepl('R', var))
  rngb <- range(c(toplo2b$Synthetic, toplo2b$val)) 
  p2b <- ggplot(toplo2b, aes(x = Synthetic, y = val, color = EBASE)) + 
    geom_point(show.legend = F, size = 0.5) +
    geom_smooth(method = 'lm', se = F, formula = y~x, show.legend = F) +
    scale_color_manual(values = c("#00BA38", "#619CFF")) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(
      xlim = c(rngb[1], rngb[2]), 
      ylim = c(rngb[1], rngb[2])
    ) +
    theme(legend.position = 'none') +
    facet_wrap(~var, scales = 'free', labeller = label_parsed, strip.position = 'top') + 
    labs(
      color = NULL,
      y = 'EBASE'
    )
  
  p <- p1 + (p2a + p2b + plot_layout(ncol = 1, guides = 'collect') & 
               theme_minimal() + 
               theme(
                 strip.background = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 strip.placement = 'outside', 
               )) + plot_layout(ncol = 2, guides = 'collect', widths = c(1, 0.4)) & theme(legend.position = 'bottom')
  
  return(p)
  
}

# create subplots for apalachicola obs v dtd metab comparisons
apacmp_plo <- function(dat, alph = 0.5, col = 'black'){
  
  toplo <- dat %>% 
    filter(var != 'D') %>% 
    mutate(
      var = factor(var, 
                   levels = c('P', 'R', 'NEM')
      ), 
      dotyp = factor(dotyp, 
                     levels = c('observed', 'detided'), 
                     labels = c('Observed', 'Detided'))
    ) %>% 
    pivot_wider(names_from = 'typ', values_from = 'val') %>% 
    na.omit()
  
  thm <- theme_minimal() + 
    theme(
      strip.background = element_blank(), 
      legend.position = 'top', 
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10)
    )
  
  subfigs <- list(
    P = '(a)~P~(mmol~m^{-2}~d^{-1})', 
    R = '(b)~R~(mmol~m^{-2}~d^{-1})', 
    NEM = '(c)~NEM~(mmol~m^{-2}~d^{-1})'
  )
  
  lvs <- levels(toplo$var)
  
  for(vr in seq_along(lvs)){
    
    var <- lvs[vr]
    
    toplotmp <- toplo %>% 
      filter(var == !!var)
    
    lims <- range(toplotmp[, c('Odum', 'EBASE')], na.rm = T)
    
    ttl <- subfigs[[var]]
    
    ptmp <- ggplot(toplotmp, aes(x = EBASE, y = Odum)) + 
      geom_point(color = col, alpha = alph, stroke = 0, size = 1.5) + 
      geom_abline(intercept = 0, slope = 1) +
      geom_smooth(formula = y ~ x ,method = 'lm', se = F, colour = 'tomato1', linewidth = 0.7) +
      facet_wrap(~dotyp, ncol = 2) +
      scale_colour_manual(values = cols) + 
      scale_x_continuous(limits = lims) + 
      scale_y_continuous(limits = lims) + 
      thm + 
      labs(
        title = parse(text = ttl)
      )
    
    assign(paste0('p', vr), ptmp)
    
  }
  
  out <- p1 + p2 + p3 + plot_layout(ncol = 1)
  
  return(out)
  
}