

# calculate paper outputs -------------------------------------------------

.print.t2 <- function( df ){
  df %<>% select(
    records, strata, paramcd, thr, events, median, q25, q75, censoring 
  )
  df %>% 
    filter( censoring == 'lcn') %>% 
    rename( records.lcn = records, events.lcn = events, median.lcn = median, q25.lcn = q25, q75.lcn = q75) %>% 
    select( -censoring ) %>% 
    left_join(
      df %>% 
        filter( censoring == 'ooe' )%>% 
        select(-censoring)
    ) %>% 
    mutate( events.obs = events , trunc.bias = median-median.lcn ) %>% 
    mutate( lcn.pct = round ( 100*(records.lcn-records)/records.lcn ) )  %>% 
    select( strata, par = paramcd, thr, N = records.lcn, events.obs, median = median.lcn, q25.lcn, q75.lcn, lcn.pct, trunc.bias)
  
    }

# Turnbull Coding ---------------------------------------------------------

.prepare_tte_raw <- function( 
    
  data, par, thr
  
  ) {
  
  data %>%
    group_by( sjid ) %>% 
    filter ( paramcd == par ) %>%
    arrange( sjid, tx. ) %>%
    mutate(
      event = ifelse(aval >= thr, 1, 0)
      # ,
      # viscount = dense_rank(tx.)
    ) %>%
    group_by(sjid) %>%
    mutate(
      first_event_time = ifelse( any(event == 1), min( tx.[event == 1] ), Inf)
    ) %>%
    filter( tx. <= first_event_time ) %>%
    ungroup() %>%
    group_by( sjid ) %>%
    summarise(
      
      t = if (any(event == 1) & any(event == 0)) max(tx.[event == 0])
      else if (all(event == 1)) NA_real_
      else max( tx. ),
      
      t2 = if (any(event == 1)) min(tx.[event == 1]) else NA_real_,
      
      pevent  = max(event),
      study   = first( study ),
      grp.    = first( grp. ),
      paramcd = par,
      thr = thr,
      .groups = "drop"
    ) %>%
    mutate(
      # Turn event codes into survival types:
      #  0 = right-censored
      #  1 = exact event
      #  2 = left-censored
      #  3 = interval-censored
      event = case_when(
        is.na(t) & !is.na(t2) ~ 2,
        pevent == 0           ~ 0,
        t == t2               ~ 1,
        is.na(t2)             ~ 0,
        TRUE                  ~ 3
      ),
      
      # Build interval bounds based on event type
      L = case_when(
        event == 2 ~ 0,         # left-censored
        TRUE       ~ t          # all others
      ),
      R = case_when(
        event == 0 ~ Inf,       # right-censored
        TRUE       ~ t2         # all others
      ),
      
      # CRITICAL: Turnbull() expects a reversed logic here
      # - Interval-censored → 0 (not censored)
      # - All others        → 1 (treated as censored by Turnbull)
      # This is unintuitive and reversed from most survival logic
      censor = ifelse(event == 3, 0, 1)
    )
}

# main loop ---------------------------------------------------------------

.run_tte_analysis <- function(
    
    param.thresholds.list,
    tm.  = "dur",                       # which time axis to use
    grp. = "sev.o",                     # which subgroup to use
    save.outputs = TRUE,
    save.dir = "DATA derived",
    verbose  = TRUE
    
) {
  
  require(dplyr)
  require(glue)
  require(ReIns)
  
  # === Load preprocessed data
  
  dm.     <- readRDS(file.path(save.dir, "demo.UNIFAI.rds")) %>% 
    mutate ( grp. = .data[[grp.]] ) %>% 
    select ( study, sjid, grp. )
  
  scores. <- readRDS(file.path(save.dir, "scores.UNIFAI.rds"))
  
  param.list <- names( param.thresholds.list )
  out.tte    <- tibble()
  mod.tte    <- list()
  
  # --- Filter scores and define time axis
  
  scores. <- scores. %>%
    filter( paramcd %in% param.list ) %>%
    inner_join( dm. ) %>%
    mutate(
      tx.   = .data[[ tm. ]],
    ) %>%
    select( study, sjid, grp., tx., avisitn, paramcd, aval)
  
  # === Main Loop over parameters and thresholds
  
  for (par in param.list) {
    
    thresholds <- param.thresholds.list[[par]]
    
    for (thr in thresholds) {
      
      tte.raw <- .prepare_tte_raw( scores., par, thr)
      
      # Optionally save raw file
      safe.name <- gsub("[^a-zA-Z0-9._-]", "_", paste(par, thr, sep = "."))
      if (save.outputs) {
        saveRDS(tte.raw, file = file.path(save.dir, paste0("tte.raw.", safe.name, ".rds")))
      }
      
      out.grp <- tibble()
      labels <- levels(droplevels(tte.raw$grp.))
      
      for (bx in labels) {
        tte.grp <- tte.raw %>% filter( grp. == bx, L <= R)
        
        if (3 %in% tte.grp$event || (0 %in% tte.grp$event && 1 %in% tte.grp$event)) {
          key.lcn <- sprintf("%s.%s.%s.lcn", par, thr, bx)
          key.ooe <- sprintf("%s.%s.%s.ooe", par, thr, bx)
          
          if (verbose) message(glue("Running {key.lcn}"))
          fit1 <- Turnbull(
            x = seq(0, 1),
            L = tte.grp$L,
            R = tte.grp$R,
            censored = tte.grp$censor,
            trunclower = -Inf,
            truncupper = Inf
          )
          mod.tte[[key.lcn]] <- fit1$fit
          
          q1 <- quantile(fit1$fit, conf.int = FALSE)
          out1 <- summary(fit1$fit)$table %>%
            t() %>% as_tibble() %>%
            mutate(strata = bx, censoring = "lcn", q25 = q1[1], q75 = q1[3])
          
          if (verbose) message(glue("Running {key.ooe}"))
          tte.grp.ooe <- tte.grp %>% filter(L != 0)
          
          fit2 <- Turnbull(
            x = seq(0, 1),
            L = tte.grp.ooe$L,
            R = tte.grp.ooe$R,
            censored = tte.grp.ooe$censor,
            trunclower = -Inf,
            truncupper = Inf
          )
          
          mod.tte[[key.ooe]] <- fit2$fit
          
          q2 <- quantile(fit2$fit, conf.int = FALSE)
          
          out2 <- summary(fit2$fit)$table %>%
            t() %>% as_tibble() %>%
            mutate(strata = bx, censoring = "ooe", q25 = q2[1], q75 = q2[3])
          
          out.grp <- bind_rows(out.grp, out1, out2)
          
        } else {
          
          warning(glue("⚠ No usable events for {par}, threshold {thr}, group {bx} — skipping."))
          out.grp <- bind_rows(out.grp, tibble(
            records = NA, n.max = NA, n.start = NA, events = NA,
            rmean = NA, se_rmean = NA, median = NA,
            LCL = NA, UCL = NA, strata = bx,
            censoring = "lcn", q25 = NA, q75 = NA
            
          ))
        }
      }
      
      out.tte <- bind_rows(out.tte, out.grp %>% mutate(paramcd = par, thr = thr))
      
    }
  }
  
  # === Save Output
  if (save.outputs) {
    saveRDS(out.tte, file.path(save.dir, "tte.summary.rds"))
    saveRDS(mod.tte, file.path(save.dir, "fits.tte.rds"))
  }
}


# plotting ----------------------------------------------------------------


.plot_overlay_surv <- function(
    plot.set,
    mod.tte,
    labs = "pstc",
    time.breaks = pretty(plot.set$median, n = 10),
    risk.table = TRUE,
    plot.mods = NULL,
    censoring = TRUE,
    surv.median.line = "hv"
) {
  
  require(dplyr)
  require(ggplot2)
  require(purrr)
  require(tibble)
  require(ggpubr)
  
  label_vars <- list(
    p = "paramcd",
    s = "strata",
    t = "thr",
    c = "censoring"
  )
  
  label.by <- unname(unlist(label_vars[strsplit(labs, "")[[1]]]))
  plot.set <- plot.set %>% mutate(.order = row_number())
  
  for (lab in label.by) {
    if (lab %in% names(plot.set)) {
      lab.order <- plot.set %>% pull(!!sym(lab)) %>% unique()
      plot.set <- plot.set %>% mutate(!!lab := factor(!!sym(lab), levels = lab.order))
    }
  }
  
  keys.df <- plot.set %>%
    mutate(
      key = paste(paramcd, thr, strata, censoring, sep = "."),
      label = if (length(label.by) > 0) {
        purrr::pmap_chr(across(all_of(label.by)), ~ paste(..., sep = " | "))
      } else {
        "Group"
      }
    ) %>%
    select(key, label, .order, all_of(label.by)) %>%
    distinct() %>%
    arrange(.order) %>%
    mutate(label = factor(label, levels = unique(label[order(.order)])))
  
  keys.df <- keys.df %>% filter(key %in% names(mod.tte))
  if (nrow(keys.df) == 0) {
    warning("\u26a0 No curves could be plotted.")
    return(ggplot() + labs(title = "No curves found"))
  }
  
  curves <- purrr::pmap_dfr(keys.df, function(key, label, ...) {
    fit <- mod.tte[[key]]
    tibble(
      time = fit$time,
      surv = fit$surv,
      n.censor = fit$n.censor,
      label = label
    )
  })
  
  if (is.null(time.breaks)) {
    time.breaks <- pretty(range(curves$time, na.rm = TRUE), n = 10)
  }
  
  risk.data <- purrr::pmap_dfr(keys.df, function(key, label, ...) {
    fit <- mod.tte[[key]]
    s <- summary(fit, times = time.breaks)
    tibble(time = s$time, n.risk = round(s$n.risk), label = label)
  })
  
  surv.plot <- ggplot(curves, aes(x = time, y = surv, color = label)) +
    geom_step(linewidth = .75)
    # geom_step(linewidth = .5, alpha = 0.9)
  
  if (censoring) {
    surv.plot <- surv.plot +
      geom_point(
        data = curves %>% filter(n.censor > 0),
        shape = 3, size = 2,
        aes(x = time, y = surv, color = label)
      )
  }
  
  if (surv.median.line != "none") {
    median.data <- purrr::map_dfr(keys.df$key, function(key) {
      med <- tryCatch({
        survminer::surv_median(mod.tte[[key]])
      }, error = function(e) NULL)
      if (!is.null(med)) tibble(median = med$median, label = keys.df$label[keys.df$key == key])
    })
    
    if (nrow(median.data) > 0) {
      if (grepl("v", surv.median.line)) {
        surv.plot <- surv.plot +
          geom_segment(
            data = median.data,
            aes(x = median, xend = median, y = 0, yend = 0.5, color = label),
            linetype = "dashed"
          )
      }
      if (grepl("h", surv.median.line)) {
        surv.plot <- surv.plot +
          geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40")
      }
    }
  }
  
  surv.plot <- surv.plot +
    scale_x_continuous(
      limits = range(time.breaks, na.rm = TRUE),
      breaks = time.breaks
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Time",
      y = "Survival Probability",
      color = if (labs == "") "Group" else paste(label.by, collapse = " + ")
    ) +
    theme_minimal()
  
  if (!is.null(plot.mods)) {
    for (m in plot.mods) surv.plot <- surv.plot + m
  }
  
  risk.plot <- ggplot(risk.data, aes(x = time, y = label, label = n.risk)) +
    geom_text(size = 3, hjust = 0.5, vjust = 0.5) +
    scale_y_discrete(limits = rev(levels(keys.df$label))) +
    scale_x_continuous(breaks = time.breaks) +
    labs(x = "Time", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (risk.table) {
    ggpubr::ggarrange(
      surv.plot,
      risk.plot,
      ncol = 1,
      heights = c(3, 1),
      align = "v"
    )
  } else {
    surv.plot
  }
}


# slide -------------------------------------------------------------------

.tteslide <- function(
    
    ti       = "TTE results",
    l        = "TTE_results",
    i        = 2,
    template = .ppt.template.file,
    m        = "CR",
    image    = FALSE
    
) {
  # Retrieve the plot
  if (!exists("pp")) {
    pp <- if (exists("p")) p else last_plot()
  }
  
  # Fix minus signs
  .fix_plot_minuses <- function(p) {
    p + labs(
      title    = gsub("-", "\u2212", p$labels$title),
      subtitle = gsub("-", "\u2212", p$labels$subtitle),
      caption  = gsub("-", "\u2212", p$labels$caption),
      x        = gsub("-", "\u2212", p$labels$x),
      y        = gsub("-", "\u2212", p$labels$y)
    )
  }
  
  pp <- .fix_plot_minuses(pp)
  
  if (is.null(pp)) stop("No plot found. Please create a plot before calling .tteslide().")
  if (!exists("plot.set")) stop("plot.set not found in workspace.")
  
  # Table from plot.set
  tab <- plot.set %>% .print.t2() %>% 
    select(-starts_with('q')) 
  
  # Output file
  target_file <- .create_output_filepath(ti)
  
  # Create slide
  ppt <- read_pptx(template) %>% 
    add_slide(layout = l, master = m)
  
  # Add plot to placeholder index i
  ppt <- if (!image) {
    ppt %>%
      ph_with(
        dml(print(pp, newpage = FALSE)),
        location = ph_location_type(type = "body", type_idx = i)
      )
  } else {
    ppt %>%
      ph_with(
        print(pp, newpage = FALSE),
        location = ph_location_type(type = "body", type_idx = i)
      )
  }
  
  # Add table to placeholder index (i - 1)
  ppt <- ppt %>%
    ph_with(
      flextable::flextable(tab) %>% autofit(),
      location = ph_location_type(type = "body", type_idx = 3)
    )
  
  # Add title and notes
  ppt <- ppt %>%
    ph_with(ti, location = ph_location_type(type = "title")) %>%
    set_notes(
      value = paste("Created at", format(Sys.time(), "%Y-%m-%d %H-%M-%S")),
      location = notes_location_type("body")
    ) %>%
    print(target = target_file)
  
  rm(pp)
}
