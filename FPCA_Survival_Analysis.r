library(mxfda)
library(ggplot2)
library(openxlsx)
library(survival)
library(dplyr)
library(metafor)
library(survminer)

analysis_dat_list1 = readRDS("data_part1.rds")
analysis_dat_list2 = readRDS("data_part2.rds")

metrics = names(analysis_dat_list1)

analysis_dat_list_combined = lapply(metrics, function(m){
  dat = dplyr::bind_rows(
    analysis_dat_list1[[m]],
    analysis_dat_list2[[m]] %>%
      dplyr::mutate(dplyr::across(c(block, histotype),
                                  ~ as.character(.x)),
                    stage_cat = ifelse(stage_cat == 1, "high stage", "low state"))
  )
  #abundance is not cohort dependent
  #func_observed is not cohort dependent
  cohort_dat = split(dat, dat$cohort)
  return(cohort_dat)
})

### Summary Table Info
AAS_summary_dat_list <- analysis_dat_list_combined[[1]]$AAS |> filter(r ==25)
NCO_summary_dat_list <- analysis_dat_list_combined[[1]]$NCO |> filter(r ==25)
NECC_summary_dat_list <- analysis_dat_list_combined[[1]]$NECC |> filter(r ==25)
NHSI_summary_dat_list <- analysis_dat_list_combined[[1]]$NHS |> filter(r ==25)
NHSII_summary_dat_list <- analysis_dat_list_combined[[1]]$NHSII |> filter(r ==25)

AAS_summary_dat_list_cyto <- analysis_dat_list_combined[[2]]$AAS |> filter(r ==25)
NCO_summary_dat_list_cyto <- analysis_dat_list_combined[[2]]$NCO |> filter(r ==25)
NECC_summary_dat_list_cyto <- analysis_dat_list_combined[[2]]$NECC |> filter(r ==25)
NHSI_summary_dat_list_cyto <- analysis_dat_list_combined[[2]]$NHS |> filter(r ==25)
NHSII_summary_dat_list_cyto <- analysis_dat_list_combined[[2]]$NHSII |> filter(r ==25)

# Age at Diagnosis
mean(AAS_summary_dat_list$ageatdx)
sd(AAS_summary_dat_list$ageatdx)
mean(NCO_summary_dat_list$ageatdx)
sd(NCO_summary_dat_list$ageatdx)
mean(NECC_summary_dat_list$ageatdx)
sd(NECC_summary_dat_list$ageatdx)
mean(NHSI_summary_dat_list$ageatdx)
sd(NHSI_summary_dat_list$ageatdx)
mean(NHSII_summary_dat_list$ageatdx)
sd(NHSII_summary_dat_list$ageatdx)

# Stage
table(AAS_summary_dat_list$stage_cat)
table(NCO_summary_dat_list$stage_cat)
table(NECC_summary_dat_list$stage_cat)
table(NHSI_summary_dat_list$stage_cat)
table(NHSII_summary_dat_list$stage_cat)

# Vital Status
table(AAS_summary_dat_list$death)
table(NCO_summary_dat_list$death)
table(NECC_summary_dat_list$death)
table(NHSI_summary_dat_list$death)
table(NHSII_summary_dat_list$death)

# Survival Stats
median(AAS_summary_dat_list$mossurv)
c(quantile(AAS_summary_dat_list$mossurv, probs = 0.10),quantile(AAS_summary_dat_list$mossurv, probs = 0.90))
median(NCO_summary_dat_list$mossurv)
c(quantile(NCO_summary_dat_list$mossurv, probs = 0.10),quantile(NCO_summary_dat_list$mossurv, probs = 0.90))
median(NECC_summary_dat_list$mossurv)
c(quantile(NECC_summary_dat_list$mossurv, probs = 0.10),quantile(NECC_summary_dat_list$mossurv, probs = 0.90))
median(NHSI_summary_dat_list$mossurv)
c(quantile(NHSI_summary_dat_list$mossurv, probs = 0.10),quantile(NHSI_summary_dat_list$mossurv, probs = 0.90))
median(NHSII_summary_dat_list$mossurv)
c(quantile(NHSII_summary_dat_list$mossurv, probs = 0.10),quantile(NHSII_summary_dat_list$mossurv, probs = 0.90))


#create the mxfda objects

mxfda_objects = lapply(analysis_dat_list_combined, function(co){
  co_mxfda_objects = lapply(co, function(l){
    l = l %>%
      dplyr::mutate(sample = block) #because mxfda expects there to be image/sample IDs and subject IDs, duplicate
    clin = l %>%
      dplyr::select(-c(func_weighted, r)) %>%
      dplyr::distinct()
    func = l %>%
      dplyr::select(sample, r, func_weighted)
    
    mxfda_obj = make_mxfda(metadata = clin,
                           spatial = NULL,
                           subject_key = "block",
                           sample_key = "sample")
    mxfda_obj = add_summary_function(mxfda_obj,
                                     summary_function_data = func,
                                     metric = 'uni g') #dummy spot
    message(unique(clin$Marker), "\t", unique(clin$Metric), "\t", unique(clin$cohort))
    return(mxfda_obj)
  })
  return(co_mxfda_objects)
})
names(mxfda_objects) = sapply(analysis_dat_list_combined, function(l) paste0(unique(l[[1]]$Marker), " ", unique(l[[1]]$Metric)))

#run FPCA
mxfda_objects = lapply(mxfda_objects, function(co){
  co_mxfda_objects = lapply(co, function(obj){
    obj = run_fpca(obj,
                   metric = 'uni g',
                   r = 'r',
                   value = 'func_weighted',
                   pve = 0.99)
    summary(obj)
    return(obj)
  })
  return(co_mxfda_objects)
})

#plot the fpcs
pls = lapply(names(mxfda_objects), function(nam){
  co_plots = lapply(names(mxfda_objects[[nam]]), function(co){
    obj = mxfda_objects[[nam]][[co]]
    p1 = plot(obj, what = 'uni g fpca', pc_choice = 1) + lims(x = c(0, 100)) #show that not using 0
    p2 = plot(obj, what = 'uni g fpca', pc_choice = 2) + lims(x = c(0, 100))
    plots = ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2) %>%
      ggpubr::annotate_figure(left = ggpubr::text_grob(co, rot = 90))
    data = dplyr::bind_rows(p1$data %>% dplyr::mutate(FPC = "fpc1"),
                           p2$data %>% dplyr::mutate(FPC = "fpc2")) %>%
      dplyr::mutate(Cohort = co)
    #want to return the data for flipping directions later
    list(plot = plots,
         data = data)
  })
  plist = lapply(co_plots, function(x) x$plot)
  data = lapply(co_plots, function(x) x$data) %>% do.call(dplyr::bind_rows, .)
  pl_panel = ggpubr::ggarrange(plotlist = plist, nrow = 5, ncol = 1) %>%
    ggpubr::annotate_figure(top = ggpubr::text_grob(nam))
  return(list(plot = pl_panel,
              data = data))
})
names(pls) = names(mxfda_objects)


# fix FPC direction -------------------------------------------------------
#because the FPCs don't seem to all be in the same direction

fpca_pl_dat = lapply(pls, function(met){
  fpca_dat_list = split(met$data, f = list(met$data$Cohort, met$data$FPC))
  fixed_fpca_dat = lapply(fpca_dat_list, function(d){
    if(d$plus[d$index == 25] < d$mu[d$index == 25]) {
      t = d$plus
      d$plus = d$minus
      d$minus = t
      d$adjust = -1
    } else {
      d$adjust = 1
    }
    return(d)
  }) %>%
    do.call(dplyr::bind_rows, .)
  pl = fixed_fpca_dat %>%
    ggplot() + 
    geom_point(aes(x = index, y = mu)) + 
    facet_grid(Cohort ~ FPC) +
    geom_point(aes(x = index, y = plus), color = "blue", shape = "+", size = 2) +
    geom_point(aes(x = index, y = minus), color = "red", shape = "-", size = 2) +
    theme_minimal() +
    labs(y = "Observed - CSR", x = "Radius (px)")
  return(list(plot = pl,
              data = fixed_fpca_dat))
})

step_wise_models = lapply(mxfda_objects, function(m){
  lapply(m, function(obj){
    meta = obj@Metadata
    
    #model 1 - no pcs
    mod1 = coxph(Surv(mossurv, death) ~ ageatdx + stage_cat + abund_binary + func_observed, data = meta) #removed cohort as covariate because at cohort level already
    
    #subset to func observed only
    obj2 = obj
    obj2@Metadata = obj2@Metadata %>%
      dplyr::filter(func_observed == 1)
    obj2@univariate_summaries$Gest = obj2@univariate_summaries$Gest %>%
      dplyr::filter(sample %in% obj2@Metadata$sample)
    
    #remove those that are low abundance 
    obj2 = run_fpca(obj2,
                    metric = 'uni g',
                    r = 'r',
                    value = 'func_weighted',
                    pve = 0.99)
    
    #now since we removed low abundance, have to redo plots and things
    p1 = plot(obj2, what = 'uni g fpca', pc_choice = 1) + lims(x = c(0, 100)) #show that not using 0
    p2 = plot(obj2, what = 'uni g fpca', pc_choice = 2) + lims(x = c(0, 100))
    data = dplyr::bind_rows(p1$data %>% dplyr::mutate(FPC = "fpc1"),
                            p2$data %>% dplyr::mutate(FPC = "fpc2")) %>%
      dplyr::mutate(Cohort = unique(meta$cohort))
    sp_dat = split(data, data$FPC)
    data_fixed = lapply(sp_dat, function(d){
      if(d$plus[d$index == 25] < d$mu[d$index == 25]) {
        t = d$plus
        d$plus = d$minus
        d$minus = t
        d$adjust = -1
      } else {
        d$adjust = 1
      }
      return(d)
    }) %>%
      do.call(dplyr::bind_rows, .)
    
    #get the adjustment value
    adj_val = data_fixed %>%
      dplyr::select(Cohort, FPC, adjust) %>%
      dplyr::distinct()
    
    fpc_dat = dplyr::inner_join(
      obj2@Metadata,
      obj2@functional_pca$Gest$score_df
    ) 
    
    for(fpc in adj_val$FPC){
      fpc_dat[[fpc]] = fpc_dat[[fpc]] * adj_val$adjust[adj_val$FPC == fpc]
    }
    
    fpc_dat2 = fpc_dat %>%
      dplyr::mutate(dplyr::across(dplyr::contains("fpc"),
                                  ~ factor(ifelse(.x >= 0, "high", "low"),
                                           levels = c('high', 'low')), #set reference
                                  .names = "{col}_cat"))
      
    fpc_dat2_summ = fpc_dat2 %>%
      dplyr::summarise(dplyr::across(dplyr::contains("fpc"),
                                     .fns = list(high = function(x){sum(x == "high")},
                                                 low = function(x){sum(x == "low")}),
                                     .names = "{col}_{fn}")) %>%
      dplyr::select(dplyr::contains("cat")) %>%
      tidyr::pivot_longer(cols = dplyr::everything(),
                          names_to = c("fpc", "w", "cat"),
                          names_pattern = "(.*)_(.*)_(.*)") %>%
      dplyr::select(-w) %>%
      tidyr::pivot_wider(names_from = "cat", values_from = "value")
    
    mod2 = coxph(Surv(mossurv, death) ~ 
                   ageatdx + stage_cat +
                   abund_binary * fpc1_cat +
                   abund_binary * fpc2_cat, data = fpc_dat2)
    mod3 = coxph(Surv(mossurv, death) ~ 
                   ageatdx + stage_cat +
                   abund_binary + fpc1_cat + fpc2_cat, data = fpc_dat2)
    
    return(list(
      models = list(
        `Abundance and Func Observed` = mod1,
        `FPC Score Interactions` = mod2,
        `FPC Score No Interactions` = mod3
      ),
      data = list(
        `High Abundance Plot` = data_fixed,
        `High Abundance Scores` = fpc_dat,
        `High Abundance Score Cats` = fpc_dat2_summ
      )
    ))
  })
})

#high abundance fpca plot data
fpca_plots = lapply(step_wise_models, function(m){
  pl_dat = lapply(m, function(x) x$data$`High Abundance Plot`) %>% 
    do.call(dplyr::bind_rows, .)
  pl_dat %>%
    ggplot() + 
    geom_point(aes(x = index, y = mu)) + 
    facet_grid(Cohort ~ FPC) +
    geom_point(aes(x = index, y = plus), color = "blue", shape = "+", size = 2) +
    geom_point(aes(x = index, y = minus), color = "red", shape = "-", size = 2) +
    theme_minimal() +
    labs(y = "Observed - CSR", x = "Radius (px)")
})

#extract model information
outs = lapply(step_wise_models, function(met){
  mod_summaries = lapply(met, function(co){
    tmp = lapply(names(co$models), function(n){
      m = co$models[[n]]
      res = summary(m)
      
      covars = data.frame(res$coefficients, check.names = FALSE) %>%
        tibble::rownames_to_column("factor") %>%
        dplyr::inner_join(
          data.frame(res$conf.int, check.names = FALSE) %>%
            tibble::rownames_to_column("factor")
        )
      model = data.frame(t(res$waldtest), check.names = FALSE) %>%
        dplyr::mutate(info = 'waldtest', .before = 1)
      model = model[rep(1, nrow(covars)),]
      row.names(model) = NULL
      dplyr::bind_cols(covars, model) %>%
        dplyr::mutate(model = n, .before = 1)
    })
    names(tmp) = names(co$models)
    return(tmp)
  })
  abund_func_mods = lapply(names(mod_summaries), function(co){
    mod_summaries[[co]]$`Abundance and Func Observed` %>%
      dplyr::mutate(cohort = co)
  }) %>%
    do.call(dplyr::bind_rows, .)
  func_int_mods = lapply(names(mod_summaries), function(co){
    mod_summaries[[co]]$`FPC Score Interactions` %>%
      dplyr::mutate(cohort = co)
  }) %>%
    do.call(dplyr::bind_rows, .)
  func_no_int_mods = lapply(names(mod_summaries), function(co){
    mod_summaries[[co]]$`FPC Score No Interactions` %>%
      dplyr::mutate(cohort = co)
  }) %>%
    do.call(dplyr::bind_rows, .)
  return(list(`Abundance and Func Observed` = abund_func_mods,
              `FPC Score Interactions` = func_int_mods,
              `FPC Score No Interactions` = func_no_int_mods))
})

#getting Ns
marker_ns = lapply(step_wise_models, function(m){ #metric
  lapply(m, function(s){ #study
    lapply(s$models, function(mod){
      data.frame(Samples = mod$n)
    }) %>%
      dplyr::bind_rows(.id = "Model") %>%
      tidyr::pivot_wider(names_from = "Model", values_from = "Samples")
  }) %>%
    dplyr::bind_rows(.id = "Cohort")
})
marker_ns



#################### Predicted Survival Curves ####################

# refit a fpc1-only interaction model using the data already inside a saved coxph model
refit_fpc1_only <- function(cox_model_with_data) {
  mf <- model.frame(cox_model_with_data)
  
  # survival outcome is stored as a single Surv object
  y <- model.response(mf)
  mf$mossurv <- y[, 1]
  mf$death   <- y[, 2]
  
  # refit Model 4 (FPC1 only, with interaction)
  coxph(Surv(mossurv, death) ~ ageatdx + stage_cat + abund_binary * fpc1_cat, data = mf)
}



# build the contrast-based group logHRs & SEs from a fitted model
extract_group_contrasts <- function(mod) {
  b <- coef(mod)
  V <- vcov(mod)
  
  nm_ab  <- "abund_binary"
  nm_f1  <- grep("^fpc1_cat", names(b), value = TRUE)           # e.g., "fpc1_catlow"
  nm_int <- grep("^abund_binary:fpc1_cat", names(b), value = TRUE)
  
  if (length(nm_f1) != 1) stop("Expected exactly 1 fpc1_cat coefficient; found: ", paste(nm_f1, collapse=", "))
  if (length(nm_int) != 1) stop("Expected exactly 1 abundance:fpc1_cat interaction; found: ", paste(nm_int, collapse=", "))
  
  L_LA_LC <- setNames(rep(0, length(b)), names(b)); L_LA_LC[nm_f1] <- 1
  L_HA_HC <- setNames(rep(0, length(b)), names(b)); L_HA_HC[nm_ab] <- 1
  L_HA_LC <- setNames(rep(0, length(b)), names(b)); L_HA_LC[c(nm_ab, nm_f1, nm_int)] <- 1
  
  get_est <- function(L) {
    est <- sum(L * b)
    se  <- sqrt(as.numeric(t(L) %*% V %*% L))
    c(logHR = est, SE = se)
  }
  
  out <- rbind(
    `Low abund + Low clustering`    = get_est(L_LA_LC),
    `High abund + High clustering`  = get_est(L_HA_HC),
    `High abund + Low clustering`   = get_est(L_HA_LC)
  )
  
  out <- as.data.frame(out)
  out$Group <- rownames(out)
  rownames(out) <- NULL
  out
}


meta_analyze_groups <- function(df) {
  # df must contain: Group, cohort, logHR, SE
  df %>%
    dplyr::group_by(Group) %>%
    dplyr::group_modify(~{
      fit <- metafor::rma(yi = .x$logHR, sei = .x$SE, method = "REML")
      
      data.frame(
        # pooled effect
        logHR = as.numeric(fit$b),
        SE    = as.numeric(fit$se),
        z     = as.numeric(fit$zval),
        p     = as.numeric(fit$pval),
        HR    = exp(as.numeric(fit$b)),
        LCL   = exp(as.numeric(fit$b) - 1.96 * as.numeric(fit$se)),
        UCL   = exp(as.numeric(fit$b) + 1.96 * as.numeric(fit$se)),
        
        # heterogeneity
        tau2  = as.numeric(fit$tau2),
        I2    = as.numeric(fit$I2),
        Q     = as.numeric(fit$QE),
        Q_p   = as.numeric(fit$QEp),
        
        # counts
        k     = fit$k
      )
    }) %>%
    dplyr::ungroup()
}



# utility: pick the "high stage" level robustly
pick_high_stage <- function(stage_factor) {
  lev <- levels(stage_factor)
  hi  <- lev[grepl("high", lev, ignore.case = TRUE)][1]
  if (is.na(hi)) stop("Could not find a 'high' level in stage_cat.")
  hi
}

make_meta_reference_curve <- function(mod_list_by_cohort, time_step = 1) {
  # mod_list_by_cohort: named list of coxph models (already fpc1-only)
  # returns: data.frame(time, Sref, Href)
  
  cohort_names <- names(mod_list_by_cohort)
  
  # overall mean age across all cohorts (weighted by cohort sample size)
  all_age <- lapply(cohort_names, function(co) {
    mf <- model.frame(mod_list_by_cohort[[co]])
    mf$ageatdx
  }) |> unlist()
  
  overall_age_mean <- mean(all_age, na.rm = TRUE)
  
  # build per-cohort reference curves
  ref_curves <- lapply(cohort_names, function(co) {
    mod <- mod_list_by_cohort[[co]]
    mf  <- model.frame(mod)
    
    stage_hi <- pick_high_stage(mf$stage_cat)
    
    nd <- data.frame(
      ageatdx = overall_age_mean,  # <-- changed from cohort mean to overall mean
      stage_cat = factor(stage_hi, levels = levels(mf$stage_cat)),
      abund_binary = 0,
      fpc1_cat = factor("high", levels = levels(mf$fpc1_cat))
    )
    
    sf <- survfit(mod, newdata = nd)
    
    data.frame(
      cohort = co,
      time = sf$time,
      Sref = sf$surv
    )
  }) %>% dplyr::bind_rows()
  
  # common time grid
  tmax <- max(ref_curves$time, na.rm = TRUE)
  grid <- seq(0, tmax, by = time_step)
  
  # convert to hazards and interpolate onto grid, then average hazards (weighted by N)
  hazards <- lapply(cohort_names, function(co) {
    mod <- mod_list_by_cohort[[co]]
    mf  <- model.frame(mod)
    nco <- nrow(mf)
    
    d <- ref_curves %>%
      dplyr::filter(cohort == co) %>%
      dplyr::mutate(Href = -log(pmax(Sref, 1e-12)))  # guard against log(0)
    
    Hgrid <- approx(x = d$time, y = d$Href, xout = grid, rule = 2)$y
    data.frame(cohort = co, n = nco, time = grid, Href = Hgrid)
  }) %>% dplyr::bind_rows()
  
  meta_href <- hazards %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      Href = weighted.mean(Href, w = n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Sref = exp(-Href))
  
  # attach overall mean age so it can be reported easily
  attr(meta_href, "overall_age_mean") <- overall_age_mean
  
  meta_href
}

make_four_group_curves <- function(Sref_df, pooled_group_HRs) {
  # pooled_group_HRs: output of meta_analyze_groups() with columns Group, HR
  # add reference group explicitly
  pooled <- bind_rows(
    data.frame(Group = "Low abund + High clustering (ref)", HR = 1, LCL = 1, UCL = 1),
    pooled_group_HRs
  )
  
  curves <- lapply(seq_len(nrow(pooled)), function(i) {
    g  <- pooled$Group[i]
    hr <- pooled$HR[i]
    lcl <- pooled$LCL[i]
    ucl <- pooled$UCL[i]
    
    # point estimate
    surv_hat <- Sref_df$Sref ^ hr
    
    # 95% CI via HR limits
    surv_lo  <- Sref_df$Sref ^ ucl
    surv_hi  <- Sref_df$Sref ^ lcl
    
    data.frame(
      time = Sref_df$time,
      Group = g,
      Surv  = surv_hat,
      Surv_lo = surv_lo,
      Surv_hi = surv_hi
      
    )
  }) %>% bind_rows()
  
  list(curves = curves, hr_table = pooled)
}

extract_anova_pvalues <- function(mod) {
  a <- anova(mod)
  a_df <- as.data.frame(a)
  a_df$term <- rownames(a_df)
  rownames(a_df) <- NULL
  
  # try to identify a p-value column more flexibly
  pcol <- grep("Pr\\(|P\\(", names(a_df), value = TRUE)
  
  # if still not found, try common alternatives
  if (length(pcol) == 0) {
    pcol <- grep("Pr|p", names(a_df), value = TRUE, ignore.case = TRUE)
  }
  
  if (length(pcol) == 0) {
    stop(
      "No p-value column found in anova(mod). Column names were: ",
      paste(names(a_df), collapse = ", ")
    )
  }
  
  # if multiple possible matches, use the first one
  pcol <- pcol[1]
  
  data.frame(
    term = a_df$term,
    p_value = a_df[[pcol]]
  )
}

extract_ph_tests <- function(mod_list) {
  ph_out <- lapply(names(mod_list), function(co) {
    zph <- cox.zph(mod_list[[co]], transform = "km")
    out <- as.data.frame(zph$table)
    out$term <- rownames(out)
    rownames(out) <- NULL
    out$cohort <- co
    out
  }) %>% dplyr::bind_rows()
  
  ph_out
}


save_ph_plots <- function(mod_list, file_name) {
  pdf(file_name, width = 8, height = 10)
  for (co in names(mod_list)) {
    zph <- cox.zph(mod_list[[co]], transform = "km")
    plot(zph, main = paste("PH diagnostics:", co))
  }
  dev.off()
}

build_km_data <- function(mod_list_by_cohort) {
  km_dat <- lapply(names(mod_list_by_cohort), function(co) {
    mod <- mod_list_by_cohort[[co]]
    mf  <- model.frame(mod)
    
    y <- model.response(mf)
    mf$mossurv <- y[, 1]
    mf$death   <- y[, 2]
    mf$cohort  <- co
    
    mf <- mf %>%
      dplyr::mutate(
        Group = dplyr::case_when(
          abund_binary == 0 & fpc1_cat == "high" ~ "Low abund + High clustering (ref)",
          abund_binary == 0 & fpc1_cat == "low"  ~ "Low abund + Low clustering",
          abund_binary == 1 & fpc1_cat == "high" ~ "High abund + High clustering",
          abund_binary == 1 & fpc1_cat == "low"  ~ "High abund + Low clustering",
          TRUE ~ NA_character_
        )
      )
    
    mf
  }) %>% dplyr::bind_rows()
  
  km_dat$Group <- factor(
    km_dat$Group,
    levels = c(
      "Low abund + High clustering (ref)",
      "Low abund + Low clustering",
      "High abund + High clustering",
      "High abund + Low clustering"
    )
  )
  
  km_dat
}

group_colors <- c(
  "Low abund + High clustering (ref)" = "#e7298a",
  "Low abund + Low clustering"  = "#d95f02",
  "High abund + High clustering" = "#7570b3",
  "High abund + Low clustering"  = "#1b9e77"
)

run_metric_survival_curves <- function(step_wise_models, metric_name) {
  
  # 1) pull the per-cohort model object
  # use the saved "FPC Score Interactions" model (then refit fpc1-only)
  base_models <- lapply(step_wise_models[[metric_name]], function(x) x$models$`FPC Score Interactions`)
  base_models <- base_models[!sapply(base_models, is.null)]
  
  # 2) refit fpc1-only model in each cohort
  fpc1_models <- lapply(base_models, refit_fpc1_only)
  
  anova_pvals <- lapply(names(fpc1_models), function(co) {
    out <- extract_anova_pvalues(fpc1_models[[co]])
    out$cohort <- co
    out
  }) %>% bind_rows()
  
  ph_tests <- extract_ph_tests(fpc1_models)
  
  
  # 3) cohort-specific group contrasts
  cohort_contrasts <- lapply(names(fpc1_models), function(co) {
    out <- extract_group_contrasts(fpc1_models[[co]])
    out$cohort <- co
    out
  }) %>% bind_rows()
  
  # 4) meta-analyze each group contrast
  pooled <- meta_analyze_groups(cohort_contrasts)
  
  # 5) build meta reference curve and transform to 4 curves
  Sref <- make_meta_reference_curve(fpc1_models, time_step = 1)
  overall_age_mean <- attr(Sref, "overall_age_mean")
  out  <- make_four_group_curves(Sref, pooled)
  
  # 6) plot
  p <- ggplot(out$curves, aes(x = time, y = Surv, color = Group)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = group_colors) +
    theme_minimal() +
    labs(
      title = paste0(metric_name, " (Model 4: abundance × FPC1.b)"),
      x = "Time",
      y = "Predicted survival probability"
    )
  
  p2 <- ggplot(out$curves, aes(x = time, y = Surv, color = Group, fill = Group)) +
    geom_ribbon(aes(ymin = Surv_lo, ymax = Surv_hi), alpha = 0.18, color = NA) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    theme_minimal() +
    labs(
      title = paste0(metric_name, " (Model 4: abundance × FPC1.b)"),
      x = "Time",
      y = "Predicted survival probability"
    )
  
  
  list(
    plot = p, plot2 = p2, hr_table = out$hr_table, pooled_contrasts = pooled,
    cohort_contrasts = cohort_contrasts, anova_pvals = anova_pvals,
    overall_age_mean = overall_age_mean,
    ph_tests = ph_tests
  )
}

summarize_censoring <- function(dat, cohort_name) {
  data.frame(
    cohort = cohort_name,
    n = nrow(dat),
    events = sum(dat$death == 1, na.rm = TRUE),
    censored = sum(dat$death == 0, na.rm = TRUE),
    censoring_pct = mean(dat$death == 0, na.rm = TRUE) * 100,
    median_followup = median(dat$mossurv, na.rm = TRUE)
  )
}

censor_summary <- bind_rows(
  summarize_censoring(AAS_summary_dat_list, "AACES"),
  summarize_censoring(NCO_summary_dat_list, "NCOCS"),
  summarize_censoring(NECC_summary_dat_list, "NECC"),
  summarize_censoring(NHSI_summary_dat_list, "NHSI"),
  summarize_censoring(NHSII_summary_dat_list, "NHSII")
)

censor_summary



res_cd3   <- run_metric_survival_curves(step_wise_models, metric_name = "CD3+ Gest")
res_cd3cd8   <- run_metric_survival_curves(step_wise_models, metric_name = "CD3+ CD8+ Gest")

res_cd3$plot <- res_cd3$plot +
  labs(title = "CD3+ G-est (Model 4: FPC1-Only (binary))")
res_cd3cd8$plot <- res_cd3cd8$plot +
  labs(title = "CD3+ CD8+ G-est (Model 4: FPC1-Only (binary))")

res_cd3$plot2 <- res_cd3$plot2 +
  labs(title = "CD3+ G-est (Model 4: FPC1-Only (binary))")
res_cd3cd8$plot2 <- res_cd3cd8$plot2 +
  labs(title = "CD3+ CD8+ G-est (Model 4: FPC1-Only (binary))")

res_cd3$plot
res_cd3cd8$plot

res_cd3$plot2
res_cd3cd8$plot2

res_cd3$hr_table
res_cd3cd8$hr_table

res_cd3$overall_age_mean
res_cd3cd8$overall_age_mean

res_cd3$ph_tests
res_cd3cd8$ph_tests


ggsave("model4_cd3_predict_survival_curves.tif", res_cd3$plot, width = 8, height = 5, dpi = 300)
ggsave("model4_cd3cd8_predict_survival_curves.tif", res_cd3cd8$plot, width = 8, height = 5, dpi = 300)

ggsave("model4_cd3_predict_survival_curves_wCI.tif", res_cd3$plot2, width = 8, height = 5, dpi = 300)
ggsave("model4_cd3cd8_predict_survival_curves_wCI.tif", res_cd3cd8$plot2, width = 8, height = 5, dpi = 300)


fpc1_models_cd3 <- lapply(step_wise_models[["CD3+ Gest"]], function(x) 
  refit_fpc1_only(x$models$`FPC Score Interactions`))

fpc1_models_cd3cd8 <- lapply(step_wise_models[["CD3+ CD8+ Gest"]], function(x) 
  refit_fpc1_only(x$models$`FPC Score Interactions`))

save_ph_plots(fpc1_models_cd3, "PH_diagnostics_CD3_Gest.pdf")
save_ph_plots(fpc1_models_cd3cd8, "PH_diagnostics_CD3CD8_Gest.pdf")

km_dat_cd3 <- build_km_data(fpc1_models_cd3)
km_dat_cd3cd8 <- build_km_data(fpc1_models_cd3cd8)

km_fit_cd3 <- survfit(Surv(mossurv, death) ~ Group, data = km_dat_cd3)
km_fit_cd3cd8 <- survfit(Surv(mossurv, death) ~ Group, data = km_dat_cd3cd8)

cols <- group_colors[levels(km_dat_cd3$Group)]
plot(km_fit_cd3,
     col = cols, lwd = 2, mark.time = TRUE,
     xlab = "Time (months)", ylab = "Survival probability",
     main = "CD3+ empirical Kaplan–Meier curves")
legend("topright", legend = levels(km_dat_cd3$Group), col = cols, lwd = 2, bty = "n")

plot(km_fit_cd3cd8,
     col = cols, lwd = 2, mark.time = TRUE,
     xlab = "Time (months)", ylab = "Survival probability",
     main = "CD3+CD8+ empirical Kaplan–Meier curves")
legend("topright", legend = levels(km_dat_cd3cd8$Group), col = cols, lwd = 2, bty = "n")


tiff("KM_CD3_groups.tiff", width = 8, height = 5, units = "in", res = 300, compression = "lzw")
plot(km_fit_cd3,
     col = cols, lwd = 2, mark.time = TRUE,
     xlab = "Time (months)", ylab = "Survival probability",
     main = "CD3+ empirical Kaplan–Meier curves")
legend("topright", legend = levels(km_dat_cd3$Group), col = cols, lwd = 2, bty = "n")
dev.off()

tiff("KM_CD3CD8_groups.tiff", width = 8, height = 5, units = "in", res = 300, compression = "lzw")
plot(km_fit_cd3cd8,
     col = cols, lwd = 2, mark.time = TRUE,
     xlab = "Time (months)", ylab = "Survival probability",
     main = "CD3+CD8+ empirical Kaplan–Meier curves")
legend("topright", legend = levels(km_dat_cd3cd8$Group), col = cols, lwd = 2, bty = "n")
dev.off()








