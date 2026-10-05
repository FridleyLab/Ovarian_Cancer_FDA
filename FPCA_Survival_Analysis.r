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
    analysis_dat_list2[[m]] |>
      dplyr::mutate(dplyr::across(c(block, histotype),
                                  ~ as.character(.x)),
                    stage_cat = ifelse(stage_cat == 1, "high stage", "low state"))
  )
  #abundance is not cohort dependent
  #func_observed is not cohort dependent
  cohort_dat = split(dat, dat$cohort)
  return(cohort_dat)
})

names(analysis_dat_list_combined) <- metrics

build_fixed_radius_five_group_data <- function(
    cohort_dat,
    radius_value,
    cutoff = c("median", "zero")
) {
  
  cutoff <- match.arg(cutoff)
  
  required_variables <- c(
    "block",
    "cohort",
    "death",
    "mossurv",
    "ageatdx",
    "stage_cat",
    "abund_binary",
    "func_observed",
    "r",
    "func_weighted"
  )
  
  missing_variables <- setdiff(
    required_variables,
    names(cohort_dat)
  )
  
  if (length(missing_variables) > 0) {
    stop(
      "Missing required variables: ",
      paste(missing_variables, collapse = ", ")
    )
  }
  
  # One metadata row per subject, including subjects without
  # an estimable spatial trajectory
  all_subject_metadata <- cohort_dat |>
    dplyr::select(
      block,
      cohort,
      death,
      mossurv,
      ageatdx,
      stage_cat,
      abund_binary,
      func_observed
    ) |>
    dplyr::distinct(block, .keep_all = TRUE)
  
  # Fixed-radius spatial value among subjects with an estimable trajectory
  radius_data <- cohort_dat |>
    dplyr::filter(func_observed == 1, r == radius_value) |>
    dplyr::select(block, func_weighted) |>
    dplyr::distinct()
  
  duplicate_check <- radius_data |>
    dplyr::count(block) |>
    dplyr::filter(n > 1)
  
  if (nrow(duplicate_check) > 0) {
    stop("More than one fixed-radius value was found per subject at radius ", radius_value, ".")
  }
  
  if (nrow(radius_data) == 0) {
    stop("No estimable fixed-radius values were found at radius ", radius_value, ".")
  }
  
  cutoff_value <- switch(
    cutoff,
    median = stats::median(radius_data$func_weighted, na.rm = TRUE),
    zero = 0
  )
  
  radius_data <- radius_data |>
    dplyr::mutate(
      fixed_cluster_cat = factor(
        ifelse(func_weighted >= cutoff_value, "high", "low"),
        levels = c("high", "low")
      )
    )
  
  five_group_data <- all_subject_metadata |>
    dplyr::left_join(radius_data, by = "block") |>
    dplyr::mutate(
      stage_cat = factor(as.character(stage_cat), levels = c("high stage", "low stage")),
      abund_binary = as.numeric(as.character(abund_binary)),
      
      Group5 = dplyr::case_when(
        func_observed == 0 | is.na(func_weighted) ~ "Low Immune Activity (ref)",
        abund_binary == 0 & fixed_cluster_cat == "high" ~ "Low abundance + high clustering",
        abund_binary == 0 & fixed_cluster_cat == "low" ~ "Low abundance + low clustering",
        abund_binary == 1 & fixed_cluster_cat == "high" ~ "High abundance + high clustering",
        abund_binary == 1 & fixed_cluster_cat == "low" ~ "High abundance + low clustering",
        TRUE ~ NA_character_
      ),
      
      Group5 = factor(
        Group5,
        levels = c(
          "Low Immune Activity (ref)",
          "Low abundance + high clustering",
          "Low abundance + low clustering",
          "High abundance + high clustering",
          "High abundance + low clustering"
        )
      )
    )
  
  attr(five_group_data, "cutoff_value") <- cutoff_value
  attr(five_group_data, "radius") <- radius_value
  attr(five_group_data, "cutoff_method") <- cutoff
  
  five_group_data
}

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
    l = l |>
      dplyr::mutate(sample = block) #because mxfda expects there to be image/sample IDs and subject IDs, duplicate
    clin = l |>
      dplyr::select(-c(func_weighted, r)) |>
      dplyr::distinct()
    func = l |>
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
    obj = run_fpca(obj, metric = 'uni g', r = 'r', value = 'func_weighted', pve = 0.99)
    summary(obj)
    return(obj)
  })
  return(co_mxfda_objects)
})

#plot the fpcs
pls = lapply(names(mxfda_objects), function(nam){
  co_plots = lapply(names(mxfda_objects[[nam]]), function(co){
    obj = mxfda_objects[[nam]][[co]]
    eigvals = obj@functional_pca$Gest$fpc_object$evalues
    pve = eigvals / sum(eigvals) * 100
    p1 = plot(obj, what = 'uni g fpca', pc_choice = 1) + lims(x = c(0, 100)) #show that not using 0
    p2 = plot(obj, what = 'uni g fpca', pc_choice = 2) + lims(x = c(0, 100))
    plots = ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2) |>
      ggpubr::annotate_figure(left = ggpubr::text_grob(co, rot = 90))
    data = dplyr::bind_rows(
      p1$data |>
        dplyr::mutate(FPC = "fpc1", VarExp = pve[1]),
      p2$data |>
        dplyr::mutate(FPC = "fpc2", VarExp = pve[2])
    ) |>
      dplyr::mutate(Cohort = co)
    #want to return the data for flipping directions later
    list(plot = plots,
         data = data)
  })
  plist = lapply(co_plots, function(x) x$plot)
  data = lapply(co_plots, function(x) x$data) |> do.call(dplyr::bind_rows, .)
  pl_panel = ggpubr::ggarrange(plotlist = plist, nrow = 5, ncol = 1) |>
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
  }) |>
    do.call(dplyr::bind_rows, .)
  label_dat <- fixed_fpca_dat |>
    group_by(Cohort, FPC) |>
    summarise(VarExp = first(VarExp), .groups = "drop") |>
    mutate(label = paste0("Variance = ", round(VarExp, 1), "%"))
  pl = fixed_fpca_dat |>
    ggplot() + 
    geom_point(aes(x = index, y = mu)) + 
    facet_grid(Cohort ~ FPC) +
    geom_point(aes(x = index, y = plus), color = "blue", shape = "+", size = 2) +
    geom_point(aes(x = index, y = minus), color = "red", shape = "-", size = 2) +
    theme_minimal() +
    geom_text(data = label_dat, aes(x = 5, y = Inf, label = label),
      hjust = 0, vjust = -0.3, size = 3.2, inherit.aes = FALSE) +
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5),
      panel.spacing = unit(1.5, "lines"),
      strip.text.x = element_text(
        margin = margin(b = 10), # move labels upward
        size = 12
      )
    ) +
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
    obj2@Metadata = obj2@Metadata |>
      dplyr::filter(func_observed == 1)
    obj2@univariate_summaries$Gest = obj2@univariate_summaries$Gest |>
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
    data = dplyr::bind_rows(p1$data |> dplyr::mutate(FPC = "fpc1"),
                            p2$data |> dplyr::mutate(FPC = "fpc2")) |>
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
    }) |>
      do.call(dplyr::bind_rows, .)
    
    #get the adjustment value
    adj_val = data_fixed |>
      dplyr::select(Cohort, FPC, adjust) |>
      dplyr::distinct()
    
    fpc_dat = dplyr::inner_join(obj2@Metadata, obj2@functional_pca$Gest$score_df) 
    
    for(fpc in adj_val$FPC){
      fpc_dat[[fpc]] = fpc_dat[[fpc]] * adj_val$adjust[adj_val$FPC == fpc]
    }
    
    fpc_dat2 = fpc_dat |>
      dplyr::mutate(dplyr::across(dplyr::contains("fpc"),
                                  ~ factor(ifelse(.x >= 0, "high", "low"),
                                           levels = c('high', 'low')), #set reference
                                  .names = "{col}_cat"))
    
    # ------------------------------------------------------------
    # Create five-group data, including subjects with <6 cells
    # ------------------------------------------------------------
    
    # One metadata row per subject/sample, including subjects for whom
    # the spatial function could not be estimated
    all_subject_metadata <- meta |> dplyr::distinct(block, .keep_all = TRUE)
    
    # Retain only the identifiers and sign-aligned FPC1 category from
    # subjects with an estimable spatial trajectory
    fpc1_group_data <- fpc_dat2 |>
      dplyr::select(dplyr::any_of(c("block", "sample")), fpc1, fpc1_cat) |>
      dplyr::distinct()
    
    # Join FPC1 information back to the complete cohort metadata
    five_group_dat <- all_subject_metadata |>
      dplyr::left_join(
        fpc1_group_data,
        by = intersect(
          c("block", "sample"),
          intersect(names(all_subject_metadata), names(fpc1_group_data))
        )
      ) |>
      dplyr::mutate(
        Group5 = dplyr::case_when(
          func_observed == 0 ~ "Low Immune Activity (ref)",
          func_observed == 1 & abund_binary == 0 & fpc1_cat == "high" ~ "Low abundance + high clustering",
          func_observed == 1 & abund_binary == 0 & fpc1_cat == "low" ~ "Low abundance + low clustering",
          func_observed == 1 & abund_binary == 1 & fpc1_cat == "high" ~ "High abundance + high clustering",
          func_observed == 1 & abund_binary == 1 & fpc1_cat == "low" ~ "High abundance + low clustering",
          TRUE ~ NA_character_
        ),
        Group5 = factor(
          Group5,
          levels = c(
            "Low Immune Activity (ref)",
            "Low abundance + high clustering",
            "Low abundance + low clustering",
            "High abundance + high clustering",
            "High abundance + low clustering"
          )
        )
      )
    
    # Fit the five-category model with the <6-cell group as reference
    mod4_five_group <- coxph(
      survival::Surv(mossurv, death) ~ ageatdx + stage_cat + Group5,
      data = five_group_dat,
      x = TRUE,
      model = TRUE
    )
    
    fpc_dat2_summ = fpc_dat2 |>
      dplyr::summarise(dplyr::across(dplyr::contains("fpc"),
                                     .fns = list(high = function(x){sum(x == "high")},
                                                 low = function(x){sum(x == "low")}),
                                     .names = "{col}_{fn}")) |>
      dplyr::select(dplyr::contains("cat")) |>
      tidyr::pivot_longer(cols = dplyr::everything(),
                          names_to = c("fpc", "w", "cat"),
                          names_pattern = "(.*)_(.*)_(.*)") |>
      dplyr::select(-w) |>
      tidyr::pivot_wider(names_from = "cat", values_from = "value")
    
    mod2 = coxph(
      Surv(mossurv, death) ~ ageatdx + stage_cat +
                   abund_binary * fpc1_cat + abund_binary * fpc2_cat,
      data = fpc_dat2)
    mod3 = coxph(
      Surv(mossurv, death) ~ ageatdx + stage_cat +
                   abund_binary + fpc1_cat + fpc2_cat,
      data = fpc_dat2)
    
    return(list(
      models = list(
        `Abundance and Func Observed` = mod1,
        `FPC Score Interactions` = mod2,
        `FPC Score No Interactions` = mod3,
        `Five Group Model` = mod4_five_group
      ),
      data = list(
        `High Abundance Plot` = data_fixed,
        `High Abundance Scores` = fpc_dat,
        `High Abundance Score Cats` = fpc_dat2_summ,
        `Five Group Data` = five_group_dat
      )
    ))
  })
})

#high abundance fpca plot data
fpca_plots = lapply(step_wise_models, function(m){
  pl_dat = lapply(m, function(x) x$data$`High Abundance Plot`) |> 
    do.call(dplyr::bind_rows, .)
  pl_dat |>
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
      
      covars = data.frame(res$coefficients, check.names = FALSE) |>
        tibble::rownames_to_column("factor") |>
        dplyr::inner_join(data.frame(res$conf.int, check.names = FALSE) |>
                            tibble::rownames_to_column("factor")
        )
      model = data.frame(t(res$waldtest), check.names = FALSE) |>
        dplyr::mutate(info = 'waldtest', .before = 1)
      model = model[rep(1, nrow(covars)),]
      row.names(model) = NULL
      dplyr::bind_cols(covars, model) |>
        dplyr::mutate(model = n, .before = 1)
    })
    names(tmp) = names(co$models)
    return(tmp)
  })
  abund_func_mods = lapply(names(mod_summaries), function(co){
    mod_summaries[[co]]$`Abundance and Func Observed` |>
      dplyr::mutate(cohort = co)
  }) |>
    do.call(dplyr::bind_rows, .)
  func_int_mods = lapply(names(mod_summaries), function(co){
    mod_summaries[[co]]$`FPC Score Interactions` |>
      dplyr::mutate(cohort = co)
  }) |>
    do.call(dplyr::bind_rows, .)
  func_no_int_mods = lapply(names(mod_summaries), function(co){
    mod_summaries[[co]]$`FPC Score No Interactions` |>
      dplyr::mutate(cohort = co)
  }) |>
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
    }) |>
      dplyr::bind_rows(.id = "Model") |>
      tidyr::pivot_wider(names_from = "Model", values_from = "Samples")
  }) |>
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
  df |>
    dplyr::group_by(Group) |>
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
    }) |>
    dplyr::ungroup()
}

meta_analyze_five_groups <- function(df) {
  
  df |>
    dplyr::filter(is.finite(logHR), is.finite(SE), SE > 0) |>
    dplyr::group_by(Group) |>
    dplyr::group_modify(~{
      
      dat_g <- .x
      
      if (nrow(dat_g) < 2) {
        return(
          data.frame(
            logHR = NA_real_,
            SE = NA_real_,
            z = NA_real_,
            p = NA_real_,
            HR = NA_real_,
            LCL = NA_real_,
            UCL = NA_real_,
            tau2 = NA_real_,
            I2 = NA_real_,
            Q = NA_real_,
            Q_p = NA_real_,
            k = nrow(dat_g)
          )
        )
      }
      
      fit <- metafor::rma(yi = dat_g$logHR, sei = dat_g$SE, method = "REML")
      
      data.frame(
        logHR = as.numeric(fit$b),
        SE = as.numeric(fit$se),
        z = as.numeric(fit$zval),
        p = as.numeric(fit$pval),
        HR = exp(as.numeric(fit$b)),
        LCL = exp(as.numeric(fit$ci.lb)),
        UCL = exp(as.numeric(fit$ci.ub)),
        tau2 = as.numeric(fit$tau2),
        I2 = as.numeric(fit$I2),
        Q = as.numeric(fit$QE),
        Q_p = as.numeric(fit$QEp),
        k = fit$k
      )
    }) |>
    dplyr::ungroup()
}

extract_five_group_effects <- function(mod) {
  
  sm <- summary(mod)
  
  as.data.frame(sm$coefficients, check.names = FALSE) |>
    tibble::rownames_to_column("term") |>
    dplyr::filter(grepl("^Group5", term)) |>
    dplyr::transmute(
      Group = sub("^Group5", "", term),
      logHR = coef,
      SE = `se(coef)`
    )
}


# utility: pick the "high stage" level robustly
pick_high_stage <- function(stage_factor) {
  
  lev <- unique(as.character(stage_factor))
  
  hi <- lev[grepl("high", lev, ignore.case = TRUE)][1]
  
  if(length(hi) == 0 || is.na(hi)) {
    stop("Could not find a 'high' level in stage_cat. Values were: ", paste(lev, collapse = ", "))
  }
  
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
      ageatdx = overall_age_mean,
      stage_cat = stage_hi,
      abund_binary = 0,
      fpc1_cat = factor("high", levels = levels(mf$fpc1_cat))
    )
    
    sf <- survfit(mod, newdata = nd)
    
    data.frame(cohort = co, time = sf$time, Sref = sf$surv)
  }) |> dplyr::bind_rows()
  
  # common time grid
  tmax <- max(ref_curves$time, na.rm = TRUE)
  grid <- seq(0, tmax, by = time_step)
  
  # convert to hazards and interpolate onto grid, then average hazards (weighted by N)
  hazards <- lapply(cohort_names, function(co) {
    mod <- mod_list_by_cohort[[co]]
    mf  <- model.frame(mod)
    nco <- nrow(mf)
    
    d <- ref_curves |>
      dplyr::filter(cohort == co) |>
      dplyr::mutate(Href = -log(pmax(Sref, 1e-12)))  # guard against log(0)
    
    Hgrid <- approx(x = d$time, y = d$Href, xout = grid, rule = 2)$y
    data.frame(cohort = co, n = nco, time = grid, Href = Hgrid)
  }) |> dplyr::bind_rows()
  
  meta_href <- hazards |>
    dplyr::group_by(time) |>
    dplyr::summarise(Href = weighted.mean(Href, w = n, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(Sref = exp(-Href))
  
  # attach overall mean age so it can be reported easily
  attr(meta_href, "overall_age_mean") <- overall_age_mean
  
  meta_href
}

make_five_group_reference_curve <- function(
    mod_list_by_cohort,
    time_step = 1
) {
  
  cohort_names <- names(mod_list_by_cohort)
  
  # Overall mean age among subjects included in the five-group models
  all_age <- lapply(cohort_names, function(co) {
    mf <- model.frame(mod_list_by_cohort[[co]])
    mf$ageatdx
  }) |>
    unlist()
  
  overall_age_mean <- mean(all_age, na.rm = TRUE)
  
  reference_label <- "Low Immune Activity (ref)"
  
  # Cohort-specific adjusted reference curves
  ref_curves <- lapply(cohort_names, function(co) {
    
    mod <- mod_list_by_cohort[[co]]
    mf <- model.frame(mod)
    
    stage_hi <- pick_high_stage(mf$stage_cat)
    
    if (!reference_label %in% levels(mf$Group5)) {
      stop("Reference group is absent from Group5 levels for cohort: ", co)
    }
    
    nd <- data.frame(
      ageatdx = overall_age_mean,
      stage_cat = stage_hi,
      Group5 = factor(reference_label, levels = levels(mf$Group5))
    )
    
    sf <- survival::survfit(mod, newdata = nd)
    
    data.frame(cohort = co, time = sf$time, Sref = sf$surv)
  }) |>
    dplyr::bind_rows()
  
  tmax <- max(ref_curves$time, na.rm = TRUE)
  grid <- seq(0, tmax, by = time_step)
  
  hazards <- lapply(cohort_names, function(co) {
    
    mod <- mod_list_by_cohort[[co]]
    mf <- model.frame(mod)
    nco <- nrow(mf)
    
    d <- ref_curves |>
      dplyr::filter(cohort == co) |>
      dplyr::mutate(Href = -log(pmax(Sref, 1e-12)))
    
    Hgrid <- approx(x = d$time, y = d$Href, xout = grid,rule = 2)$y
    
    data.frame(
      cohort = co,
      n = nco,
      time = grid,
      Href = Hgrid
    )
  }) |>
    dplyr::bind_rows()
  
  meta_href <- hazards |>
    dplyr::group_by(time) |>
    dplyr::summarise(Href = weighted.mean(Href, w = n, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(Sref = exp(-Href))
  
  attr(meta_href, "overall_age_mean") <- overall_age_mean
  
  meta_href
}

make_five_group_curves <- function(
    Sref_df,
    pooled_group_HRs
) {
  
  reference_label <- "Low Immune Activity (ref)"
  
  pooled <- dplyr::bind_rows(
    data.frame(
      Group = reference_label,
      logHR = 0,
      SE = 0,
      z = NA_real_,
      p = NA_real_,
      HR = 1,
      LCL = 1,
      UCL = 1,
      tau2 = NA_real_,
      I2 = NA_real_,
      Q = NA_real_,
      Q_p = NA_real_,
      k = NA_integer_
    ),
    pooled_group_HRs
  )
  
  curves <- lapply(seq_len(nrow(pooled)), function(i) {
    
    g <- pooled$Group[i]
    hr <- pooled$HR[i]
    lcl <- pooled$LCL[i]
    ucl <- pooled$UCL[i]
    
    if (!is.finite(hr)) {
      return(NULL)
    }
    
    surv_hat <- Sref_df$Sref ^ hr
    
    if (is.finite(lcl) && is.finite(ucl)) {
      surv_lo <- Sref_df$Sref ^ ucl
      surv_hi <- Sref_df$Sref ^ lcl
    } else {
      surv_lo <- NA_real_
      surv_hi <- NA_real_
    }
    
    data.frame(
      time = Sref_df$time,
      Group = g,
      Surv = surv_hat,
      Surv_lo = surv_lo,
      Surv_hi = surv_hi
    )
  }) |>
    dplyr::bind_rows()
  
  curves$Group <- factor(
    curves$Group,
    levels = c(
      reference_label,
      "Low abundance + high clustering",
      "Low abundance + low clustering",
      "High abundance + high clustering",
      "High abundance + low clustering"
    )
  )
  
  list(
    curves = curves,
    hr_table = pooled
  )
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
  }) |> bind_rows()
  
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
    stop("No p-value column found in anova(mod). Column names were: ", paste(names(a_df), collapse = ", "))
  }
  
  # if multiple possible matches, use the first one
  pcol <- pcol[1]
  
  data.frame(term = a_df$term, p_value = a_df[[pcol]])
}

extract_ph_tests <- function(mod_list) {
  ph_out <- lapply(names(mod_list), function(co) {
    zph <- cox.zph(mod_list[[co]], transform = "km")
    out <- as.data.frame(zph$table)
    out$term <- rownames(out)
    rownames(out) <- NULL
    out$cohort <- co
    out
  }) |> dplyr::bind_rows()
  
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
    
    mf <- mf |>
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
  }) |> dplyr::bind_rows()
  
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

five_group_colors <- c(
  "Low Immune Activity (ref)" = "#666666",
  "Low abundance + high clustering" = "#e7298a",
  "Low abundance + low clustering" = "#d95f02",
  "High abundance + high clustering" = "#7570b3",
  "High abundance + low clustering" = "#1b9e77"
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
  }) |> bind_rows()
  
  ph_tests <- extract_ph_tests(fpc1_models)
  
  
  # 3) cohort-specific group contrasts
  cohort_contrasts <- lapply(names(fpc1_models), function(co) {
    out <- extract_group_contrasts(fpc1_models[[co]])
    out$cohort <- co
    out
  }) |> bind_rows()
  
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

run_metric_five_group_analysis <- function(
    step_wise_models,
    metric_name,
    time_step = 1
) {
  
  metric_results <- step_wise_models[[metric_name]]
  
  # Pull the newly fitted five-group Cox model from each cohort
  five_group_models <- lapply(metric_results, function(x) x$models$`Five Group Model`)
  
  five_group_models <- five_group_models[!vapply(five_group_models, is.null, logical(1))]
  
  if (length(five_group_models) == 0) {
    stop("No five-group models were found for metric: ", metric_name,
      ". Rerun step_wise_models after adding the Five Group Model."
    )
  }
  
  # Group sample and event counts by cohort
  group_counts <- lapply(names(five_group_models), function(co) {
    
    mf <- model.frame(five_group_models[[co]])
    y <- model.response(mf)
    
    data.frame(Group = mf$Group5, death = y[, 2]) |>
      dplyr::group_by(Group, .drop = FALSE) |>
      dplyr::summarise(
        n = dplyr::n(),
        events = sum(death == 1, na.rm = TRUE),
        censored = sum(death == 0, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(cohort = co, .before = 1)
  }) |>
    dplyr::bind_rows()
  
  # Cohort-specific HRs versus the <6-cell reference
  cohort_effects <- lapply(names(five_group_models), function(co) {
    
    out <- extract_five_group_effects(five_group_models[[co]])
    out$cohort <- co
    out
  }) |>
    dplyr::bind_rows()
  
  # Random-effects meta-analysis for each of the four contrasts
  pooled_effects <- meta_analyze_five_groups(cohort_effects)
  
  # Construct the illustrative reference survival curve
  Sref <- make_five_group_reference_curve(five_group_models, time_step = time_step)
  
  overall_age_mean <- attr(Sref, "overall_age_mean")
  
  # Apply pooled HRs to the <6-cell reference curve
  curve_output <- make_five_group_curves(Sref_df = Sref, pooled_group_HRs = pooled_effects)
  
  # PH tests
  ph_tests <- extract_ph_tests(five_group_models)
  
  # Overall tests for the five-level group factor
  # Overall likelihood-ratio test for the five-level group factor
  group_tests <- lapply(names(five_group_models), function(co) {
    
    original_mod <- five_group_models[[co]]
    mf <- model.frame(original_mod)
    
    # Recover survival time and event status from the stored Surv response
    y <- model.response(mf)
    
    analysis_dat <- data.frame(
      mossurv = y[, 1],
      death = y[, 2],
      ageatdx = mf$ageatdx,
      stage_cat = mf$stage_cat,
      Group5 = mf$Group5
    )
    
    # Reduced model without the five-level group variable
    reduced_mod <- survival::coxph(
      survival::Surv(mossurv, death) ~ ageatdx + stage_cat,
      data = analysis_dat,
      ties = original_mod$method,
      x = TRUE,
      model = TRUE
    )
    
    # Refit the full model from the same stored analysis data
    full_mod <- survival::coxph(
      survival::Surv(mossurv, death) ~ ageatdx + stage_cat + Group5,
      data = analysis_dat,
      ties = original_mod$method,
      x = TRUE,
      model = TRUE
    )
    
    lrt <- anova(reduced_mod, full_mod, test = "LRT")
    
    data.frame(
      cohort = co,
      df = as.numeric(lrt[2, "Df"]),
      chisq = as.numeric(lrt[2, "Chisq"]),
      p = as.numeric(lrt[2, "Pr(>|Chi|)"])
    )
  }) |>
    dplyr::bind_rows()
  
  p <- ggplot2::ggplot(
    curve_output$curves,
    ggplot2::aes(
      x = time,
      y = Surv,
      color = Group
    )
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::scale_color_manual(values = five_group_colors, drop = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0(metric_name, ": five-group standardized survival curves"),
      x = "Time (months)",
      y = "Predicted survival probability",
      color = "Group"
    )
  
  p_ci <- ggplot2::ggplot(
    curve_output$curves,
    ggplot2::aes(
      x = time,
      y = Surv,
      color = Group,
      fill = Group
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = Surv_lo,
        ymax = Surv_hi
      ),
      alpha = 0.12,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::scale_color_manual(values = five_group_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = five_group_colors, drop = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0(metric_name, ": five-group standardized survival curves"),
      x = "Time (months)",
      y = "Predicted survival probability",
      color = "Group",
      fill = "Group"
    )
  
  list(
    models = five_group_models,
    plot = p,
    plot_ci = p_ci,
    curves = curve_output$curves,
    hr_table = curve_output$hr_table,
    pooled_effects = pooled_effects,
    cohort_effects = cohort_effects,
    group_counts = group_counts,
    group_tests = group_tests,
    overall_age_mean = overall_age_mean,
    ph_tests = ph_tests,
    reference_curve = Sref
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

res5_cd3 <- run_metric_five_group_analysis(
  step_wise_models = step_wise_models,
  metric_name = "CD3+ Gest",
  time_step = 1
)

res5_cd3cd8 <- run_metric_five_group_analysis(
  step_wise_models = step_wise_models,
  metric_name = "CD3+ CD8+ Gest",
  time_step = 1
)


fit_fixed_radius_five_group_model <- function(
    cohort_dat,
    radius_value,
    cutoff = "median"
) {
  
  radius_dat <-
    build_fixed_radius_five_group_data(
      cohort_dat = cohort_dat,
      radius_value = radius_value,
      cutoff = cutoff
    ) |>
    dplyr::filter(
      !is.na(mossurv),
      !is.na(death),
      !is.na(ageatdx),
      !is.na(stage_cat),
      !is.na(Group5)
    )
  
  reference_group <- "Low Immune Activity (ref)"
  
  if (!reference_group %in% as.character(radius_dat$Group5)) {
    stop("Reference group was not present for cohort ", unique(radius_dat$cohort),
         " at radius ", radius_value, ".")
  }
  
  radius_dat$Group5 <- factor(
    radius_dat$Group5,
    levels = c(
      "Low Immune Activity (ref)",
      "Low abundance + high clustering",
      "Low abundance + low clustering",
      "High abundance + high clustering",
      "High abundance + low clustering"
    )
  )
  
  coxph(
    Surv(mossurv, death) ~ ageatdx + stage_cat + Group5,
    data = radius_dat,
    ties = "efron",
    x = TRUE,
    model = TRUE
  )
}


############################################################
# Fixed-Radius Five-Group Survival Curves
############################################################

fixed25_cd3 <- lapply(
  analysis_dat_list_combined[["CD3+ G"]],
  fit_fixed_radius_five_group_model,
  radius_value = 25
)

fixed50_cd3 <- lapply(
  analysis_dat_list_combined[["CD3+ G"]],
  fit_fixed_radius_five_group_model,
  radius_value = 50
)

fixed75_cd3 <- lapply(
  analysis_dat_list_combined[["CD3+ G"]],
  fit_fixed_radius_five_group_model,
  radius_value = 75
)

fixed25_cd3cd8 <- lapply(
  analysis_dat_list_combined[["CD3+ CD8+ G"]],
  fit_fixed_radius_five_group_model,
  radius_value = 25
)

fixed50_cd3cd8 <- lapply(
  analysis_dat_list_combined[["CD3+ CD8+ G"]],
  fit_fixed_radius_five_group_model,
  radius_value = 50
)

fixed75_cd3cd8 <- lapply(
  analysis_dat_list_combined[["CD3+ CD8+ G"]],
  fit_fixed_radius_five_group_model,
  radius_value = 75
)

cohort_effects_25_cd3 <- lapply(
  names(fixed25_cd3),
  function(co){
    out <- extract_five_group_effects(fixed25_cd3[[co]])
    out$cohort <- co
    out
  }
) |>
  bind_rows()

pooled_25_cd3 <- meta_analyze_five_groups(cohort_effects_25_cd3)

Sref_25_cd3 <- make_five_group_reference_curve(fixed25_cd3)

curve25_cd3 <- make_five_group_curves(Sref_25_cd3, pooled_25_cd3)

cohort_effects_50_cd3 <- lapply(
  names(fixed50_cd3),
  function(co){
    out <- extract_five_group_effects(fixed50_cd3[[co]])
    out$cohort <- co
    out
  }
) |>
  bind_rows()

pooled_50_cd3 <- meta_analyze_five_groups(cohort_effects_50_cd3)

Sref_50_cd3 <- make_five_group_reference_curve(fixed50_cd3)

curve50_cd3 <- make_five_group_curves(Sref_50_cd3, pooled_50_cd3)

cohort_effects_75_cd3 <- lapply(
  names(fixed75_cd3),
  function(co){
    out <- extract_five_group_effects(fixed75_cd3[[co]])
    out$cohort <- co
    out
  }
) |>
  bind_rows()

pooled_75_cd3 <- meta_analyze_five_groups(cohort_effects_75_cd3)

Sref_75_cd3 <- make_five_group_reference_curve(fixed75_cd3)

curve75_cd3 <- make_five_group_curves(Sref_75_cd3, pooled_75_cd3)

cohort_effects_25_cd3cd8 <- lapply(
  names(fixed25_cd3cd8),
  function(co){
    out <- extract_five_group_effects(fixed25_cd3cd8[[co]])
    out$cohort <- co
    out
  }
) |>
  bind_rows()

pooled_25_cd3cd8 <- meta_analyze_five_groups(cohort_effects_25_cd3cd8)

Sref_25_cd3cd8 <- make_five_group_reference_curve(fixed25_cd3cd8)

curve25_cd3cd8 <-make_five_group_curves(Sref_25_cd3cd8, pooled_25_cd3cd8)

cohort_effects_50_cd3cd8 <- lapply(
  names(fixed50_cd3cd8),
  function(co){
    out <- extract_five_group_effects(fixed50_cd3cd8[[co]])
    out$cohort <- co
    out
  }
) |>
  bind_rows()

pooled_50_cd3cd8 <- meta_analyze_five_groups(cohort_effects_50_cd3cd8)

Sref_50_cd3cd8 <- make_five_group_reference_curve(fixed50_cd3cd8)

curve50_cd3cd8 <- make_five_group_curves(Sref_50_cd3cd8, pooled_50_cd3cd8)

cohort_effects_75_cd3cd8 <- lapply(
  names(fixed75_cd3cd8),
  function(co){
    out <- extract_five_group_effects(fixed75_cd3cd8[[co]])
    out$cohort <- co
    out
  }
) |>
  bind_rows()

pooled_75_cd3cd8 <- meta_analyze_five_groups(cohort_effects_75_cd3cd8)

Sref_75_cd3cd8 <- make_five_group_reference_curve(fixed75_cd3cd8)

curve75_cd3cd8 <- make_five_group_curves(Sref_75_cd3cd8, pooled_75_cd3cd8)

plot25_cd3 <- ggplot(curve25_cd3$curves,
                     aes(x = time, y = Surv, color = Group)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = five_group_colors) +
  theme_minimal() +
  labs(
    title = "CD3+ G(25)",
    x = "Time (months)",
    y = "Predicted survival probability"
  )

plot50_cd3 <- ggplot(curve50_cd3$curves,
                     aes(x = time, y = Surv, color = Group)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = five_group_colors) +
  theme_minimal() +
  labs(
    title = "CD3+ G(50)",
    x = "Time (months)",
    y = "Predicted survival probability"
  )

plot75_cd3 <- ggplot(curve75_cd3$curves,
                     aes(x = time, y = Surv, color = Group)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = five_group_colors) +
  theme_minimal() +
  labs(
    title = "CD3+ G(75)",
    x = "Time (months)",
    y = "Predicted survival probability"
  )

plot25_cd3cd8 <- ggplot(curve25_cd3cd8$curves,
                        aes(x = time,y = Surv,color = Group)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = five_group_colors) +
  theme_minimal() +
  labs(
    title = "CD3+CD8+ G(25)",
    x = "Time (months)",
    y = "Predicted survival probability"
  )

plot50_cd3cd8 <- ggplot(curve50_cd3cd8$curves,
                        aes(x = time, y = Surv, color = Group)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = five_group_colors) +
  theme_minimal() +
  labs(
    title = "CD3+CD8+ G(50)",
    x = "Time (months)",
    y = "Predicted survival probability"
  )

plot75_cd3cd8 <- ggplot(curve75_cd3cd8$curves,
                        aes(x = time, y = Surv, color = Group)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = five_group_colors) +
  theme_minimal() +
  labs(
    title = "CD3+CD8+ G(75)",
    x = "Time (months)",
    y = "Predicted survival probability"
  )

res25_cd3 <- list(
  plot = plot25_cd3,
  hr_table = curve25_cd3$hr_table,
  curves = curve25_cd3$curves,
  pooled_effects = pooled_25_cd3,
  cohort_effects = cohort_effects_25_cd3
)

res50_cd3 <- list(
  plot = plot50_cd3,
  hr_table = curve50_cd3$hr_table,
  curves = curve50_cd3$curves,
  pooled_effects = pooled_50_cd3,
  cohort_effects = cohort_effects_50_cd3
)

res75_cd3 <- list(
  plot = plot75_cd3,
  hr_table = curve75_cd3$hr_table,
  curves = curve75_cd3$curves,
  pooled_effects = pooled_75_cd3,
  cohort_effects = cohort_effects_75_cd3
)

res25_cd3cd8 <- list(
  plot = plot25_cd3cd8,
  hr_table = curve25_cd3cd8$hr_table,
  curves = curve25_cd3cd8$curves,
  pooled_effects = pooled_25_cd3cd8,
  cohort_effects = cohort_effects_25_cd3cd8
)

res50_cd3cd8 <- list(
  plot = plot50_cd3cd8,
  hr_table = curve50_cd3cd8$hr_table,
  curves = curve50_cd3cd8$curves,
  pooled_effects = pooled_50_cd3cd8,
  cohort_effects = cohort_effects_50_cd3cd8
)

res75_cd3cd8 <- list(
  plot = plot75_cd3cd8,
  hr_table = curve75_cd3cd8$hr_table,
  curves = curve75_cd3cd8$curves,
  pooled_effects = pooled_75_cd3cd8,
  cohort_effects = cohort_effects_75_cd3cd8
)

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

res25_cd3$plot <- res25_cd3$plot +
  labs(title = "CD3+ G(25)")

res50_cd3$plot <- res50_cd3$plot +
  labs(title = "CD3+ G(50)")

res75_cd3$plot <- res75_cd3$plot +
  labs(title = "CD3+ G(75)")

res25_cd3cd8$plot <- res25_cd3cd8$plot +
  labs(title = "CD3+CD8+ G(25)")

res50_cd3cd8$plot <- res50_cd3cd8$plot +
  labs(title = "CD3+CD8+ G(50)")

res75_cd3cd8$plot <- res75_cd3cd8$plot +
  labs(title = "CD3+CD8+ G(75)")

res5_cd3$plot
res5_cd3cd8$plot

res5_cd3$plot_ci
res5_cd3cd8$plot_ci

res5_cd3$hr_table
res5_cd3cd8$hr_table

res5_cd3$cohort_effects
res5_cd3cd8$cohort_effects

res5_cd3$pooled_effects
res5_cd3cd8$pooled_effects

res5_cd3$group_counts
res5_cd3cd8$group_counts

res5_cd3$group_tests
res5_cd3cd8$group_tests

res5_cd3$ph_tests
res5_cd3cd8$ph_tests

res25_cd3$plot
res25_cd3cd8$plot

res50_cd3$plot
res50_cd3cd8$plot

res75_cd3$plot
res75_cd3cd8$plot


# ggsave("model4_cd3_predict_survival_curves.tif", res_cd3$plot, width = 8, height = 5, dpi = 300)
# ggsave("model4_cd3cd8_predict_survival_curves.tif", res_cd3cd8$plot, width = 8, height = 5, dpi = 300)

# ggsave("model4_cd3_predict_survival_curves_wCI.tif", res_cd3$plot2, width = 8, height = 5, dpi = 300)
# ggsave("model4_cd3cd8_predict_survival_curves_wCI.tif", res_cd3cd8$plot2, width = 8, height = 5, dpi = 300)

# ggplot2::ggsave(
#   "five_group_cd3_standardized_survival_curves.tif",
#   res5_cd3$plot,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "five_group_cd3cd8_standardized_survival_curves.tif",
#   res5_cd3cd8$plot,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "five_group_cd3_standardized_survival_curves_G25.tif",
#   res25_cd3$plot,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "five_group_cd3cd8_standardized_survival_curves_G25.tif",
#   res25_cd3cd8$plot,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "five_group_cd3_standardized_survival_curves_G50.tif",
#   res50_cd3$plot,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "five_group_cd3cd8_standardized_survival_curves_G50.tif",
#   res50_cd3cd8$plot,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "five_group_cd3_standardized_survival_curves_G75.tif",
#   res75_cd3$plot,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "five_group_cd3cd8_standardized_survival_curves_G75.tif",
#   res75_cd3cd8$plot,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "fpca_plot3_cd3_Gest_wVariances.tif",
#   fpca_pl_dat$`CD3+ Gest`$plot,
#   width = 6,
#   height = 8,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "fpca_plot3_cd3cd8_Gest_wVariances.tif",
#   fpca_pl_dat$`CD3+ CD8+ Gest`$plot,
#   width = 6,
#   height = 8,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "fpca_plot3_cd3_Kest_wVariances.tif",
#   fpca_pl_dat$`CD3+ Kest`$plot,
#   width = 6,
#   height = 8,
#   dpi = 300,
#   compression = "lzw"
# )

# ggplot2::ggsave(
#   "fpca_plot3_cd3cd8_Kest_wVariances.tif",
#   fpca_pl_dat$`CD3+ CD8+ Kest`$plot,
#   width = 6,
#   height = 8,
#   dpi = 300,
#   compression = "lzw"
# )


fpc1_models_cd3 <- lapply(step_wise_models[["CD3+ Gest"]], function(x) 
  refit_fpc1_only(x$models$`FPC Score Interactions`))

fpc1_models_cd3cd8 <- lapply(step_wise_models[["CD3+ CD8+ Gest"]], function(x) 
  refit_fpc1_only(x$models$`FPC Score Interactions`))

# save_ph_plots(fpc1_models_cd3, "PH_diagnostics_CD3_Gest.pdf")
# save_ph_plots(fpc1_models_cd3cd8, "PH_diagnostics_CD3CD8_Gest.pdf")

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

