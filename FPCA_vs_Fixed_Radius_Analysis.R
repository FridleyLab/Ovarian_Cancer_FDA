library(mxfda)
library(ggplot2)
library(openxlsx)
library(survival)
library(dplyr)
library(metafor)
library(survminer)
library(tidyverse)
library(readxl)
library(knitr)
library(meta)
library(grid)

fixed_radii <- c(25, 50, 75)

analysis_dat_list1 <- readRDS("data_part1.rds")
analysis_dat_list2 <- readRDS("data_part2.rds")

metrics <- names(analysis_dat_list1)

analysis_dat_list_combined <- lapply(metrics, function(m) {
  
  dat <- dplyr::bind_rows(
    analysis_dat_list1[[m]],
    analysis_dat_list2[[m]] |>
      dplyr::mutate(
        dplyr::across(
          c(block, histotype),
          as.character
        ),
        stage_cat = ifelse(
          stage_cat == 1,
          "high stage",
          "low stage"
        )
      )
  ) |>
    dplyr::mutate(
      block = as.character(block),
      stage_cat = factor(
        stage_cat,
        levels = c("high stage", "low stage")
      ),
      abund_binary = as.numeric(abund_binary)
    )
  
  split(dat, dat$cohort)
})

names(analysis_dat_list_combined) <- metrics

build_fixed_radius_data <- function(
    cohort_dat,
    radii = c(25, 50, 75)
) {
  
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
  
  missing_variables <- setdiff(required_variables, names(cohort_dat))
  
  if (length(missing_variables) > 0) {
    stop("Missing required variables: ", paste(missing_variables, collapse = ", "))
  }
  
  dat <- cohort_dat |>
    dplyr::filter(func_observed == 1, r %in% radii) |>
    dplyr::select(dplyr::all_of(required_variables)) |>
    dplyr::distinct()
  
  duplicate_check <- dat |>
    dplyr::count(block, r) |>
    dplyr::filter(n > 1)
  
  if (nrow(duplicate_check) > 0) {
    stop("More than one subject-level trajectory value was found for at least one block-radius combination.")
  }
  
  complete_radius_ids <- dat |>
    dplyr::group_by(block) |>
    dplyr::summarise(n_radii = dplyr::n_distinct(r), .groups = "drop") |>
    dplyr::filter(n_radii == length(radii)) |>
    dplyr::pull(block)
  
  dat <- dat |>
    dplyr::filter(block %in% complete_radius_ids) |>
    dplyr::group_by(r) |>
    dplyr::mutate(
      fixed_value_z = as.numeric(scale(func_weighted)),
      fixed_cluster_cat = factor(
        ifelse(
          func_weighted >= median(
            func_weighted,
            na.rm = TRUE
          ),
          "high",
          "low"
        ),
        levels = c("high","low")
      )
    ) |>
    dplyr::ungroup()
  
  dat
}

fit_fixed_radius_models <- function(
    cohort_dat,
    radii = c(25, 50, 75)
) {
  
  fixed_dat <- build_fixed_radius_data(cohort_dat = cohort_dat, radii = radii)
  
  split(fixed_dat, fixed_dat$r) |>
    lapply(function(dat_r) {
      
      if (nlevels(droplevels(dat_r$fixed_cluster_cat)) < 2) {
        stop(
          "Only one fixed-radius clustering category was present ",
          "at radius ", unique(dat_r$r), " in cohort ", unique(dat_r$cohort), "."
        )
      }
      
      coxph(
        survival::Surv(mossurv, death) ~ ageatdx + stage_cat + abund_binary * fixed_cluster_cat,
        data = dat_r,
        ties = "efron",
        x = TRUE,
        model = TRUE
      )
    })
}

fit_fixed_radius_models_cont <- function(
    cohort_dat,
    radii = c(25, 50, 75)
) {
  
  fixed_dat <- build_fixed_radius_data(cohort_dat = cohort_dat, radii = radii)
  
  split(fixed_dat, fixed_dat$r) |>
    lapply(function(dat_r) {
      
      if (
        all(is.na(dat_r$fixed_value_z)) ||
        stats::sd(dat_r$fixed_value_z, na.rm = TRUE) == 0
      ) {
        stop(
          "The standardized fixed-radius value has no variation ",
          "at radius ", unique(dat_r$r), " in cohort ", unique(dat_r$cohort), "."
        )
      }
      
      coxph(survival::Surv(mossurv,death) ~ ageatdx + stage_cat + abund_binary * fixed_value_z,
        data = dat_r,
        ties = "efron",
        x = TRUE,
        model = TRUE
      )
    })
}

fixed_radius_models <- list(
  `CD3+` = lapply(
    analysis_dat_list_combined[["CD3+ G"]],
    fit_fixed_radius_models,
    radii = fixed_radii
  ),
  `CD3+CD8+` = lapply(
    analysis_dat_list_combined[["CD3+ CD8+ G"]],
    fit_fixed_radius_models,
    radii = fixed_radii
  )
)

fixed_radius_models_cont <- list(
  `CD3+` = lapply(
    analysis_dat_list_combined[["CD3+ G"]],
    fit_fixed_radius_models_cont,
    radii = fixed_radii
  ),
  
  `CD3+CD8+` = lapply(
    analysis_dat_list_combined[["CD3+ CD8+ G"]],
    fit_fixed_radius_models_cont,
    radii = fixed_radii
  )
)

lapply(
  fixed_radius_models_cont,
  function(phenotype_models) {
    sapply(phenotype_models, length)
  }
)

summary(
  fixed_radius_models_cont[["CD3+"]][["AAS"]][["25"]]
)

extract_model_coefficients <- function(
    mod,
    phenotype,
    cohort,
    radius
) {
  
  sm <- summary(mod)
  
  coef_dat <- as.data.frame(sm$coefficients, check.names = FALSE) |>
    tibble::rownames_to_column("factor")
  
  ci_dat <- as.data.frame(sm$conf.int, check.names = FALSE) |>
    tibble::rownames_to_column("factor") |>
    dplyr::select(factor, `exp(coef)`, `lower .95`, `upper .95`)
  
  coef_dat |>
    dplyr::left_join(ci_dat, by = "factor") |>
    dplyr::mutate(
      phenotype = phenotype,
      cohort = cohort,
      radius = as.numeric(radius),
      .before = 1
    )
}

fixed_radius_coefficients <- purrr::imap_dfr(
  fixed_radius_models,
  function(phenotype_models, phenotype_name) {
    
    purrr::imap_dfr(
      phenotype_models,
      function(cohort_models, cohort_name) {
        
        purrr::imap_dfr(
          cohort_models,
          function(mod, radius_name) {
            
            extract_model_coefficients(
              mod = mod,
              phenotype = phenotype_name,
              cohort = cohort_name,
              radius = radius_name
            )
          }
        )
      }
    )
  }
)

fixed_radius_coefficients_cont <- purrr::imap_dfr(
  fixed_radius_models_cont,
  function(
    phenotype_models,
    phenotype_name
  ) {
    
    purrr::imap_dfr(
      phenotype_models,
      function(
    cohort_models,
    cohort_name
      ) {
        
        purrr::imap_dfr(
          cohort_models,
          function(
    mod,
    radius_name
          ) {
            
            extract_model_coefficients(
              mod = mod,
              phenotype = phenotype_name,
              cohort = cohort_name,
              radius = radius_name
            )
          }
        )
      }
    )
  }
)

extract_fixed_radius_contrasts <- function(
    mod,
    phenotype,
    cohort,
    radius
) {
  
  b <- stats::coef(mod)
  V <- stats::vcov(mod)
  
  abundance_term <- "abund_binary"
  
  clustering_term <- grep("^fixed_cluster_cat", names(b), value = TRUE)
  
  interaction_term <- grep(
    "abund_binary:fixed_cluster_cat|fixed_cluster_cat.*:abund_binary",
    names(b),
    value = TRUE
  )
  
  if (length(clustering_term) != 1) {
    stop("Expected one fixed clustering coefficient. Found: ", paste(clustering_term, collapse = ", "))
  }
  
  if (length(interaction_term) != 1) {
    stop("Expected one abundance-by-clustering interaction coefficient. Found: ",
         paste(interaction_term, collapse = ", "))
  }
  
  make_contrast <- function(terms_to_add) {
    
    L <- stats::setNames(rep(0, length(b)), names(b))
    L[terms_to_add] <- 1
    used_terms <- names(L)[L != 0]
    
    if (
      anyNA(b[used_terms]) ||
      any(!is.finite(b[used_terms]))
    ) {
      return(c(logHR = NA_real_, SE = NA_real_)
      )
    }
    
    logHR <- sum(L[used_terms] * b[used_terms])
    V_used <- V[used_terms, used_terms, drop = FALSE]
    L_used <- L[used_terms]
    SE <- sqrt(as.numeric(t(L_used) %*% V_used %*% L_used))
    
    c(logHR = logHR, SE = SE)
  }
  
  contrast_results <- rbind(
    `Low abundance + low clustering` = make_contrast(clustering_term),
    `High abundance + high clustering` = make_contrast(abundance_term),
    `High abundance + low clustering` = make_contrast(c(abundance_term, clustering_term, interaction_term))
  ) |>
    as.data.frame()
  
  contrast_results |>
    tibble::rownames_to_column("Group") |>
    dplyr::mutate(
      phenotype = phenotype,
      cohort = cohort,
      radius = as.numeric(radius),
      HR = exp(logHR),
      LCL = exp(logHR - 1.96 * SE),
      UCL = exp(logHR + 1.96 * SE),
      .before = 1
    )
}

fixed_radius_contrasts <- purrr::imap_dfr(
  fixed_radius_models,
  function(phenotype_models, phenotype_name) {
    
    purrr::imap_dfr(
      phenotype_models,
      function(cohort_models, cohort_name) {
        
        purrr::imap_dfr(
          cohort_models,
          function(mod, radius_name) {
            
            extract_fixed_radius_contrasts(
              mod = mod,
              phenotype = phenotype_name,
              cohort = cohort_name,
              radius = radius_name
            )
          }
        )
      }
    )
  }
)

extract_fixed_radius_model_metrics <- function(
    mod,
    phenotype,
    cohort,
    radius
) {
  
  mf <- model.frame(mod)
  y <- model.response(mf)
  
  analysis_dat <- data.frame(
    mossurv = y[, 1],
    death = y[, 2],
    ageatdx = mf$ageatdx,
    stage_cat = mf$stage_cat,
    abund_binary = mf$abund_binary,
    fixed_cluster_cat = mf$fixed_cluster_cat
  )
  
  # Make categorical reference levels explicit
  analysis_dat$stage_cat <- factor(
    analysis_dat$stage_cat,
    levels = c("high stage", "low stage")
  )
  
  analysis_dat$fixed_cluster_cat <- factor(
    analysis_dat$fixed_cluster_cat,
    levels = c("high", "low")
  )
  
  # Reduced model:
  # age + stage + abundance only
  reduced_mod <- coxph(survival::Surv(mossurv, death) ~ ageatdx + stage_cat + abund_binary,
    data = analysis_dat,
    ties = mod$method,
    x = TRUE,
    model = TRUE
  )
  
  # Full fixed-radius model:
  # age + stage + abundance + clustering + abundance x clustering
  full_mod <- coxph(survival::Surv(mossurv, death) ~ ageatdx + stage_cat + abund_binary * fixed_cluster_cat,
    data = analysis_dat,
    ties = mod$method,
    x = TRUE,
    model = TRUE
  )
  
  # Joint likelihood-ratio test of the clustering main effect
  # and abundance-by-clustering interaction
  lrt <- anova(reduced_mod, full_mod, test = "LRT")
  
  concordance_result <- summary(full_mod)$concordance
  
  data.frame(
    phenotype = phenotype,
    cohort = cohort,
    radius = as.numeric(radius),
    n = full_mod$n,
    events = full_mod$nevent,
    logLik = as.numeric(logLik(full_mod)),
    AIC = AIC(full_mod),
    BIC = BIC(full_mod),
    C_index = unname(concordance_result[1]),
    C_index_SE = unname(concordance_result[2]),
    spatial_terms_LRT_df = as.numeric(lrt[2, "Df"]),
    spatial_terms_LRT_chisq = as.numeric(lrt[2, "Chisq"]),
    spatial_terms_LRT_p = as.numeric(lrt[2, "Pr(>|Chi|)"])
  )
}

extract_fixed_radius_model_metrics_cont <- function(
    mod,
    phenotype,
    cohort,
    radius
) {
  
  mf <- model.frame(mod)
  y <- model.response(mf)
  
  analysis_dat <- data.frame(
    mossurv = y[, 1],
    death = y[, 2],
    ageatdx = mf$ageatdx,
    stage_cat = mf$stage_cat,
    abund_binary = mf$abund_binary,
    fixed_value_z = mf$fixed_value_z
  )
  
  analysis_dat$stage_cat <- factor(
    analysis_dat$stage_cat,
    levels = c("high stage", "low stage")
  )
  
  # Reduced model:
  # age + stage + abundance
  reduced_mod <- coxph(survival::Surv(mossurv,death) ~ ageatdx + stage_cat + abund_binary,
    data = analysis_dat,
    ties = mod$method,
    x = TRUE,
    model = TRUE
  )
  
  # Full continuous fixed-radius Model 3:
  # age + stage + abundance + G(r) + abundance x G(r)
  full_mod <- coxph(survival::Surv(mossurv,death) ~ ageatdx + stage_cat + abund_binary * fixed_value_z,
    data = analysis_dat,
    ties = mod$method,
    x = TRUE,
    model = TRUE
  )
  
  # Joint 2-df test of:
  # fixed_value_z
  # abundance x fixed_value_z
  lrt <- anova(reduced_mod, full_mod, test = "LRT")
  concordance_result <- summary(full_mod)$concordance
  coefficient_table <- summary(full_mod)$coefficients
  spatial_main_term <- "fixed_value_z"
  spatial_interaction_term <- "abund_binary:fixed_value_z"
  data.frame(
    phenotype = phenotype,
    cohort = cohort,
    radius = as.numeric(radius),
    n = full_mod$n,
    events = full_mod$nevent,
    logLik = as.numeric(logLik(full_mod)),
    AIC = stats::AIC(full_mod),
    BIC = stats::BIC(full_mod),
    C_index = unname(concordance_result[1]),
    C_index_SE = unname(concordance_result[2]),
    spatial_main_logHR = coefficient_table[spatial_main_term, "coef"],
    spatial_main_SE = coefficient_table[spatial_main_term, "se(coef)"],
    spatial_main_p = coefficient_table[spatial_main_term, "Pr(>|z|)"],
    interaction_logHR = coefficient_table[spatial_interaction_term, "coef"],
    interaction_SE = coefficient_table[spatial_interaction_term, "se(coef)"],
    interaction_p = coefficient_table[spatial_interaction_term, "Pr(>|z|)"],
    spatial_terms_LRT_df = as.numeric(lrt[2, "Df"]),
    spatial_terms_LRT_chisq = as.numeric(lrt[2, "Chisq"]),
    spatial_terms_LRT_p = as.numeric(lrt[2, "Pr(>|Chi|)"])
  )
}

fixed_radius_model_metrics <- purrr::imap_dfr(
  fixed_radius_models,
  function(phenotype_models, phenotype_name) {
    
    purrr::imap_dfr(
      phenotype_models,
      function(cohort_models, cohort_name) {
        
        purrr::imap_dfr(
          cohort_models,
          function(mod, radius_name) {
            
            extract_fixed_radius_model_metrics(
              mod = mod,
              phenotype = phenotype_name,
              cohort = cohort_name,
              radius = radius_name
            )
          }
        )
      }
    )
  }
)

fixed_radius_model_metrics_cont <- purrr::imap_dfr(
  fixed_radius_models_cont,
  function(
    phenotype_models,
    phenotype_name
  ) {
    
    purrr::imap_dfr(
      phenotype_models,
      function(
    cohort_models,
    cohort_name
      ) {
        
        purrr::imap_dfr(
          cohort_models,
          function(
    mod,
    radius_name
          ) {
            
            extract_fixed_radius_model_metrics_cont(
              mod = mod,
              phenotype = phenotype_name,
              cohort = cohort_name,
              radius = radius_name
            )
          }
        )
      }
    )
  }
)

fixed_radius_group_counts <- purrr::imap_dfr(
  fixed_radius_models,
  function(phenotype_models, phenotype_name) {
    
    purrr::imap_dfr(
      phenotype_models,
      function(cohort_models, cohort_name) {
        
        purrr::imap_dfr(
          cohort_models,
          function(mod, radius_name) {
            
            mf <- model.frame(mod)
            y <- model.response(mf)
            
            mf |>
              dplyr::mutate(
                death = y[, 2],
                Group = dplyr::case_when(
                  abund_binary == 0 & fixed_cluster_cat == "high" ~ "Low abundance + high clustering",
                  abund_binary == 0 & fixed_cluster_cat == "low" ~ "Low abundance + low clustering",
                  abund_binary == 1 & fixed_cluster_cat == "high" ~ "High abundance + high clustering",
                  abund_binary == 1 & fixed_cluster_cat == "low" ~ "High abundance + low clustering",
                  TRUE ~ NA_character_
                )
              ) |>
              dplyr::count(Group, death, name = "n") |>
              dplyr::mutate(
                phenotype = phenotype_name,
                cohort = cohort_name,
                radius = as.numeric(radius_name),
                .before = 1
              )
          }
        )
      }
    )
  }
)

fixed_radius_ph_tests_cont <- purrr::imap_dfr(
  fixed_radius_models_cont,
  function(
    phenotype_models,
    phenotype_name
  ) {
    
    purrr::imap_dfr(
      phenotype_models,
      function(
    cohort_models,
    cohort_name
      ) {
        
        purrr::imap_dfr(
          cohort_models,
          function(
    mod,
    radius_name
          ) {
            
            zph <- survival::cox.zph(mod, transform = "km")
            
            as.data.frame(zph$table) |>
              tibble::rownames_to_column("term") |>
              dplyr::mutate(
                phenotype = phenotype_name,
                cohort = cohort_name,
                radius = as.numeric(
                  radius_name
                ),
                .before = 1
              )
          }
        )
      }
    )
  }
)

fixed_radius_convergence_check_cont <-
  fixed_radius_coefficients_cont |>
  dplyr::mutate(
    unstable = dplyr::case_when(
      !is.finite(coef) ~ TRUE,
      !is.finite(`se(coef)`) ~ TRUE,
      abs(coef) > 10 ~ TRUE,
      `se(coef)` > 10 ~ TRUE,
      TRUE ~ FALSE
    )
  )


fixed_radius_convergence_check_cont |>
  dplyr::filter(unstable) |>
  dplyr::select(
    phenotype,
    cohort,
    radius,
    factor,
    coef,
    `se(coef)`
  )

extract_continuous_clustering_slopes <- function(
    mod,
    phenotype,
    cohort,
    radius
) {
  b <- stats::coef(mod)
  V <- stats::vcov(mod)
  main_term <- "fixed_value_z"
  interaction_term <- "abund_binary:fixed_value_z"
  required_terms <- c(main_term, interaction_term)
  
  if (any(!required_terms %in% names(b))) {
    stop("Expected continuous spatial terms were not found.")
  }
  
  make_estimate <- function(
    terms_to_add
  ) {
    
    if (
      anyNA(b[terms_to_add]) ||
      any(!is.finite(b[terms_to_add]))
    ) {
      return(
        c(logHR = NA_real_, SE = NA_real_)
      )
    }
    
    L <- stats::setNames(rep(0, length(b)), names(b))
    L[terms_to_add] <- 1
    L_used <- L[terms_to_add]
    V_used <- V[terms_to_add, terms_to_add, drop = FALSE]
    logHR <- sum(L_used * b[terms_to_add])
    SE <- sqrt(as.numeric(t(L_used) %*% V_used %*% L_used))
    c(logHR = logHR, SE = SE)
  }
  
  results <- rbind(
    `Low abundance: clustering slope` = make_estimate(main_term),
    `High abundance: clustering slope` = make_estimate(c(main_term,interaction_term))
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("Contrast")
  
  results |>
    dplyr::mutate(
      phenotype = phenotype,
      cohort = cohort,
      radius = as.numeric(radius),
      HR = exp(logHR),
      LCL = exp(logHR - 1.96 * SE),
      UCL = exp(logHR + 1.96 * SE),
      .before = 1
    )
}

fixed_radius_slopes_cont <- purrr::imap_dfr(
  fixed_radius_models_cont,
  function(
    phenotype_models,
    phenotype_name
  ) {
    
    purrr::imap_dfr(
      phenotype_models,
      function(
    cohort_models,
    cohort_name
      ) {
        
        purrr::imap_dfr(
          cohort_models,
          function(
    mod,
    radius_name
          ) {
            
            extract_continuous_clustering_slopes(
              mod = mod,
              phenotype = phenotype_name,
              cohort = cohort_name,
              radius = radius_name
            )
          }
        )
      }
    )
  }
)

dim(fixed_radius_slopes_cont)

fixed_radius_slopes_cont |>
  dplyr::arrange(
    phenotype,
    radius,
    cohort,
    Contrast
  )

fixed_radius_output_file <- paste0("fixed_radius_25_50_75_cox_results.xlsx")

openxlsx::write.xlsx(
  list(
    `Model Coefficients` = fixed_radius_coefficients,
    `Group Contrasts` = fixed_radius_contrasts,
    `Model Metrics` = fixed_radius_model_metrics,
    `Group Counts` = fixed_radius_group_counts
  ),
  file = fixed_radius_output_file,
  overwrite = TRUE
)

fixed_radius_cont_output_file <- paste0("fixed_radius_25_50_75_continuous_cox_results.xlsx")

continuous_output_objects <- list(
  `Model Coefficients` = fixed_radius_coefficients_cont,
  `Model Metrics` = fixed_radius_model_metrics_cont,
  `PH Tests` = fixed_radius_ph_tests_cont,
  `Convergence Check` = fixed_radius_convergence_check_cont
)

continuous_output_row_counts <- vapply(continuous_output_objects, nrow, integer(1))

print(continuous_output_row_counts)

if (any(continuous_output_row_counts == 0)) {
  stop(
    "At least one continuous-model output is empty: ",
    paste(
      names(continuous_output_row_counts)[continuous_output_row_counts == 0],
      collapse = ", "
    )
  )
}

# openxlsx::write.xlsx(
#   continuous_output_objects,
#   file = fixed_radius_cont_output_file,
#   overwrite = TRUE
# )
