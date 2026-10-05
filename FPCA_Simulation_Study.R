library(scSpatialSIM)
library(spatialTIME)
library(mxfda)
library(survival)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

set.seed(9)

simulation_settings <- list(
  n_replicates = 100,
  n_subjects = 700,
  mean_positive_proportion = 0.03,
  abundance_cutoff = 0.01,
  radii = 5:100,
  fixed_radii = c(25, 75),
  min_positive_cells = 6,
  lambda = 50,
  internal_window = spatstat.geom::owin(
    xrange = c(0, 10),
    yrange = c(0, 10)
  ),
  coordinate_multiplier = 200,
  probability_multiplier = c(
    "1" = 10,
    "6" = 2
  ),
  censoring_rate = 0.003,
  baseline_rate = 0.01
)

make_subject_design <- function(n_subjects,seed) {
  
  set.seed(seed)
  
  if (n_subjects %% 4 != 0) {
    stop("n_subjects must be divisible by 4 for a balanced abundance-by-spatial design.")
  }
  
  design <- tidyr::expand_grid(
    abund_binary = c(0L, 1L),
    spatial_binary = c(0L, 1L),
    replicate_within_group = seq_len(n_subjects / 4)
  ) |>
    dplyr::slice_sample(prop = 1) |>
    dplyr::mutate(subject_id = sprintf("S%03d", dplyr::row_number()),
                  target_proportion = dplyr::if_else(
                    abund_binary == 0L,
                    0.005,
                    0.055
                  ),
                  k = dplyr::if_else(
                    spatial_binary == 0L,
                    6L,
                    1L
                  )
    ) |>
    dplyr::select(subject_id,abund_binary,target_proportion,spatial_binary,k)
  
  design
}

simulate_spatial_replicate <- function(replicate_id,settings) {
  
  design <- make_subject_design(
    n_subjects = settings$n_subjects,
    seed = 10000 + replicate_id
  )
  
  spatial_files <- vector("list",settings$n_subjects)
  
  for (i in seq_len(settings$n_subjects)) {
    
    set.seed(replicate_id * 100000 + i)
    
    current_multiplier <- settings$probability_multiplier[as.character(design$k[i])]
    positivity_parameter <- min(design$target_proportion[i] * current_multiplier, 0.95)
    
    sim_object <-
      suppressMessages(
        suppressWarnings(scSpatialSIM::CreateSimulationObject(
          window = settings$internal_window,
          sims = 1,
          cell_types = 1
        ) |>
          scSpatialSIM::GenerateSpatialPattern(lambda = settings$lambda) |>
          scSpatialSIM::GenerateCellPositivity(
            probs = c(0, positivity_parameter),
            k = design$k[i],
            xmin = 0,
            xmax = 10,
            ymin = 0,
            ymax = 10,
            Force = TRUE,
            use_window = FALSE
          )
        )
      )
    
    
    one_spatial_file <- scSpatialSIM::CreateSpatialList(sim_object, single_df = FALSE)[[1]] |>
      dplyr::mutate(
        x = x * settings$coordinate_multiplier,
        y = y * settings$coordinate_multiplier,
        subject_id = design$subject_id[i],
        sample_id = design$subject_id[i]
      )
    
    spatial_files[[i]] <- one_spatial_file
    
    if(i %% 100 == 0){
      cat("\n", format(Sys.time(), "%H:%M:%S"), "- Replicate", replicate_id,
          "- completed", i, "/", settings$n_subjects, "subjects\n")
    }
  }
  
  names(spatial_files) <- design$subject_id
  
  list(design = design, spatial_files = spatial_files)
}

summarize_simulated_abundance <- function(
    spatial_files,
    design,
    marker_column,
    abundance_cutoff = 0.01,
    min_positive_cells = 6
) {
  
  observed <- purrr::imap_dfr(spatial_files, function(dat, subject_name) {
    
    positive_count <- sum(dat[[marker_column]] == 1, na.rm = TRUE)
    total_count <- nrow(dat)
    tibble::tibble(
      subject_id = subject_name,
      total_cells = total_count,
      positive_cells = positive_count,
      observed_proportion = positive_count / total_count,
      abund_binary_observed = as.integer(observed_proportion >= abundance_cutoff),
      func_observed = as.integer(positive_count >= min_positive_cells)
    )
  }
  )
  
  design |> dplyr::left_join(observed, by = "subject_id")
}

simulate_survival <- function(
    subject_data,
    beta_abundance,
    beta_spatial,
    beta_interaction,
    baseline_rate,
    censoring_rate,
    seed = NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  linear_predictor <-
    beta_abundance * subject_data$abund_binary_observed +
    beta_spatial * subject_data$spatial_binary +
    beta_interaction * subject_data$abund_binary_observed * subject_data$spatial_binary
  
  event_time <- rexp(nrow(subject_data), rate = baseline_rate * exp(linear_predictor))
  
  censor_time <- rexp(nrow(subject_data), rate = censoring_rate)
  
  subject_data |>
    dplyr::mutate(
      mossurv = pmin(event_time, censor_time),
      death = as.integer(event_time <= censor_time),
      true_linear_predictor = linear_predictor
    )
}

make_spatialtime_object <- function(spatial_files, subject_data) {
  
  clinical_data <- subject_data |>
    dplyr::transmute(
      subject_id,
      mossurv,
      death,
      abund_binary,
      abund_binary_observed,
      spatial_binary,
      func_observed
    )
  
  sample_data <- subject_data |>
    dplyr::transmute(
      subject_id,
      sample_id = subject_id,
      total_cells,
      positive_cells,
      observed_proportion
    )
  
  spatialTIME::create_mif(
    clinical_data = clinical_data,
    sample_data = sample_data,
    spatial_list = spatial_files,
    patient_id = "subject_id",
    sample_id = "sample_id"
  )
}

extract_interaction_result <- function(
    mod,
    method_name,
    replicate_id,
    interaction_pattern
) {
  
  coef_table <- summary(mod)$coefficients
  interaction_row <- grep(interaction_pattern, rownames(coef_table), value = TRUE)
  
  if (length(interaction_row) != 1) {
    return(
      tibble::tibble(
        replicate = replicate_id,
        method = method_name,
        estimate = NA_real_,
        SE = NA_real_,
        HR = NA_real_,
        LCL = NA_real_,
        UCL = NA_real_,
        significant = NA,
        estimable = FALSE
      )
    )
  }
  
  estimate <- coef_table[interaction_row, "coef"]
  standard_error <- coef_table[interaction_row, "se(coef)"]
  
  estimable <-
    is.finite(estimate) &&
    is.finite(standard_error) &&
    standard_error > 0 &&
    standard_error <= 10
  
  tibble::tibble(
    replicate = replicate_id,
    method = method_name,
    estimate = estimate,
    SE = standard_error,
    HR = exp(estimate),
    LCL = exp(estimate - 1.96 * standard_error),
    UCL = exp(estimate + 1.96 * standard_error),
    significant = estimable && (LCL > 1 || UCL < 1),
    estimable = estimable
  )
}

run_single_replicate <- function(
    replicate_id,
    settings,
    beta_abundance,
    beta_spatial,
    beta_interaction,
    num_permutations = 50
) {
  
  message("Starting replicate ", replicate_id)
  
  # Simulate spatial point patterns
  
  sim_data <- simulate_spatial_replicate(
    replicate_id = replicate_id,
    settings = settings
  )
  
  # Calculate realized abundance and spatial availability
  
  subject_data <- summarize_simulated_abundance(
    spatial_files = sim_data$spatial_files,
    design = sim_data$design,
    marker_column = "Cell 1 Assignment",
    abundance_cutoff = settings$abundance_cutoff,
    min_positive_cells = settings$min_positive_cells
  )
  
  # Generate survival outcomes
  
  subject_data <- simulate_survival(
    subject_data = subject_data,
    beta_abundance = beta_abundance,
    beta_spatial = beta_spatial,
    beta_interaction = beta_interaction,
    baseline_rate = settings$baseline_rate,
    censoring_rate = settings$censoring_rate,
    seed = 20000 + replicate_id
  )
  
  abundance_diagnostic <- subject_data |>
    dplyr::summarise(
      latent_observed_agreement = mean(abund_binary == abund_binary_observed, na.rm = TRUE),
      observed_high_rate = mean(abund_binary_observed == 1, na.rm = TRUE)
    )
  
  print(abundance_diagnostic)
  
  print(
    table(
      latent = subject_data$abund_binary,
      observed = subject_data$abund_binary_observed
    )
  )
  
  print(
    table(
      observed_abundance = subject_data$abund_binary_observed,
      spatial_truth = subject_data$spatial_binary
    )
  )
  
  cat("\nLatent abundance by observed abundance and spatial truth:\n")
  
  print(
    with(
      subject_data,
      table(
        latent_abundance = abund_binary,
        observed_abundance = abund_binary_observed,
        spatial_truth = spatial_binary
      )
    )
  )
  
  cat("\nObserved abundance summaries by latent design and spatial truth:\n")
  
  print(
    subject_data |>
      dplyr::group_by(abund_binary, spatial_binary) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_positive_cells = mean(positive_cells, na.rm = TRUE),
        median_positive_cells = median(positive_cells, na.rm = TRUE),
        mean_observed_proportion = mean(observed_proportion, na.rm = TRUE),
        median_observed_proportion = median(observed_proportion, na.rm = TRUE),
        observed_high_rate = mean(abund_binary_observed == 1, na.rm = TRUE),
        eligible_rate = mean(func_observed == 1, na.rm = TRUE),
        .groups = "drop"
      )
  )
  
  # Create spatialTIME object
  
  mif_object <- make_spatialtime_object(
    spatial_files = sim_data$spatial_files,
    subject_data = subject_data
  )
  
  # Calculate nearest-neighbor G trajectories
  
  mif_object <- spatialTIME::NN_G(
    mif = mif_object,
    mnames = "Cell 1 Assignment",
    r_range = settings$radii,
    num_permutations = num_permutations,
    edge_correction = "rs",
    workers = 1,
    overwrite = TRUE,
    xloc = "x",
    yloc = "y"
  )
  
  # Extract the centered G trajectories
  
  g_traj <- mif_object$derived$univariate_NN |>
    dplyr::ungroup() |>
    dplyr::filter(Marker == "Cell 1 Assignment") |>
    dplyr::transmute(sample_id, r, clustering = `Degree of Clustering Permutation`) |>
    dplyr::left_join(
      subject_data |>
        dplyr::select(subject_id, func_observed),
      by = c("sample_id" = "subject_id")
    )
  
  # Prepare eligible trajectories for FPCA
  
  eligible_subject_data <- subject_data |>
    dplyr::filter(func_observed == 1)
  
  mxfda_metadata <- eligible_subject_data |>
    dplyr::transmute(
      block = subject_id,
      sample = subject_id,
      mossurv,
      death,
      abund_binary,
      abund_binary_observed,
      spatial_binary,
      k,
      positive_cells,
      observed_proportion
    )
  
  mxfda_function_data <- g_traj |>
    dplyr::filter(
      func_observed == 1,
      r %in% settings$radii,
      is.finite(clustering)
    ) |>
    dplyr::transmute(sample = sample_id, r, func_weighted = clustering)
  
  # Stop cleanly if too few subjects have eligible trajectories
  eligible_trajectory_count <- dplyr::n_distinct(mxfda_function_data$sample)
  
  if (eligible_trajectory_count < 8) {
    stop("Replicate ", replicate_id," has fewer than 8 eligible trajectories.")
  }
  
  # Create mxfda object and run FPCA
  
  sim_mxfda_object <- mxfda::make_mxfda(
    metadata = mxfda_metadata,
    spatial = NULL,
    subject_key = "block",
    sample_key = "sample"
  )
  
  sim_mxfda_object <- mxfda::add_summary_function(
    sim_mxfda_object,
    summary_function_data = mxfda_function_data,
    metric = "uni g"
  )
  
  sim_mxfda_object <- mxfda::run_fpca(
    sim_mxfda_object,
    metric = "uni g",
    r = "r",
    value = "func_weighted",
    pve = 0.99
  )
  
  # Extract and orient FPC1
  
  fpca_scores <-
    sim_mxfda_object@functional_pca$Gest$score_df |>
    dplyr::select(sample, fpc1)
  
  fpca_scores <- fpca_scores |>
    left_join(
      subject_data |>
        select(subject_id, spatial_binary),
      by = c("sample" = "subject_id")
    )
  
  fpc1_direction <- cor(
    fpca_scores$fpc1,
    fpca_scores$spatial_binary,
    use = "complete.obs"
  )
  
  if(is.finite(fpc1_direction) && fpc1_direction > 0){
    fpca_scores$fpc1 <- -fpca_scores$fpc1
  }
  
  fpca_scores <- fpca_scores |>
    dplyr::select(sample, fpc1)
  
  # Extract G(25) and G(75)
  
  fixed_radius_data <- g_traj |>
    dplyr::filter(r %in% settings$fixed_radii) |>
    dplyr::select(sample_id, r, clustering) |>
    tidyr::pivot_wider(
      names_from = r,
      values_from = clustering,
      names_prefix = "g"
    )
  
  expected_fixed_columns <- paste0("g", settings$fixed_radii)
  
  if (!all(expected_fixed_columns %in% names(fixed_radius_data))) {
    stop("Fixed-radius columns were not generated in replicate ", replicate_id, ".")
  }
  
  # Construct analysis dataset
  
  analysis_data <- subject_data |>
    dplyr::left_join(fpca_scores, by = c("subject_id" = "sample")) |>
    dplyr::left_join(fixed_radius_data, by = c("subject_id" = "sample_id"))
  
  # This implementation assumes the two fixed radii are 25 and 75
  if (!all(c("g25", "g75") %in% names(analysis_data))) {
    stop("Expected g25 and g75 were not found in replicate ", replicate_id, ".")
  }
  
  # Create comparable high/low spatial categories
  
  g25_cutoff <- stats::median(analysis_data$g25[analysis_data$func_observed == 1], na.rm = TRUE)
  g75_cutoff <- stats::median(analysis_data$g75[analysis_data$func_observed == 1], na.rm = TRUE)
  
  analysis_data <- analysis_data |>
    dplyr::mutate(
      fpc1_cat = factor(
        dplyr::case_when(
          !is.finite(fpc1) ~ NA_character_,
          fpc1 >= 0 ~ "high",
          TRUE ~ "low"
        ),
        levels = c("high", "low")
      ),
      
      g25_cat = factor(
        dplyr::case_when(
          !is.finite(g25) ~ NA_character_,
          g25 < g25_cutoff ~ "high",
          TRUE ~ "low"
        ),
        levels = c("high", "low")
      ),
      
      g75_cat = factor(
        dplyr::case_when(
          !is.finite(g75) ~ NA_character_,
          g75 < g75_cutoff ~ "high",
          TRUE ~ "low"
        ),
        levels = c("high", "low")
      )
    )
  
  # Restrict every method to the same complete-case sample
  
  analysis_data_model <- analysis_data |>
    dplyr::filter(
      func_observed == 1,
      !is.na(fpc1_cat),
      !is.na(g25_cat),
      !is.na(g75_cat),
      !is.na(mossurv),
      !is.na(death),
      !is.na(abund_binary_observed)
    )
  
  abundance_agreement <- mean(
    subject_data$abund_binary == subject_data$abund_binary_observed,
    na.rm = TRUE
  )
  
  observed_high_rate <- mean(subject_data$abund_binary_observed == 1, na.rm = TRUE)
  
  fpc1_median_cut <- median(analysis_data_model$fpc1, na.rm = TRUE)
  
  analysis_data_model$fpc1_cat_median <-
    factor(
      ifelse(analysis_data_model$fpc1 >= fpc1_median_cut, "high", "low"),
      levels = c("high","low")
    )
  
  print(table(analysis_data_model$spatial_binary, analysis_data_model$fpc1_cat_median))
  
  cat("\nFPC1 summary:\n")
  print(summary(analysis_data_model$fpc1))
  
  cat("\nFPC1 quantiles:\n")
  print(
    quantile(
      analysis_data_model$fpc1,
      probs = c(
        0,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        1
      ),
      na.rm = TRUE
    )
  )
  
  fda_tab <- table(
    analysis_data_model$abund_binary_observed,
    analysis_data_model$fpc1_cat
  )
  
  min_fda_cell <- min(fda_tab)
  
  fda_event_tab <- with(
    analysis_data_model,
    table(
      abund_binary_observed,
      fpc1_cat,
      death
    )
  )
  
  min_fda_event_cell <- min(fda_event_tab[, , "1"])
  
  if (nrow(analysis_data_model) < 20) {
    stop("Replicate ", replicate_id, " has fewer than 20 common complete cases.")
  }
  
  # Check that all four combinations exist for every method
  category_check <- function(abundance, spatial_category) {
    tab <- table(abundance, spatial_category)
    all(dim(tab) == c(2, 2)) && all(tab > 0)
  }
  
  if (
    !category_check(
      analysis_data_model$abund_binary_observed,
      analysis_data_model$fpc1_cat
    ) ||
    !category_check(
      analysis_data_model$abund_binary_observed,
      analysis_data_model$g25_cat
    ) ||
    !category_check(
      analysis_data_model$abund_binary_observed,
      analysis_data_model$g75_cat
    )
  ) {
    stop("At least one abundance-by-spatial cell is empty in replicate ", replicate_id, ".")
  }
  
  # Fit Cox interaction models
  
  fda_warning <- FALSE
  
  fda_model <- withCallingHandlers(
    coxph(survival::Surv(mossurv,death) ~ abund_binary_observed * fpc1_cat,
                    data = analysis_data_model,
                    ties = "efron"
    ),
    
    warning = function(w) {
      
      if(grepl("coefficient may be infinite", conditionMessage(w))){
        fda_warning <<- TRUE
      }
      
      invokeRestart("muffleWarning")
    }
  )
  
  g25_model <- coxph(survival::Surv(mossurv,death) ~ abund_binary_observed * g25_cat,
                        data = analysis_data_model,
                        ties = "efron"
  )
  
  g75_model <- coxph(survival::Surv(mossurv,death) ~ abund_binary_observed * g75_cat,
                        data = analysis_data_model,
                        ties = "efron"
  )
  
  oracle_model <- coxph(survival::Surv(mossurv,death) ~ abund_binary_observed * spatial_binary,
                        data = analysis_data_model,
                        ties = "efron"
  )
  
  # Extract the interaction results
  
  fda_result <- extract_interaction_result(
    mod = fda_model,
    method_name = "FDA",
    replicate_id = replicate_id,
    interaction_pattern = "abund_binary_observed:fpc1_cat|fpc1_cat.*:abund_binary_observed"
  )
  
  g25_result <- extract_interaction_result(
    mod = g25_model,
    method_name = "G25",
    replicate_id = replicate_id,
    interaction_pattern = "abund_binary_observed:g25_cat|g25_cat.*:abund_binary_observed"
  )
  
  g75_result <- extract_interaction_result(
    mod = g75_model,
    method_name = "G75",
    replicate_id = replicate_id,
    interaction_pattern = "abund_binary_observed:g75_cat|g75_cat.*:abund_binary_observed"
  )
  
  oracle_result <- extract_interaction_result(
    mod = oracle_model,
    method_name = "Oracle",
    replicate_id = replicate_id,
    interaction_pattern = "abund_binary_observed:spatial_binary|spatial_binary.*:abund_binary_observed"
  )
  
  # Add replicate-level diagnostics
  
  n_common <- nrow(analysis_data_model)
  n_events <- sum(analysis_data_model$death == 1, na.rm = TRUE)
  censoring_percent <- mean(analysis_data_model$death == 0, na.rm = TRUE) * 100
  mean_observed_proportion <- mean(subject_data$observed_proportion, na.rm = TRUE)
  n_spatial_unavailable <- sum(subject_data$func_observed == 0, na.rm = TRUE)
  
  classification_accuracy <- function(truth,predicted) {
    
    keep <- !is.na(truth) & !is.na(predicted)
    
    if (!any(keep)) {
      return(NA_real_)
    }
    
    truth_label <- ifelse(truth[keep] == 0, "high", "low")
    mean(as.character(predicted[keep]) == truth_label)
  }
  
  fda_accuracy <- classification_accuracy(
    truth = analysis_data_model$spatial_binary,
    predicted = analysis_data_model$fpc1_cat
  )
  
  g25_accuracy <- classification_accuracy(
    truth = analysis_data_model$spatial_binary,
    predicted = analysis_data_model$g25_cat
  )
  
  g75_accuracy <- classification_accuracy(
    truth = analysis_data_model$spatial_binary,
    predicted = analysis_data_model$g75_cat
  )
  
  results <- dplyr::bind_rows(
    fda_result |>
      dplyr::mutate(classification_accuracy = fda_accuracy),
    
    g25_result |>
      dplyr::mutate(classification_accuracy = g25_accuracy),
    
    g75_result |>
      dplyr::mutate(classification_accuracy = g75_accuracy),
    
    oracle_result |>
      dplyr::mutate(classification_accuracy = 1)
    
  ) |>
    dplyr::mutate(
      n_used = n_common,
      n_events = n_events,
      censoring_percent = censoring_percent,
      mean_observed_proportion = mean_observed_proportion,
      n_spatial_unavailable = n_spatial_unavailable,
      abundance_agreement = abundance_agreement,
      observed_high_rate = observed_high_rate,
      min_fda_cell = min_fda_cell,
      min_fda_event_cell = min_fda_event_cell,
      fda_warning = fda_warning,
      g25_cutoff = g25_cutoff,
      g75_cutoff = g75_cutoff
    )
  
  message("Completed replicate ", replicate_id," using ", n_common, " common subjects.")
  
  results
}

summarize_simulation_results <- function(results, true_logHR) {
  
  results |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      n_attempted = dplyr::n(),
      n_estimable = sum(estimable, na.rm = TRUE),
      estimable_rate = mean(estimable, na.rm = TRUE),
      mean_classification_accuracy = mean(classification_accuracy, na.rm = TRUE),
      mean_estimate = mean(estimate[estimable], na.rm = TRUE),
      bias = mean(estimate[estimable] - true_logHR, na.rm = TRUE),
      empirical_SE = stats::sd(estimate[estimable], na.rm = TRUE),
      mean_model_SE = mean(SE[estimable], na.rm = TRUE),
      coverage =
        mean(
          LCL[estimable] <= exp(true_logHR) &
            UCL[estimable] >= exp(true_logHR),
          na.rm = TRUE
        ),
      detection_rate = mean(significant[estimable], na.rm = TRUE),
      .groups = "drop"
    )
}

replicate_results <- vector("list", simulation_settings$n_replicates)

for (replicate_id in seq_len(simulation_settings$n_replicates)) {
  
  cat("\nRunning replicate", replicate_id, "of",simulation_settings$n_replicates, "\n")
  
  replicate_results[[replicate_id]] <- tryCatch(
    run_single_replicate(
      replicate_id = replicate_id,
      settings = simulation_settings,
      beta_abundance = log(0.75),
      beta_spatial = log(1.25),
      beta_interaction = 0,
      num_permutations = 50
    ),
    
    error = function(e) {
      
      message("Replicate ", replicate_id, " failed: ", conditionMessage(e))
      
      tibble::tibble(
        replicate = replicate_id,
        method = c("FDA", "G25", "G75", "Oracle"),
        estimate = NA_real_,
        SE = NA_real_,
        HR = NA_real_,
        LCL = NA_real_,
        UCL = NA_real_,
        significant = NA,
        estimable = FALSE,
        failure_reason = conditionMessage(e)
      )
    }
  )
}

simulation_results <-  dplyr::bind_rows(replicate_results)

all_summary <- summarize_simulation_results(simulation_results, true_logHR = 0)

as.data.frame(all_summary)

performance_summary <-
  simulation_results |>
  dplyr::filter(estimable) |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    n_estimable = dplyr::n(),
    mean_estimate = mean(estimate),
    bias = mean(estimate - 0),
    empirical_SE = sd(estimate),
    mean_model_SE = mean(SE),
    coverage = mean(LCL <= 1 & UCL >= 1),
    detection_rate = mean(significant),
    .groups = "drop"
  )

performance_summary

simulation_results |>
  distinct(replicate, n_used, n_events) |>
  as.data.frame()

simulation_results |>
  filter(method == "FDA") |>
  select(replicate, estimate, SE, significant) |>
  as.data.frame()

simulation_results |>
  filter(method == "Oracle") |>
  select(replicate, estimate, SE, significant) |>
  as.data.frame()

paired_detection <-
  simulation_results |>
  dplyr::select(replicate, method, significant) |>
  tidyr::pivot_wider(
    names_from = method,
    values_from = significant
  )

with(paired_detection, table(FDA = FDA, G25 = G25))

with(paired_detection, table(FDA = FDA, G75 = G75))

mcnemar.test(table(paired_detection$FDA, paired_detection$G25), correct = TRUE)

mcnemar.test(table(paired_detection$FDA, paired_detection$G75), correct = TRUE)

# save(replicate_results, file = "FDA_Simulation_Study_Data_Null.RData")

### Single Images

example_design <- tidyr::expand_grid(
  abundance_group = c("Low abundance", "High abundance"),
  clustering_group = c("High clustering", "Low clustering"),
  example_id = 1:20
) |>
  dplyr::mutate(
    target_proportion = dplyr::if_else(abundance_group == "Low abundance", 0.005, 0.055),
    k = dplyr::if_else(clustering_group == "High clustering", 6L, 1L),
    group_label = paste(abundance_group, clustering_group, sep = " + ")
  )

simulate_example_subject <- function(target_proportion, k, settings, seed){
  
  set.seed(seed)
  
  current_multiplier <- settings$probability_multiplier[as.character(k)]
  positivity_parameter <- min(target_proportion * current_multiplier, 0.95)
  
  sim_object <-
    suppressMessages(
      suppressWarnings(
        scSpatialSIM::CreateSimulationObject(
          window = settings$internal_window,
          sims = 1,
          cell_types = 1
        ) |>
          scSpatialSIM::GenerateSpatialPattern(lambda = settings$lambda) |>
          scSpatialSIM::GenerateCellPositivity(
            probs = c(0, positivity_parameter),
            k = k,
            xmin = 0,
            xmax = 10,
            ymin = 0,
            ymax = 10,
            Force = TRUE,
            use_window = FALSE
          )
      )
    )
  
  dat <- scSpatialSIM::CreateSpatialList(sim_object, single_df = FALSE)[[1]]
  
  dat
}

example_spatial <- purrr::pmap(
  
  list(
    target_proportion = example_design$target_proportion,
    k = example_design$k,
    seed = 1:nrow(example_design)
  ),
  
  ~ simulate_example_subject(
    target_proportion = ..1,
    k = ..2,
    settings = simulation_settings,
    seed = ..3
  )
)

example_design$spatial_data <- example_spatial

example_summary <-
  purrr::map2_dfr(
    
    example_design$spatial_data,
    seq_len(nrow(example_design)),
    
    function(dat, i){
      
      tibble::tibble(
        row_id = i,
        group = example_design$group_label[i],
        positive_cells = sum(dat$`Cell 1 Assignment` == 1),
        total_cells = nrow(dat),
        observed_proportion = mean(dat$`Cell 1 Assignment` == 1)
      )
    }
  )

representative_examples <-
  example_summary |>
  dplyr::group_by(group) |>
  dplyr::slice_min(abs(observed_proportion - median(observed_proportion)), n = 1) |>
  dplyr::ungroup()

plot_data <-
  purrr::map2_dfr(
    representative_examples$row_id,
    representative_examples$group,
    
    function(idx, grp){
      example_design$spatial_data[[idx]] |>
        dplyr::mutate(group = grp)
    }
  )

plot_single_group <- function(idx, group_name){
  
  dat <- example_design$spatial_data[[idx]]
  
  p <- ggplot(dat, aes(x, y)) +
    
    geom_point(data = subset(dat, `Cell 1 Assignment` == 0), color = "grey80", size = 0.4) +
    geom_point(data = subset(dat, `Cell 1 Assignment` == 1), color = "red", size = 1.2) +
    coord_equal() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  file_name <- paste0(gsub("[[:space:]]+", "_", tolower(group_name)), ".tiff")
  
  ggsave(
    filename = file_name,
    plot = p,
    width = 3.5,
    height = 3.5,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  
  invisible(p)
}

purrr::walk2(
  representative_examples$row_id,
  representative_examples$group,
  plot_single_group
)
