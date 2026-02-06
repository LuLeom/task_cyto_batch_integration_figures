# Create custom funky heatmap for the paper.
# Code is highly based on openproblems report template qmd file.

params <- list(
  task_results_json = "/Users/putri.g/Documents/cytobenchmark/run_2026-02-06_14-25-48/report/task_cyto_batch_integration/local_20260206_231526/combined_output.json",
  functions = "render_report/report-functions.R",
  outdir = "output/2026-02-06_14-25-48"
)
source(params$functions)
dir.create(params$outdir, recursive = TRUE)

task_results <- jsonlite::read_json(
  params$task_results_json,
  simplifyVector = FALSE,
  simplifyDataFrame = FALSE
)

# load dataset information

dataset_info <- task_results$dataset_info

dataset_summary <- purrr::map_dfr(dataset_info, function(.dataset) {
  data.frame(
    dataset = .dataset$name,
    label = .dataset$label,
    summary = .dataset$summary
  )
})

dataset_details <- purrr::map_dfr(dataset_info, function(.dataset) {
  data.frame(
    description = .dataset$description,
    modalities = paste(.dataset$modalities, collapse = ", "),
    organisms = paste(.dataset$organisms, collapse = ", "),
    file_size_mb = .dataset$file_size_mb %||% NA_real_,
    commit = .dataset$commit,
    source_urls = .dataset$source_urls |>
      as.character() |>
      format_html_link() |>
      paste(collapse = "<br>"),
    common_dataset_names = paste(.dataset$common_dataset_names, collapse = ", "),
    date_created = .dataset$date_created %||% NA_character_
  )
})

# load method information
method_info <- task_results$method_info

method_summary <- purrr::map_dfr(method_info, function(.method) {
  data.frame(
    method = .method$name,
    label = .method$label,
    type = .method$type,
    summary = .method$summary
  )
})

method_details <- purrr::map_dfr(method_info, function(.method) {
  method_data <- purrr::map(.method, \(.x) {ifelse(is.null(.x), "", .x)})

  data.frame(
    description = method_data$description,
    commit = method_data$commit,
    version = method_data$version,
    link_code = method_data$link_code,
    link_documentation = method_data$link_documentation,
    link_implementation = method_data$link_implementation,
    link_container_image = method_data$link_container_image
  )
})

# load metric information
metric_info <- task_results$metric_info

metric_summary <- purrr::map_dfr(metric_info, function(.metric) {
  data.frame(
    metric = .metric$name,
    label = .metric$label,
    summary = .metric$summary
  )
})

metric_details <- purrr::map_dfr(metric_info, function(.metric) {
  metric_data <- purrr::map(.metric, \(.x) {ifelse(is.null(.x), "", .x)})

  data.frame(
    description = metric_data$description,
    component_name = metric_data$component_name,
    commit = metric_data$commit,
    version = metric_data$version,
    maximize = metric_data$maximize,
    link_implementation = metric_data$link_implementation,
    link_container_image = metric_data$link_container_image
  )
})

# normalise scores

dataset_names <- purrr::map_chr(dataset_info, "name")
method_names <- purrr::map_chr(method_info, "name")
metric_names <- purrr::map_chr(metric_info, "name")

control_method_names <- method_summary$method[method_summary$type == "control_method"]

n_controls <- purrr::map_dfr(task_results$results, function(.result) {
  data.frame(
    dataset = .result$dataset_name,
    method = .result$method,
    succeeded = .result$succeeded
  )
}) |>
  dplyr::filter(
    method %in% control_method_names,
    succeeded
  ) |>
  dplyr::group_by(dataset) |>
  dplyr::count(name = "n_controls") |>
  dplyr::ungroup() |>
  dplyr::mutate(dataset = factor(dataset, levels = dataset_names)) |>
  tidyr::complete(dataset, fill = list(n_controls = 0))

has_controls <- all(n_controls$n_controls >= 2)

dataset_details <- purrr::map_dfr(dataset_info, function(.dataset) {
  data.frame(
    dataset = .dataset$name,
    dataset_label = .dataset$label
  )
}) |>
  dplyr::arrange(dataset)

# TODO: custom naming for now. but real fix should fix the actual h5ad file
dataset_details <- dataset_details |>
  dplyr::mutate(
    dataset_label = dplyr::case_when(
      dataset == "lille_spectral_flow_cytometry" ~ "Mouse Spleen Spectral Flow Cytometry",
      TRUE ~ dataset_label
    )
  )

method_details <- purrr::map_dfr(method_info, function(.method) {
  data.frame(
    method = .method$name,
    method_label = .method$label,
    method_type = .method$type
  )
}) |>
  dplyr::mutate(
    method_label = stringr::str_replace_all(
      method_label,
      stringr::regex("control", ignore_case = TRUE),
      "reference"
    )
  ) |>
  dplyr::arrange(method)

# rename "Shuffle Integration — within cell type" to Shuffle Integration Within cell type"
method_details <- method_details |>
  dplyr::mutate(
    method_label = stringr::str_replace_all(
      method_label,
      "Shuffle Integration — within cell type",
      "Shuffle Integration Within cell type"
    )
  )

# read in the metric_types to get the metric type and join
metric_types <- read.csv("data/metric_types.csv") |>
  dplyr::select(metric, metric_type)

metric_details <- purrr::map_dfr(metric_info, function(.metric) {
  data.frame(
    metric = .metric$name,
    metric_label = .metric$label,
    metric_maximize = .metric$maximize
  )
}) |> 
  dplyr::left_join(metric_types, by = "metric") |>
  dplyr::mutate(metric_label = ifelse(metric_maximize, metric_label, paste0(metric_label, "*"))) |>
  dplyr::arrange(metric)

metric_reverse <- metric_details |>
  dplyr::filter(metric_maximize == FALSE) |>
  dplyr::pull(metric)

scores <- purrr::map_dfr(task_results$results, function(.result) {
  if (!.result$succeeded) {
    return(NULL)
  }

  if (all(unlist(purrr::map(.result$metric_components, "succeeded")) == FALSE)) {
    return(NULL)
  }

  if (length(.result$metric_values) == 0) {
    warning(
      "At least one metric component succeeded but there are no metric values ",
      "for method '", .result$method_name, "' on dataset '",
      .result$dataset_name, "'",
      call. = FALSE
    )
    return(NULL)
  }

  tibble::tibble(
    dataset = .result$dataset_name,
    method = .result$method,
    metric = unlist(.result$metric_names),
    value = unlist(.result$metric_values)
  )
})

control_ranges <- scores |>
  dplyr::left_join(method_details, by = "method") |>
  dplyr::filter(method_type == "control_method") |>
  dplyr::group_by(dataset, metric) |>
  dplyr::summarise(
    control_min = min(value, na.rm = TRUE),
    control_max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

scaled_scores <- scores |>
  dplyr::left_join(control_ranges, by = c("dataset", "metric")) |>
  dplyr::mutate(
    scaled_value = (value - control_min) / (control_max - control_min),
    scaled_value = dplyr::if_else(
      metric %in% metric_reverse,
      1 - scaled_value,
      scaled_value
    )
  )

complete_scores <- tidyr::expand_grid(
  dataset = dataset_names,
  method = method_names,
  metric = metric_names
) |>
  dplyr::left_join(dataset_details, by = "dataset") |>
  dplyr::relocate(method, metric, .after = dplyr::last_col()) |>
  dplyr::left_join(method_details, by = "method") |>
  dplyr::relocate(metric, .after = dplyr::last_col()) |>
  dplyr::left_join(metric_details, by = "metric") |>
  dplyr::left_join(scaled_scores, by = c("dataset", "method", "metric")) |>
#   tidyr::replace_na(list(scaled_value = 0)) |>
  dplyr::arrange(dataset, method, metric)



# we don't want emd_mean_ct_vert here for human blood mass cytometry data
# as it could not be calculated!
complete_scores <- complete_scores |>
  dplyr::filter(
    !(dataset == "human_blood_mass_cytometry" & metric == "emd_mean_ct_vert")
  )

# results table

mean_scores <- complete_scores |>
  dplyr::group_by(dataset, method, metric_type) |>
  dplyr::summarise(
    mean_score = aggregate_scores(scaled_value),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = metric_type,
    values_from = mean_score,
    names_prefix = "mean_"
  ) |> 
  dplyr::rowwise() |>
  dplyr::mutate(
    mean_score = mean(c(mean_batch_mixing, mean_bio_conservation), na.rm = TRUE)
  ) |>
  dplyr::ungroup()

# export raw values separately
dataset_scores <- complete_scores |>
  dplyr::select(dataset, method, metric, value) |>
  tidyr::pivot_wider(
    names_from = metric,
    values_from = value
  )
write.csv(dataset_scores, file.path(params$outdir, "raw_dataset_scores.csv"), row.names = FALSE)

# by dataset scores, the NA for emd mean vertical for human blood mass cytometry will be visible.
# because of the nature of the group by.
# hence we need to disable the aggregate scores so NA is not replaced by 0..
# as we just didn't calculate it. It does not mean it existed!
dataset_scores <- complete_scores |>
  dplyr::select(dataset, method, metric, scaled_value) |>
  tidyr::pivot_wider(
    names_from = metric,
    values_from = scaled_value
  ) |>
  dplyr::left_join(mean_scores, by = c("dataset", "method"))


overall_scores <- dataset_scores |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    dataset = "overall",
    dplyr::across(
      tidyselect::where(is.numeric),
      aggregate_scores
    ),
    .groups = "drop"
  )


resources <- purrr::map_dfr(task_results$results, function(.result) {
  data.frame(
    dataset = .result$dataset_name,
    method = .result$method
  )
}) |>
  dplyr::mutate(
    run_duration_secs = purrr::map(task_results$results, "run_duration_secs"),
    run_cpu_pct = purrr::map(task_results$results, "run_cpu_pct"),
    run_peak_memory_mb = purrr::map(task_results$results, "run_peak_memory_mb"),
    run_disk_read_mb = purrr::map(task_results$results, "run_disk_read_mb"),
    run_disk_write_mb = purrr::map(task_results$results, "run_disk_write_mb")
  ) |>
  # Summarise per task
  dplyr::mutate(
    run_cpu_pct = purrr::map_dbl(run_cpu_pct, function(.values) {
      if (length(.values) == 0) {
        return(NA_real_)
      }

      mean(unlist(.values), na.rm = TRUE)
    }),
    run_peak_memory_mb = purrr::map_dbl(run_peak_memory_mb, function(.values) {
      if (length(.values) == 0) {
        return(NA_real_)
      }

      max(unlist(.values), na.rm = TRUE)
    }),
    run_disk_read_mb = purrr::map_dbl(run_disk_read_mb, function(.values) {
      if (length(.values) == 0) {
        return(NA_real_)
      }

      sum(unlist(.values), na.rm = TRUE)
    }),
    run_disk_write_mb = purrr::map_dbl(run_disk_write_mb, function(.values) {
      if (length(.values) == 0) {
        return(NA_real_)
      }

      sum(unlist(.values), na.rm = TRUE)
    }),
    run_duration_secs = purrr::map_dbl(run_duration_secs, function(.values) {
      if (length(.values) == 0) {
        return(NA_real_)
      }

      sum(unlist(.values), na.rm = TRUE)
    })
  ) |>
  # Summarise by method
  dplyr::group_by(method) |>
  dplyr::summarise(
    mean_cpu_pct = mean(run_cpu_pct, na.rm = TRUE),
    mean_peak_memory_mb = mean(run_peak_memory_mb, na.rm = TRUE),
    mean_disk_read_mb = mean(run_disk_read_mb, na.rm = TRUE),
    mean_disk_write_mb = mean(run_disk_write_mb, na.rm = TRUE),
    mean_duration_secs = mean(run_duration_secs, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    mean_peak_memory_mb_log = -log10(mean_peak_memory_mb),
    mean_peak_memory_label = paste0(" ", label_memory(mean_peak_memory_mb), " "),
    mean_disk_read_mb_log = -log10(mean_disk_read_mb),
    mean_disk_read_label = paste0(" ", label_memory(mean_disk_read_mb), " "),
    mean_disk_write_mb_log = -log10(mean_disk_write_mb),
    mean_disk_write_label = paste0(" ", label_memory(mean_disk_write_mb), " "),
    mean_duration_secs_log = -log10(mean_duration_secs),
    mean_duration_label = paste0(" ", label_time(mean_duration_secs), " ")
  )

# collate figure data

figure_data <- overall_scores |>
  dplyr::select(-dataset) |>
  dplyr::relocate(mean_score, mean_batch_mixing, mean_bio_conservation, .after = method) |>
  dplyr::left_join(
    resources |>
      dplyr::select(
        method,
        mean_peak_memory_mb_log,
        mean_peak_memory_label,
        mean_duration_secs_log,
        mean_duration_label
      ),
    by = "method"
  ) |>
  # Resources are not 0-1 so need to be rescaled
  dplyr::mutate(
    mean_peak_memory_mb_log = scales::rescale(mean_peak_memory_mb_log),
    mean_duration_secs_log = scales::rescale(mean_duration_secs_log)
  ) |>
  dplyr::arrange(dplyr::desc(mean_score)) |>
  dplyr::mutate(
    method = factor(
      method,
      levels = method_details$method,
      labels = method_details$method_label
    )
  ) |>
  dplyr::rename(id = method)

# relocate the metrics so we show batch mixing then bio
batch_mix_metrics <- metric_details |>
  dplyr::filter(metric_type == "batch_mixing") |>
  dplyr::select(metric, metric_label) |>
  dplyr::arrange(metric)

bio_conserve_metrics <- metric_details |>
  dplyr::filter(metric_type == "bio_conservation") |>
  dplyr::select(metric, metric_label) |>
  dplyr::arrange(metric)

figure_data <- figure_data |>
  dplyr::relocate(
    tidyselect::all_of(batch_mix_metrics$metric),
    .after = mean_bio_conservation
  ) |>
  dplyr::relocate(
    tidyselect::all_of(bio_conserve_metrics$metric),
    .after = tidyselect::all_of(batch_mix_metrics$metric)
  )

column_info <- tibble::tibble(
  id = colnames(figure_data),
  name = c(
    "Method",
    "Overall mean score",
    "Batch mixing (mean)",
    "Bio conservation (mean)",
    batch_mix_metrics$metric_label,
    bio_conserve_metrics$metric_label,
    "Peak memory",
    "",
    "Duration",
    ""
  ),
  geom = c(
    "text",
    rep("funkyrect", 3),
    rep("funkyrect", length(metric_names)),
    rep(c("rect", "text"), 2))
  ,
  group = c(
    NA,
    rep("overall", 3),
    rep("batch mixing metrics", length(batch_mix_metrics$metric)),
    rep("bio conservation metrics", length(bio_conserve_metrics$metric)),
    rep("resources", 4)
  ),
  palette = c(
    NA,
    rep("overall_palette", 3),
    rep("batch_mix_palette", length(batch_mix_metrics$metric)),
    rep("bio_conserve_palette", length(bio_conserve_metrics$metric)),
    rep(c("resources_palette", "black"), 2)
  ),
  width = c(
    10,
    rep(1, 3),
    rep(1, length(metric_names)),
    rep(1, 4)
  ),
  overlay = c(
    FALSE,
    rep(FALSE, 3),
    rep(FALSE, length(metric_names)),
    rep(c(FALSE, TRUE), 2)
  ),
  hjust = c(
    0,
    rep(0.5, 3),
    rep(0.5, length(metric_names)),
    rep(0.5, 4)
  )
)

column_groups <- tibble::tibble(
  group = c("overall", "batch mixing metrics", "bio conservation metrics", "resources"),
  category = c("Overall", "Batch mixing metrics", "Bio conservation metrics", "Resources"),
  palette = c("overall_palette", "batch_mix_palette", "bio_conserve_palette", "resources_palette"),
)

palettes <- list(
  overall_palette = "Blues",
  batch_mix_palette = "Reds",
  bio_conserve_palette = "Greens",
  resources_palette = "YlOrBr",
  black = c("black", "black")
)

legends <- list(
  list(
    geom = "funkyrect",
    title = "Score",
    colour = "white",
    label_hjust = 0
  ),
  list(
    palette = "overall_palette",
    enabled = FALSE
  ),
  list(
    palette = "batch_mix_palette",
    enabled = FALSE
  ),
  list(
    palette = "bio_conserve_palette",
    enabled = FALSE
  ),
  list(
    palette = "resources_palette",
    enabled = FALSE
  )
)

# TODO the footnote is not working yet..
funkyheatmap::funky_heatmap(
  figure_data,
  column_info = column_info,
  column_groups = column_groups,
  palettes = palettes,
  legends = legends,
  scale_column = FALSE,
  position_args = funkyheatmap::position_arguments(
    col_space = 0.5,
    col_bigspace = 1.5,
    col_annot_offset = 4
  )
) + ggplot2::labs(caption = "Note: Metrics with an asterisk (*) have been reversed so higher values indicate better performance.")
ggplot2::ggsave(file.path(params$outdir, "funky_heatmap.png"), width = 20, height = 16, limitsize = FALSE)











