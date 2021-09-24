#' Title Create calibration plot based on observed and predicted outcomes.
#'
#' @param data Data include observed and predicted outcomes.
#' @param obs Name of observed outcome in the input data.
#' @param follow_up Name of follow-up time (if applicable) in the input data.
#' @param pred_1 Name of first predicted outcome in the input data.
#' @param pred_2 Name of second (if applicable) predicted outcome in the input data.
#' @param nTiles Number of tiles (e.g., 10 for deciles) in the calibration plot.
#' @param legendPosition Legend position on the calibration plot.
#' @param title Title on the calibration plot.
#' @param x_lim Limits of x-axis on the calibration plot.
#' @param y_lim Limits of y-axis on the calibration plot.
#' @param model_label_1 Label of the first prediction set.
#' @param model_label_2 Label of the first prediction set.
#' @param xlab Label of x-axis on the calibration plot.
#' @param ylab Label of y-axis on the calibration plot.
#' @param points_col_list Points' color on the calibration plot.
#' @param data_summary Logical indicates whether a summary of the predicted and observed outcomes.
#' needs to be included in the output.
#' @examples 
#' library(predtools)
#' library(dplyr)
#' x <- rnorm(100, 10, 2)
#' y <- x + rnorm(100,0, 1)
#' data <- tibble(x,y)
#' calibration_plot(data, obs = x, pred_1 = y)
#' @return Returns calibration plot (a ggplot object) and a dataset including summary statistics of
#' the predicted and observed outcomes (if data_summary set to be TRUE).
#' @export
#'
calibration_plot <- function(data,
                             obs,
                             follow_up = NULL,
                             pred_1,
                             pred_2 = NULL,
                             nTiles = 10,
                             legendPosition = "right",
                             title = NULL,
                             x_lim = NULL,
                             y_lim = NULL,
                             model_label_1 = NULL,
                             model_label_2 = NULL,
                             xlab = "Prediction",
                             ylab = "Observation",
                             points_col_list = NULL,
                             data_summary = FALSE) {
  
  if (! exists("obs") | ! exists("pred_1")) stop("obs and pred_1 can not be null.")
  
  data %<>%
    mutate(decile = ntile(!!sym(pred_1), nTiles))
  
  if (is.null(follow_up)) data$follow_up <- 1

  if (! is.null(pred_2)) {
    data %>%
      group_by(decile) %>%
      summarise(obsRate = mean(!!sym(obs) / follow_up, na.rm = T),
                obsRate_SE = sd(!!sym(obs) / follow_up, na.rm = T) / sqrt(n()),
                obsNo = n(),
                model = ifelse(is.null(model_label_1), "Model 1", model_label_1),
                predRate = mean(!!sym(pred_1), na.rm = T)) -> dataDec_mod_1
    data %>%
      group_by(decile) %>%
      summarise(obsRate = mean(!!sym(obs) / follow_up, na.rm = T),
                obsRate_SE = sd(!!sym(obs) / follow_up, na.rm = T) / sqrt(n()),
                obsNo = n(),
                model = ifelse(is.null(model_label_1), "Model 2", model_label_2),
                predRate = mean(!!sym(pred_2), na.rm = T)) -> dataDec_mod_2
    dataDec_mods <- rbind(dataDec_mod_1, dataDec_mod_2)
  }
  else {
    data %>%
      group_by(decile) %>%
      summarise(obsRate = mean(!!sym(obs) / follow_up, na.rm = T),
                obsRate_SE = sd(!!sym(obs) / follow_up, na.rm = T) / sqrt(n()),
                obsNo = n(),
                model = model_label_1,
                predRate = mean(!!sym(pred_1), na.rm = T)) -> dataDec_mods
  }

  dataDec_mods$obsRate_UCL <- dataDec_mods$obsRate + 1.96 * dataDec_mods$obsRate_SE
  dataDec_mods$obsRate_LCL <- dataDec_mods$obsRate - 1.96 * dataDec_mods$obsRate_SE
  
  
  if (! is.null(pred_2)) {
    calibPlot_obj <-
      ggplot(data = dataDec_mods, aes(y = obsRate, x = predRate, group = model, color = model)) +
      geom_point() +
      lims(x = ifelse(c(is.null(x_lim), is.null(x_lim)), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim),
           y = ifelse(c(is.null(y_lim), is.null(y_lim)),
                      c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)) +
      geom_errorbar(aes(ymax = obsRate_UCL, ymin = obsRate_LCL)) +
      geom_abline(intercept = 0, slope = 1) +
      scale_color_manual(values = ifelse(c(is.null(points_col_list), is.null(points_col_list)),
                                         c(ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[4],
                                           ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5]),
                                         points_col_list)) +
      labs(x = ifelse(is.null(xlab), pred_1, xlab),
           y = ifelse(is.null(ylab), obs, ylab),
           title = title) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(fill = "white"),
            axis.text = element_text(colour = "black", size = 12),
            legend.position = legendPosition)
  }
  else {
    calibPlot_obj <-
      ggplot(data = dataDec_mods, aes(y = obsRate, x = predRate)) +
      geom_point(color = ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5]) +
      lims(x = ifelse(c(is.null(x_lim), is.null(x_lim)), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim),
           y = ifelse(c(is.null(y_lim),  is.null(y_lim)),
                      c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)) +
      geom_errorbar(aes(ymax = obsRate_UCL, ymin = obsRate_LCL),
                    col = ifelse(is.null(points_col_list),
                                 ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5],
                                 points_col_list)) +
      geom_abline(intercept = 0, slope = 1) +
      scale_color_manual(values = ifelse(is.null(points_col_list),
                                         ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5],
                                         points_col_list)) +
      labs(x = ifelse(is.null(xlab), pred_1, xlab),
           y = ifelse(is.null(ylab), obs, ylab),
           title = title) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text = element_text(colour = "black", size = 12),
            legend.position = legendPosition)
  }
  
  res_list <- list(calibration_plot = calibPlot_obj)
  if (data_summary) res_list$data_summary <- dataDec_mods

  return(res_list)
}
