#' Title Create calibration plot based on observed and predicted outcomes.
#'
#' @param data Data include observed and predicted outcomes.
#' @param obs Name of observed outcome in the input data.
#' @param follow_up Name of follow-up time (if applicable) in the input data.
#' @param pred Name of first predicted outcome in the input data.
#' @param group Name of grouping column (if  applicable) in the input data.
#' @param nTiles Number of tiles (e.g., 10 for deciles) in the calibration plot.
#' @param legendPosition Legend position on the calibration plot.
#' @param title Title on the calibration plot.
#' @param x_lim Limits of x-axis on the calibration plot.
#' @param y_lim Limits of y-axis on the calibration plot.
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
#' data <- data.frame(x, y)
#' calibration_plot(data, obs = "x", pred = "y")
#' @return Returns calibration plot (a ggplot object) and a dataset including summary statistics of
#' the predicted and observed outcomes (if data_summary set to be TRUE).
#' @export
#'
calibration_plot <- function(data,
                             obs,
                             follow_up = NULL,
                             pred,
                             group = NULL,
                             nTiles = 10,
                             legendPosition = "right",
                             title = NULL,
                             x_lim = NULL,
                             y_lim = NULL,
                             xlab = "Prediction",
                             ylab = "Observation",
                             points_col_list = NULL,
                             data_summary = FALSE) {
  
  if (! exists("obs") | ! exists("pred")) stop("obs and pred can not be null.")
  
  n_groups <- length(unique(data[ , group]))
  
  data %<>%
    mutate(decile = ntile(!!sym(pred), nTiles))
  
  if (is.null(follow_up)) data$follow_up <- 1

  if (! is.null(group)) {
    data %>%
      group_by(.data$decile, !!sym(group)) %>%
      summarise(obsRate = mean(!!sym(obs) / follow_up, na.rm = T),
                obsRate_SE = sd(!!sym(obs) / follow_up, na.rm = T) / sqrt(n()),
                obsNo = n(),
                predRate = mean(!!sym(pred), na.rm = T)) -> dataDec_mods
    colnames(dataDec_mods)[colnames(dataDec_mods) == "group"] <- group
  }
  else {
    data %>%
      group_by(.data$decile) %>%
      summarise(obsRate = mean(!!sym(obs) / follow_up, na.rm = T),
                obsRate_SE = sd(!!sym(obs) / follow_up, na.rm = T) / sqrt(n()),
                obsNo = n(),
                predRate = mean(!!sym(pred), na.rm = T)) -> dataDec_mods
  }

  dataDec_mods$obsRate_UCL <- dataDec_mods$obsRate + 1.96 * dataDec_mods$obsRate_SE
  dataDec_mods$obsRate_LCL <- dataDec_mods$obsRate - 1.96 * dataDec_mods$obsRate_SE
  
  dataDec_mods <- as.data.frame(dataDec_mods)
  
  if (! is.null(group)) {
    dataDec_mods[ , group] <- factor(dataDec_mods[ , group])
    calibPlot_obj <-
      ggplot(data = dataDec_mods, aes(y = .data$obsRate, x = .data$predRate, group = !!sym(group), color = !!sym(group))) +
      geom_point() +
      lims(x = ifelse(rep(is.null(x_lim), 2), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim),
           y = ifelse(rep(is.null(y_lim), 2), c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)) +
      geom_errorbar(aes(ymax = .data$obsRate_UCL, ymin = .data$obsRate_LCL)) +
      geom_abline(intercept = 0, slope = 1) +
      scale_color_manual(values = ifelse(rep(is.null(points_col_list), n_groups),
                                         (ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[c(4 : 8, 1 : 3)])[c(1 : n_groups)],
                                         points_col_list)) +
      labs(x = ifelse(is.null(xlab), pred, xlab),
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
      ggplot(data = dataDec_mods, aes(y = .data$obsRate, x = .data$predRate)) +
      geom_point(color = ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5]) +
      lims(x = ifelse(rep(is.null(x_lim), 2), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim),
           y = ifelse(rep(is.null(y_lim), 2), c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)) +
      geom_errorbar(aes(ymax = .data$obsRate_UCL, ymin = .data$obsRate_LCL),
                    col = ifelse(is.null(points_col_list),
                                 ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5],
                                 points_col_list)) +
      geom_abline(intercept = 0, slope = 1) +
      scale_color_manual(values = ifelse(is.null(points_col_list),
                                         ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5],
                                         points_col_list)) +
      labs(x = ifelse(is.null(xlab), pred, xlab),
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
