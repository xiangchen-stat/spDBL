#' Quick raster heatmap
#'
#' Creates a \code{ggplot2} raster heatmap using the \code{col_bgr} colour
#' palette with a squished colour scale capped at \code{max_y}.
#'
#' @param dt Data frame with columns \code{x} (horizontal position),
#'   \code{y} (vertical position), and \code{fill} (value to display).
#' @param max_y Numeric scalar. Upper limit of the colour scale. Values above
#'   this are squished to the maximum colour.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{quick_save}}
#' @export
quick_heat <- function(dt, max_y){
  p <- ggplot(dt, aes(x = x, y = y, fill = fill)) +
    geom_raster() +
    scale_fill_gradientn(colours = col_bgr,
                         limits = c(0, max_y),
                         oob = scales::squish) +
    labs(x = "x", y = "y", fill = "Value")+
    theme(text = element_text(size=28),
          legend.text = element_text(size = 28),
          legend.key.size = unit(1.5, "cm"))
  return(p)
}

#' Save a ggplot to a timestamped PNG file
#'
#' Saves \code{plot} to a PNG file in \code{path_fig} whose name is
#' \code{filename} followed by the current Unix timestamp, ensuring unique
#' file names across repeated calls.
#'
#' @param filename Character. Base name for the output file (no extension).
#' @param path_fig Character. Directory path where the file is saved.
#' @param plot A \code{ggplot} object to save.
#'
#' @return Invisibly returns \code{NULL} (called for its side effect).
#'
#' @seealso \code{\link{quick_heat}}
#' @export
quick_save <- function(filename, path_fig, plot){
  ggsave(filename = paste(filename, "_", as.numeric(Sys.time()), ".png", sep = ""),
         path = path_fig,
         plot = plot,
         device = "png",
         width = 42,
         height = 35,
         units = "cm",
         dpi = 100
  )
}

# Emulation ----

#' Plot a 3-by-3 panel of heatmaps across selected time stamps
#'
#' Produces nine raster heatmaps arranged in a 3-by-3 grid, one for each time
#' stamp in \code{tstamp}. Each panel shows the spatial field for a given input
#' (row) index and time step.
#'
#' @param dat List of length \eqn{\geq \max(\code{tstamp})}. Each element is a
#'   matrix of dimension \code{c(n_inputs, Nx * Ny)}.
#' @param input_num Integer. Row index within each time-step matrix to plot.
#' @param tstamp Integer vector of length 9. Time step indices to display.
#' @param max_y Numeric scalar. Upper limit of the colour scale.
#' @param Nx Integer. Number of grid points in the x-direction.
#' @param Ny Integer. Number of grid points in the y-direction.
#' @param nT Integer. Total number of time steps in \code{dat}.
#' @param filename Character. Base name for the saved PNG (used when
#'   \code{savei = TRUE}). Defaults to \code{"plot_panel"}.
#' @param savei Logical. Whether to save the panel to disk. Defaults to
#'   \code{TRUE}.
#'
#' @return A list of nine \code{ggplot} objects (one per panel).
#'
#' @export
plot_panel_heatmap_9 <- function(dat, input_num, tstamp, max_y, Nx, Ny, nT, filename = "plot_panel", savei = T){
  plot_ls <- list()
  ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
  ind_plot <- 1

  for (i in tstamp) {
    temp <- dat[[i]][input_num,]
    rownames(temp) <- NULL
    colnames(temp) <- NULL
    dt <- data.frame(row = ind_sp$row, col = ind_sp$col, sol = temp)%>%
      as.data.frame()

    p <- ggplot(dt, aes(x = col, y = row, fill = sol)) +
      geom_raster() +
      scale_fill_gradientn(colours = col_bgr,
                           limits = c(0, max_y),
                           oob = scales::squish) +
      labs(x = "x", y = "y", fill = "Value")+
      theme(text = element_text(size=28),
            legend.text = element_text(size = 28),
            legend.key.size = unit(1.5, "cm"))

    plot_ls[[ind_plot]] <- p
    ind_plot <- ind_plot + 1
  }
  if(savei){
    ggsave(filename = paste(filename, ".png", sep = ""),
           path = path_fig,
           plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
                            plot_ls[[4]], plot_ls[[5]], plot_ls[[6]],
                            plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
                            ncol = 3, nrow = 3,
                            labels = c(paste("PDE: t =", tstamp[1]-1),
                                       paste("PDE: t =", tstamp[2]-1),
                                       paste("PDE: t =", tstamp[3]-1),
                                       paste("PDE: t =", tstamp[4]-1),
                                       paste("PDE: t =", tstamp[5]-1),
                                       paste("PDE: t =", tstamp[6]-1),
                                       paste("PDE: t =", tstamp[7]-1),
                                       paste("PDE: t =", tstamp[8]-1),
                                       paste("PDE: t =", tstamp[9]-1)),
                            font.label = list(size = 28),
                            vjust = 1.2,
                            align = "hv",
                            common.legend = T,
                            legend = "right"

           ),
           device = "png",
           width = 60,
           height = 50,
           units = "cm",
           dpi = 100
    )
  }
  return(plot_ls)
}


#' Compute median and 95% credible interval across rows
#'
#' For each row of the matrix \code{X}, computes the median and the 2.5th and
#' 97.5th percentiles across columns (samples).
#'
#' @param X Numeric matrix. Rows correspond to spatial locations or parameters;
#'   columns correspond to posterior samples.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{med}{Row-wise median.}
#'     \item{lower}{Row-wise 2.5th percentile.}
#'     \item{upper}{Row-wise 97.5th percentile.}
#'   }
#'
#' @seealso \code{\link{cal_errorbar_mean}}
#' @export
cal_errorbar <- function(X){
  out <- data.frame(med = apply(X = X, MARGIN = 1, FUN = median),
                    lower = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.025),
                    upper = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.975))
  return(out)
}


#' Compute mean and 95% credible interval across rows
#'
#' For each row of the matrix \code{X}, computes the mean and the 2.5th and
#' 97.5th percentiles across columns (samples).
#'
#' @param X Numeric matrix. Rows correspond to spatial locations or parameters;
#'   columns correspond to posterior samples.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{med}{Row-wise mean.}
#'     \item{lower}{Row-wise 2.5th percentile.}
#'     \item{upper}{Row-wise 97.5th percentile.}
#'   }
#'
#' @seealso \code{\link{cal_errorbar}}
#' @export
cal_errorbar_mean <- function(X){
  out <- data.frame(med = apply(X = X, MARGIN = 1, FUN = mean),
                    lower = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.025),
                    upper = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.975))
  return(out)
}

# calibration ----

#' Plot a 3-by-3 panel of calibration heatmaps
#'
#' Produces nine raster heatmaps for calibration data arranged in a 3-by-3
#' grid, one for each time stamp in \code{tstamp}. Locations are specified via
#' \code{loc_cal} rather than a regular grid.
#'
#' @param dat Numeric matrix of dimension \code{c(nT, n_locations)}. Rows are
#'   time steps; columns are calibration locations.
#' @param tstamp Integer vector of length 9. Time-step indices (rows of
#'   \code{dat}) to display.
#' @param max_y Numeric scalar. Upper limit of the colour scale.
#' @param loc_cal Data frame with columns \code{row} and \code{col} giving the
#'   spatial coordinates of calibration locations.
#' @param Nx Integer. Number of grid points in the x-direction (used for
#'   reference only).
#' @param Ny Integer. Number of grid points in the y-direction (used for
#'   reference only).
#' @param filename Character. Base name for the saved PNG. Defaults to
#'   \code{"plot_panel"}.
#' @param savei Logical. Whether to save the panel to disk. Defaults to
#'   \code{TRUE}.
#'
#' @return A list of nine \code{ggplot} objects.
#'
#' @seealso \code{\link{plot_panel_heatmap_9_cal_nolab}}
#' @export
plot_panel_heatmap_9_cal <- function(dat, tstamp, max_y, loc_cal, Nx, Ny, filename = "plot_panel", savei = T){
  plot_ls <- list()
  ind_plot <- 1 # start counting

  for (i in tstamp) {
    temp <- dat[i,]
    rownames(temp) <- NULL
    colnames(temp) <- NULL
    dt <- data.frame(row = loc_cal$row, col = loc_cal$col, sol = temp)%>%
      as.data.frame()

    p <- ggplot(dt, aes(x = col, y = row, fill = sol)) +
      geom_raster() +
      scale_fill_gradientn(colours = col_bgr,
                           limits = c(0, max_y),
                           oob = scales::squish) +
      labs(x = "x", y = "y", fill = "Value")+
      theme(text = element_text(size=28),
            legend.text = element_text(size = 28),
            legend.key.size = unit(1.5, "cm"))

    plot_ls[[ind_plot]] <- p
    ind_plot <- ind_plot + 1
  }

  if(savei){
    ggsave(filename = paste(filename, ".png", sep = ""),
           path = path_fig,
           plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
                            plot_ls[[4]], plot_ls[[5]], plot_ls[[6]],
                            plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
                            ncol = 3, nrow = 3,
                            labels = c(paste("t =", tstamp[1]-1),
                                       paste("t =", tstamp[2]-1),
                                       paste("t =", tstamp[3]-1),
                                       paste("t =", tstamp[4]-1),
                                       paste("t =", tstamp[5]-1),
                                       paste("t =", tstamp[6]-1),
                                       paste("t =", tstamp[7]-1),
                                       paste("t =", tstamp[8]-1),
                                       paste("t =", tstamp[9]-1)),
                            font.label = list(size = 28),
                            vjust = 1.2,
                            align = "hv",
                            common.legend = T,
                            legend = "right"

           ),
           device = "png",
           width = 60,
           height = 50,
           units = "cm",
           dpi = 100
    )
  }

  return(plot_ls)
}

#' Plot a 3-by-3 panel of calibration heatmaps without axis labels
#'
#' Like \code{\link{plot_panel_heatmap_9_cal}} but suppresses x- and y-axis
#' labels and supports a custom lower colour-scale limit (\code{min_y}).
#'
#' @param dat Numeric matrix of dimension \code{c(nT, n_locations)}.
#' @param tstamp Integer vector of length 9. Time-step indices to display.
#' @param min_y Numeric scalar. Lower limit of the colour scale. Defaults to
#'   \code{0}.
#' @param max_y Numeric scalar. Upper limit of the colour scale.
#' @param loc_cal Data frame with columns \code{row} and \code{col}.
#' @param Nx Integer. Number of grid points in the x-direction.
#' @param Ny Integer. Number of grid points in the y-direction.
#' @param filename Character. Base name for the saved PNG. Defaults to
#'   \code{"plot_panel"}.
#' @param savei Logical. Whether to save the panel to disk. Defaults to
#'   \code{TRUE}.
#'
#' @return A list of nine \code{ggplot} objects.
#'
#' @seealso \code{\link{plot_panel_heatmap_9_cal}}
#' @export
plot_panel_heatmap_9_cal_nolab <- function(dat, tstamp, min_y = 0, max_y, loc_cal, Nx, Ny, filename = "plot_panel", savei = T){
  plot_ls <- list()
  ind_plot <- 1 # start counting

  for (i in tstamp) {
    temp <- dat[i,]
    rownames(temp) <- NULL
    colnames(temp) <- NULL
    dt <- data.frame(row = loc_cal$row, col = loc_cal$col, sol = temp)%>%
      as.data.frame()

    p <- ggplot(dt, aes(x = col, y = row, fill = sol)) +
      geom_raster() +
      scale_fill_gradientn(colours = col_bgr,
                           limits = c(min_y, max_y),
                           oob = scales::squish) +
      labs(x = NULL, y = NULL, fill = "Value")+
      theme(text = element_text(size=28),
            legend.text = element_text(size = 28),
            legend.key.size = unit(1.5, "cm"))

    plot_ls[[ind_plot]] <- p
    ind_plot <- ind_plot + 1
  }

  if(savei){
    ggsave(filename = paste(filename, ".png", sep = ""),
           path = path_fig,
           plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
                            plot_ls[[4]], plot_ls[[5]], plot_ls[[6]],
                            plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
                            ncol = 3, nrow = 3,
                            labels = c(paste("t =", tstamp[1]-1),
                                       paste("t =", tstamp[2]-1),
                                       paste("t =", tstamp[3]-1),
                                       paste("t =", tstamp[4]-1),
                                       paste("t =", tstamp[5]-1),
                                       paste("t =", tstamp[6]-1),
                                       paste("t =", tstamp[7]-1),
                                       paste("t =", tstamp[8]-1),
                                       paste("t =", tstamp[9]-1)),
                            font.label = list(size = 28),
                            vjust = 1.2,
                            align = "hv",
                            common.legend = T,
                            legend = "right"

           ),
           device = "png",
           width = 60,
           height = 50,
           units = "cm",
           dpi = 100
    )
  }

  return(plot_ls)
}
