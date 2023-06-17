newggplot.ggdend <- function (data, segments = TRUE, labels = TRUE, nodes = TRUE, 
                              horiz = FALSE, theme = theme_dendro(), offset_labels = 0, ...) {
  data <- prepare.ggdend(data)
  #angle <- ifelse(horiz, 0, 90)
  #hjust <- ifelse(horiz, 0, 1)
  p <- ggplot()
  if (segments) {
    p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend", colour = "col", linetype = "lty", size = "lwd"), lineend = "square") + 
      guides(linetype = FALSE, col = FALSE) + scale_colour_identity() + 
      scale_size_identity() + scale_linetype_identity()
  }
  if (nodes) {
    p <- p + geom_point(data = data$nodes, aes_string(x = "x", y = "y", colour = "col", shape = "pch", size = "cex")) + 
      guides(shape = FALSE, col = FALSE, size = FALSE) + 
      scale_shape_identity()
  }
  if (labels) {
    data$labels$cex <- 5 * data$labels$cex
    data$labels$y <- data$labels$y + offset_labels
    p <- p + geom_text(data = data$labels, aes_string(x = "x", y = "y", label = "label", colour = "col", size = "cex", angle = "angle", hjust = "hjust", vjust = "vjust"))#edited
  }
  if (horiz) {
    p <- p + coord_flip() + scale_y_reverse(expand = c(0.2, 0))
  }
  if (!is.null(theme)) {
    p <- p + theme
  }
  p
}