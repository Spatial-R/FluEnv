# Draw a graph with two y-axis, one on the left and one on the right,
# from two ggplot2 elements.
#
# Usage:
# dual_axis_graph(graph_1, graph_2, label_1, label_2)
# with graph_1 and graph_2 two ggplot2 elements, label_1 and label_2
# two strings corresponding to the y-axis label
#
# Exemples:
# see dual_axis_test.R for usage exemples
#
# Author:
# gbaquiast@quantmetry.com

library(gdata)
library(ggplot2)
library(reshape2)
#library(plyr)
library(gtable)
library(grid)

dual_axis_graph <- function(graph_1, graph_2, label_1, label_2)
{
  grid.newpage()
  
  graph_2 <- graph_2 +
    theme(panel.background = element_rect(fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # extract gtable
  grob_1 <- ggplot_gtable(ggplot_build(graph_1))
  grob_2 <- ggplot_gtable(ggplot_build(graph_2))
  
  # overlap the panel of 2nd plot on that of 1st plot
  pp <- c(subset(grob_1$layout, name == "panel", se = t:r))
  dual_grob <- gtable_add_grob(grob_1, grob_2$grobs[[which(grob_2$layout$name == "panel")]], pp$t,
                               pp$l, pp$b, pp$l)
  
  # axis tweaks
  ia <- which(grob_2$layout$name == "axis-l")
  ga <- grob_2$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  dual_grob <- gtable_add_cols(dual_grob, grob_2$widths[grob_2$layout[ia, ]$l], length(dual_grob$widths) - 1)
  dual_grob <- gtable_add_grob(dual_grob, ax, pp$t, length(dual_grob$widths) - 1, pp$b)
  
  # get colour labels
  g1 <- ggplot_build(graph_1)
  colour_left_label <- NA
  if ("fill" %in% colnames(g1$data[[1]])) {
    colour_left_label <- unique(g1$data[[1]]["fill"])[1,]
  }
  
  if (is.na(colour_left_label) & ("colour" %in% colnames(g1$data[[1]]))) {
    colour_left_label <- unique(g1$data[[1]]["colour"])[1,]
  }
  
  if (is.na(colour_left_label)){
    colour_left_label <- "black"
  }
  
  g2 <- ggplot_build(graph_2)
  colour_right_label <- NA
  if ("fill" %in% colnames(g2$data[[1]])) {
    colour_right_label <- unique(g2$data[[1]]["fill"])[1,]
  }
  
  if (is.na(colour_right_label) & ("colour" %in% colnames(g2$data[[1]]))) {
    colour_right_label <- unique(g2$data[[1]]["colour"])[1,]
  }
  
  if (is.na(colour_left_label)){
    colour_right_label <- "black"
  }
  
  left_label = textGrob(label_1, x = 0, y = 1, just = c("left", "top"),
                        gp = gpar(col =  colour_left_label, fontsize = 13))
  right_label =  textGrob(label_2, x = 1, y = 1, just = c("right", "top"),
                          gp = gpar(col =  colour_right_label, fontsize = 13))
  wrap_labels = gTree("Labs", children = gList(left_label, right_label))
  
  height <- unit(1.5, "grobheight", left_label)
  dual_grob <- gtable_add_rows(dual_grob, height, 2)
  dual_grob <- gtable_add_grob(dual_grob, wrap_labels, t=3, l=3, r=5, z = 10)
  
  height <- unit(1, "grobheight", left_label)
  dual_grob <- gtable_add_rows(dual_grob, height, 3)
  dual_grob <- gtable_add_grob(dual_grob, rectGrob(gp=gpar(col="white"), height = height),
                               t=4, l=3, r=5, z = 2)
  
  # grid.draw(dual_grob)
  return(dual_grob)
}

