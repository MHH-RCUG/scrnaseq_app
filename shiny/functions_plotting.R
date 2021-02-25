#' Out plotting style.
#'
#'@param title The plot title.
#'@param col A vector of colours to use.
#'@param fill A vector of fill colours to use.
#'@param legend_title The legend title.
#'@param legend_position The legend position.
#'@param xlab The title of the x-axis.
#'@param ylab The title of the y-axis.
#'@return None, add as theme.
AddStyle = function(title=NULL, col=NULL, fill=NULL, legend_title=NULL, legend_position=NULL, xlab=NULL, ylab=NULL) {
  list(
    theme_light() + theme(panel.border = element_blank()), 
    if (!is.null(title)) ggtitle(title), 
    if (length(col) > 0) scale_colour_manual(values=col),
    if (length(fill) > 0)  scale_fill_manual(values=fill),
    if (!is.null(legend_title)) {
      labs(color=legend_title, fill=legend_title)
    } else {
      theme(legend.title = element_blank()) 
    },
    if (!is.null(legend_position)) theme(legend.position=legend_position),
    if (!is.null(xlab)) xlab(xlab),
    if (!is.null(ylab)) ylab(ylab)
  )
}