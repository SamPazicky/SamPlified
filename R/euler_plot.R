#' euler_plot
#'
#' This function generates an Euler diagram to visualize the relationships between multiple sets.
#'
#' @param sets A data frame with two columns: 'id' (set identifiers) and 'condition'.
#' @param set.shape A character string specifying the shape of the Euler diagram. Can be either "circle" or "ellipse".
#' @param plot.name A character string representing the title of the plot.
#' @param plot.label.size A numeric value indicating the font size of the labels in the plot.
#' @param quantities.type A character string indicating the type of quantities to display, typically "counts" or other quantification types.
#' @param quantities.font.type A character string specifying the font style for quantities, such as "plain", "bold", or "italic".
#' @param quantities.font.size A numeric value specifying the font size for the quantities in the plot.
#' @param color.scale A character string specifying the color scale to use. Can be "rainbow" or "distinct".
#' @param custom.scale A vector of colors for custom color scaling.
#' @param plot.euler A logical value indicating whether to generate the Euler plot. Default is TRUE.
#'
#' @import eulerr
#' @import ggplot2
#' @import ggplotify
#' @import gridExtra
#' @import gtools
#'
#' @return A list containing:
#' \item{counts}{A named list of the number of elements in the intersections of each combination of sets.}
#' \item{combinations}{A named list of the actual sets of elements for each combination of sets.}
#' \item{union_counts}{A named list of the number of elements in the union of each combination of sets.}
#' \item{unions}{A named list of the actual sets of elements in the union for each combination of sets.}
#' \item{euler_diagram}{A plot object containing the generated Euler diagram (or NULL if `plot.euler` is FALSE).}
#' \item{legend}{A grob (graphical object) representing the legend for the plot.}
#' \item{plot}{A ggplot object representing the generated Euler plot.}
#'
#' @examples
#' euler_plot(sets)
#' @export

euler_plot <- function(
    sets=NULL,                  # Data frame with columns 'id' and 'condition'.
    set.shape="circle",         # Shape of Euler diagram: "circle" or "ellipse"
    plot.name="",               # Title of the plot
    plot.label.size=0,          # Font size of plot labels
    quantities.type="counts",   # Type of quantities displayed (e.g., counts)
    quantities.font.type="plain",  # Font type for quantities: "plain", "bold", "italic", etc.
    quantities.font.size=1,     # Font size for quantities
    color.scale="rainbow",      # Color scale: "rainbow" or "distinct"
    custom.scale=NULL,          # Custom color scale (vector of colors)
    plot.euler=TRUE             # Boolean: whether to generate an Euler plot
) {

  # Validate input: quantities.font.size should be numeric
  if (!is.numeric(quantities.font.size)) {
    stop("quantities.font.size must be a number")
  }

  # Validate input: sets should be a list
  sets.prep <- prepare.data(sets,cols=c("id","condition"),no.cols=TRUE)
  if(sets.prep$message!="ok") {
    cat("Problem with sets parameter:\n")
    cat(sets.prep$message)
    if(sets.prep$status=="error") {
      stop()
    }
  }
  sets <- sets.prep$data
  rm(sets.prep)

  setlist <- sets %>%
    dplyr::select(id,condition) %>%
    unstack()


  # Validate input: each element in setlist should be a vector
  for (i in seq_along(setlist)) {
    if (!is.vector(setlist[[i]])) {  # Fix: was checking `i` instead of `setlist[[i]]`
      stop("Elements of list 'setlist' must be vectors")
    }
  }

  # Validate input: plot.label.size should be numeric
  if (!is.numeric(plot.label.size)) {
    stop("plot.label.size must be a number")
  }

  # Get unique set names sorted to ensure color consistency
  unique_set_names <- mixedsort(unique(names(setlist)))

  Pal25 <- c(  # A predefined color palette with 25 colors
    "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
    "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
    "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
    "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
    "green1", "yellow4", "yellow3", "darkorange4", "brown"
  )

  if (color.scale == "rainbow") {
    plot.scale <- rainbow(length(setlist))
  } else if (color.scale == "distinct") {
    plot.scale <- rep(Pal25, length.out=length(setlist))  # Fix: use length.out instead of ceiling(...)
  } else {
    if (is.null(custom.scale)) {
      stop("Either set plot.scale to 'rainbow' or 'distinct', or provide a custom color scale")
    } else if (!is.vector(custom.scale)) {
      stop("custom.scale must be a vector of colors")
    } else {
      plot.scale <- custom.scale
    }
  }

  # Ensure enough colors for all sets
  if (length(unique_set_names) > length(plot.scale)) {
    stop("Not enough colors in fixed palette. Add more colors.")
  }

  # Assign colors consistently
  color_map <- setNames(plot.scale[seq_along(unique_set_names)], unique_set_names)
  plot.scale <- color_map[names(setlist)]

  # Function to extract legend from a ggplot object
  get_legend <- function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) > 0) {
      return(tmp$grobs[[leg]])
    } else {
      return(NULL)
    }
  }

  # Design settings for plot appearance
  customPlot <- list(
    theme_bw(base_size = 12),
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  )

  # Validate font type for quantities
  qtable <- data.frame(name = c("plain", "bold", "italic", "bold italic"),
                       number = c(1, 2, 3, 4))
  qtype <- qtable$number[match(quantities.font.type, qtable$name)]

  if (is.na(qtype)) {
    stop("quantities.font.type must be one of: plain, bold, italic, or bold italic")
  }

  # Create a fake legend for sets
  fakedata <- data.frame(x = rep(1, length(setlist)),
                         y = rep(2, length(setlist)),
                         Set = names(setlist))
  fakedata$Set <- factor(fakedata$Set, levels = mixedsort(unique(fakedata$Set)))

  legend <- get_legend(
    ggplot(fakedata, aes(x=x, y=y, color=Set)) +
      geom_point() +
      scale_color_manual(values=plot.scale) +
      customPlot +
      theme(legend.position="bottom", legend.title=element_blank()) +
      guides(color=guide_legend(ncol=length(setlist), nrow=1))
  )

  # Generate all possible combinations of sets
  combs <- unlist(lapply(
    seq_along(names(setlist)),
    function(y) combn(names(setlist), y, simplify=FALSE)
  ), recursive=FALSE)

  # Prepare storage for combinations
  combs_numbers <- vector()
  combs_items <- list()
  union_numbers <- vector()
  union_items <- list()

  # Progress bar for intersection calculations
  pb <- txtProgressBar(min = 0, max = length(combs), style = 3)
  for (i in seq_along(combs)) {
    item_indeces <- paste(combs[[i]], collapse="&")

    combs_items[[item_indeces]] <- Reduce(intersect, setlist[combs[[i]]])
    combs_numbers[[item_indeces]] <- length(combs_items[[item_indeces]])

    union_items[[item_indeces]] <- Reduce(union, setlist[combs[[i]]])
    union_numbers[[item_indeces]] <- length(union_items[[item_indeces]])

    setTxtProgressBar(pb, i)
  }
  close(pb)  # Close the progress bar

  # Generate Euler plot if requested
  if (plot.euler) {
    cat("Producing Euler diagram... might take a few minutes...\n")

    # Set Euler diagram parameters
    eulerr_options(labels = list(fontsize = plot.label.size), padding = grid::unit(1, "cm"))

    # Compute Euler diagram
    venn_euler <- euler(combs_numbers, shape = set.shape, input = "union")

    # Generate the plot
    plot <- plot(venn_euler,
                 edges = plot.scale,
                 fills = list(fill = plot.scale, alpha = rep(0.4, length(setlist))),
                 quantities = list(type = quantities.type, cex = quantities.font.size, font = qtype),
                 main = plot.name)

    # Convert plot to ggplot object
    plot <- as.ggplot(plot)

    # Arrange Euler diagram with legend
    euler_grid <- grid.arrange(grobs = list(plot, legend), heights = c(5, 1))
  } else {
    euler_grid <- NULL
  }

  print("Plotting complete")

  # Return results
  return(list(
    counts = combs_numbers,
    combinations = combs_items,
    union_counts = union_numbers,
    unions = union_items,
    euler_diagram = euler_grid,
    legend = as.grob(legend),
    plot = plot
  ))
}
