# ---------------------------------------------------------------------------- #
#                                                                              #
#                             Themes and styling                               #
#                                                                              #
# ---------------------------------------------------------------------------- #

theme_nice <- function(size = 12) {
    baseSize = size
    largeSize = size + 2
    smallSize = size - 2
    theme_bw() +
        theme(
            ## axis
            axis.text = element_text(size = smallSize, color = "black"),
            axis.title = element_text(size = baseSize, color = "black"),
            axis.ticks=element_line(color="black"),
            axis.ticks.length=unit(.15, "cm"),
            strip.text.x = element_text(size = baseSize, color = "black"),
            ## legend
            legend.text = element_text(size = baseSize, color = "black"),
            legend.title = element_text(size = baseSize, color = "black"),
            ## title
            title = element_text(size = largeSize, color = "black"),
            plot.title = element_text(size = largeSize, color = "black"),
            plot.subtitle = element_text(size = smallSize, color = "black"),
            plot.title.position = 'plot',
            ## caption
            plot.caption = element_text(hjust = 0,
                                        face= "italic", size = smallSize),
            plot.caption.position = "plot",
            # panel background
            panel.border=element_rect(color="black", fill = NA),
            panel.background = element_blank(),
            plot.background = element_blank(),
        )
}


theme_nice_pie <- function(size = 12) {
    baseSize = size
    largeSize = size + 2
    smallSize = size - 2
    theme_void() +
        theme(## legend
            legend.text = element_text(size = baseSize),
            legend.title = element_text(size = baseSize),
            ## title
            title = element_text(size = largeSize),
            plot.title = element_text(size = largeSize),
            plot.subtitle = element_text(size = smallSize),
            plot.title.position = 'plot',
            ## caption
            plot.caption = element_text(hjust = 0,
                                        face= "italic", size = smallSize),
            plot.caption.position = "plot",
        )
}

theme_pub <- function(){
    theme_nice() + 
        theme(axis.text = element_text(size = 6),
              axis.title = element_text(size = 8),
              title = element_text(size = 10),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              plot.title = element_text(size = 10),
              strip.text.x = element_text(size = 8))
} 

theme_md <- function() {
    theme(              ## nicer text formatting
        plot.title = ggtext::element_markdown(),
        plot.caption = ggtext::element_markdown(),
        plot.subtitle = ggtext::element_markdown(),
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown())
}

colors.maps.color <- function(){
    col = c("Canonical" = "#F2D388", "Alt. putative" = "#808080",
            "Alt. used" = "#B4CDE6", "Alt. used + Sig." = "#DD5353")
    return(col)
}

colors.simple.fill <- function(){
    col = c("Alt. putative" = "#808080",
            "Alt. used" = "#B4CDE6",
            "Alt. used + Sig." = "#DD5353")
    return(col)
}

colors.simple.color <- function(){
    col = c("Alt. putative" = "#404040",
            "Alt. used" = "#204060",
            "Alt. used + Sig." = "#6a1515")
    return(col)
} 

colors.maps.enhanced.color <- function(){
    col = c("Canonical" = "#F2D388",
            "Alt. putative" = "#808080",
            "Not alt. used + Sig." = "#628E90",
            "Alt. used + Not sig." = "#B4CDE6",
            "Alt. used + Sig." = "#DD5353")
    return(col)
}


colors.sig.fill <- function(){
    col = c("FALSE" = "#808080", "TRUE" = "#DD5353")
    return(col)
}
colors.sig.color <- function(){
    col = c("FALSE" = "#404040", "TRUE" = "#6a1515")
    return(col)
}


colors.enhanced.fill <- function(){
    col = c("Event not sig." = "#808080",
            "Event sig." = "#628E90",
            "Alt. used and not sig." = "#B4CDE6",
            "Alt. used and sig" = "#DD5353")
    return(col)
}
colors.enhanced.color <- function(){
    col = c("Event not sig." = "#404040",
            "Event sig." = "#344b4c",
            "Alt. used and not sig." = "#264d73",
            "Alt. used and sig" = "#6a1515")
    return(col)
}

nice_colors = list(
    bright_full = c("#808080", "#628E90", "#B4CDE6", "#DD5353", "#F2D388"),
    dark_full = c("#404040", "#344b4c", "#264d73", "#6a1515", "#73540d"),
    
    bright_map = c("#808080", "#B4CDE6", "#DD5353", "#F2D388"),
    dark_map = c("#404040", "#264d73", "#6a1515", "#73540d"),
    
    bright_grey_red = c("#808080", "#DD5353"),
    dark_grey_red = c("#404040", "#6a1515"),
    
    bright_blue_red = c("#B4CDE6", "#DD5353"),
    dark_blue_red = c("#264d73", "#6a1515"),
    
    bright_grey_green = c("#808080", "#628E90"),
    dark_grey_green = c("#404040", "#344b4c")
)


nice_palettes = function(name, n, all_palettes = nice_colors, type = c("discrete", "continous")) {
    palette = all_palettes[[name]]
    if (missing(n)) {
        n = length(palette)
    }
    type = match.arg(type)
    out = switch(type,
                 continuous = grDevices::colorRampPalette(palette)(n),
                 discrete = palette[1:n]
    )
    structure(out, name = name, class = "palette")
}

# ggplot functions for color
scale_colour_nice_full_b <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("bright_full", type = type))
}
scale_colour_nice_full_d <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("dark_full", type = type))
}

scale_colour_nice_map_b <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("bright_map", type = type))
}
scale_colour_nice_map_d <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("dark_map", type = type))
}

scale_colour_nice_gr_b <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("bright_grey_red", type = type))
}
scale_colour_nice_gr_d <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("dark_grey_red", type = type))
}
scale_colour_nice_br_b <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("bright_blue_red", type = type))
}
scale_colour_nice_br_d <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("dark_blue_red", type = type))
}
scale_colour_nice_gg_b <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("bright_grey_green", type = type))
}
scale_colour_nice_gg_d <- function(type = "discrete") {
    ggplot2::scale_colour_manual(values = nice_palettes("dark_grey_green", type = type))
}
# fix spelling
scale_color_nice_full_b = scale_colour_nice_full_b
scale_color_nice_full_d = scale_colour_nice_full_d
scale_color_nice_map_b = scale_colour_nice_map_b
scale_color_nice_map_d = scale_colour_nice_map_d
scale_color_nice_gr_b = scale_colour_nice_gr_b
scale_color_nice_gr_d = scale_colour_nice_gr_d
scale_color_nice_br_b = scale_colour_nice_br_b
scale_color_nice_br_d = scale_colour_nice_br_d
scale_color_nice_gg_b = scale_colour_nice_gg_b
scale_color_nice_gg_d = scale_colour_nice_gg_d


# ggplot functions for fill
scale_fill_nice_full_b <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("bright_full", type = type))
}
scale_fill_nice_full_d <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("dark_full", type = type))
}

scale_fill_nice_map_b <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("bright_map", type = type))
}
scale_fill_nice_map_d <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("dark_map", type = type))
}

scale_fill_nice_gr_b <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("bright_grey_red", type = type))
}
scale_fill_nice_gr_d <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("dark_grey_red", type = type))
}
scale_fill_nice_br_b <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("bright_blue_red", type = type))
}
scale_fill_nice_br_d <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("dark_blue_red", type = type))
}
scale_fill_nice_gg_b <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("bright_grey_green", type = type))
}
scale_fill_nice_gg_d <- function(type = "discrete") {
    ggplot2::scale_fill_manual(values = nice_palettes("dark_grey_green", type = type))
}