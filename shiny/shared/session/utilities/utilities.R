# add transparency to a vector of colors 
addAlphaToColors <- Vectorize(function(color, alpha = 0.25) {
    x <- col2rgb(color)
    rgb(x[1], x[2], x[3], max = 255, alpha = alpha * 255)
})

# convert a matrix of color names/HEX to a imager::cimg object
matrixToCImg <- function(m, width, height) {
    (col2rgb(m) / 255) %>%
    t() %>% 
    array(dim = c(width, height, 1, 3)) %>% # input and output width and height are not necessarily the same
    imager::as.cimg()
}

# common plot settings
paInteractiveFrame <- list(
    Plot_Frame = list(
        Width_Inches = 4,
        Height_Inches = 3,
        Font_Size = 9
    )
)
paRainbow <- function(n, alpha = NULL) {# eliminates magenta colors that are visually too close to red
    colors <- rainbow(n / 0.9)[1:n]
    if(!is.null(alpha)) colors <- addAlphaToColors(colors, alpha)
    colors
}
titledMar <- c(4, 4, 1.7, 0.2) + 0.1
addStageXAxis <- function(stages, ylim, side = 1) {
    axis(side = side, labels = FALSE)
    text(1:length(stages), ylim[1] - diff(ylim) / 8, stages, xpd = TRUE, srt = 45, adj = 1)
}
