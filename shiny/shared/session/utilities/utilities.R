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
