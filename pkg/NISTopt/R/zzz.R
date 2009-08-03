# Ensure that the stats package is available.
.onLoad <- function(lib, pkg) {
    require("stats", character = TRUE, quietly = TRUE)
}
