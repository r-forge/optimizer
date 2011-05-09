#step functions


#no relaxation steps (pure classic fixed point algorithms)
purestep <- function(k) 1


#constant-decreasing steps
decrstep <- function(k, param) ifelse(k <= param, 1/2, 1/2/(k-param))

decrstep10 <- function(k) decrstep(k, 10)

decrstep20 <- function(k) decrstep(k, 20)

decrstep30 <- function(k) decrstep(k, 30)
