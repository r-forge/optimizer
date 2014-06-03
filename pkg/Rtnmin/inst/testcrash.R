source("crash.R")
# testcrash

low <- c(0, 1, 2)
up  <- c(4, 5, 6)
x <- c(2,2,2)
print(crash(x, low, up))
x <- c(5, 3, 3)
print(crash(x, low, up))
low <- c(4, 1, 2)
print(crash(x, low, up))
