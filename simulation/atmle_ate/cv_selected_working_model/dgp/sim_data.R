library(purrr)
`%+%` <- function(a, b) paste0(a, b)

# source DGP scripts
dgp_num <- seq(2)
f_names <- "dgp/sim_data_" %+% dgp_num %+% ".R"
walk(f_names, source)
