source("run.R")
gamma <- 5
res_df <- run(gamma = gamma)
write.csv(res_df, file = "out/res_gamma_" %+% gamma %+% "_" %+% timestamp %+% ".csv", row.names = FALSE)
