source("run.R")
dgp <- 1
source("dgps/sim_data_" %+% dgp %+% ".R")
res_df <- run(sim_data = sim_data)
write.csv(res_df, file = "out/dgp_" %+% dgp %+% "_" %+% timestamp %+% ".csv", row.names = FALSE)
