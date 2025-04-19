source("run.R")
dgp <- 2
gamma <- 2
source("dgps/sim_data_" %+% dgp %+% ".R")
res_df <- run(sim_data = sim_data, gamma = gamma)
write.csv(res_df, file = "out/dgp_" %+% dgp %+% "_" %+% gamma %+% "_" %+% timestamp %+% ".csv", row.names = FALSE)
