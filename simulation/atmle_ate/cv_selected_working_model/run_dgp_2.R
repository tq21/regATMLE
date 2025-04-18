source("run.R")
dgp_num <- 2
res_df <- run(dgp_num = dgp_num)
write.csv(res_df, file = "out/res_dgp_" %+% dgp_num %+% "_" %+% timestamp %+% ".csv", row.names = FALSE)
