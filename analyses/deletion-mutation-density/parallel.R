library(parallel)
num_cores = detectCores() - 1

my_cluster = makeCluster(num_cores, type="FORK")

stopCluster(my_cluster)


z = parRapply(my_cluster, hits, function(x){countOverlaps(all_tiles[[x[1]]], snv_ranges[[x[2]]]) / snv_counts[x[2]];})