library(parallel)
num_cores = detectCores() - 1

my_cluster = makeCluster(num_cores, type="FORK")

stopCluster(my_cluster)


z = parRapply(my_cluster, hits, function(x){countOverlaps(all_tiles[[x[1]]], snv_ranges[[x[2]]]) / snv_counts[x[2]];})

get_binned_densities <- function(deletion_filter, snv_filter, cluster){
  
  if(!is.null(deletion_filter)){
    filtered_deletion_ranges = rowRanges(deletions)[-deletion_filter]
    filtered_deletion_genotypes = deletion_genotypes[-deletion_filter,]
    filtered_deletion_carrier_mask = deletion_carrier_mask[-deletion_filter,]
  }else{
    filtered_deletion_ranges = rowRanges(deletions)
    filtered_deletion_genotypes = deletion_genotypes
    filtered_deletion_carrier_mask = deletion_carrier_mask
  }
  
  if(!is.null(snv_filter)){
    filtered_snv_ranges = snv_ranges[-snv_filter]
    filtered_snv_counts = snv_counts[-snv_filter]
    filtered_deletion_carrier_mask = filtered_deletion_carrier_mask[,-snv_filter]
  }else{
    filtered_snv_ranges = snv_ranges
    filtered_snv_counts = snv_counts
  }
  
  pre_deletion_ranges = trim(flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = T))
  post_deletion_ranges = trim(flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = F))
  
  deletion_tiles = tile(filtered_deletion_ranges, n = 10)
  pre_deletion_tiles = tile(pre_deletion_ranges, n = 10)
  post_deletion_tiles = tile(post_deletion_ranges, n = 10)
  
  all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)
  
  hits = which(filtered_deletion_carrier_mask[,] > 0, arr.ind = T)
  pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])
  
  binned_densities = parRapply(cluster, hits, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
  
  return(binned_densities)
  
} 

bins_by_cancer_project <- function(cluster){
  #dir.create(path="/home/centos/deletion_analysis_data/del_density/by_project")
  project_codes = unique(donor_meta[,dcc_project_code])
  deletion_filter = NULL
  
  sapply(project_codes, 
         function(x){print(x); get_binned_densities(NULL, which(donor_meta$dcc_project_code == x), cluster)},
         simplify = F,
         USE.NAMES = T)
  
}

project_based_bins = bins_by_cancer_project(my_cluster)