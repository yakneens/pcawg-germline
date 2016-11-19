get_binned_densities <- function(deletion_filter, snv_filter){
  
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
  
  binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
  
  return(binned_densities)
  
} 
  #z4 = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], snv_ranges[[x[2]]][seqnames(snv_ranges[[x[2]]]) == as.character(seqnames(filtered_deletion_ranges[x[1]]))])  / snv_counts[x[2]];})
  #snv_intervals = lapply(snv_ranges, function(x){ temp <- as.data.table(ranges(x[seqnames(x) == "21"]))[,.(start, end)]; if(!isEmpty(temp)){Intervals(temp)}else{Intervals()}})
  #deletion_intervals = Intervals(as.data.table(ranges(filtered_deletion_ranges))[,.(start, end)])
  #z4 = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles_intervals[[x[1]]], snv_intervals[[x[2]]]) / snv_counts[x[2]];})

# res = list()
# system.time(for(i in 1:10){
#   carriers = GRangesList(snv_ranges[which(filtered_deletion_carrier_mask[i,] == 1)])
#   over = findOverlaps(all_tiles[[i]], carriers)
#   over_t = as.table(over)
#   over_t[queryHits(over)] = over_t[queryHits(over)] * (1 / snv_counts[subjectHits(over)])
#   res[[i]] = over_t
# })
# 
# snv_ranges_glist = GRangesList(snv_ranges)
# system.time(for(i in 1:10){
#   #carriers = GRangesList(snv_ranges[which(filtered_deletion_carrier_mask[i,] == 1)])
#   over = findOverlaps(all_tiles[[i]], snv_ranges_glist[which(filtered_deletion_carrier_mask[i,] == 1)])
#   over_t = as.table(over)
#   over_t[queryHits(over)] = over_t[queryHits(over)] * (1 / snv_counts[subjectHits(over)])
#   res[[i]] = over_t
# })
# 
# get_binned_densities2(hit_vec){}
# 
# apply(filtered_deletion_carrier_mask, 1, function(x){x})
