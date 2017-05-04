my_del = "DEL00105731"
print(deletion_ranges[my_del])
print(width(deletion_ranges[my_del]))

print(donor_meta[which(deletion_carrier_mask[my_del,] > 0)]$dcc_project_code)
print(reduced_melted_binned_densities_carriers_dt[del == my_del])

donor_hits = as.vector(colnames(deletion_carrier_mask[which(deletion_carrier_mask[my_del,] > 0)])[unique(ceiling(which(melted_binned_densities_carriers_dt[del == my_del]$snv_density > 0) / 30))])
print(donor_hits)
print(donor_meta[donor_unique_id %in% donor_hits])

pre_deletion_ranges = trim(flank(deletion_ranges[my_del], width(deletion_ranges[my_del]), start = T))
post_deletion_ranges = trim(flank(deletion_ranges[my_del], width(deletion_ranges[my_del]), start = F))

deletion_tiles = tile(deletion_ranges[my_del], n = as.integer(10))
pre_deletion_tiles = tile(pre_deletion_ranges, n = as.integer(10))
post_deletion_tiles = tile(post_deletion_ranges, n = as.integer(10))

all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)

# for(my_donor in donor_hits){
#   overlaps = countOverlaps(all_tiles[[1]], filtered_snv_ranges[my_donor][[1]])
#   print(overlaps)
#   print(filtered_snv_ranges[my_donor][[1]][which(as.character(seqnames(filtered_snv_ranges[my_donor][[1]])) == seqnames(deletion_ranges[my_del]))])
# }

hit_counts = as.data.table(do.call(rbind, donor_hits %>% map(~countOverlaps(all_tiles[[1]],filtered_snv_ranges[.][[1]]))))
hit_counts[, donor_id := donor_hits]
