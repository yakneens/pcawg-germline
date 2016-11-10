pre_deletion_ranges = trim(flank(deletion_ranges, width(deletion_ranges), start = T))
post_deletion_ranges = trim(flank(deletion_ranges, width(deletion_ranges), start = F))

deletion_tiles = tile(deletion_ranges, n = 10)
pre_deletion_tiles = tile(pre_deletion_ranges, n = 10)
post_deletion_tiles = tile(post_deletion_ranges, n = 10)

#all_tiles = cbind(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)

lapply(snv_samples[1:2], function(x) do.call(rbind, lapply(deletion_tiles, function(y) countOverlaps(y, x))))

i = 1
hits = which(deletion_carrier_mask[,] > 0, arr.ind = T)
apply(hits, 1, function(x){ do.call(rbind, lapply(deletion_tiles[x[1]], function(y) {countOverlaps(y, rowRanges(snv_samples[[x[2]]])) / snv_counts[x[2]];})); i <<- i + 1; print(i);})