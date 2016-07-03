count_dt = data.table(count_df)
count_dt[, `:=`(chr=as.vector(seqnames(deletion_ranges)), width_bin=findInterval(widths, c(500,1000,2000,5000,10000,100000,250000,1000000)), carrier_count_bin=findInterval(carrier_counts, seq(0,300,25)))]
ggplot(count_dt, aes(x=-log10(carriers), y=-log10(non_carriers),colour=factor(carrier_count_bin))) + geom_point()

ggplot(melt(count_dt, measure.vars = c("carriers", "non_carriers")), aes(-log10(value), y=..scaled..,fill=variable)) + geom_density(alpha=0.6) + xlab("-log10(normalized_snv_count)") + ylab("Density")

flank_dt = data.table(flank_df)
flank_dt[, `:=`(chr=as.vector(seqnames(deletion_ranges)), width_bin=findInterval(widths, c(500,1000,2000,5000,10000,100000,250000,1000000)))]
ggplot(flank_dt, aes(x=-log10(dels), y=-log10(flanks))) + geom_point()
ggplot(melt(flank_dt, measure.vars = c("dels", "flanks")), aes(-log10(value), y=..scaled..,fill=variable)) + geom_density(alpha=0.6) + xlab("-log10(normalized_snv_count)") + ylab("Density")



