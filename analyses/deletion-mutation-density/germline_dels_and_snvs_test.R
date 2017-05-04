library(VariantAnnotation)
library(data.table)
library(devtools)
devtools::load_all(pkg="~/Documents/workspace/pcawg-germline/analyses/deletion-mutation-density/pcawg.common")
library(pcawg.common)

num_bins = 10

df <- data.frame(seqnames=seqnames(chr20_del_ranges),
                 starts=start(chr20_del_ranges)-1,
                 ends=end(chr20_del_ranges),
                 names=c(rep(".", length(chr20_del_ranges))),
                 strands=strand(chr20_del_ranges))
write.table(df, file="~/Downloads/pcawg_data/del_density/germline_dels_chr_20.bed", quote=F, sep="\t", row.names=F, col.names=F)

germline_snv_chr_20 = read_germline_snvs_from_vcf("~/Downloads/pcawg_data/del_density/germline_snv_chr_20_deletions_subset.vcf.gz")
germline_snv_chr_20_ranges = rowRanges(germline_snv_chr_20)
germline_snv_genotypes = geno(germline_snv_chr_20)$GT

load("~/Downloads/pcawg_data/del_density/input_data/deletion_carrier_mask.RData")
load("~/Downloads/pcawg_data/del_density/input_data/donor_meta.RData")
load("~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData")

germline_snv_genotypes = germline_snv_genotypes[,match(donor_meta$normal_wgs_aliquot_id, colnames(germline_snv_genotypes))] 
missing_samples = which(is.na(colnames(germline_snv_genotypes)))
#missing_samples = 1836
germline_snv_genotypes = germline_snv_genotypes[,-missing_samples]
save(germline_snv_genotypes, file="~/Downloads/pcawg_data/del_density/input_data/germline_snv_genotypes_chr_20.RData")
donor_meta = donor_meta[-missing_samples,]
deletion_carrier_mask = deletion_carrier_mask[, -missing_samples]

#snv_counts = unlist(lapply(snv_ranges, length))
load(file="~/Downloads/pcawg_data/del_density/input_data/germline_snv_genotypes_chr_20.RData", verbose=T)

deletion_filter = which(as.character(seqnames(deletion_ranges)) != "20")
filtered_deletion_carrier_mask = deletion_carrier_mask
filtered_snv_counts = snv_counts

if(!is.null(deletion_filter)){
  filtered_deletion_ranges = deletion_ranges[-deletion_filter]
  filtered_deletion_carrier_mask = filtered_deletion_carrier_mask[-deletion_filter,]
}else{
  filtered_deletion_ranges = deletion_ranges
  filtered_deletion_carrier_mask = deletion_carrier_mask
}

pre_deletion_ranges = trim(flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = T))
post_deletion_ranges = trim(flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = F))

deletion_tiles = tile(filtered_deletion_ranges, n = as.integer(num_bins))
pre_deletion_tiles = tile(pre_deletion_ranges, n = as.integer(num_bins))
post_deletion_tiles = tile(post_deletion_ranges, n = as.integer(num_bins))

all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)

hits = which(filtered_deletion_carrier_mask > 0, arr.ind = T)
library(logging)
loginfo("Analyzing %s hits", dim(hits))
library(progress)

dels_by_donor = as.data.table(hits)[, list(lapply(.SD, c)), by=col]
pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(dels_by_donor)[1])

get_overlap_counts<- function(donor_index, del_list){
  tiles_by_donor = unlist(all_tiles[unlist(del_list)])
  variants_by_donor = germline_snv_chr_20_ranges[which(!germline_snv_genotypes[,donor_index] %in% c(".", "0/0"))]
  overlaps = as.data.table(findOverlaps(tiles_by_donor, variants_by_donor))
  overlap_counts = overlaps[, .N, queryHits]
  x = vector(mode="integer", length = 30 * length(dels_by_donor$V1[[donor_index]]))
  #x = matrix(0, nrow=30, ncol=length(dels_by_donor$V1[[donor_index]]))
  x[overlap_counts$queryHits] = overlap_counts$N / length(variants_by_donor)
  return(x)
}

binned_densities_germline_snv_carriers = do.call(cbind, lapply(apply(dels_by_donor, 1, function(x){pb$tick(); get_overlap_counts(x[[1]], x[[2]])}), function(x){matrix(x, nrow=30);}))


save(binned_densities_germline_snv_carriers, file="binned_densities_germline_snv_carriers.RData")
load("binned_densities_germline_snv_carriers.RData", verbose=T)

binned_densities_germline_snv_carriers_dt = as.data.table(binned_densities_germline_snv_carriers)
del_list = copy(colnames(binned_densities_germline_snv_carriers_dt))
binned_densities_germline_snv_carriers_dt[, bin_index := seq(30)]
melted_binned_densities_germline_snv_carriers_dt = melt(binned_densities_germline_snv_carriers_dt, id.vars="bin_index",variable.name="del", value.name="snv_density")
del_list = unlist(lapply(del_list, function(x) rep(x, 30)))
melted_binned_densities_germline_snv_carriers_dt$del = del_list
save(melted_binned_densities_germline_snv_carriers_dt, file="~/Downloads/temp.RData")

reduced_melted_binned_densities_germline_snv_carriers_dt = melted_binned_densities_germline_snv_carriers_dt[snv_density > 0][order(bin_index)]
reduced_melted_binned_densities_germline_snv_carriers_dt[,bin_index := as.factor(bin_index)]

bin_hits_counts = reduced_melted_binned_densities_germline_snv_carriers_dt[, .N, by=bin_index]
bin_hits_counts[bin_index %in% (10 + seq(10)), mean(N)]
bin_hits_counts[!(bin_index %in% (10 + seq(10))), mean(N)]

ggplot(bin_hits_counts, aes(x=bin_index, y=N)) + 
  geom_bar(stat="identity")  + 
  xlab("Bin") + 
  ylab("# Bins that overlap >=1 SNV") + 
  theme(text=element_text(size=20))

binned_densities_germline_snv_carriers_for_plot = data.table(rowSums(binned_densities_germline_snv_carriers), seq(30))
setnames(binned_densities_germline_snv_carriers_for_plot, c("snv_density", "bin_index"))

pre_del = binned_densities_germline_snv_carriers_for_plot[bin_index %in% seq(10), snv_density]
pre_del_t = t.test(pre_del)
in_del = binned_densities_germline_snv_carriers_for_plot[bin_index %in% (10 + seq(10)), snv_density]
in_del_t = t.test(in_del)
post_del = binned_densities_germline_snv_carriers_for_plot[bin_index %in% (20 + seq(10)), snv_density]
post_del_t = t.test(post_del)

max_density = max(binned_densities_germline_snv_carriers_for_plot$snv_density) 

ggplot(binned_densities_germline_snv_carriers_for_plot, aes(x=bin_index, y=snv_density)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")

z = lapply(seq_along(filtered_deletion_ranges), function(x){as.list(table(germline_snv_genotypes[subjectHits(findOverlaps(filtered_deletion_ranges[x],germline_snv_chr_20_ranges))]))})
z2 = rbindlist(z, fill=T)
setnames(z2, c("ref", "hom", "nocall", "het"))
z2[which(is.na(z2), arr.ind = T)] = 0

ggplot(z2, aes(x=het, y=hom)) + 
  geom_point(size=0.5) + 
  geom_abline(slope=1, intercept=0) + 
  coord_fixed(ratio=1, xlim=c(0,3000), ylim=c(0,3000))


ggplot(z2, aes(x=het, y=hom)) + 
  geom_point(size=0.5) + 
  geom_abline(slope=1, intercept=0) + 
  coord_fixed(ratio=1, xlim=c(0,50), ylim=c(0,50))



z2[order(-het)][1:20]
summary(z2$het)
summary(z2$hom)
boxplot(z2[,c("hom","het")])
boxplot(log(z2[,c("hom","het")]))
boxplot(z2[,c("hom","het")])
boxplot(z2[c("hom","het")])
z2
boxplot(log(z2))
boxplot(z2)
boxplot(z2[,hom/het])
hist(z2[,hom/het])
z2[]
z2[,hom/het]
hist(log(z2[,hom/het][1:482]))
hist(z2[,hom/het][1:482])
hist(z2[,hom/het][1:482])
hist(z2[,hom/het][1:482])
z2[,hom/het]
z2 = z2[order(-het)]
z2[order(-het)][1:100]
z2[order(-het)]
z2
hist(log(z2[,het]))
hist(log(z2[,hom]))
hist(z2[,hom])
hist(z2[,hom])
z2
hist(log(z2[hom>0,het/hom]))
hist(z2[hom>0,het/hom])
z2[hom>0,het/hom]
z2[,het/hom]
z2[,hom/het]
summary(z2[,hom/het])
hist(log(z2[,hom/het]))
hist(z2[,hom/het])
z2[,hom/het]
hist(log(na.omit(z2[,hom/het])))
z2
which(is.na(z2))
z2[which(is.na(z2), arr.ind = T)]

z2[which(is.na(z2))] = 0
z2[which(is.na(z2))]
which(is.na(z2))
which(z2 == "N/A")
z2[which(z2 == "N/A")]
hist(log(na.omit(z2[,hom/het])))
hist(na.omit(z2[,hom/het]))
na.omit(z2[,hom/het])
z2[,hom/het]
