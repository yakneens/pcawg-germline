library(ggplot2)
library(plyr)
library(data.table)

analysis_runs = read.table("/Users/siakhnin/Documents/workspace/pcawg-germline/analyses/freebayes-job-durations/freebayes_discovery_regenotype.csv", header=T, sep="\t", stringsAsFactors = F)

analysis_runs$created_date = strptime(analysis_runs$created_date, "%Y-%m-%d %H:%M:%S")
analysis_runs$run_start_date = strptime(analysis_runs$run_start_date, "%Y-%m-%d %H:%M:%S")
analysis_runs$run_end_date = strptime(analysis_runs$run_end_date, "%Y-%m-%d %H:%M:%S")
difftime(analysis_runs$run_end_date, analysis_runs$run_start_date,unit="mins")

chroms = fread("/Users/siakhnin/Documents/workspace/pcawg-germline/analyses/freebayes-job-durations/chroms.csv")
colnames(chroms)[3] = "basepairs"
chroms$basepairs = as.numeric(gsub(",", "", chroms$basepairs))

task_instances = fread("/Users/siakhnin/Documents/workspace/pcawg-germline/analyses/freebayes-job-durations/freebayes-discovery-regenotype-job-durations.tsv")
task_instances$duration = task_instances$duration/60
task_instances$chr = factor(gsub("freebayes_(.*)", "\\1", task_instances$task_id), levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
task_instances = task_instances[order(chr)]

ggplot(task_instances, aes(x=reorder(task_id, duration, FUN=median), y=duration, fill=task_id)) + geom_boxplot() +  coord_cartesian(ylim = c(0, 350)) + xlab("Chromosome") + ylab("Duration (mins)") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

sum_stats_discovery_regenotype = ddply(task_instances,~task_id,summarise,mean=mean(duration),median=median(duration),sd=sd(duration), min=min(duration), max=max(duration))
sum_stats_discovery_regenotype = sum_stats_discovery_regenotype[match(chroms$Chromosome, gsub("freebayes_(.*)", "\\1", sum_stats_discovery_regenotype$task_id)),]

cor(as.numeric(sum_stats_discovery_regenotype$mean), chroms$basepairs)

ggplot(task_instances, aes(duration, colour = task_id)) + 
  geom_density() + facet_wrap(~chr, ncol=3, scales="free") + coord_cartesian(xlim=c(0,350)) + guides(color = guide_legend(ncol=1))


