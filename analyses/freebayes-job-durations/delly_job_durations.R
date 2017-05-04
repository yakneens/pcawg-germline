library(ggplot2)
library(plyr)
library(data.table)
library(xtable)


chroms = fread("/Users/siakhnin/Documents/workspace/pcawg-germline/analyses/freebayes-job-durations/chroms.csv")
colnames(chroms)[3] = "basepairs"
chroms$basepairs = as.numeric(gsub(",", "", chroms$basepairs))

del_task_instances = fread("/Users/siakhnin/Documents/workspace/pcawg-germline/analyses/freebayes-job-durations/delly_deletion_regenotype.tsv")
del_task_instances$duration = del_task_instances$duration/3600


ggplot(del_task_instances, aes(duration)) + 
  geom_histogram(bins=50) + 
  xlab("Duration (hours)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

sum_stats_discovery_regenotype = ddply(del_task_instances,~task_id,summarise,mean=mean(duration),median=median(duration),sd=sd(duration), min=min(duration), max=max(duration))


xtable(sum_stats_discovery_regenotype)

cor(as.numeric(sum_stats_discovery_regenotype$mean), chroms$basepairs)

ggplot(del_task_instances, aes(duration, colour = task_id)) + 
  geom_density()

dup_task_instances = fread("/Users/siakhnin/Documents/workspace/pcawg-germline/analyses/freebayes-job-durations/delly_dup_regenotype.tsv")
dup_task_instances$duration = dup_task_instances$duration/3600

sum_stats_discovery_regenotype = ddply(dup_task_instances,~task_id,summarise,mean=mean(duration),median=median(duration),sd=sd(duration), min=min(duration), max=max(duration))
xtable(sum_stats_discovery_regenotype)

ggplot(dup_task_instances, aes(duration)) + 
  geom_histogram(bins=50) + 
  xlab("Duration (hours)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
