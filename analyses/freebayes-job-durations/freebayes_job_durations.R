jobs = read.table("~/Downloads/freebayes_job_durations_genotyping", header=T, sep="\t")

ggplot(jobs, aes(x=reorder(task_id, duration, FUN=median), y=duration/60, fill=task_id)) + geom_boxplot() + coord_cartesian(ylim = c(0, 250)) + xlab("Chromosome") + ylab("Duration (mins)")
library(plyr)
sum_stats = ddply(jobs,~task_id,summarise,mean=mean(duration),meadian=median(duration),sd=sd(duration))