#### Filtering donors by SNV count
filter_curves <- function(data, eval_seq, filter_col, measure_col, y_label){
  result = as.data.table(
              do.call(rbind,
                      lapply(eval_seq, function(i){
                        c(sum(data[get(filter_col) > i, get(measure_col)] == 0),
                          sum(data[get(filter_col) > i, get(measure_col)] != 0),
                          sum(data[get(filter_col) < i, get(measure_col)] == 0),
                          sum(data[get(filter_col) < i, get(measure_col)] != 0),
                          i)
                      })
              )
            )

  setnames(result, c(paste("mean SNV load in ", measure_col," = 0 retained"),
                        paste("mean SNV load in ",measure_col," > 0 retained"),
                        paste("mean SNV load in ",measure_col," = 0 filtered"),
                        paste("mean SNV load in ",measure_col," > 0 filtered"),
                        filter_col))
  result = melt(result[get(paste("mean SNV load in ", measure_col," = 0 filtered")) > 0 & get(paste("mean SNV load in ", measure_col," = 0 retained")) > 0], id.vars = filter_col)
  result[, value := log10(value)][!is.finite(value), value := 0]
  localenv <- environment()
  print(
    ggplot(data=result, aes(x=log10(get(filter_col)), y=value, group=variable, col=variable)) +
      stat_smooth(se=F, span=0.1) +
      scale_x_continuous() +
      theme(legend.position="bottom") +
      xlab(paste("log10(", filter_col, ")")) +
      ylab(y_label) +
      ylim(0,4)
  )

  return(result)
}

filter_curves(results_by_donor, seq(10, 6000,20), "snv_counts", "dels", "log10(# Donors)")
filter_curves(results_by_donor, seq(10, 6000,20), "snv_counts", "flanks", "log10(# Donors)")
filter_curves(results_by_del, seq(2,3000), "carrier_counts", "carriers", "log10(# Deletions)")
filter_curves(results_by_del, seq(2,3000), "carrier_counts", "non_carriers", "log10(# Deletions)")
filter_curves(results_by_del, seq(2,3000), "carrier_counts", "flanks", "log10(# Deletions)")
