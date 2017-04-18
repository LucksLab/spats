#setwd("/Users/cole/sr_mapping/SHAPE-Seq/data/run_0/spats_out")
#setwd("/Volumes/lpdata/SHAPE/run_0/no_barcode")
setwd("/Volumes/Abel/SHAPE/spats_out/")

# lib_length_file = "treated_library_length.hist"
# 
# lib_length = read.table (lib_length_file, header = TRUE)
# 
# plot(lib_length$length, 
#      lib_length$num_frags, 
#      xlab="Fragment length", 
#      ylab="Number of fragments",
#      type="h",
#      pch=19,
#      cex=0.5)
     
cexText = 0.7

plot_adducts<-function(target_name, adducts, normalize=FALSE)
{
    treated_total_mods = sum(adducts$treated_mods)
    untreated_total_mods = sum(adducts$untreated_mods)
    if (normalize)
    {
        plot_treated = adducts$treated_mods / treated_total_mods
        plot_untreated = adducts$untreated_mods / untreated_total_mods
        treated_yl = "Treated fragment ends (%)"
        untreated_yl = "Unreated fragment ends (%)"
    }
    else
    {
        plot_treated = adducts$treated_mods
        plot_untreated = adducts$untreated_mods
        treated_yl = "Treated fragment ends"
        untreated_yl = "Untreated fragment ends"
    }
    
    xl = "5' offset"
    
    def.par<- par()
    nf <- layout(matrix(c(1,2), 2, 1),heights=c(2,2))
    par(mar=c(3,5,1,1))
    
    max_mods = max(c(max(plot_treated),max(plot_untreated)))
    
    ylimits = c(0, 1.1 * max_mods)
    xlimits = c(min(adducts$five_prime_offset), 
                max(adducts$five_prime_offset))

    plot.new()
    plot.window(xlim=xlimits, ylim=ylimits)
    title(main=list(target_name), ylab=list(treated_yl, cex=cexText + 0.1))
    axis(1, cex.axis=cexText)
    axis(2, cex.axis=cexText)
    
    points(adducts$five_prime_offset, 
           plot_treated,
           type="h")
    
    plot.new()
    plot.window(xlim=xlimits, ylim=ylimits)
    title(ylab=list(untreated_yl, cex=cexText + 0.1), xlab=list(xl, cex=cexText + 0.1))
    axis(1, cex.axis=cexText)
    axis(2, cex.axis=cexText)

    points(adducts$five_prime_offset, 
           plot_untreated,
           type="h")
    
    par(def.par)
    dev.off()
 }

adducts_file = "target_WT.adducts"
adducts = read.table (adducts_file, header = TRUE) 
#adducts$five_prime_offset = adducts$five_prime_offset + 86
plot_adducts("no barcode", adducts, normalize=TRUE)

stats_file = "mapping_stats.txt"
stats = read.table (stats_file, header = TRUE)

# barplot(stats$treated_fragments)
# title(main="Total fragments (treated)")
# barplot(stats$untreated_fragments)
# title(main="Total fragments (untreated)")
