
circular_CN_plot = function(pretty_CN_heatmap_output,
                            label_these_genes = c("CDK14","BCL6","EZH2","HIST1H1E","REL","NOL9","TNFRSF14","TOX","TP53","RB1","TCF4","CREBBP","FOXO1","HNRNPD","BCL2","NFKBIZ","TNFAIP3","PRDM1","CD70","MYC")){
  
  bin_chroms = pretty_CN_heatmap_output$chromosome_columns

  chroms_u = unique(bin_chroms)
  xlims = matrix(ncol=2,nrow=length(chroms_u))
  rownames(xlims)=chroms_u
   chrom_bins = list()
  for(chrom in chroms_u){
    these_i = bin_chroms==chrom
    these_gains = pretty_CN_heatmap_output$cumulative_gain[these_i]
    chrom_bins[[chrom]] = names(pretty_CN_heatmap_output$cumulative_gain[these_i])
    
    xlims[chrom,1]=1
    xlims[chrom,2]=length(these_gains)
  }

  circos.initialize(chroms_u, xlim=xlims)
  
  #labels
  label_list = pretty_CN_heatmap_output$labels
  
  bin_names = names(pretty_CN_heatmap_output$cumulative_gain)
  
  bins_to_label = names(label_list)
  labels = unname(label_list)
  
  get_chrom_and_index = function(bin_name){
    chrom = str_remove(bin_name,":.+")
    this_i=which(chrom_bins[[chrom]] %in% bin_name)
    return(list(chrom,this_i))
  }
  sectors = c()
  all_labels = c()
  x = c()
  for(i in c(1:length(bins_to_label))){
    
    if(!bins_to_label[i] %in% bin_names){
      next;
    }
    label = labels[i]
    if(!label %in% label_these_genes){
      next;
    }
    ci = get_chrom_and_index(bins_to_label[i])
    
    print(paste(i,ci[[1]],ci[[2]],label))
    
    if(!is.null(ci[[2]])){
      x = c(x,ci[[2]])
      all_labels= c(all_labels,label)
      sectors = c(sectors,ci[[1]])
    }
    
  }
  circos.labels(sectors = sectors,x=x,labels=all_labels,side="outside")
  
  
  circos.trackPlotRegion(chroms_u, ylim = c(0, max(pretty_CN_heatmap_output$cumulative_gain)), 
                         track.height = 0.1)

  for(chrom in chroms_u){
    these_i = bin_chroms==chrom
    these_gains = pretty_CN_heatmap_output$cumulative_gain[these_i]
    these_losses = pretty_CN_heatmap_output$cumulative_loss[these_i]
    
    circos.lines(c(1:length(these_gains)),these_gains,type="l",area=T,
               sector.index = chrom,col ="#FF000080" )
    
  }
  #circos.trackPlotRegion(chroms_u, ylim = c(0, max(pretty_CN_heatmap_output$cumulative_loss)), 

  #                       track.height = 0.1)
  for(chrom in chroms_u){
    these_i = bin_chroms==chrom
    these_gains = pretty_CN_heatmap_output$cumulative_gain[these_i]
    these_losses = pretty_CN_heatmap_output$cumulative_loss[these_i]
    
    
    circos.lines(c(1:length(these_losses)),these_losses,type="l",area=T,
                 sector.index = chrom,col ="#0000FF80")
  }
 mid_points = c()
 for(chrom in chroms_u){
   chrom_l = length(chrom_bins[[chrom]])
   mid_points = c(mid_points,round(chrom_l/2))
   
 }
 circos.labels(sectors = chroms_u,x=mid_points,labels=chroms_u,side="inside")
}

