
#' Circular overview of copy number data across a cohort
#'
#' @param pretty_CN_heatmap_output Output from the [GAMBLR.results::pretty_CN_heatmap] call.
#' @param ideogram Logical value indicating whether to plot ideogram. Default is TRUE.
#' @param track_height Change this to increase/decrease the height of the tracks. (0.1)
#' @param labelTheseGenes Specify a vector of gene names to label in the plot
#' @param del_col Optionally specify a different colour to use for the CNV deletion track
#' @param gain_col Optionally specify a different colour to use for the CNV gain track
#' @param calculate_correlations Experimental! Calculate the correlation between CNVs between different chromosomes and link highly correlated regions
#' @param min_correlation Minimum correlation to consider when plotting links
#' @param max_neg_correlation Maximum negative value for correlations <1 to consider when plotting links
#' @param link_transparency Specify a different alpha to increase or decrease the transparency of links
#' 
#'
#' @return Nothing
#' @import GAMBLR.helpers
#' @export
#'
#' @examples
#' 
#' CN_out = pretty_CN_heatmap(cn_state_matrix=all_states_binned,
#'                            these_samples_metadata = dlbcl_genome_meta,
#'                            return_data = T,
#'                            labelTheseGenes = labelTheseGenes)
#' 
#' circular_CN_plot(CN_out)
#' 
circular_CN_plot = function(pretty_CN_heatmap_output,
                            ideogram=TRUE,
                            track_height=0.1,
                            min_correlation = 0.35,
                            max_neg_correlation = -0.06,
                            del_col="#0000FF80",
                            gain_col="#FF000080",
                            calculate_correlations = FALSE,
                            link_transparency=0.8,
                            labelTheseGenes = c("CD58","TLR2",
                                                "MCL1","CDKN2A",
                                                "TMEM30A","RHOA",
                                                "B2M","PTEN","FAS",
                                                "ETV6","GRB2","FCGR2B",
                                                "CCND3","CUX1",
                                                "MIR17HG",
                                                "TFPT","CD274","JAK2",
                                                "CDK14","BCL6","EZH2",
                                                "HIST1H1E","REL","NOL9",
                                                "TNFRSF14","TOX","TP53",
                                                "RB1","TCF4",
                                                "HNRNPD","BCL2",
                                                "NFKBIZ","TNFAIP3",
                                                "PRDM1","CD70","MYC")
                            ){
  
  label_list = pretty_CN_heatmap_output$labels
  
  bin_names = names(pretty_CN_heatmap_output$cumulative_gain)
  
  bins_to_label = names(label_list)
  labels = unname(label_list)
  
  circos.par("track.height" = track_height)
  bin_chroms = pretty_CN_heatmap_output$chromosome_columns
  orig_data = pretty_CN_heatmap_output$data
  flipped_data = orig_data
  
  
  #erase (set to diploid) any value opposing the most common event in that bin
  for(i in c(1:length(pretty_CN_heatmap_output$cumulative_gain))){
    bin_name = names(pretty_CN_heatmap_output$cumulative_gain)[i]
    if(pretty_CN_heatmap_output$cumulative_gain[i] > pretty_CN_heatmap_output$cumulative_loss[i]){
      flipped_data[orig_data[,i]<2,i]=2
    }else{
      flipped_data[orig_data[,i]>2,i]=2
    }
  }
  bin_chroms = factor(bin_chroms,levels=unique(bin_chroms))
  flipped_data = flipped_data -2
  if(calculate_correlations){
    correlations = cor(abs(flipped_data))
    link_1_chrom = c()
    link_1_start = c()
    link_1_end = c()
    link_1_value = c()
    link_2_chrom = c()
    link_2_start =c()
    link_2_end = c()
    link_2_value = c()
    link_1_col=c()
    link_1_val = c()
    for(i in c(1:length(bin_chroms))){
      for(j in c(1:length(bin_chroms))){
        if(bin_names[i] %in% bins_to_label & bin_names[j] %in% bins_to_label){
          if(bin_chroms[i]==bin_chroms[j]){
            correlations[i,j]=0
            
          }
        }
        
      }
    }
    colfun2 = colorRamp2(c(min(correlations),max(correlations)), c("blue", "red"),transparency=link_transparency)
    for(i in c(1:length(bin_chroms))){
      for(j in c(1:length(bin_chroms))){
        
        
        if(bin_names[i] %in% bins_to_label & bin_names[j] %in% bins_to_label){
          if(correlations[i,j]==0){
            next;
          }
          if(correlations[i,j] > 0 & correlations[i,j] <min_correlation){
            next;
          }
          if(correlations[i,j]< 0 & correlations[i,j] > max_neg_correlation){
            next;
          }
          region1 = names(pretty_CN_heatmap_output$cumulative_gain)[i]
          region2 = names(pretty_CN_heatmap_output$cumulative_gain)[j]
          chunks1 = region_to_chunks(region1)
          chunks2 = region_to_chunks(region2)
          if(chunks1$chromosome == chunks2$chromosome){
            print(paste("skipping",region1,region2,correlations[i,j]))
            next;
          }
          link_1_chrom = c(link_1_chrom,chunks1$chromosome)
          link_2_chrom = c(link_2_chrom,chunks2$chromosome)
          link_1_start = c(link_1_start,chunks1$start)
          link_1_end = c(link_1_end,chunks1$end)
          link_2_start = c(link_2_start,chunks2$start)
          link_2_end = c(link_2_end,chunks2$end)
          link_1_col = c(link_1_col,colfun2(correlations[i,j]))
          link_1_val = c(link_1_val,correlations[i,j])
        }
        
      }
    }
    bed1 = data.frame(chrom=link_1_chrom,start=as.integer(link_1_start),end=as.integer(link_1_end),correlation = link_1_val,colour=link_1_col)
    bed2 = data.frame(chrom=link_2_chrom,start=as.integer(link_2_start),end=as.integer(link_2_end))
    print(bed1)
    print("=-=")
    print(bed2)
  }
  
  
  

  
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
   if(ideogram){
     circos.initializeWithIdeogram(chromosome.index=unique(chroms_u))
     
   }else{
     circos.initialize(chroms_u, xlim=xlims)
     
   }
  
  #labels
  
  
  get_chrom_and_index = function(bin_name){
    chrom = str_remove(bin_name,":.+")
    this_i=which(chrom_bins[[chrom]] %in% bin_name)
    return(list(chrom,this_i))
  }
  sectors = c()
  all_labels = c()
  x = c()
  label_chrom = c()
  label_start = c()
  label_end = c()
  label_text = c()
  for(i in c(1:length(bins_to_label))){
    chunks = region_to_chunks(bins_to_label[i])
    
    if(!bins_to_label[i] %in% bin_names){
      message(paste("NOT FOUND:",bins_to_label[i],labels[i]))
      next;
    }
    label = labels[i]
    if(!label %in% labelTheseGenes){
      message(paste("SKIP:",bins_to_label[i],labels[i]))
      next;

    }
    
    ci = get_chrom_and_index(bins_to_label[i])
    
    #print(paste(bins_to_label[i],i,ci[[1]],ci[[2]],label))
    
    if(!is.null(ci[[2]])){
      x = c(x,ci[[2]])
      all_labels= c(all_labels,label)
      sectors = c(sectors,ci[[1]])
      label_chrom = c(label_chrom,chunks$chromosome)
      label_start = c(label_start,chunks$start)
      label_end = c(label_end,chunks$end)
      label_text = c(label_text,label[[1]])
    }
    
  }
  #print(label_chrom)
  #print(label_start)
  #print(label_text)
  label_bed = data.frame(chr=label_chrom,
                         start=as.integer(label_start),
                         end=as.integer(label_end),
                         value1=label_text)
  
  
  if(ideogram){
    #print(tail(label_bed))
    
  }else{
    circos.trackPlotRegion(chroms_u_nochr, ylim = c(0, max(pretty_CN_heatmap_output$cumulative_gain)), 
                           track.height = 0.1)
  }


  for(chrom in chroms_u){
    these_i = bin_chroms==chrom
    these_gains = pretty_CN_heatmap_output$cumulative_gain[these_i]
    these_losses = pretty_CN_heatmap_output$cumulative_loss[these_i]
    if(!ideogram){
      circos.lines(c(1:length(these_gains)),these_gains,type="l",area=T,
               sector.index = chrom,col =gain_col )
    }
    
  }
  if(ideogram){
    these_gains = pretty_CN_heatmap_output$cumulative_gain
    gain_df = data.frame(region=names(these_gains),value1=unname(these_gains))
    #print(head(gain_df))
    gain_df = data.frame(region=names(these_gains),value1=unname(these_gains)) %>%
      separate(region,into=c("chr","coords"),sep=":") %>%
      separate(coords,into=c("start","end"),sep="-") %>% 
      mutate(start=as.integer(start),end=as.integer(end))
    #print(head(gain_df))
    circos.genomicTrack(gain_df, panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "l",area=T, col =gain_col,  ...)
    })
  }
  #circos.trackPlotRegion(chroms_u, ylim = c(0, max(pretty_CN_heatmap_output$cumulative_loss)), 

  #                       track.height = 0.1)
  
  for(chrom in chroms_u){
    these_i = bin_chroms==chrom
    these_gains = pretty_CN_heatmap_output$cumulative_gain[these_i]
    these_losses = pretty_CN_heatmap_output$cumulative_loss[these_i]
    
    if(!ideogram){
      circos.lines(c(1:length(these_losses)),these_losses,type="l",area=T,
                   sector.index = chrom,col =del_col)
    }
    
  }
  if(ideogram){
    these_losses = pretty_CN_heatmap_output$cumulative_loss
    
    loss_df = data.frame(region=names(these_losses),value1=unname(these_losses)) %>%
      separate(region,into=c("chr","coords"),sep=":") %>%
      separate(coords,into=c("start","end"),sep="-") %>% 
      mutate(start=as.integer(start),end=as.integer(end))
    #print(head(loss_df))
    circos.genomicTrack(loss_df, panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "l",area=T, col = del_col,  ...)
    })
    if(calculate_correlations){
      bed1 = mutate(bed1,start=start-8000000,start=ifelse(start<1,1,start))
      
      bed2 = mutate(bed2,start=start-8000000,start=ifelse(start<1,1,start))
      circos.genomicLink(bed1, bed2, col = bed1$colour, border = NA)
    }
    
  }
  
 mid_points = c()
 for(chrom in chroms_u){
   chrom_l = length(chrom_bins[[chrom]])
   mid_points = c(mid_points,round(chrom_l/2))
   
 }
 if(ideogram){
   circos.genomicLabels(label_bed, labels.column = 4, side = "inside")
 }else{
   circos.labels(sectors = sectors,x=x,labels=all_labels,side="inside")
 }
 
}


#' Categorize arm-level and chromosomal CNV events
#'
#' @param pretty_CN_heatmap_output The output from running the pretty_CN_heatmap function
#'
#' @return List of data frames
#' @export
#'
#' @examples
#' 
#' cn_out = pretty_CN_heatmap(cn_state_matrix=all_states_binned,
#'     scale_by_sample = T,
#'     these_samples_metadata = all_genome_meta,
#'     metadataColumns = c("pathology","seq_type"),
#'     return_data = T)
#'     
#' aneuploidies = categorize_CN_events(cn_out)
#' 
#' 
#' 
categorize_CN_events = function(pretty_CN_heatmap_output){
  CN_mat = pretty_CN_heatmap_output$data
  labels = pretty_CN_heatmap_output$labels
  
  cytoband_df = circlize::read.cytoband()$df
  #whole chromosome, arm-level or focal
  chromosome_cols = pretty_CN_heatmap_output$chromosome_columns
  chroms_u = unique(chromosome_cols)
  if(is.null(chromosome_cols)){
    stop("problem with input")
  }
  regions_df = data.frame(regions=colnames(CN_mat)) %>% 
    mutate(name=regions) %>%
    separate(regions,into=c("chromosome","coords"),sep=":") %>%
    separate(coords,into=c("start","end"),sep="-") %>%
    mutate(start=as.integer(start),end=as.integer(end)) %>%
    mutate(arm=NA)
  chrom_arm_ranges = list()
  arm_df = expand.grid(chromosome=chroms_u,arm=c("p","q")) %>% mutate(name=paste0(chromosome,arm)) %>%
    column_to_rownames("name") %>%
    mutate(start=0,end=0) 
  for(chrom in chroms_u){
    these_col= which(chromosome_cols %in% chrom)
    p_coords1 = filter(cytoband_df,V1==chrom,grepl("p",V4)) %>% 
      pull(V2) 
    p_coords2 = filter(cytoband_df,V1==chrom,grepl("p",V4)) %>% 
      pull(V3) 
    p_coords = range(c(p_coords1,p_coords2))
    q_coords1 = filter(cytoband_df,V1==chrom,grepl("q",V4)) %>% 
      pull(V2)
    q_coords2 = filter(cytoband_df,V1==chrom,grepl("q",V4)) %>% 
      pull(V3)
    q_coords = range(c(q_coords1,q_coords2))
    
    chrom_arm_ranges[[chrom]] = list(p=p_coords,q=q_coords)
    arm_df[paste0(chrom,"p"),"start"] = p_coords[[1]] 
    arm_df[paste0(chrom,"p"),"end"] = p_coords[[2]] 
    arm_df[paste0(chrom,"q"),"start"] = q_coords[[1]] 
    arm_df[paste0(chrom,"q"),"end"] = q_coords[[2]] 
    regions_df = mutate(regions_df,
                        arm=case_when(chromosome==chrom & start < p_coords[2] & start > p_coords[1] ~"p",
                                      chromosome==chrom & start < q_coords[2] & start > q_coords[1] ~"q",
                                      TRUE ~ arm))
    
  }
  chrom_events = list()
  chrom_arm_events = list()
  skip_arms = c("chr20p","chr15p","chr14p","chr21p","chr22p","chr13p","chrXp","chrXq")
  unique_chrom = unique(chromosome_cols)
  for(chrom in unique_chrom){
    
    events = c()
    p_events = c()
    q_events = c()
    for(sample in rownames(CN_mat)){
      p_col = filter(regions_df,chromosome==chrom,arm=="p") %>% pull(name)
      q_col = filter(regions_df,chromosome==chrom,arm=="q") %>% pull(name)
      if(length(p_col)>0){
        p_mean = mean(as.numeric(CN_mat[sample,p_col])) - 2
      }else{
        p_mean = 0
      }
      
      if(length(q_col)>0){
        q_mean = mean(as.numeric(CN_mat[sample,q_col])) - 2
      }else{
        q_mean = 0
      }
      chrom_mean = mean(as.numeric(CN_mat[sample,c(p_col,q_col)])) -2
      event = NA
      if(p_mean >0.8 & q_mean > 0.8){
        event = paste0("chrom","_","gain")
        p_event = p_mean
        q_event = q_mean
        #print(paste("whole chromosome gain",chrom,sample,p_mean,q_mean,event))
 
      }else if(p_mean < -0.8 & q_mean < -0.8){
        event = paste0("chrom","_","loss")
        p_event = p_mean
        q_event = q_mean
        #print(paste("whole chromosome loss",chrom,sample,p_mean,q_mean,event))
        
        
      }else if(q_mean < -0.8 & p_mean > 0.8){
        p_event = p_mean
        q_event = q_mean
        event = paste0("iso-","qp_","lossgain")
        #print(paste("loss p, gain q",chrom,sample,p_mean,q_mean,event))
      }else if(p_mean < -0.8 & q_mean > 0.8){
        p_event = p_mean
        q_event = q_mean
        event = paste0("iso-","pq_","lossgain")
        #print(paste("loss p, gain q",chrom,sample,p_mean,q_mean,event))
      }else if(p_mean < -0.8  | q_mean < -0.8){
        if(p_mean < -0.8){
          p_event = p_mean
          q_event = 0
          event = paste0("arm-","p_","loss")
          #print(paste("loss p",chrom,sample,p_mean,q_mean,event))
        }
        if(q_mean < -0.8){
          p_event = 0
          q_event = q_mean
          event = paste0("arm-","q_","loss")
          #print(paste("loss q",chrom,sample,p_mean,q_mean,event))
        }
      }else if(p_mean > 0.8  | q_mean > 0.8){
        if(p_mean > 0.8){
          p_event = p_mean
          q_event = 0
          event = paste0("arm-","p_","gain")
          #print(paste("gain p",chrom,sample,p_mean,q_mean,event))
        }
        if(q_mean > 0.8){
          p_event = 0
          q_event = q_mean
          event = paste0("arm-","q_","gain")
          #print(paste("gain q",chrom,sample,p_mean,q_mean,event))
        }
      }
      #if(!is.na(event)){
        events = c(events,event)
      if(is.na(event)){
        p_event = 0
        q_event = 0
      }
      if(!paste0(chrom,"p") %in% skip_arms){
        p_events = c(p_events,p_event)
      }
      if(!paste0(chrom,"q") %in% skip_arms){
          q_events = c(q_events,q_event)
      }
       
      #}
      
    }
    chrom_events[[chrom]] = events
    chrom_arm_events[[paste0(chrom,"p")]] = p_events
    chrom_arm_events[[paste0(chrom,"q")]] = q_events
    
  }
  chrom_events_df= do.call("bind_cols",chrom_events) %>% as.data.frame()
  rownames(chrom_events_df)= rownames(CN_mat)
  arm_events_df = do.call("bind_cols",chrom_arm_events) %>% as.data.frame()
  rownames(arm_events_df)= rownames(CN_mat)
  # give arm-level events meaningful names and only report the most common event (drop the others) per arm
  arm_events_simplified = arm_events_df
  arm_means = colMeans(arm_events_df)
  for(i in c(1:ncol(arm_events_simplified))){
    arm = colnames(arm_events_simplified)[i]
    if(arm_means[i]<0){
      #deletion
      new_name = paste0(arm,"_LOSS")
      arm_events_simplified[,i]=ifelse(arm_events_simplified[,i]<0,1,0)
    }else if(arm_means[i]>0){
      #gain 
      new_name = paste0(arm,"_GAIN")
      arm_events_simplified[,i]=ifelse(arm_events_simplified[,i]>0,1,0)
    }
    colnames(arm_events_simplified)[i] = new_name
    
  }
  return(list(chrom_events=chrom_events_df,arm_events=arm_events_df,arm_events_simplified=arm_events_simplified))
}

