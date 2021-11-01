library("ggplot2")
library("mgcv")
library("dplyr")
library("permute")



##############################################
#             input parameters
##############################################

#setwd("/path/to/find_hotspot/dir")

#for control
path_data_control = "input_example/_merged_0h_Arc_Nr4a1_.csv"

#for treatment
path_data_treat = "input_example/_merged_1h_Arc_Nr4a1_.csv"

list_genes = c("Arc","Nr4a1","Egr2")

output_file_name = "output_example/hotspot_res.pdf"

shuffle_data = TRUE

shuffle_file_name = "output_example/hotspot_shuffled.pdf"

##############################################
#             functions
##############################################

make_bin = function(d_0,gene_col){
  percentage = 0.33333
  q = as.numeric(quantile(d_0[d_0[,gene_col]>1,gene_col],probs = seq(0,1,percentage)))
  bin = c(-1,1,q[2],q[3],10000)
  return(bin)
}

find_hotspot_all_rep = function(all_data,gene_col, gene_cutoff, level_cutoff, rotateData =TRUE){
  data_split = split(all_data,all_data$rep_num)
  for (i in 1:length(data_split)){
    
    print(paste0("analysis for replica number ",i))
    
    d= data_split[[i]]
    if(rotateData){
      d= rotate_data(d)
    }
    # find the hotspot for this replica
    d_hot=find_hotspot(d, gene_col, gene_cutoff, level_cutoff)
    d_hot$rep_num = i
    
    # find the bounderis of polygot for this replica
    all_poly_bnd = find_poly_bnd(d, gene_col, gene_cutoff, level_cutoff) 
    all_poly_bnd$num_replica = i
    
    
    # merge polygon bounderis for all replicas
    if (i == 1){
      all_poly_bnd_replicas = all_poly_bnd
    }else{
      all_poly_bnd_replicas = rbind(all_poly_bnd_replicas,all_poly_bnd)
    }
    
    
    # merge hotspot data for all
    if (i == 1){
      all_hotspot = d_hot
    }else{
      all_hotspot = rbind(all_hotspot,d_hot)
    }
  }
  
  # This function returns: 
  # 1. data with in/out hot spot  
  # 2. polygones for overlap
  results_list <- list(all_hotspot, all_poly_bnd_replicas)
  
}

find_hotspot = function(d, gene_col, gene_cutoff, level_cutoff){
  d_only_gene_cutoff = d [ d[,gene_col] > gene_cutoff,]
  gg <- ggplot(d_only_gene_cutoff, aes(x = center_x, y = center_y)) + 
    geom_point() + 
    geom_density2d() + 
    stat_density2d(aes(fill = ..level..), geom = "polygon")
  
  gb <- ggplot_build(gg)
  
  data = gb$data[[3]]
  
  q = quantile(data$level, probs = seq(0,1,0.1))
  cutoff = q[level_cutoff]
  
  Poly = data[as.numeric(data$level) == cutoff,]
  all_poly = split(Poly, as.character(Poly$group))
  
  
  gr = cbind(d$center_x, d$center_y)
  in_hot_spot = d
  in_hot_spot$hotSpot = "no"
  
  for (i in 1:length(all_poly)){
    myPoly = all_poly[[i]]
    
    max_x = max(d$center_x)
    max_y = max(d$center_y)
    
    bnd = cbind(myPoly$x, myPoly$y)
    inside <- in.out(bnd,gr)
    in_hot_spot$hotSpot[inside] = "yes"
  }

  return(in_hot_spot)
}

rotate_data = function(d){
  d$center_y = -1*(d$center_y)
  to_add= min(d$center_y)
  d$center_y = d$center_y+ (-1*to_add)
  return(d)
}

find_poly_bnd = function(d, gene_col, gene_cutoff, level_cutoff){
  d_only_gene_cutoff = d [ d[,gene_col] > gene_cutoff,]
  gg <- ggplot(d_only_gene_cutoff, aes(x = center_x, y = center_y)) + 
    geom_point() + 
    geom_density2d() + 
    stat_density2d(aes(fill = ..level..), geom = "polygon")
  
  gb <- ggplot_build(gg)
  
  data = gb$data[[3]]
  
  q = quantile(data$level, probs = seq(0,1,0.1))
  cutoff = q[level_cutoff]
  
  Poly = data[as.numeric(data$level) == cutoff,]
  all_poly = split(Poly, as.character(Poly$group))
  
  for (i in 1:length(all_poly)){
    myPoly = all_poly[[i]]
    bnd = cbind(myPoly$x, myPoly$y)
    bnd = cbind(bnd,i)
    bnd = as.data.frame(bnd)
    names(bnd) = c("x","y","mum_poly")
    if (i==1){
      all_poly_bnd = bnd
    } else {
      all_poly_bnd = rbind(all_poly_bnd,bnd)
    }
  }
  return(all_poly_bnd)
}

plot_poly_all_genes = function(all_poly_genes, gene, all_data, condition){

  if (gene != "all"){
    all_poly_bnd_replicas = all_poly_genes[all_poly_genes$gene == gene,]
    poly_rep = split(all_poly_bnd_replicas,all_poly_bnd_replicas$num_replica)
    x_max = max(all_data$center_x)
    y_max = max(all_data$center_y)
    plot(NULL, xlim=c(0,x_max), ylim=c(0,y_max), ylab="y", xlab="x", main = paste(gene,condition))
    par(new=TRUE)
    list_replica = vector()
    
    for (j in 1:length(poly_rep)){
      par(new=TRUE)
      rep = poly_rep[[j]]
      list_replica[j] = paste0("rep_",unique(rep$num_replica))
      poly_in_rep = split(rep, rep$mum_poly)
      for (k in 1:length(poly_in_rep)){
        poly_to_plot = poly_in_rep[[k]]
        lines(poly_to_plot$x,poly_to_plot$y, col =j)
      }
    }
    legend("topleft",legend =list_replica, fill = c(1:j))
  }
  
  else{
    poly_gene = split(all_poly_genes,all_poly_genes$gene)
    x_max = max(all_data$center_x)
    y_max = max(all_data$center_y)
    plot(NULL, xlim=c(0,x_max), ylim=c(0,y_max), ylab="y", xlab="x", main = paste(gene,condition))
    
    list_genes = vector()
    
    for (n in 1:length(poly_gene)){
      gene_data = poly_gene[[n]]
      list_genes[n] = unique(gene_data$gene)
      poly_rep = split(gene_data,gene_data$num_replica)
      for (j in 1:length(poly_rep)){
        par(new=TRUE)
        rep = poly_rep[[j]]
        poly_in_rep = split(rep, rep$mum_poly)
        for (k in 1:length(poly_in_rep)){
          poly_to_plot = poly_in_rep[[k]]
          lines(poly_to_plot$x,poly_to_plot$y, col =n)
        }
      }
    }
    legend("topleft",list_genes, fill = c(1:n))
  }
  
}

get_data_one_gene = function(data1,data_for_cutoff,current_gene, cond){
  print(cond)
  print(current_gene)
  
  
  # get cutoff level of the current gene
  gene_col = grep(current_gene,names(data_for_cutoff))
  bin = make_bin(data_for_cutoff,gene_col)
  gene_cutoff = bin[length(bin)-1] -1
  
  gene_col = grep(current_gene,names(data1))
  all_data_current_gene = find_hotspot_all_rep(data1,gene_col, gene_cutoff, level_cutoff = 9)
  
  
  gene1_poly= all_data_current_gene[[2]]
  gene1_poly$gene = current_gene
  gene1_poly$cond = cond
  #assign(paste(current_gene, cond,sep ="_"),gene1_poly)
  
  plot_poly_all_genes(gene1_poly, current_gene, data1, cond )
  
  in_out_gene1 = all_data_current_gene[[1]]
  in_out_gene1$gene = current_gene
  in_out_gene1$cond = cond
  
  return(list(gene1_poly,in_out_gene1))
}

suf_data = function(data, list_genes){
  for (gene in list_genes){
    d = data[,gene]
    new_d = d[shuffle(d)]
    data[,gene] = new_d
  }
  return(data)
}

plot_data = function(data_trea,data_con){
  data_treat$condition = "treat"
  data_con$condition = "con"
  data_to_plot = rbind(data_con,data_treat)
  data_to_plot = rotate_data(data_to_plot)
  p = ggplot(data_to_plot, aes(x=center_x, y = center_y)) +
    geom_point() + facet_grid(cols = vars(condition), rows = vars(rep_num)) +
    theme_bw(base_size = 20)
  plot(p)
}

run_analysis = function(all_data, data_for_bin, output_file) {
  pdf(output_file)
  
  # view of the data
  plot_data(all_data, data_for_bin)
  
  polygons = list()
  in_out_expression = list()
  i=1
  
  # plot polygons of the treatment condition
  for (current_gene in list_genes){
    res = get_data_one_gene(all_data,data_for_bin,current_gene, "treat")
    polygons[[i]] =res[[1]]
    in_out_expression[[i]] = res[[2]]
    i= i+1
  }
  
  
  # plot polygons of the control condition
  for (current_gene in list_genes){
    res = get_data_one_gene(data_for_bin,data_for_bin,current_gene, "con")
    polygons[[i]] = res[[1]]
    in_out_expression[[i]] = res[[2]]
    i= i+1
  }
  
  
  all_poly = bind_rows(polygons)
  all_expression =  bind_rows(in_out_expression)
  
  
  # calculate mean expression levels within polygons
  mean_exp = aggregate(all_expression$Egr2, by = list(all_expression$rep_num,
                                                      all_expression$cond,
                                                      all_expression$gene), mean)
  names(mean_exp) = c("rep_num","cond","gene","mean_exp")
  
  all_poly$id = paste(all_poly$mum_poly, all_poly$num_replica, all_poly$gene, all_poly$cond)
  
  # add expression information to polygon info
  i = 1
  for(i in 1:nrow(mean_exp)){
    ind = all_poly$cond == mean_exp$cond[i] & all_poly$num_replica == mean_exp$rep_num[i] & all_poly$gene == mean_exp$gene[i]
    all_poly$mean_expression[ind]= mean_exp$mean_exp[i]
  }
  
  # Plot global view of hotpots: polygons and expression levels
  p <- ggplot(all_poly, aes(x = x, y = y)) +
    geom_polygon(aes(group = id , fill = gene, alpha= mean_expression)) +
    scale_alpha("mean_expression", range = c(0.1,0.7)) +
    facet_grid(cols = vars(factor(cond)), rows = vars(gene)) +
    theme_bw() +
    theme(axis.text=element_text(size=15),
          legend.text=element_text(size=15),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    ylab("") + xlab("")
  
  plot(p)
  
  
  dev.off()
}



##############################################
#             analysis
##############################################


# data of the treatment condition
data_treat = read.csv(path_data_treat)


# data of control condition
data_con = read.csv(path_data_control)



# run analysis to find hotspots in data
run_analysis(data_treat, data_con, output_file_name)

# run analysis on shuffled data
if (shuffle_data == TRUE){
  data_treat = suf_data(data_treat, list_genes)
  data_con = suf_data(data_con, list_genes)
}

run_analysis(data_treat, data_con, shuffle_file_name)



