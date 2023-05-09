




comparison_time_points_2d_limma_bandwidth = function(d1_data,
                                                     d2_data,
                                                     s1_col_name,
                                                     s2_col_name,
                                                     nna_cutoff,
                                                     permute_time,
                                                     working_dir,
                                                     compare_name,
                                                     col_annot_file,
                                                     bandF = 1,
                                                     is_purity_adjusted)

{
  
  # d1_data = d1_data
  # d2_data = d2_data
  # s1_col_name = s1_col_name
  # s2_col_name = s2_col_name[slices[[k]]]
  # nna_cutoff = nna_cutoff
  # permute_time = permute_time
  # working_dir = working_dir
  # compare_name = cn
  # col_annot_file = annot_file
  # bandF = bandwidth_factor
  # is_purity_adjusted = adjust_purity_flag
  # 
  
  
  d1_data = as.data.frame(d1_data)
  d2_data = as.data.frame(d2_data)
  
  
  # s1_col_name = intersect(colnames(d1_data))
  
  
  s1_d1 = d1_data[,s1_col_name]
  s1_d2 = d2_data[,s1_col_name]
  
  s2_d1 = d1_data[,s2_col_name]
  s2_d2 = d2_data[,s2_col_name]
  
  #### filter data
  
  s1_s2_use = s1_s2_filter_unpair_2d(s1_d1 = s1_d1,
                                     s1_d2 = s1_d2,
                                     s2_d1 = s2_d1,
                                     s2_d2 = s2_d2,
                                     names = d1_data[,1],
                                     nna_cutoff = nna_cutoff)
  
  
  
  cat("data filtered.", "\n")
  #### generate Z
  s1_d1_use = s1_s2_use[[1]]
  s1_d2_use = s1_s2_use[[2]]
  s2_d1_use = s1_s2_use[[3]]
  s2_d2_use = s1_s2_use[[4]]
  
  
  

  
  
  
  s1_s2_Z_use_d1 = generate_twoSample_tstat_1d_limma(s1_df1 = s1_d1_use,
                                                     s2_df1 = s2_d1_use,
                                                     rname = s1_s2_use[[5]],
                                                     col_annot_file = col_annot_file,
                                                     adjust_purity = is_purity_adjusted) 
  
  
  s1_s2_Z_use_d2 = generate_twoSample_tstat_1d_limma(s1_df1 = s1_d2_use,
                                                     s2_df1 = s2_d2_use,
                                                     rname = s1_s2_use[[5]],
                                                     col_annot_file = col_annot_file,
                                                     adjust_purity = is_purity_adjusted) 
  

  s1_mean_d1 = apply(s1_d1_use, 1, mean, na.rm = T)
  s2_mean_d1 = apply(s2_d1_use, 1, mean, na.rm = T)
  
  s1_s2_lfc_d1 = s1_mean_d1 - s2_mean_d1
  
  s1_mean_d2 = apply(s1_d2_use, 1, mean, na.rm = T)
  s2_mean_d2 = apply(s2_d2_use, 1, mean, na.rm = T)
  
  s1_s2_lfc_d2 = s1_mean_d2 - s2_mean_d2
  
  
  sn = min(c(ncol(s1_d1_use), ncol(s2_d1_use)))
  cn = sn/2
  
  #### 0508  check this part again 
  
  null_permute  = lapply(1:permute_time,function(x) {

    sset = (x-1)*4 + c(1:4)
    
    
    set.seed(sset[1])
    s1_d1_sel = sample(c(1:ncol(s1_d1_use)), sn)
    set.seed(sset[2])
    s2_d1_sel = sample(c(1:ncol(s2_d1_use)), sn)
    
    
    if(x%%2 ==0)
    {
      set.seed(sset[3])
      s1_change = sample(s1_d1_sel,floor(cn))
      set.seed(sset[4])
      s2_change = sample(s2_d1_sel,floor(cn))
      
      
    }else{
      set.seed(sset[3])
      s1_change = sample(s1_d1_sel,ceiling(cn))
      set.seed(sset[4])
      s2_change = sample(s2_d1_sel,ceiling(cn))
      
    }
    
    this_list = list(setdiff(s1_d1_sel,s1_change),
                     s1_change, 
                     setdiff(s2_d1_sel,s2_change),
                     s2_change)
    
    return(this_list)
    
    
  })
  
  

  
  
  s1_s2_z_null_use = lapply(1:permute_time, function(x)
  {
    # x = 1
    this_assign = null_permute[[x]]
    
    #first = which(this_assign ==1)
    
    this_c1_d1 = cbind(s1_d1_use[,this_assign[[1]]], s2_d1_use[,this_assign[[4]]])
    this_c2_d1 = cbind(s1_d1_use[,this_assign[[2]]], s2_d1_use[,this_assign[[3]]])
    
    this_c1_d2 = cbind(s1_d2_use[,this_assign[[1]]], s2_d2_use[,this_assign[[4]]])
    this_c2_d2 = cbind(s1_d2_use[,this_assign[[2]]], s2_d2_use[,this_assign[[3]]])
    
    permute_use = s1_s2_filter_unpair_2d(s1_d1 = this_c1_d1,
                                         s2_d1 = this_c2_d1,
                                         s1_d2 = this_c1_d2,
                                         s2_d2 = this_c2_d2,
                                         names = s1_s2_use[[5]],
                                         nna_cutoff = nna_cutoff)
    
    
    
    

    
    this_z_d1 = generate_twoSample_tstat_1d_limma(s1_df1 = permute_use[[1]],
                                                  s2_df1 = permute_use[[3]],
                                                  rname = permute_use[[5]],
                                                  col_annot_file = col_annot_file,
                                                  adjust_purity = is_purity_adjusted) 
    
    
    this_z_d2 = generate_twoSample_tstat_1d_limma(s1_df1 = permute_use[[2]],
                                                  s2_df1 = permute_use[[4]],
                                                  rname = permute_use[[5]],
                                                  col_annot_file = col_annot_file,
                                                  adjust_purity = is_purity_adjusted) 
    
    
    
    
    this_z = cbind(this_z_d1, this_z_d2)
    
     if(x%%100 == 0)
    cat(x, "\n")
    
    return(this_z)
    
  })
  
  
  
  s1_s2_z_null_use_cat =  Reduce(rbind, s1_s2_z_null_use)
  
  
  
  
  xlim = c(min(c(min(s1_s2_Z_use_d1), min(s1_s2_z_null_use_cat[,1]))),
           max(c(max(s1_s2_Z_use_d1), max(s1_s2_z_null_use_cat[,1]))))
  
  ylim = c(min(c(min(s1_s2_Z_use_d2), min(s1_s2_z_null_use_cat[,2]))),
           max(c(max(s1_s2_Z_use_d2), max(s1_s2_z_null_use_cat[,2]))))
  
  
  h1 = geth(s1_s2_Z_use_d1)
  h2 = geth(s1_s2_Z_use_d2)
  
  
  
  
  
  Z_dens = kde2d(s1_s2_Z_use_d1,
                 s1_s2_Z_use_d2,
                 n = 1000,
                 h = c(bandF*h1, bandF*h2),  
                 lims = c(xlim, ylim))
  
  Z_x_dens = density(s1_s2_Z_use_d1)
  Z_y_dens = density(s1_s2_Z_use_d2)
  
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_Z_dens_function.pdf"),
                 x = Z_dens$x,
                 y = Z_dens$y,
                 z = Z_dens$z,
                 x_margin_dens_x = Z_x_dens$x,
                 x_margin_dens_y = Z_x_dens$y,
                 y_margin_dens_x = Z_y_dens$x,
                 y_margin_dens_y = Z_y_dens$y,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "Z",
                 my_theta = 30,
                 my_phi = 15,
                 my_shade = 0.5)
  
  

  set.seed(123)
  sel_row = sample(c(1:nrow(s1_s2_z_null_use_cat)), 10000)
  
  
  
  h1 = geth(s1_s2_z_null_use_cat[sel_row,1])
  h2 = geth(s1_s2_z_null_use_cat[sel_row,2])
  
  z_null_dens = kde2d(s1_s2_z_null_use_cat[sel_row,1],
                      s1_s2_z_null_use_cat[sel_row,2],
                      n = 1000, 
                      h = c(h1*(max((bandF-3.5),1)), h2*(max((bandF-3.5),1))), ### understand why 3.5 here  
                      lims = c(xlim, ylim)) 
  
  z_null_x_dens =  density(s1_s2_z_null_use_cat[sel_row,1])
  z_null_y_dens =  density(s1_s2_z_null_use_cat[sel_row,2])
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_z_null_dens_function.pdf"),
                 x = z_null_dens$x,
                 y = z_null_dens$y,
                 z = z_null_dens$z,
                 x_margin_dens_x = z_null_x_dens$x,
                 x_margin_dens_y = z_null_x_dens$y,
                 y_margin_dens_x= z_null_y_dens$x,
                 y_margin_dens_y = z_null_y_dens$y,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "z",
                 my_theta = 30,
                 my_phi = 15,
                 my_shade = 0.5)
  
  
  
  #############################################
  #############################################
  #############################################
  
  
  
  
  
  
  Z_zero_which_x = which(abs(Z_dens$x-0) == min(abs(Z_dens$x-0)))[1]
  Z_zero_which_y = which(abs(Z_dens$y-0) == min(abs(Z_dens$y-0)))[1]
  
  if(Z_dens$x[Z_zero_which_x]>0)
  {
    Z_zero_x2 = Z_zero_which_x-1
  }else{
    Z_zero_x2 = Z_zero_which_x+1
  }
  
  
  if(Z_dens$y[Z_zero_which_y]>0)
  {
    Z_zero_y2 = Z_zero_which_y-1
  }else{
    Z_zero_y2 = Z_zero_which_y+1
  }
  
  Z_zero_dens1 = Z_dens$z[Z_zero_which_x, Z_zero_which_y]
  Z_zero_dens2 = Z_dens$z[Z_zero_which_x, Z_zero_y2]
  Z_zero_dens3 = Z_dens$z[Z_zero_x2, Z_zero_which_y]
  Z_zero_dens4 = Z_dens$z[Z_zero_x2, Z_zero_y2]
  
  
  z_null_zero_dens1 = z_null_dens$z[Z_zero_which_x, Z_zero_which_y]
  z_null_zero_dens2 = z_null_dens$z[Z_zero_which_x, Z_zero_y2]
  z_null_zero_dens3 = z_null_dens$z[Z_zero_x2, Z_zero_which_y]
  z_null_zero_dens4 = z_null_dens$z[Z_zero_x2, Z_zero_y2]
  
  s1_s2_p0_use_vector = c(Z_zero_dens1, Z_zero_dens2, Z_zero_dens3, Z_zero_dens4)/c(z_null_zero_dens1, z_null_zero_dens2, z_null_zero_dens3, z_null_zero_dens4)
  s1_s2_p0_use = mean(s1_s2_p0_use_vector)
  
  
  
  f02f = rep(0, length(s1_s2_Z_use_d1))
  
  for(i in 1:length(s1_s2_Z_use_d1))
  {
    this_value_x = s1_s2_Z_use_d1[i]
    this_value_y = s1_s2_Z_use_d2[i]
    
    dis_this_x = abs(this_value_x-Z_dens$x)
    nearest_x = which(dis_this_x == min(dis_this_x))[1]
    dis_this_y = abs(this_value_y-Z_dens$y)
    nearest_y = which(dis_this_y == min(dis_this_y))[1]
    
    this_Z_f = Z_dens$z[nearest_x, nearest_y]
    
    this_Z_f0 = z_null_dens$z[nearest_x, nearest_y]
    
    
    f02f[i] = this_Z_f0/this_Z_f
    
  }
  
  Z_p0 = s1_s2_p0_use*f02f
  Z_p1 = 1-Z_p0
  
  p1_corrected = Z_p1
  p1_corrected[which(Z_p1<0)]=0
  lfdr = 1-p1_corrected
  
  
  
  
  posterior_pdf_name = paste0(working_dir,compare_name, "_2d_Z_posterior.pdf")
  
  ### when plot remove the outlier ones 
  pdf(posterior_pdf_name,useDingbats = F)
  
  get_col = rep(rgb(0.1,0,0,0.5), length(p1_corrected))
  get_col[which(p1_corrected>0.9)] = rgb(0.5,0,0,0.5)
  get_col[which(p1_corrected>0.99)] = rgb(1,0.5,0,0.5)
  
  
  
  scatterplot3d(x = s1_s2_Z_use_d1,
                y = s1_s2_Z_use_d2,
                z = p1_corrected,
                # type = "h",
                color = get_col,
                # box = F,
                pch = 16,
                xlab = "kinase",
                ylab = "substrate",
                zlab = "posterior")
  
  dev.off()
  
  
  
  
  result_df = data.frame(name =s1_s2_use[[5]],
                         Z_kinase = s1_s2_Z_use_d1,
                         Z_substrate = s1_s2_Z_use_d2,
                         fc_kinase = s1_s2_lfc_d1,
                         fc_substrate = s1_s2_lfc_d2,
                         p1 = Z_p1,
                         p1_corrected = p1_corrected,
                         fdr = 1-p1_corrected,
                         abs_fc_kinase = abs(s1_s2_lfc_d1),
                         abs_fc_substrate = abs(s1_s2_lfc_d2),
                         stringsAsFactors = F)%>%
    dplyr::arrange(desc(p1))
  
  result_df_name = paste0(working_dir, compare_name,"_Z_result.tsv")
  
  write.table(result_df, result_df_name,
              quote = F, row.names = F, sep = "\t")
  
  f1_dens  = Z_dens$z
  
  for(i in 1:length(Z_dens$x))
  {
    for(j in 1:length(Z_dens$y))
    {
      
      this_dens = Z_dens$z[i,j]
      
      this_null_dens = z_null_dens$z[i, j]
      
      f1_dens[i,j] = (this_dens-s1_s2_p0_use*this_null_dens)/(1-s1_s2_p0_use)
    }
    
    
  }
  
  f1_dens[which(f1_dens<0)] = 0
  
  
  f1_x_marginal = rowSums(f1_dens)/sum( rowSums(f1_dens))
  f1_y_marginal = colSums(f1_dens)/sum(colSums(f1_dens))
  
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_f1_function.pdf"),
                 x = Z_dens$x,
                 y = Z_dens$y,
                 z = f1_dens,
                 x_margin_dens_x = Z_dens$x,
                 x_margin_dens_y = f1_x_marginal,
                 y_margin_dens_x= Z_dens$y,
                 y_margin_dens_y = f1_y_marginal,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "f1",
                 my_theta = 30,
                 my_phi = 15,
                 my_shade = 0.5)
  
  
  
  cat("posterior calculated.", "\n")
  
  
  return(s1_s2_p0_use)
  
  
}











