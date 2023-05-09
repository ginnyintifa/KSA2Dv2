


geth = function(x)
{
  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2] - r[1])/1.34
 t =  4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
 
 
 return(t)
}



noNA = function(x)
{
  return(length(which(!is.na(x))))
}


sum_not_na = function(x)
{
  it = which(!is.na(x))
  return(length(it))
}


get_na_df = function(v1, 
                     v2,
                     merge_df,
                     site_df)
{
  
  v1 = col_ncc_wgii_high
  v2 = col_ncc_wgii_low
  merge_df = filter_merge_r
  site_df = fread("/Users/ginny/My Drive/nonCCRCC20211020/dataFreeze/phosphosite/phosphosite_MD_allGene.tsv",
                  stringsAsFactors = F, data.table = F)
  
  na_eva = rbindlist(lapply(1:nrow(merge_df), function(x) {
   # x = 724
    tp = merge_df$name[x]
    site = gsub("_prot","", tp)
    
    t = paste(unlist(strsplit(site, split = "_"))[2:3], collapse = "_")
    w = which(site_df$geneSite == t)
    
    v1d = site_df[w, v1]
    v2d = site_df[w,v2]
    
    v1_t = sum(is.na(v1d)==F) #### number of non-missing values 
    v2_t = sum(is.na(v2d)==F)
    
    df = data.frame(name = tp, notNA1 = v1_t, notNA2 = v2_t, stringsAsFactors = F)
    
    return(df)
    
  }))
  
  
  return(na_eva)
}





get_slice = function(x,num_slice, slice_seeds)
{
  
  sx = c(1:x)
  set.seed(slice_seeds)
  s = sample(sx)
  slice_n  = floor(x/num_slice)
  left_n = x%%num_slice
  
  t = vector(mode = "list", length = num_slice)
  
  for(i in 1:num_slice)
  {
    startn = (i-1)*slice_n +1
    endn = i*slice_n 
    t[[i]] = s[startn:endn]
  }
  
  if(left_n >0)
  {
    
    for(m in 1:left_n)
      t[[m]] = c(t[[m]], s[x-m+1])
    
  }
  
  return(t)
}



balance_sample = function(d1_data,
                          d2_data,
                          s1_col_name,
                          s2_col_name,
                          nna_cutoff,
                          seeds,
                          permute_time,
                          working_dir,
                          compare_name)
{
  
   
  s1_len = length(s1_col_name)
  s2_len = length(s2_col_name)
  
  if(s1_len <= s2_len)
  {
    m = floor(s2_len/s1_len)
    
    slices = get_slice(s2_len, m)
    
    
    for(k in 1:m)
    {
      
      #k = 1
      cn = paste0("slice", k, "_", compare_name)
      
      slice_result = comparison_time_points_2d(d1_data = d1_data,
                                               d2_data = d2_data,
                                               s1_col_name = s1_col_name,
                                               s2_col_name = s2_col_name[slices[[k]]],
                                               nna_cutoff = nna_cutoff,
                                               seeds = seeds,
                                               permute_time = permute_time,
                                               working_dir = working_dir,
                                               compare_name = cn)
      
      cat(k, slice_result)    
      
    }
    
  }
  
  
  
  ############## the other scenario 
  
  if(s1_len > s2_len)
  {
    m = floor(s1_len/s2_len)
    
    
    slices = get_slice(s1_len, m)
    
    
    for(k in 1:m)
    {
      cn = paste0("slice", k, "_", compare_name)
      
      slice_result = comparison_time_points_2d(d1_data = d1_data,
                                               d2_data = d2_data,
                                               s1_col_name = s1_col_name[slices[[k]]],
                                               s2_col_name = s2_col_name,
                                               nna_cutoff = nna_cutoff,
                                               seeds = seeds,
                                               permute_time = permute_time,
                                               working_dir = working_dir,
                                               compare_name = cn)
      
      cat(k, slice_result)    
      
    }
    
  }
  
  
  
}





s1_s2_filter_unpair_2d = function(s1_d1,
                                  s1_d2,
                                  s2_d1,
                                  s2_d2,
                                  names,
                                  nna_cutoff)
{

  
  s1_d1_na = apply(s1_d1, 1, sum_not_na)
  s2_d1_na = apply(s2_d1, 1, sum_not_na)
  
  #### filter out genes with variance 0 in either
  s1_d1_zv = apply(s1_d1, 1, var, na.rm = T)
  s2_d1_zv = apply(s2_d1, 1, var, na.rm = T)
  
  
  s1_d2_na = apply(s1_d2, 1, sum_not_na)
  s2_d2_na = apply(s2_d2, 1, sum_not_na)
  
  #### filter out genes with variance 0 in either
  s1_d2_zv = apply(s1_d2, 1, var, na.rm = T)
  s2_d2_zv = apply(s2_d2, 1, var, na.rm = T)
  
  s1_d1_keep = intersect(which(s1_d1_na >= nna_cutoff), which(s1_d1_zv>0.01))
  s2_d1_keep = intersect(which(s2_d1_na >= nna_cutoff), which(s2_d1_zv>0.01))
  
  s1_d2_keep = intersect(which(s1_d2_na >= nna_cutoff), which(s1_d2_zv>0.01))
  s2_d2_keep = intersect(which(s2_d2_na >= nna_cutoff), which(s2_d2_zv>0.01))
  
  
  d1_both_keep = intersect(s1_d1_keep, s2_d1_keep)
  d2_both_keep = intersect(s1_d2_keep, s2_d2_keep)
  
  both_keep = intersect(d1_both_keep, d2_both_keep)
  
  
  s1_s2_use = list(s1_d1[both_keep,],
                   s1_d2[both_keep,],
                   s2_d1[both_keep,],
                   s2_d2[both_keep,],
                   names[both_keep])
  
  return(s1_s2_use)
  
}


density_persp = function(pdf_name,
                         x,
                         y,
                         z,
                         x_margin_dens_x,
                         x_margin_dens_y,
                         y_margin_dens_x,
                         y_margin_dens_y,
                         my_xlab,
                         my_ylab,
                         my_zlab,
                         my_main,
                         my_theta,
                         my_phi,
                         my_shade)

{
  
  
  pdf(pdf_name, useDingbats = F)
  
  
  zl= c(0, 1.5*max(z))
  
  
  col.pal = colorRampPalette(c("blue","red"))
  colors = col.pal(100)
  z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
  z.facet.range<-cut(z.facet.center, 100)
  
  
  
  trmat = persp(x,y,z, 
                theta=my_theta,
                phi = my_phi,
                box = T,
                border = NA,
                zlim = zl, 
                xlab = my_xlab,
                ylab = my_ylab,
                zlab = my_zlab,
                main = my_main,
                col = colors[z.facet.range],
                shade=my_shade, 
                lwd = 0.01,
                ticktype = "detailed")
  
  
  ### calculate a scale factor for the two margnial distributions 
  
  max_z = zl[2]
  max_xy = max(c(max(x_margin_dens_y),max(y_margin_dens_y)))
  ms = (max_z/max_xy)*0.9
  
  
  
  lines( trans3d( x_margin_dens_x, max(y), ms*x_margin_dens_y,
                  trmat), lwd=1, col=rgb(0.2,0.6,1,0.3) )
  
  
  lines( trans3d( c(0,0), c(max(y),max(y)),c(0,zl[2]),
                  trmat), lwd=0.5,lty = 2)
  
  
  lines( trans3d( min(x), y_margin_dens_x, ms*y_margin_dens_y,
                  trmat), lwd=1, col=rgb(0.2,0.6,1,0.3) )
  
  lines( trans3d( c(min(x),min(x)), c(0,0),c(0,zl[2]),
                  trmat), lwd=0.5,lty = 2)
  
  
  
  dev.off()
  
  
}



####### adjust for purity 
######## use limma to calculate the tmod after adding purity as a variable


generate_twoSample_tstat_1d_limma= function(col_annot_file, # contain purity data   
                                 s1_df1,
                                 s2_df1,
                                 rname,
                         adjust_purity)
  
{
  
  #  
  col_annot = fread(col_annot_file, stringsAsFactors = F, data.table = F)
  
  all_sample = col_annot$caseID
  p_sample = c(colnames(s1_df1), colnames(s2_df1))
  
  inter_sample = intersect(all_sample, p_sample)
  
  
  ### trim data to samples in intersection 
  
  data_table = cbind(rname,cbind(s1_df1, s2_df1))
  
  col_annot = col_annot%>%
    dplyr::filter(caseID %in% inter_sample)
  
  wi = which(p_sample %in% inter_sample)
  wi = wi+1
  
  data_table = cbind(rname,data_table[, wi])
  
  
  ### retrieve samples of interest
  
  interest_sample1 = colnames(s1_df1)
  
  s1 = col_annot%>%
    dplyr::filter(caseID %in% interest_sample1)
  
  s1_ID = intersect(s1$caseID, colnames(data_table))
  
  s1p = data.frame(caseID = s1_ID, stringsAsFactors = F)%>%
    dplyr::left_join(col_annot, by = "caseID")
  s1_purity = s1p$purity
  
  
  interest_sample2 = colnames(s2_df1)
  
  
  s2 = col_annot%>%
    dplyr::filter(caseID %in% interest_sample2)
  
  s2_ID = intersect(s2$caseID, colnames(data_table))
  
  
  
  s2p = data.frame(caseID = s2_ID, stringsAsFactors = F)%>%
    dplyr::left_join(col_annot, by = "caseID")
  s2_purity = s2p$purity
  
  
  
    
    group = as.factor(c(rep(1, length(s1_ID)),rep(2, length(s2_ID))))
    
    design = model.matrix( ~ group)
    
    if(adjust_purity == T)
    {
      purity = c(s1_purity, s2_purity)
      
      design = model.matrix( ~  group + purity)
      
    }
    
    
   contrast = suppressWarnings(makeContrasts(0-group2, levels = colnames(design)))
    
    
    colnames(contrast) = "group2"
    rownames(contrast)[1] = "(Intercept)"
    
    channels <- c(s1_ID, s2_ID)
    
    
    
    td = data_table[, channels]
    
    rownames(td) = data_table[,1]
    
    m1 = rowMeans(td[, s1_ID], na.rm = T)
    m1[is.nan(m1)] = NA
    m2 = rowMeans(td[, s2_ID], na.rm = T)
    m2[is.nan(m2)] = NA
    lfc = m1-m2
    
    de.r<- eb.fit(td, design, contrast,lfc)
    
    de.r$ID <- rownames(de.r)
    
    ###### so as long as only 1 of the samples has data it can output some signficance, but if t test can not be conducted, 
    de.r = de.r%>%
      dplyr::select(ID, log2FC, everything())
    
    
    colnames(de.r)[1:2] = c("ID","log2fc")
    
    return(de.r$t.mod)

}






eb.fit <- function(dat, design, contrast, lfc){
  
  # 
  # dat = td
  # design = design
  # contrast = contrast
  # lfc = lfc
  # 
  
  
  n <- dim(dat)[1]
  fit <- lmFit(dat, as.matrix(design))
  
  
  
  fit2 <- contrasts.fit(fit, as.matrix(contrast))
  fit.eb <- eBayes(fit2)
  log2FC <- fit.eb$coefficients[, 1]
  log2FC_ori = lfc
  
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 1]/fit.eb$sigma/fit.eb$stdev.unscaled[, 1]
  t.mod <- fit.eb$t[, 1]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 1]
  
  #q.ord <- qvalue::qvalue(p.ord)$q   ###### there is a bug with qvalue pacakge 
  
  q.ord <- p.adjust(p.ord, method = "BH")
  q.mod <- p.adjust(p.mod, method = "BH")
  
  
  results.eb <- data.frame(
    log2FC,
    log2FC_ori,
    t.ord,
    t.mod,
    p.ord,
    p.mod,
    q.ord,
    q.mod,
    df.r,
    df.0,
    s2.0,
    s2,
    s2.post
  )
  
  # results.eb <- results.eb[order(results.eb$p.mod), ]
  
  return(results.eb)
}





