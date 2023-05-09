


#' kinase_substrate mapping function
#' 
#' map the knase and substrate information
#' @param protData_filename absolute directory to the file of the proteome data, in the format of the tidied up FragPipe TMT-I output
#' @param psiteData_filename  absolute directory to the file of the single site level phospho data, in the format of the tidied up FragPipe TMT-I output, normalized or not normalized
#' @param ksNetwork_filename absolute directory to the kinase substrate network file, as an example, merge_ks_omni.tsv. 
#' @param ks_outputName the file name of the output file of this function 
#' @import dplyr data.table magrittr 
#' @export
#' 
ks_map = function(protData_filename,
                  psiteData_filename,
                  ksNetwork_filename,
                  ks_outputName)
  
  
{
  
  
  
  sel_prot= fread(protData_filename,
                  stringsAsFactors = F,
                  data.table = F)
  
  p_psite = fread(psiteData_filename,
                  stringsAsFactors = F,
                  data.table = F)
  
  ks_network = fread(ksNetwork_filename,
                     stringsAsFactors = F,
                     data.table = F)
  
  
  ks_network = ks_network%>%
    unique()%>%
    dplyr::mutate(gs = paste0(substrate, "_", substrate_site))%>%
    dplyr::filter(gs %in% p_psite$Gene_site)
  
  
  p_psite = p_psite%>%
    dplyr::filter(Gene_site %in% ks_network$gs)%>%
    dplyr::select(-SequenceWindow)%>%
    dplyr::group_by(Gene_site)%>%
    dplyr::summarise_if(is.numeric, mean, na.rm = T)%>%
    ungroup()%>%
    dplyr::mutate_all(function(x) ifelse(is.nan(x), NA, x))%>%
    as.data.frame()
  
  
  
  all_ks = unique(ks_network$kinase)
  
  inter_sample = intersect(colnames(sel_prot), colnames(p_psite))
  
  data_ks = rbindlist(lapply(1:length(all_ks), function(x) {
    
    if(x%%100 ==0)
      cat(x, "\n")
    kin = all_ks[x]
    
    kin_prot = sel_prot%>%
      dplyr::filter(Gene_name == kin)
    
    if(nrow(kin_prot)>0)
    {
      sel_ks = ks_network%>%
        dplyr::filter(kinase == kin)
      
      subsite_psite = p_psite%>%
        dplyr::filter(Gene_site %in% sel_ks$gs)
      
      subsite_psite_data = as.matrix(subsite_psite[,inter_sample])
      kin_prot_data = as.matrix(kin_prot[rep(1,nrow(subsite_psite_data)),inter_sample])
      
      fuse_mat = matrix(NA, nrow = 2*nrow(subsite_psite_data), ncol = ncol(subsite_psite_data))
      odd_col = c(1:nrow(subsite_psite_data))*2-1
      even_col = c(1:nrow(subsite_psite_data))*2
      
      fuse_mat[odd_col,] = kin_prot_data
      fuse_mat[even_col,] = subsite_psite_data
      
      name = rep("", nrow(fuse_mat))
      name[odd_col] = paste0(kin,"_", subsite_psite$Gene_site, "_", "prot")
      name[even_col] = paste0(kin,"_", subsite_psite$Gene_site, "_", "psite")
      
      fuse_df = data.frame(name, fuse_mat, stringsAsFactors = F)
      colnames(fuse_df) = c("name", inter_sample)
      
      return(fuse_df)
      
      
    }
    
  }))
  
  
  write.table(data_ks, 
              ks_outputName,
              quote = F, row.names = F, sep = "\t")
  
  
  npair = nrow(data_ks)/2
  p_row  = c(1:npair)*2-1
  s_row = c(1:npair)*2
  
  kprot_data = data_ks[p_row,]
  subsite_data = data_ks[s_row,]
  
  pair_data_list = list(kprot_data, subsite_data)
  
  return(pair_data_list)
  
}





#' main EB function
#' 
#' calculate posterior probability of each relationship being significant
#' @param s1_col_name names of the columns in the first condition 
#' @param s2_col_name names of the columns in the second condition
#' @param d1_data dataframe for the data of the first dimension, rows are genes, columns are patients in 2 conditions
#' @param d2_data dataframe for the data of the second dimension, rows are genes, columns are patients in 2 conditions
#' @param nna_cutoff genes with values in at least the number of samples are included in comparison
#' @param permute_time number of permutations 
#' @param working_dir the directory output files will be deposited in 
#' @param compare_name a label for the anlaysis 
#' @param bandwidth_factor a numeric value for smoothing the empirical distribution, the larger the smoother, default to 1
#' @param adjust_purity_flag a boolean value if set to True tumor purity is adjusted and the annot_file must be provided 
#' @param annot_file absolute directory to the sample annotation file, column name must be purity, ID column name must be caseID
#' @import dplyr data.table magrittr splines stringr MASS scatterplot3d limma
#' @keywords 2D unpaired comparison
#' @export
#' @examples
#' balanced sampling 

balance_sample_purity_bandwidth = function(d1_data,
                                           d2_data,
                                           s1_col_name,
                                           s2_col_name,
                                           nna_cutoff,
                                           permute_time,
                                           working_dir,
                                           compare_name,
                                           bandwidth_factor = 1, 
                                           adjust_purity_flag,
                                           annot_file)
{
  

    # 
  s1_col_name = intersect(colnames(d1_data), s1_col_name)
  s2_col_name = intersect(colnames(d1_data), s2_col_name)
  
  
  s1_len = length(s1_col_name)
  s2_len = length(s2_col_name)
  
  

  if(s1_len <= s2_len)
  {
    m = floor(s2_len/s1_len)
    
    slices = get_slice(s2_len, m, seeds[1])
    
    
    for(k in 1:m)
    {
      #k = 1
      
      cn = paste0("slice", k, "_", compare_name)
      
      slice_result = comparison_time_points_2d_limma_bandwidth(d1_data = d1_data,
                                                               d2_data = d2_data,
                                                               s1_col_name = s1_col_name,
                                                               s2_col_name = s2_col_name[slices[[k]]],
                                                               nna_cutoff = nna_cutoff,
                                                               permute_time = permute_time,
                                                               working_dir = working_dir,
                                                               compare_name = cn,
                                                               col_annot_file = annot_file,
                                                               bandF = bandwidth_factor,
                                                               is_purity_adjusted = adjust_purity_flag)
      
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
      #k = 1
      
      cn = paste0("slice", k, "_", compare_name)
      
      slice_result = comparison_time_points_2d_limma_bandwidth(d1_data = d1_data,
                                                               d2_data = d2_data,
                                                               s1_col_name = s1_col_name[slices[[k]]],
                                                               s2_col_name = s2_col_name,
                                                               nna_cutoff = nna_cutoff,
                                                               permute_time = permute_time,
                                                               working_dir = working_dir,
                                                               compare_name = cn,
                                                               col_annot_file = annot_file,
                                                               bandF = bandwidth_factor,
                                                               is_purity_adjusted = adjust_purity_flag)
      
      cat(k, slice_result)    
      
    }
    
  }
  
  
  
}




#### extract results and make them into edge and nodes files 



#' generate node and edge file for cytoscape visualization
#' 
#' map the knase and substrate information
#' @param working_dir the directory of the files that KSA2D main results are in  
#' @param output_dir the directory for the output files of this function 
#' @param freq_cutoff an integer indicating the number of times a kinase to substrate relationship is evaluated across all slices, default to half of the number of total slices 
#' @param fdr_cutoff the cutoff below which kinase substrate relationships are seen as significant, default to 0.05
#' @param kinase_abs_fc_cutoff a positive numerical value indicating the absolute fold change cutoff of kinase 
#' @param subsite_abs_fc_cutoff a positive numerical value indicating the absolute fold change cutoff of substrate site
#' @param ksNetwork_filename absolute directory to the kinase substrate network file, as an example, merge_ks_omni.tsv. 
#' @import dplyr data.table magrittr 
#' @export
#' 
ksa2d_result_cytoscape = function(working_dir,
                                  output_dir,
                                  freq_cutoff = NULL,
                                  fdr_cutoff = 0.05, 
                                  kinase_abs_fc_cutoff, 
                                  subsite_abs_fc_cutoff,
                                  ksNetwork_filename)
{
  
  # working_dir = "/Users/ginny/My Drive/KSA2D_v2/test_output/ncc_105_1/"
  # output_dir = "/Users/ginny/My Drive/KSA2D_v2/test_output/ncc_105_1/test1/"
  # freq_cutoff = NULL
  # fdr_cutoff = 0.1
  # kinase_abs_fc_cutoff = 0.05
  # subsite_abs_fc_cutoff = 0.5
  # ksNetwork_filename = "/Users/ginny/My Drive/nonCCRCC20211020/merge_ks_omni.tsv"
  # 
  
  
  all_files = list.files(working_dir)
  all_results_files = grep("\\.tsv", all_files, value = T)
  all_results_files = grep("slice", all_results_files, value = T)
  
rl = list()
names = c()

k = length(all_results_files)

for(i in 1: k)
{

  file_name = all_results_files[i]
  file_dir = paste0(working_dir, file_name)
  s = fread(file_dir,stringsAsFactors = F, data.table = F)
  
  names = c(names, s$name)
  
  rl[[i]] = s
  
  
}

if(is.null(freq_cutoff))
{
  freq_cutoff_use = ceiling(k/2)
  
}else{
  freq_cutoff_use = freq_cutoff
}

name_freq = as.data.frame(table(names))%>%
  dplyr::arrange(desc(Freq))%>%
  dplyr::filter(Freq >= freq_cutoff_use)


merge_r = rbindlist(lapply(1:length(rl), function(x) {
  
  it = rl[[x]]
  
  return(it)

}))

filter_merge_r = merge_r%>%
  dplyr::filter(name %in% name_freq$names)%>%
  dplyr::group_by(name)%>%
  dplyr::summarise_if(is.numeric, mean, na.rm = T)%>%
  dplyr::ungroup()%>%
  as.data.frame()


dir.create(output_dir)

write.table(filter_merge_r, paste0(output_dir,"merge_all.tsv"),
                              quote = F, row.names = F, sep = "\t")
            

merge_sig = filter_merge_r%>%
  dplyr::filter(fdr<fdr_cutoff, abs(fc_kinase)>kinase_abs_fc_cutoff, abs(fc_substrate)>subsite_abs_fc_cutoff)



write.table(merge_sig, paste0(output_dir,"merge_sig.tsv"),
            quote = F, row.names = F, sep = "\t")


omni = fread(ksNetwork_filename,
             stringsAsFactors = F, data.table = F)


sig_k = gsub("_.*","", merge_sig$name)

sig_s =rbindlist(lapply(1:nrow(merge_sig), function(x){
  
  t = merge_sig$name[x]
  st = unlist(strsplit(t, split = "_"))
  
  s = st[2]
  
  p = st[3]
  df = data.frame(substrate = s, resPos = p, stringsAsFactors = F)
  
  
  return(df)
}))

sig_df = merge_sig%>%
  dplyr::mutate(kinase = sig_k, substrate =sig_s$substrate, resPos = sig_s$resPos)%>%
  dplyr::select(kinase, substrate, resPos, everything())

all_nodes = unique(c(sig_k, sig_s$substrate))

##### build a network with these nodes 


all_nodes_network = omni%>%
  dplyr::filter(kinase %in% all_nodes, substrate %in% all_nodes)%>%
  dplyr::select(kinase, substrate)%>%
  unique()

edge = rbindlist(lapply(1:nrow(all_nodes_network), function(x) {
  
  k = all_nodes_network$kinase[x]
  s = all_nodes_network$substrate[x]

  edge_fdr = 1
  edge_found = F
  edge_num = 0
  edge_both_positive = 0
  edge_both_negative = 0
  edge_kp_sn = 0
  edge_kn_sp = 0 
  
  edge_dir = "none"
  
  sig_sub = sig_df%>%
    dplyr::filter(kinase ==k, substrate == s)
  
  if(nrow(sig_sub)>0)
  {
    
    edge_fdr = min(sig_sub$fdr)
    edge_found = T
    edge_num = nrow(sig_sub)
    
    
    signs = sig_sub$fc_kinase*sig_sub$fc_substrate
    
    
    if(sig_sub$fc_kinase[1]>0)
    {
      
      t = which(signs>0)
      
      edge_both_positive = length(t)
      
      f = which(signs<0)
      
      edge_kp_sn = length(f)
      
      if(edge_both_positive>0)
      {
        edge_dir = "bothP"
      }else{
        edge_dir = "kPsN"
      }
      
      
    }
    
    if(sig_sub$fc_kinase[1]<0)
    {
      
      t = which(signs>0)
      
      edge_both_negative = length(t)
      
      f = which(signs<0)
      
      edge_kn_sp = length(f)
      
      
      if(edge_both_negative>0)
      {
        edge_dir = "bothN"
      }else{
        edge_dir = "kNsP"
      }
      
    }
  }
  
  
  
  edge_sig = -log10(edge_fdr)
  edge_sig_bl = T
  if(edge_fdr>fdr_cutoff)
    edge_sig_bl = F
  
  
  df = data.frame(kinase = k, substrate = s, edge_found,
                  edge_num, edge_both_positive, edge_both_negative, edge_kn_sp, edge_kp_sn, edge_dir, 
                  edge_fdr, edge_sig, edge_sig_bl, stringsAsFactors = F)
  return(df)
  
}))

write.table(edge, paste0(output_dir, "edge.tsv"),
            quote = F, row.names = F, sep = "\t")


#### this is building up the egde table, I need to have a node table too 

sig_edge = edge%>%
  dplyr::filter(edge_sig_bl == T)


node = rbindlist(lapply(1:length(all_nodes), function(x) {
  
  n = all_nodes[x]
  
  kinase_bl= F
  kinase_fc =0 
  site_num = 0
  substrate_av_fc =0
  
  nsig_edge = 0

  sub_k = sig_df%>%
    dplyr::filter(kinase == n)
  
  
  if(nrow(sub_k)>0)
  {
    kinase_bl  = T
    kinase_fc = sub_k$fc_kinase[1]
  }
  
  sub_s = sig_df%>%
    dplyr::filter(substrate == n)
  
  if(nrow(sub_s)>0)
  {
    site_num = nrow(sub_s)
    substrate_av_fc = mean(sub_s$fc_substrate)
  }
  
  sig_edge1 = sig_edge%>%
    dplyr::mutate(name = paste(kinase, substrate, sep = "_"))
  nsig_edge = length(grep(n, sig_edge1$name))
  
  df = data.frame(node = n, kinase_bl, kinase_fc, site_num, nsig_edge,substrate_av_fc, stringsAsFactors = F)
  
  # nsig_edge: number of significant edges this node involves in 
}))


write.table(node, paste0(output_dir, "node.tsv"),
            quote = F, row.names = F, sep = "\t")
}







