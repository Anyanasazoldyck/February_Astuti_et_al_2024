#create seurat object####
#set a list
create_seurat_list <- function(sample_sheet, data_path){
  seurat_list <- list()
  
  for (sample_name in sample_sheet$sample_id){
    
    #record sample path
    message (sample_name)
    feat = paste0(sample_name,"_features.tsv.gz")
    bc = paste0(sample_name,"_barcodes.tsv.gz")
    mtx = paste0(sample_name,"_matrix.mtx.gz")
    
    feat_file = file.path(data_path,feat)
    bc_file = file.path(data_path,bc)
    mtx_file = file.path(data_path,mtx)
    
    if (!file.exists(mtx_file)) {
      print("mtx does not exist")
    }
    
    if (!file.exists(feat_file)) {
      print("feat_file does not exist")
    }
    
    if (!file.exists(bc_file)) {
      print("bc_file does not exist")
    }
    
    #read data 
    data <- ReadMtx  (
      mtx  = mtx_file,
      features = feat_file,
      cells = bc_file
    )
    #create seurat samples while setting the name of the asample_sheetay and project
    sc_obj= CreateSeuratObject(counts = data,project = sample_name,assay = "RNA")
    
    
    
    #add metadata 
    meta_row = sample_sheet[sample_sheet$sample_id==sample_name,] #match meta to p
    
    sc_obj$group = meta_row$group
    sc_obj$timepoint = meta_row$timepoint
    
    seurat_list[[sample_name]]<-sc_obj
  }
  
  return(seurat_list)
}
