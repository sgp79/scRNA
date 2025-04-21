library(data.table)
library(Seurat)
library(ggplot2)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 10 * 1024^3) 

nn_leiden_umap <- function(s_obj, reduction, dims=1:30, leidenres=1, grph_name=NULL, clst_name=NULL, umap_name=NULL){
    # takes the s_obj, and the name of the data layer ("readuction") to work with (in the active expt)
    # gets knn and performs leiden clustering using this
    # performs umap
    # returns object
    if(is.null(grph_name)){grph_name = paste0(reduction, '.knn')}
    if(is.null(clst_name)){clst_name = paste0(reduction, '.leiden')}
    if(is.null(umap_name)){umap_name = paste0(reduction, 'umap')}

    s_obj <- FindNeighbors(object = s_obj, reduction = reduction, graph.name = grph_name, dims = dims)
    s_obj <- FindClusters(object = s_obj, graph.name = grph_name, cluster.name = clst_name, algorithm = 4, resolution = leidenres)
    s_obj <- RunUMAP(object = s_obj, dims = dims, reduction = reduction, reduction.name = umap_name)
    return(s_obj)
}

seurat_from_csv <- function(csvfile, ...){
    # takes a path to a csv and returns a seurat object
    # additional parameters supplied are passed to fread
    # assumes that the first column contains the gene names
    DT <- fread(csvfile, header=TRUE, ...)  # load (using a fast data.table function)
    mat <- as.matrix(DT[,2:ncol(DT)])  # matrix from all columns other than gene ID
    rownames(mat) <- DT[[1]]  # and add the gene IDs as row names
    invisible(CreateSeuratObject(counts = mat, min.features = 1))
}


add_mt_pct <- function(s_obj, mt_ids=c()){
    # takes a seurat object and a vector of gene IDs, and finds the % of reads which comprise those IDs.
    # The % is added to cell metadata as $percent.mt, and the modified surat object is returned
    if(length(mt_ids) > 0){
        # if any are not in the dataset, an exception is raised, so remove missing genes from the list first
        ids_in_dataset <- dimnames(s_obj$RNA@features)
        mt_ids_in_dataset <- mt_ids[mt_ids %in% ids_in_dataset]
        if(length(mt_ids_in_dataset) == 0) {
            cat("No mitochondrial genes found in dataset\n")
        }else{
            s_obj[["percent.mt"]] <- PercentageFeatureSet(s_obj, features = mt_ids_in_dataset)
        }
    }
    return(s_obj)
}


apply_filters_to_list <- function(s_obj_list, filter_list){
    # iterate through the lists, subsetting where a filter is supplied
    for(batch in names(s_obj_list)){
        cat(paste0('Filtering ', as.character(batch), '\n'))
        if(!is.null(filter_list[[batch]])){  # if a filter is supplied for this object
            parsed_expr <- rlang::parse_expr(filter_list[[batch]])  # Convert string to expression
            # Evaluate the condition in the metadata context
            mtdata <- s_obj_list[[batch]]@meta.data
            cells_to_keep <- rlang::eval_tidy(parsed_expr, data = mtdata)
            # subset using the list of cells from the evaluation
            s_obj_list[[batch]] <- subset(s_obj_list[[batch]], cells = rownames(mtdata)[cells_to_keep])
        }else{
            cat(paste0('No filter supplied for ', as.character(batch), '\n'))
        }
    }
    return(s_obj_list)
}


filter_norm_dimreduct <- function(s_obj_list, filter_list){
    # this function takes a list of Seurat objects and a list of filter criteria, then:
    # filters each object according to the supplied criteria, then normalises using
    #     the SCTransform v2 method, and performs PCA
    # filter criteria should be supplied as a string representing the format used by Seurat::subset, 
    #     or be NULL
    # It returns a seurat object with each batch in a layer
    # Example usage: 
    #      s_obj <- (s_obj, subset = nCount_RNA >= 50000 & nFeature_RNA >= 150)

    # check lists are the same length (they should be ordered)
    # this isn't quite right at the moment, as I changed other code, but is a worthwhile check as-is
    if(length(s_obj_list) != length(filter_list)){
        stop("filter_norm_dimreduct() must be given two named lists: s_obj_list and filter_list")
    }
    # iterate through the lists, subsetting where a filter is supplied
    s_obj_list <- apply_filters_to_list(s_obj_list, filter_list)

    if(length(s_obj_list) > 1){
        cat('Combining into single object\n')
        comb_s_obj <- merge(x = s_obj_list[[1]], y = s_obj_list[2:length(s_obj_list)],
                           merge.data=FALSE)
        # split into layers according to batch
        #comb_s_obj[["RNA"]] <- split(comb_s_obj[["RNA"]], f = comb_s_obj$batch) # above does so
    }else{
        cat('Only one object supplied\n')
        comb_s_obj <- s_obj_list[[1]]
    }
    cat('Normalising and finding highly variable features\n')
    comb_s_obj <- SCTransform(comb_s_obj)  # v2 by default
    cat('Running PCA\n')
    comb_s_obj <- RunPCA(comb_s_obj, verbose=FALSE)
    # split SCT assay into layers by batch (it combines them, for some reason)
    #comb_s_obj[["SCT_layered"]] <- split(comb_s_obj[["SCT"]], f = comb_s_obj$batch)
    #DefaultAssay(comb_s_obj) <- 'SCT_layered'

    #cat('Complete, returning combined object, with each layer normalised individually\n')
    return(comb_s_obj)
}


add_metacsv <- function(fn, s_obj, cell.id.col=1, ident.col=NULL, drop=NULL, slotnames=c(), ...){
    # a function to read a metadata csv and add it to a seurat object
    # just flexible enough for the datasets we're dealing with here
    # parameters:
    #     fn: a csv file containing metadata
    #     s_obj: the seurat object to modify
    #     cell_id_col: the column name or index to use as the row index (defaults to the first column)
    #     drop: a vector giving columns in the csv to exclude
    #     slotnames: a named vector to instruct column renaming - the names to give colnames in the csv
    #                                                             the values to give the desired slot names
    #     ident.col: set the cell identity property to this value (default of NULL leaves it unchanged)
    #     additional parameters will be passed to fread
    #
    # order of operations: 1. record cell.id.column seperately and drop it
    #                      2. drop columns named in drop
    #                      3. rename columns as specified by slotnames
    #                      4. set the cell identity property to the value of the column given in ident.col
    #     -> therefore if a column that is to be renamed is to be used for identity, ident.col
    #        should contain its final name
    #
    DT <- fread(fn, header=TRUE, ...)

    DT <- DT[get(cell.id.col) %in% colnames(s_obj)]  # drop meta data with no matching cell in seurat object
    
    cell_ids <- DT[[cell.id.col]]  # record the cell IDs
    DT[,c(cell.id.col):=NULL]  # and drop the cell IDs from the table

    # drop specified columns
    if(!is.null(drop)){
        DT[, (drop) := NULL]
    }
    
    # rename columns
    if(length(slotnames) > 0){
        if(is.null(names(slotnames))){
            stop('If slotnames are supplied to read_metacsv(), they must be given as a named
                  vector, with the names being the csv column names and the values the names
                  added to the seurat object metadata.  csv columns which are not specified
                  in slotnames OR in drop will be added using the unmodified csv column name')
        }
        setnames(DT, names(slotnames), slotnames)
    }
    DF <- as.data.frame(DT)  # convert to a dataframe
    rownames(DF) <- cell_ids  # and add back the cell IDs as row names

    # now add the metadata to the object
    s_obj <- SeuratObject::AddMetaData(s_obj, DF)

    # if ident.col is set:
    if(!is.null(ident.col)){
        Idents(s_obj) <- s_obj@meta.data[, ident.col]
    }
    invisible(s_obj)  # return
}
