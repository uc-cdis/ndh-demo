# Install and require dependent library
source("http://www.bioconductor.org/biocLite.R")
load.lib<-c("httr","jsonlite","dplyr","gplots","RColorBrewer","xml2","repr","VennDiagram")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
load.s3<-c("limma")
install.s3<-load.s3[!load.s3 %in% installed.packages()]
for(s3 in install.s3) biocLite(s3)
sapply(load.lib,require,character=TRUE)
sapply(load.s3,require,character=TRUE)

# Input "credential.json" downloaded from data common profile page to generate access token for data query
add_keys <- function(filename){
    json_data <- readChar(filename, file.info(filename)$size)
    keys <- fromJSON(json_data)
    variable <- jsonlite::toJSON(list(api_key = keys$api_key), auto_unbox = TRUE)
    auth <- POST('https://niaid.bionimbus.org/user/credentials/cdis/access_token', add_headers("Content-Type" = "application/json"), body = variable)
    return(auth)
}

# Create Graphql query text based on inputs study_id(in study node), aliquote_type(in aliquote node) and file_type (file node name).
query_cell_file <- function(study_id,aliquot_type,file_type){
    query_txt = paste('{
  study(submitter_id:"',study_id,'",first:0){
    subjects(first:0){
      samples(first:0,composition:"Cell",with_path_to:{type:"aliquot",analyte_type:"',aliquot_type,'"}){
        hours_to_collection
        virus_infections(first:0){
          strain
          mutation
        }
        aliquots(first:0,with_links:"',file_type,'"){',
          file_type,'{
            file_name
            id
          }
        }
      }
    }
  }
}',sep="")
    return(query_txt)
}

# query data from API using the query text created by create_query_txt function. The content of the response is stored in an array
query_api <- function(query_txt){
    auth = add_keys("credentials.json")
    query_txt = query_txt
    query <- jsonlite::toJSON(list(query = query_txt), auto_unbox = TRUE)
    token <- paste('bearer', content(auth)$access_token, sep=" ")
    response <- POST('https://niaid.bionimbus.org/api/v0/submission/graphql',add_headers("Authorization" = token, "Content-Type" = "application/json"), body = query)
    return(content(response)$data)
}

# Organize response content and coverted to a table that has four columns "virus", "time_point","FileName" and "uuid"
parse_cell_file<- function(study_id,aliquot_type,file_type){
    query_txt = query_cell_file(study_id,aliquot_type,file_type)
    response = query_api(query_txt)
    data = response$study[[1]]$subjects[[1]][[1]]
    datalist = list()
    i = 1
    for (instance in data){
        hours_to_collection = instance$hours_to_collection
        strain = instance$virus_infections[[1]]$strain
        if (length(strain) !=0){
            strain = strain
        }else{
        strain = "mock"}
        mutation = instance$virus_infection[[1]]$mutation
        if (length(mutation)!=0){
            mutation = mutation
        }else{
        mutation = "NA"}
        entities = sapply(instance$aliquots,function(x) x[[file_type]])
        files = vector('character')
        ids = vector('character')
        j = 0
        for (entity in entities){
            if(!(is.null(entity$file))){
                j = j + 1
                files[j] = entity$file
                ids[j] = entity$id

            }else{
                m = 1
                for (ind in entity){
                    n = j + m
                    files[n] = ind$file
                    ids[n] = ind$id
                    m = m + 1
                }
                j = n
            }

        }
        if(length(files)!=0){
            for (index in 1:length(files)){
            each_aliquot = data.frame(virus = paste(strain,mutation,sep="_"), time_point = hours_to_collection, FileName = files[index],uuid = ids[index])
            datalist[[i]] = each_aliquot
            i = i+1
            }
        }
    }
    target_matrix = do.call(rbind,datalist)
    return(target_matrix)
}

# Verify if a file has been download already.  Create folder if it hasn't and proceed to download the file.  Store the files under the folder.
download_data <- function(study_id,aliquot_type,file_type){
    if(!dir.exists(file.path(paste(study_id,file_type,sep="_")))){
        dir.create(file.path(paste(study_id,file_type,sep="_")))
        metadata = parse_cell_file(study_id,aliquot_type,file_type)
        download_matrix = unique(data.frame(uuid=metadata$uuid,file_name=metadata$FileName))
        auth = add_keys("credentials.json")
        token <- paste('bearer', content(auth)$access_token, sep=" ")
        # Based on file_type, download file from aws s3 or ftp_url
        if(file_type %in% c("mRNA_microarrays","protein_expressions")){
            for (i in 1:length(download_matrix$uuid)){
            response <- GET(paste('https://niaid.bionimbus.org/user/data/download/',download_matrix$uuid[i],sep=""),add_headers("Authorization" = token))
            response_file <- GET(content(response)$url,write_disk(file.path(paste(paste(study_id,file_type,sep="_"),download_matrix$file_name[i],sep="/")),overwrite=TRUE))
        }
        }else{
           for(i in 1:length(download_matrix$uuid)){
               response <- GET(paste('https://niaid.bionimbus.org/index/index/',download_matrix$uuid[i],sep=""),add_headers("Authorization" = token))
system(paste("wget","-O",paste(paste(study_id,file_type,sep="_"),download_matrix$file_name[i],sep="/"),content(response)$urls[[1]],sep=" "),intern=TRUE)
           }
        }
        return("Finished Downloading")
    }else{return("Data Already Exist")}
}

# Use function parse_cell_file to create design matrix. Use design matrix and downloaded files as input to create the aggregated normalized gene expression data frame
array_normalization <- function(study_id){
    # Load arrays
    array_Dir <- paste(study_id,"mRNA_microarrays",sep="_")
    meta = parse_cell_file(study_id,"RNA","mRNA_microarrays")
    meta$SampleName = sub(".txt","",meta$FileName)
    meta$uuid <- NULL
    study_array <- list()
    study_array$targets <- meta
    otherColumns <- c("gIsFound", "gIsWellAboveBG", "gIsSaturated", "gIsFeatNonUnifOL", "gIsFeatPopnOL")
    study_array$rawData <- read.maimages(study_array$targets, path=array_Dir, source="agilent", green.only=TRUE, other.columns=otherColumns)
    # BG correct and normalize
    study_array$bgcData50 <- backgroundCorrect(study_array$rawData, method="normexp", normexp.method="mle", offset=50)
    array_norm <- normalizeBetweenArrays(study_array$bgcData50, method="quantile")
    probeInfo <- array_norm$genes
    array_norm <- array_norm[probeInfo$ControlType == 0,]
    study_array$normData <- array_norm
    array_norm <- avereps(array_norm,ID=array_norm$genes[, "ProbeName"])
    study_array$normData_collapsed <- array_norm
    ProbeID = study_array$normData_collapsed$genes$ProbeID
    GeneName = study_array$normData_collapsed$genes$GeneName
    norm_exp = study_array$normData_collapsed$E
    if(all(rownames(norm_exp)== ProbeID)){
        rownames(norm_exp) = GeneName
    }
    return(norm_exp)
}

# Create gene expression by averaging probes
aggregate_gene_exp <- function(study_id){
    probe_exp = array_normalization(study_id)
    gene_exp <- tbl_df(probe_exp) %>% group_by(rownames(probe_exp))%>% summarise_all(funs(mean))
    gene_exp = data.frame(gene_exp)
    rownames(gene_exp) = gene_exp[,1]
    gene_exp <-gene_exp[,-1]
    return(gene_exp)
}

# Implement limma package to perform differential gene expression analysis. The return data frame have log2FC at each timepoint. The overall p-value and adjusted p-value using Benjermin multiple test correction
DE_gene <-function(study_id,virus,exclude_time_points=NULL){
    # Create design matrix
    targets = parse_cell_file(study_id,"RNA","mRNA_microarrays")
    targets$virus = gsub(" ","_",targets$virus)
    target = paste(targets$virus,targets$time_point,sep=".")
    times = sort(unique(targets$time_point),decreasing = FALSE)
    gene_exp = aggregate_gene_exp(study_id)
    lev = unique(target)
    f = factor(target,levels=lev)
    design <- model.matrix(~0+f)
    colnames(design)<-lev
    # fit linear model
    fit<-lmFit(gene_exp,design)
    times = times[!times%in%exclude_time_points]
    # identify contrast groups
    contrast_groups = rep(NA,length(times))
    for (i in 1:length(times)){
        contr1 = paste(virus[1],times[i],sep=".")
        contr2 = paste(virus[2],times[i],sep=".")
        contrast_groups[i] = paste(contr2,contr1,sep="-")
    }
    # Calculate Log2FC and p value for contrast groups
    cont.mu <- makeContrasts(contrasts = contrast_groups,levels = design)
    fit2 = contrasts.fit(fit,cont.mu)
    fit2 <- eBayes(fit2)
    DE_gene = topTable(fit2,n=Inf,adjust="BH")
    virus_string = paste(virus,collapse = "vs")
    write.table(DE_gene,file = paste(study_id,virus_string,"DE.txt",sep="."),quote=F,sep="\t",col.names=T,row.names=T)
    return(DE_gene)
}

# For each study, select gene expression of type I interferon signature
Select_ISG <- function(ISG,dataset){
    ISG_gene = read.table(ISG)
    ISG_gene = ISG_gene$V1
    DE_gene = read.table(dataset,header=T,row.names=1)
    ISG_DE = DE_gene[rownames(DE_gene)%in%ISG_gene,]
    ISG_DE_order = ISG_DE[order(rev(ISG_DE)[5],decreasing = FALSE),]
    write.table(ISG_DE_order,paste(gsub(".txt","",dataset),"ISG.txt",sep="."),sep="\t",quote=F,col.names=T,row.names=T)
}

# Plot heatmap for each individual study ISG genes.
headmap_plot <- function(dataset){
    par(cex.main=0.4)
    ISG_DE = read.table(dataset,header=T,row.names=1)
    ISG_DE = subset(ISG_DE,select = -c(AveExpr,F,P.Value,adj.P.Val))
    ISG_DE_matrix = data.matrix(ISG_DE)
    my_palette <- colorRampPalette(c("green", "black", "red"))
    col_breaks = c(seq(-3,-0.75,length=10),seq(-0.74,0.74,length=10),seq(0.75,3,length=10))
    heatmap.2(ISG_DE_matrix,density.info="none",trace="none",margins =c(22,22),col=my_palette,Colv="NA",breaks=col_breaks,dendrogram="none",Rowv="NA",labRow=FALSE, colRow = FALSE,main=gsub(".DE.ISG.txt","",dataset),key=T, cexRow = 0.5, cexCol=0.5,keysize = 1.2,key.par = list(cex=0.5),lhei=c(1,6), lwid=c(1,3))
}

# Select ISG signatures that are common across different studies.  The input parameter- "datasets" is a vector that contains all the study ids.
common_ISG <- function(ISG,datasets){
    ISG_gene = read.table(ISG)
    ISG_gene = ISG_gene$V1
    for (i in 1:length(datasets)){
        DE_gene = read.table(datasets[i],header=T,row.names=1)
        ISG_DE = rownames(DE_gene[rownames(DE_gene)%in%ISG_gene,])
        ISG_gene = ISG_gene[ISG_gene%in%ISG_DE]
    }
    for (i in 1:length(datasets)){
        DE_gene = read.table(datasets[i],header=T,row.names=1)
        ISG_DE = DE_gene[rownames(DE_gene)%in%ISG_gene,]
        ISG_DE_order = ISG_DE[order(rev(ISG_DE)[5],decreasing = FALSE),]
        write.table(ISG_DE_order,paste(gsub(".txt","",datasets[i]),"ISG_common.txt",sep="."),sep="\t",quote=F,col.names=T,row.names=T)
    }
    return(ISG_gene)
}

# Rank ISG signatures. All the signatures were catagorized into three groups: commonly down-regulated in dataset1 and dataset2, down-regulated in dataset1 only and the rest group. Within each group, genes are ranked from smallest to largest at the last timepoint in dataset1.
order_ISG <- function(dataset1,dataset2){
    ISG_DE1 = read.table(dataset1,header=T,row.names=1)
    ISG_DE2 = read.table(dataset2,header=T,row.names=1)
    ISG_DE1_down = rownames(ISG_DE1[rev(ISG_DE1)[5]< 0,])
    ISG_DE2_down = rownames(ISG_DE2[rev(ISG_DE2)[5]< 0,])
    common_down = ISG_DE1_down[ISG_DE1_down%in%ISG_DE2_down]
    ISG_DE1_com_down = ISG_DE1[rownames(ISG_DE1)%in%common_down,]
    ISG_DE1_com_down_sort = ISG_DE1_com_down[order(rev(ISG_DE1_com_down)[5],decreasing = FALSE),]
    ISG_DE1_dis_down = ISG_DE1[!rownames(ISG_DE1)%in%common_down,]
    ISG_DE1_dis_down_sort = ISG_DE1_dis_down[order(rev(ISG_DE1_dis_down)[5],decreasing = FALSE),]
    ISG_DE1_order = rbind(ISG_DE1_com_down_sort,ISG_DE1_dis_down_sort)
    return(rownames(ISG_DE1_order))
}

# Plot heatmap across studies. The ISG signatures are ordered by function order_ISG.
heatmap_plot_across <- function(datasets,dataset1,dataset2){
    gene_list = order_ISG(dataset1,dataset2)
    par(cex.main=0.4)
    for (i in 1:length(datasets)){
        ISG_DE = read.table(datasets[i],header=T,row.names=1)
        ISG_DE_order = ISG_DE[match(gene_list,rownames(ISG_DE)),]
        ISG_DE_order = subset(ISG_DE_order,select = -c(AveExpr,F,P.Value,adj.P.Val))
        ISG_DE_matrix = data.matrix(ISG_DE_order)
        my_palette <- colorRampPalette(c("green", "black", "red"))
        col_breaks = c(seq(-3,-0.75,length=10),seq(-0.74,0.74,length=10),seq(0.75,3,length=10))
        heatmap.2(ISG_DE_matrix,density.info="none",trace="none",margins =c(22,18),col=my_palette,Colv="NA",breaks=col_breaks,dendrogram="none",Rowv="NA",labRow=FALSE, colRow = FALSE,main=gsub(".DE.ISG","",datasets[i]),key=T, cexRow = 0.5, cexCol=0.5,keysize = 1.2,key.par = list(cex=0.5),lhei=c(1,6), lwid=c(1,3))
    }
}
# Plot heatmap for single study. The ISG signature are ordered by function order_ISG. The plot is saved in .png format with high resolution
heatmap_plot_single <- function(dataset, dataset1,dataset2){
    gene_list = order_ISG(dataset1,dataset2)
    ISG_DE = read.table(dataset,header=T,row.names=1)
    ISG_DE_order = ISG_DE[match(gene_list,rownames(ISG_DE)),]
    ISG_DE_order = subset(ISG_DE_order,select = -c(AveExpr,F,P.Value,adj.P.Val))
    ISG_DE_matrix = data.matrix(ISG_DE_order)
    my_palette <- colorRampPalette(c("green", "black", "red"))
    col_breaks = c(seq(-3,-0.75,length=10),seq(-0.74,0.74,length=10),seq(0.75,3,length=10))
    par(mar=c(1,1,1,1))
    png(paste(gsub(".txt","",dataset),"png",sep="."),height=1080,width=1080,res=100)
    heatmap.2(ISG_DE_matrix,density.info="none",trace="none",col=my_palette,Colv="NA",breaks=col_breaks,dendrogram="none",Rowv="NA",labRow=FALSE, colRow = FALSE)
    dev.off()
}

# Protein differential expression were performed by t test using normalized protein expression profile as input. At each timepoint, virus infected group were compared to control group by applying t test. If both of the virus infected group and control groups having expression value for more than 2 samples, the t test is performed. Otherwise the fold change is recorded as NaN.  If the t test is performed, the t test p value is <0.05 and the virus infected group has higher protein level, the fold change is recorded as 1, otherwise is recorded as -1. If the p value is > 0.05, the fold change is recorded as 0. The fold change result is saved in file. The gene list of upregulated protein is return back.
protein_DE_test <- function(study_id,time=-1){
    # Graphql meta data associated with data
    metadata = parse_cell_file(study_id,"protein","protein_expressions")
    viruses = unique(metadata$virus)
    timepoints = as.character(sort(as.numeric(unique(metadata$time_point)),decreasing = FALSE))
    timepoints = timepoints[timepoints!=time]
    filename = unique(metadata$FileName)
    # Read normalized protein quantity
    protein_expression = read.table(paste(paste(study_id,"protein_expressions",sep="_"),filename,sep="/"),sep="\t",header=T)
    # Initialize list to store pvalue for t-test, expression change direction for t-test, pvalue for G-test and expression change direction for G-test for different timepoints and virus infections. Initialize vector to store protein symbols
    tstatlist = list()
    tfoldlist = list()
    gstatlist = list()
    gfoldlist = list()
    namelist = vector('character')
    n = 1
    # Compare virus infection and control group one timepoint at a time
    for (timepoint in timepoints){
    # Extract protein quantity for control at one specific timepoint
        control = protein_expression[,grepl(paste(paste("mock_NA",timepoint,sep="_"),"_",sep=""),names(protein_expression))]
        # Extract protein quantity for one virus infection at a time and compare with control
        for (virus in viruses){
            if(!"mock_NA"%in% virus){
                test = protein_expression[,grepl(paste(paste(virus,timepoint,sep="_"),"_",sep=""),names(protein_expression))]
                # Initiate vectors to store pvalue for t-test, expression change direction for t-test, pvalue for G-test and expression change direction for G-test at specific timepoint and virus infection
                tstats = vector('numeric')
                tfolds = vector('numeric')
                gstats = vector('numeric')
                gfolds = vector('numeric')
                # Perform t-test or G-test for each protein
                for (i in 1:nrow(control)){
                    control_data = control[i,]
                    control_data = as.numeric(control_data[!is.na(control_data)])
                    test_data = test[i,]
                    test_data = as.numeric(test_data[!is.na(test_data)])
                    # If both the control and infection groups having more than 2 samples, t-test is performed. Otherwise, G-test is performed.
                    if(length(control_data)>=2 && length(test_data)>=2){
                        tstats[i] = as.numeric(t.test(control_data, test_data,var.equal = TRUE)[c("p.value")])
                        gstats[i] = "NaN"
                        gfolds[i] = "NaN"
                        if(tstats[i] <= 0.05 && mean(test_data)!=mean(control_data)){
                            tfolds[i] = ifelse(mean(test_data)-mean(control_data)>0,1,-1)
                        }else{tfolds[i]=0}
                    }else{
                        CO1 = length(control[i,!is.na(control[i,])])
                        CA1 = length(control[i,is.na(control[i,])])
                        CO2 = length(test[i,!is.na(test[i,])])
                        CA2 = length(test[i,is.na(test[i,])])
                        EO1 = length(control[i,])*(CO1 + CO2)/(CO1+CO2+CA1+CA2)
                        EA1 = length(control[i,])*(CA1 + CA2)/(CO1+CO2+CA1+CA2)
                        EO2 = length(test[i,])*(CO1 + CO2)/(CO1+CO2+CA1+CA2)
                        EA2 = length(test[i,])*(CA1 + CA2)/(CO1+CO2+CA1+CA2)
                        G = 2*(sum(CO1*log(CO1/EO1),CA1*log(CA1/EA1),CO2*log(CO2/EO2),CA2*log(CA2/EA2),na.rm=TRUE))
                        gstats[i] = pchisq(G,1,lower.tail=FALSE)
                        tstats[i] = "NaN"
                        tfolds[i] = "NaN"
                        if(gstats[i] <= 0.05 && CO1/(CO1+CA1)!= CO2/(CO2+CA2)){
                            gfolds[i] = ifelse(CO1/(CO1+CA1) < CO2/(CO2+CA2),1,-1)
                        }else{gfolds[i]=0}
                    }
                }
                tstatlist[[n]] = tstats
                tfoldlist[[n]] = tfolds
                gstatlist[[n]] = gstats
                gfoldlist[[n]] = gfolds
                namelist[n] = paste(virus,timepoint,sep="_")
                n = n + 1
            }
        }
    }
    # Create dataframes: pvalue for t-test, expression change direction for t-test, pvalue for G-test and expression change direction for G-test at all timepoints
    tstat_matrix = data.frame(do.call(cbind,tstatlist))
    tfold_matrix = data.frame(do.call(cbind,tfoldlist))
    gstat_matrix = data.frame(do.call(cbind,gstatlist))
    gfold_matrix = data.frame(do.call(cbind,gfoldlist))
    names(tstat_matrix) = sapply(namelist,function(x) paste("Tpvalue",x,sep="_"))
    names(tfold_matrix) = sapply(namelist,function(x) paste("Tdiff",x,sep="_"))
    names(gstat_matrix) = sapply(namelist,function(x) paste("Gpvalue",x,sep="_"))
    names(gfold_matrix) = sapply(namelist,function(x) paste("Gdiff",x,sep="_"))
    row.names(tstat_matrix) = protein_expression[,2]
    row.names(tfold_matrix) = protein_expression[,2]
    row.names(gstat_matrix) = protein_expression[,2]
    row.names(gfold_matrix) = protein_expression[,2]
    write.table(tfold_matrix,paste(study_id,"protein","Tfolds",sep="_"),sep="\t",quote=F,col.names=T,row.names=T)
    write.table(tstat_matrix,paste(study_id,"protein","Ttest",sep="_"),sep="\t",quote=F,col.names=T,row.names=T)
                                 write.table(gfold_matrix,paste(study_id,"protein","Gfolds",sep="_"),sep="\t",quote=F,col.names=T,row.names=T)
                                 write.table(gstat_matrix,paste(study_id,"protein","Gtest",sep="_"),sep="\t",quote=F,col.names=T,row.names=T)
    # Select up-regulated proteins tested by T-test or G-test
    tsignificant_up = tfold_matrix[apply(tfold_matrix,MARGIN=1,function(x) 1%in%x),]
    gsignificant_up = gfold_matrix[apply(gfold_matrix,MARGIN=1,function(x) 1%in%x),]
    num_significant = nrow(tsignificant_up) + nrow(gsignificant_up)
    print(paste(num_significant,"Up-regulated protein",sep=" "))
    protein_list = c(gsub("_HUMAN","",row.names(tsignificant_up)),gsub("_HUMAN","",row.names(gsignificant_up)))
    return(head(protein_list,10))
}

# Upregulated proteins and ISG among the upregulated proteins were plotted using Veen Diagram.
ISG_DE_protein <- function(study_id,signature){
    tfold_matrix = read.table(paste(study_id,"protein","Tfolds",sep="_"),sep="\t",header=T,row.names=1)
    tsignificant_up = tfold_matrix[apply(tfold_matrix,MARGIN=1,function(x) 1%in%x),]
    gfold_matrix = read.table(paste(study_id,"protein","Gfolds",sep="_"),sep="\t",header=T,row.names=1)
    gsignificant_up = gfold_matrix[apply(gfold_matrix,MARGIN=1,function(x) 1%in%x),]
    up_genes = c(gsub("_HUMAN","",row.names(tsignificant_up)),gsub("_HUMAN","",row.names(gsignificant_up)))
    print(length(up_genes))
    ISG_gene = read.table(signature)$V1
    print("Significant up-regulated ISG proteins:")
    print(intersect(up_genes,ISG_gene))
    total_number = length(up_genes)
    common_number = length(intersect(up_genes,ISG_gene))
    options(repr.plot.width=2, repr.plot.height=2)
    draw.pairwise.venn(area1 = total_number, area2 = common_number, cross.area = common_number, category = c("Significant upregulated proteins", "total ISG"),fill=c("blue","red"),cat.pos = c(360, 200),cat.cex = 0.6)
}

