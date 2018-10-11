# install.packages("devtools") # we need this to install github packages
# options(unzip = "unzip") # the default is "", we need to to set this option to unzip for installations
# ## Install igraph from github
# devtools::install_github("gaborcsardi/pkgconfig")
# devtools::install_github("igraph/rigraph")
# devtools::install_github("kassambara/ggpubr")
## Now I can install phyloseq
load.lib<-c("httr","jsonlite","xml2","repr","utils","ggplot2","reshape2","igraph","phyloseq","ggpubr")
# install.lib<-load.lib[!load.lib %in% installed.packages()]
# for(lib in install.lib) install.packages(lib,dependencies=TRUE, lib=NULL)
# load.s3<-c("phyloseq")
# install.s3<-load.s3[!load.s3 %in% installed.packages()]
sapply(load.lib,require,character=TRUE)
# sapply(load.s3,require,character=TRUE)

options(warn=-1)
#Create access key from credential file download from profile
add_keys <- function(filename){
    json_data <- readChar(filename, file.info(filename)$size)
    keys <- fromJSON(json_data)
    variable <- jsonlite::toJSON(list(api_key = keys$api_key), auto_unbox = TRUE)
    auth <- POST('https://niaid.bionimbus.org/user/credentials/cdis/access_token', add_headers("Content-Type" = "application/json"), body = variable)
    return(auth)
}

query_microbiome_info <-function(study_id,offset){
    query_txt = paste('{
  study(submitter_id:"',study_id,'",project_id: "ndh-dait-microbiome"){
    follow_ups(first:1,offset:',offset,'){
        visit_name
        subjects{
            submitter_id
        }
        samples(first:0){
            submitter_id
            aliquots(first:0){
                experiments{
                    experiment_type
                    sequencing_results{
                        object_id
                        file_name
                    }
                }
            }
        }
    }
  }
}',sep="")
    return(query_txt)
}

query_subject_counts <-function(study_id){
    query_txt = paste('{
        study(submitter_id:"',study_id,'", project_id:"ndh-dait-microbiome"){
            _follow_ups_count
            }
        }',sep="")
     return(query_txt)
}

# Query API
query_api <- function(query_txt){
    auth = add_keys("credentials.json")
    query_txt = query_txt
    query <- jsonlite::toJSON(list(query = query_txt), auto_unbox = TRUE)
    token <- paste('bearer', content(auth)$access_token, sep=" ")
    response <- POST('https://niaid.bionimbus.org/api/v0/submission/graphql',add_headers("Authorization" = token, "Content-Type" = "application/json"), body = query)
    return(content(response)$data)
}

# Parse reponse data to data frame
parse_microbiome_info <- function(study_id){
    count_query = query_subject_counts(study_id)
    counts_response = query_api(count_query)
    print(counts_response)
    counts = counts_response$study[[1]]$`_follow_ups_count`
    offset = 0
    int_data = data.frame(subject = character(), samples=character(), visit_name=character(), organ=character(), days=numeric(), uuid=character(),file_name=character())
    while(offset < counts){
        query_txt = query_microbiome_info(study_id,offset)
        response = query_api(query_txt)
        data = response$study[[1]]
        app_data = parse_microbiome_data(data)
        int_data = rbind(int_data, app_data)
        offset = offset + 1
    }
    write.table(int_data,paste(study_id,"microbiome_info.txt", sep="_"),sep="\t", col.names=T, row.names=F, quote=F)
    return (int_data)
}

parse_microbiome_data <- function(data){
    app_data = data.frame(subject=character(), samples=character(), visit_name=character(), organ=character(), days=numeric(), uuid=character(),file_name=character() )
    for (follow_up in data$follow_ups){
        visit_name = follow_up$visit_name
        if(length(visit_name) >0){
            visit_name = visit_name
        }else{
            visit_name = "NA"
        }
        for (subject in follow_up$subjects){
            subject_id = subject$submitter_id
        }
        for (sample in follow_up$samples){
            sample_id = sample$submitter_id
            for (aliquot in sample$aliquots){
                for (experiment in aliquot$experiments){
                experiment_type = experiment$experiment_type
                if(length(experiment_type) >0){
                    experiment_type = experiment_type
                    organ = strsplit(experiment_type," ")[[1]][3]
                    days = strsplit(experiment_type," ")[[1]][2]
                 }else{
                    organ = "NA"
                    days = "NA"
                }
                    for (sequence in experiment$sequencing_results){
                        file_id = sequence$object_id
                        file_name = sequence$file_name
                        each_sample = data.frame(subject = subject_id, samples = sample_id, visit_name = visit_name, organ = organ, days = days, uuid = file_id, file_name = file_name)
                        app_data = rbind(app_data,each_sample)
                    }
                }
            }
        }
    }
    return (app_data)
}

# Download files if the files are not exist
download_data <- function(study_id){
    if(!dir.exists(file.path(study_id))){
        dir.create(file.path(study_id))
        metadata = parse_microbiome_info(study_id)
        download_matrix = unique(data.frame(uuid=metadata$uuid,file_name=metadata$file_name))
        auth = add_keys("credentials.json")
        token <- paste('bearer', content(auth)$access_token, sep=" ")
        for (i in 1:length(download_matrix$uuid)){
            response <- GET(paste('https://niaid.bionimbus.org/user/data/download/',download_matrix$uuid[i],sep=""),add_headers("Authorization" = token))
            response_file <- GET(content(response)$url,write_disk(file.path(paste(study_id,download_matrix$file_name[i],sep="/")),overwrite=TRUE))
        }
        files = list.files(study_id)
        comb_file = read.table(paste(study_id,files[1],sep="/"),sep="\t",header=T,row.names=1)
        if(length(files) >1){
            for(i in 2:length(files)){
                data = read.table(paste(study_id,files[i],sep="/"),sep="\t",header=T, row.names=1)
                comb_file = cbind(comb_file,data)
            }
            names(comb_file) = metadata$samples
        }
        write.table(comb_file,paste(study_id,"microbiome_data.txt",sep="_"),col.names=T,row.names=T,quote=T,sep="\t")
        return("Finished Downloading")
    }else{return("Data Already Exist")}
}

# Construct physeq S4 object
construct_physeq <-function(study_id){
    otumat = read.table(paste(study_id,"microbiome_data.txt",sep="_"),sep="\t",header=T,row.names=1)
    otumat = otumat[rowSums(otumat)!=0,]
    samplemat = read.table(paste(study_id,"microbiome_info.txt", sep="_"), header=T, row.names=2, stringsAsFactors=F,sep="\t")
    samplemat = samplemat[order(match(rownames(samplemat),names(otumat))),]
    OTU = otu_table(otumat,taxa_are_rows = TRUE)
    Sample = sample_data(samplemat)
    physeq = phyloseq(OTU, Sample)
    return (physeq)
}

# Plot otu prevalence vs otu counts
prevalence_count <- function(study_id){
    physeq = construct_physeq(study_id)
    Count = rowSums(otu_table(physeq))
    Prevalance = rowSums(otu_table(physeq)>0)
    par(mfrow=c(1,1))
    plot(log(Prevalance)~log(Count),cex=0.5,pch=16)
}

# Plot

# alpha estimate by organ with Gestation between 10-40 weeks
alpha_estimate <-function(study_id, sites="NA"){
    divdfs <- list()
    physeq = construct_physeq(study_id)
    erDF <- estimate_richness(physeq, split = TRUE, measures = "Shannon")
    df <- data.frame(erDF, sample_data(physeq))
    if (sites!="NA"){
        for (site in sites){
            df_site <- df[df$organ == site,]
            df_site$weeks <- floor(df_site$days/7)
            df_site <- df_site[df_site$weeks >=10,]
            df_site <- df_site[df_site$weeks <=40,]
            divdfs[[site]] <- df_site
        }
        return (divdfs)
    }else{
        return (df)
    }
}

# Plot alpha diversity by organ
plot_alpha_diversity <- function(study_id, sites){
    divdfs = alpha_estimate(study_id, sites)
    par(mfrow=c(1,4),mar=c(2,4,6,2))
    for (site in sites){
        data = divdfs[[site]]
        data$weeks <- as.numeric(data$weeks)
        plot(Shannon~weeks, data, main = site, ylim=c(0,4.5), pch=16)
        fit <- lm(Shannon~weeks, data)
        newx <- seq(min(data$weeks), max(data$weeks), length.out=100)
        preds <- predict(fit, newdata = data.frame(weeks=newx),
                 interval = 'confidence')
        polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = rgb(0, 0, 0,0.5), border = NA)
        abline(fit,col="blue")
    }
}

#Plot alpha diversity across study for cell free DNA
plot_cf_alpha_diversity <- function(studies){
    shannon_visit = data.frame(Shannon = character(),visit = character(),study = character())
    for(study in studies){
        df = alpha_estimate(study)
        df = df[df$visit_name%in%c("Trimester 1","Trimester 2","Trimester 3","24h Post-Partum"),]
        df$study = rep(study,dim(df)[1])
        df = df[,c("Shannon","visit_name","study")]
        shannon_visit = rbind(shannon_visit,df)
    }
    stacked.data = melt(shannon_visit, id = c('visit_name', 'study'))
    boxplot(value ~ visit_name + study,at = c(1, 1.8, 2.6 ,3.4 ,6, 6.8, 7.6, 8.4), xaxt='n',stacked.data,col = c('red', 'blue', 'yellow','purple'),ylim=c(0,5))
    axis(side=1, at=c(1.8, 6.8), labels=c(studies[1], studies[2]), line=0.5, lwd=0)
    axis(side=2, at=2, labels="Shannon", line=1, lwd=0)
    text(c(1, 1.8, 2.6, 3.4, 6, 6.8, 7.6,8.4), c(4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5), c('T1', 'T2', 'T3','P', 'T1', 'T2', 'T3','P'))
    title('Comparing alpha diversity among \n different Trimesters in two studies')
}
# Calculate and plot Euclidean distance for samples withion one organ
beta_diversity <- function(study_id){
    Eucl_organ = data.frame(Euclidean_distance = as.numeric(),organ = as.character())
    physeq = construct_physeq(study_id)
    physeq <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
    for(site in c("Vaginal_Swab","Stool","Saliva","Tooth_Gum")){
        physeq_site <- prune_samples(sample_data(physeq)$organ==site,physeq)
        EuclDist_site = dist(t(otu_table(physeq_site)))
        data = data.frame(Euclidean_distance = as.numeric(EuclDist_site), organ = rep(site,length(EuclDist_site)))
        Eucl_organ = rbind(Eucl_organ,data)
    }
    par(cex.lab=0.5)
    my_comparison <- list(c("Vaginal_Swab","Stool"),c("Stool","Saliva"),c("Saliva","Tooth_Gum"))
    ggboxplot(Eucl_organ, x = "organ", y="Euclidean_distance", col = "organ") + stat_compare_means(comparisons = my_comparison) + stat_compare_means(label.y = 2, method = "anova")
}

beta_diversity_onesub <- function(study_id,subject){
    physeq = construct_physeq(study_id)
    physeq <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
    # Select one subject
    physeq_sub <- prune_samples(sample_data(physeq)$subject==subject,physeq)
    physeq_sub <- subset_samples(physeq_sub,organ%in%c("Stool","Vaginal_Swab","Tooth_Gum"))
    EuclDist_sub = dist(t(otu_table(physeq_sub)))
    EuclDistMatrix_sub <- as.matrix(EuclDist_sub)
    pc_euc_sub <- cmdscale(EuclDistMatrix_sub,k=4)
    write.table(pc_euc_sub,"pc_euc_sub.txt",sep="\t",quote=F, col.names=T, row.names=T)
    pc_euc_sub <- read.table("pc_euc_sub.txt",header=T, row.names=1)
    names(pc_euc_sub) <- c("PC1","PC2","PC3","PC4")
    pc_sample_sub <- sample_data(physeq_sub)
    pc_sample_sub <- cbind(pc_euc_sub,pc_sample_sub)
    par(mar=c(3,3,3,3))
    p <- ggplot(pc_sample_sub, aes(PC1,PC2))
    p + geom_point(size = 2, alpha = 0.3)
    p + geom_point(aes(colour = factor(organ)),size=2)
}

beta_diversity_organ <- function(study_id,organ,subjects){
    physeq = construct_physeq(study_id)
    physeq <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
    physeq_sub <- prune_samples(sample_data(physeq)$subject%in%subjects,physeq)
    physeq_sub <- prune_samples(sample_data(physeq_sub)$organ==organ, physeq_sub)
    EuclDist_sub = dist(t(otu_table(physeq_sub)))
    EuclDistMatrix_sub <- as.matrix(EuclDist_sub)
    pc_euc_sub <- cmdscale(EuclDistMatrix_sub,k=4)
    write.table(pc_euc_sub,"pc_euc_sub.txt",sep="\t",quote=F, col.names=T, row.names=T)
    pc_euc_sub <- read.table("pc_euc_sub.txt",header=T, row.names=1)
    names(pc_euc_sub) <- c("PC1","PC2","PC3","PC4")
    pc_sample_sub <- sample_data(physeq_sub)
    pc_sample_sub <- cbind(pc_euc_sub,pc_sample_sub)
    par(mar=c(3,3,3,3))
    p <- ggplot(pc_sample_sub, aes(PC1,PC2))
    p + geom_point(size = 2, alpha = 0.3)
    p + geom_point(aes(colour = factor(subject)),size=2)
}
