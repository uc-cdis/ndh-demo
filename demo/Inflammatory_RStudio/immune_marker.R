library(jsonlite)
library(dplyr)
source("./Gen3AuthHelper.R")
source("./Gen3Submission.R")
options(warn=-1)

#################################
#
# Functions
#
################################

# Query demographic propeties for subject using the sub instance
query_demo<-function(subject_id){
  query_txt = paste('{subject(id:\\"',subject_id,'\\"){demographics{gender,race,ethnicity},comorbidities{bshcvstat},exposures{ageatbl,baseidu},visits(visit_type:\\"Baseline Visit\\"){summary_socio_demographics{smoke_status}}}}',sep="")
  output = sub$query(query_txt)
  data = content(output)$data$subject
  datalist = list()
  if (length(data[[1]]$comorbidities[[1]]$bshcvstat)!=0){
    baseline_hcv = data[[1]]$comorbidities[[1]]$bshcvstat
  }else{
    baseline_hcv = ""
  }
  if(length(data[[1]]$exposures)!=0){
    baseline_age = data[[1]]$exposures[[1]]$ageatbl
    baseline_idu = data[[1]]$exposures[[1]]$baseidu
  }else{
    baseline_age = ""
    baseline_idu = ""
  }
  for(i in 1:length(data[[1]]$demographics[[1]]$race)){
    if(length(data[[1]]$visits)!=0 && length(data[[1]]$visits[[1]]$summary_socio_demographics)!=0){
      each_subject = data.frame(subject = subject_id, baseline_hcv = baseline_hcv,ethnicity=data[[1]]$demographics[[1]]$ethnicity,gender=data[[1]]$demographics[[1]]$gender,race=data[[1]]$demographics[[1]]$race[[i]],baseline_age = baseline_age,baseidu=baseline_idu,smoke_status = data[[1]]$visits[[1]]$summary_socio_demographics[[1]]$smoke_status)
    }else{
      each_subject = data.frame(subject = subject_id, baseline_hcv = baseline_hcv,ethnicity=data[[1]]$demographics[[1]]$ethnicity,gender=data[[1]]$demographics[[1]]$gender,race=data[[1]]$demographics[[1]]$race[[i]],baseline_age = baseline_age,baseidu=baseline_idu,smoke_status = "")
    }
    datalist[[i]]= each_subject
  }
  return(datalist)
}

# Create dataframe for each of the three categories
query_by_category<-function(subjects){
  data_by_category = list()
  for(subject in subjects){
    each_subject = query_demo(subject)
    data_by_category = append(data_by_category,each_subject)
  }
  target_matrix = do.call(rbind,data_by_category)
  return(target_matrix)
}

# Calculate summary percentage/average table for each of the three categories
summary_table<-function(target_matrixs,enum_variables,enum_values,categories,continuous_variables){
  enum_by_category<-list()
  continuous_by_category<-list()
  for (i in 1:length(enum_variables)){
    percentage_by_category<-list()
    variable = enum_variables[[i]]
    value = enum_values[[i]]
    for (j in 1:length(target_matrixs)){
      target_matrix = target_matrixs[[j]]
      target_matrix = distinct(target_matrix[,c("subject",variable)])
      summary_table = transform(as.data.frame(table(target_matrix[,variable])),percentage_column=Freq/nrow(target_matrix)*100)
      percentage = summary_table[summary_table$Var1==value,"percentage_column"]
      percentage = format(round(percentage, 2), nsmall = 2)
      if (length(percentage)==0){
        percentage = 0
      }
      df <- data.frame(percentage)
      names(df) <- categories[[j]]
      row.names(df)<-paste(variable,"_%",value,sep="")
      percentage_by_category[[j]] = df
    }
    category_matrix = do.call(cbind,percentage_by_category)
    enum_by_category[[i]] = category_matrix
  }
  for (i in 1:length(continuous_variables)){
    mean_continuous = list()
    variable = continuous_variables[[i]]
    for (j in 1:length(target_matrixs)){
      target_matrix = target_matrixs[[j]]
      if(all(is.na(unique(target_matrix[,variable][target_matrix[,variable]!= ""])))){
        df<-data.frame("<NA>")
      }else{
        vavg = mean(as.numeric(target_matrix[,variable]),na.rm = TRUE)
        vavg = format(round(vavg, 2), nsmall = 2)
        df<-data.frame(vavg)
      }
      names(df) <- categories[[j]]
      row.names(df)<- paste("average",continuous_variables[[i]],sep="_")
      mean_continuous[[j]]=df
    }
    continuous_matrix = do.call(cbind,mean_continuous)
    continuous_by_category[[i]] = continuous_matrix
  }
  variable_matrix = do.call(rbind,append(enum_by_category,continuous_by_category))
  return(variable_matrix)
}

# Get the first visit of a specific vist year
search_oneyear_after_visit<-function(subject){
  one_year_after = subject$first_year_qualify_ltnp + 1
  follow_ups = subject$follow_ups
  for (row in 1:nrow(follow_ups)){
    if (follow_ups[row,"visit_date"] ==  one_year_after){
      return (follow_ups[row,"submitter_id"])
    }
  }
  return("NA")
}

# Get the frist visit that qualify EC criterial
search_qualify_visit<-function(subject){
  i = 0
  qualify_year=0
  ec_period = subject$ecPeriod.ec_perid_1
  for (row in 1:nrow(ec_period)){
    if (!is.na(ec_period[row,"viral_load"])&& ec_period[row,"viral_load"] < as.numeric(sup_upper_bound)){
      i = i + 1
      if(i == 2){
        qualify_visit = ec_period[row,"submitter_id"]
        qualify_year = ec_period[row,"visit_date"]
      }
    }
    if (ec_period[row,"visit_date"] == qualify_year + 1){
      one_year_after_visit = ec_period[row,"submitter_id"]
      my_list = list("qualify_visit" = qualify_visit,"one_year_after_visit" = one_year_after_visit)
      return (my_list)
    }
  }
  my_list = list("qualify_visit" = qualify_visit,"one_year_after_visit" = "NA")
  return(my_list)
}

# Data availability at critical timepoint for inflammatory markers
data_availability_percentage<-function(record, record_name, case_number, markers){
  percentage_by_marker = list()
  for (i in 1:length(markers)){
    marker = markers[[i]]
    percentage = length(record[!is.na(record[,marker]),marker])/case_number
    df <- data.frame(percentage)
    names(df) <-record_name
    row.names(df)<-marker
    percentage_by_marker[[i]] = df
  }
  return(do.call(rbind,percentage_by_marker))
}

data_availability_matrix<-function(records,records_name,case_numbers, markers){
  percentage_matrix = list()
  for (i in 1:length(records)){
    percentage_dataframe = data_availability_percentage(records[[i]],records_name[[i]],case_numbers[[1]],markers)
    percentage_matrix[[i]] = percentage_dataframe
  }
  return(do.call(cbind,percentage_matrix))
}

# Box plots to compare inflammatory markers among different cohorts
box_plot_markers<-function(plot_markers,records,record_names){
  par(mfrow = c(2,4),mar=c(2,2,2,2))
  for (marker in plot_markers){
    inflammatory_marker_values = list()
    for (i in 1:length(records)){
      record = records[[i]]
      df = data.frame(value = record[,marker])
      df$visit = rep(record_names[[i]],nrow(record))
      inflammatory_marker_values[[i]] = df
    }
    combine_data = do.call(rbind,inflammatory_marker_values)
    boxplot(value~visit, data=combine_data,names=c("ptc_on","ptc_off","ec","ltnp"),cex.axis=0.5)
    title(main=marker)
    stripchart(value~visit, vertical = TRUE, data = combine_data, method = "jitter", add = TRUE, pch = 20, col = 'blue')
  }
}


#################################
#
#Input required for running functions
#
#################################

# Please specify the json file path export from HIV cohort selection App
PTC_json<-"~/Documents/NIAID/Charlie/LTNP_EC_PTC_analysis/ptc-cohort-vload-400-months-24.json"
LTNP_json<-"~/Documents/NIAID/Charlie/LTNP_EC_PTC_analysis/ltnp-cohort-CD4-500-years-5.json"
EC_json<-"~/Documents/NIAID/Charlie/LTNP_EC_PTC_analysis/ec-cohort-suppressvload-50-spikevload-1000-visits-2.json"

# Specify data common end point, program, project, nodetype, immune_markers for data availability table and box plot.
endpoint<-"https://aids.niaiddata.org"
program = "HIV"
project = "CHARLIE"
node_type = "summary_lab_result"
inflammatory_markers = list("baff","bca1","crp","eotaxin","gmcsf","gp130","ifng","il10","il12","il17","il1b","il2","il4","il6","il8","mcp1","sil2ra","tnfa")
plot_markers = list("il2","il4","il6","ifng","gp130","gmcsf","bca1","tnfa")


#################################
#
#Run functions
#
#################################

# Create sub instance for the submission Class defined in Gen3 R SDK to export or query data from data common
auth <- Gen3AuthHelper(endpoint, refresh_file="credentials.json")
sub <- Gen3Submission(endpoint, auth)

# Read json file to identify subjects for each category
PTC<-fromJSON(PTC_json)
PTC_subjects <- PTC$subjects$subject_id
LTNP<-fromJSON(LTNP_json)
LTNP_subjects <- LTNP$subjects$subject_id
EC<-fromJSON(EC_json)
EC_subjects <- EC$subjects$subject_id

# Create data frame for each category
PTC_dataframe<-query_by_category(PTC_subjects)
LTNP_dataframe<-query_by_category(LTNP_subjects)
EC_dataframe<-query_by_category(EC_subjects)

# Data availability at critical time point
# Extract visit submitter_id at critical time point for PTC cases
PTC_haart_end_visit<-PTC$subjects$consecutive_haart_treatments_end_at_followup
PTC_maintain_visit<-PTC$subjects$stop_treatments_maintain_viral_load_at_followup
PTC_critical_visit<-data.frame(subject_id = PTC_subjects, PTC_haart_end_visit = PTC_haart_end_visit, PTC_maintain_visit = PTC_maintain_visit)

# Extract visit submitter_id at critical time point for LTNP cases
LTNP_first_hiv_visit<-LTNP$subjects$first_hiv_positive_visit
LTNP_first_qualify_visit<-LTNP$subjects$first_visit_qualify_ltnp
LTNP_first_qualify_year<-LTNP$subjects$first_year_qualify_ltnp
LTNP_qualify_oneyear_visit<-apply(LTNP$subjects,1,search_oneyear_after_visit)
LTNP_critical_visit<-data.frame(LTNP_subjects = LTNP_subjects, LTNP_first_hiv_visit = LTNP_first_hiv_visit, LTNP_first_qualify_visit = LTNP_first_qualify_visit,LTNP_qualify_oneyear_visit = LTNP_qualify_oneyear_visit)

# Extract visit submitter_id at critical time point for EC cases
EC_first_hiv_visit<-EC$subjects$first_hiv_positive_visit
sup_upper_bound <- EC$viral_load_sup_upper_bound
EC_qualify = apply(EC$subjects,1,search_qualify_visit)
EC_first_qualify_visit = sapply(EC_qualify,function(l) l$qualify_visit)
EC_qualify_oneyear_visit = sapply(EC_qualify,function(l) l$one_year_after_visit)

EC_critical_visit<-data.frame(EC_subjects = EC_subjects, EC_first_hiv_visit = EC_first_hiv_visit, EC_first_qualify_visit = EC_first_qualify_visit, EC_qualify_oneyear_visit = EC_qualify_oneyear_visit)

# export summary_lab_result node
filename = paste(tolower(project),node_type,sep="_")
filename = paste(filename,"tsv",sep=".")
sub$export_node(program,project,node_type,"tsv",filename)

#subset immune marker for visits at critical time point for three cohorts (LTNP, EC, PTC)
lab_result = read.table(filename,header=T, sep="\t")
#LTNP cohort
num_ltnp_cases = dim(LTNP_critical_visit)[1]
ltnp_first_hiv_visit_records<-merge(LTNP_critical_visit,lab_result,by.x="LTNP_first_hiv_visit",by.y="visits.submitter_id")
ltnp_first_qualify_visit_records<-merge(LTNP_critical_visit,lab_result,by.x="LTNP_first_qualify_visit",by.y="visits.submitter_id")
ltnp_qualify_oneyear_visit_records<-merge(LTNP_critical_visit,lab_result,by.x="LTNP_qualify_oneyear_visit",by.y="visits.submitter_id")

#EC cohort
num_ec_cases = dim(EC_critical_visit)[1]
ec_first_hiv_visit_records<-merge(EC_critical_visit,lab_result, by.x = "EC_first_hiv_visit", by.y="visits.submitter_id")
ec_first_qualify_visit_records<-merge(EC_critical_visit,lab_result, by.x = "EC_first_qualify_visit", by.y="visits.submitter_id")
ec_qualify_oneyear_visit_records<-merge(EC_critical_visit,lab_result, by.x = "EC_qualify_oneyear_visit", by.y="visits.submitter_id")

#PTC cohort
num_ptc_cases = dim(PTC_critical_visit)[1]
ptc_haart_end_visit_records<-merge(PTC_critical_visit,lab_result, by.x = "PTC_haart_end_visit", by.y="visits.submitter_id")
ptc_maintain_visit_records<-merge(PTC_critical_visit,lab_result, by.x = "PTC_maintain_visit",by.y="visits.submitter_id")

#################################
#
#Output tables or plots
#
#################################

# Specify the argument, arguments are configurable if user prefer show different cohorts, enumerate attributes, enumerate attribute values, column names for the table and numerical attributes.

target_matrixs = list(EC_dataframe,PTC_dataframe,LTNP_dataframe)
enum_variables = list("baseline_hcv","ethnicity","race","race","smoke_status","gender","baseidu")
enum_values = list("HCV positive","Hispanic or Latino","White","Black or African American", "Not currently smoking","female",TRUE)
categories = list("EC","PTC","LTNP")
continuous_variables = list("baseline_age")

summary_table(target_matrixs,enum_variables,enum_values,categories,continuous_variables)

# Show data availability for different cohorts at critical time points. The arguments are configurable if user want to show different cohorts at critical time points.

records = list(ltnp_first_hiv_visit_records,ltnp_first_qualify_visit_records,ltnp_qualify_oneyear_visit_records,ec_first_hiv_visit_records,ec_first_qualify_visit_records,ec_qualify_oneyear_visit_records,ptc_haart_end_visit_records,ptc_maintain_visit_records)

records_name = list("ltnp_first_hiv_visit","ltnp_first_qualify_visit","ltnp_qualify_oneyear_visit","ec_first_hiv_visit","ec_first_qualify_visit","ec_qualify_oneyear_visit","ptc_haart_end_visit","ptc_maintain_visit")

case_numbers = list(num_ltnp_cases,num_ltnp_cases,num_ltnp_cases,num_ec_cases,num_ec_cases,num_ec_cases,num_ptc_cases,num_ptc_cases)

data_availability_matrix(records,records_name,case_numbers, inflammatory_markers)

# Box plots to compare inflammatory markers for different cohorts. Arguments are configurable if user want to compare different cohorts at critcal time points
plot_records = list(ptc_haart_end_visit_records,ptc_maintain_visit_records,ec_first_qualify_visit_records,ltnp_first_qualify_visit_records)
plot_names = list("ptc_haart_end_visit","ptc_maintain_visit","ec_first_qualify_visit","ltnp_first_qualify_visit")
box_plot_markers(plot_markers,plot_records,plot_names)
