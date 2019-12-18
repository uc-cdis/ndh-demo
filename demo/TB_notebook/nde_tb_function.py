from gen3.auth import Gen3Auth
from gen3.submission import Gen3Submission
from gen3.file import Gen3File
from sh import gunzip
import time
import subprocess
import flatten_json
import warnings
import os
import pandas as pd
import requests
import uuid
import json
import fnmatch
import hashlib
warnings.filterwarnings('ignore')


res2Drug={"Ethambutol":'conferring resistance to ethambutol',"Isoniazid":'conferring resistance to isoniazid',"Pyrazinamide":'conferring resistance to pyrazinamide',"Rifampicin":'conferring resistance to rifampicin'}
#res2Drug={"ethambutol":'ethambutol',"isoniazid":'isoniazid',"pyrazinamide":'pyrazinamide',"rifampicin":'rifampicin'}
# Set up Gen3 SDK
endpoint = "https://tb.niaiddata.org/"
auth = Gen3Auth(endpoint, refresh_file = "/home/jovyan/pd/credentials.json")
sub = Gen3Submission(endpoint, auth)
file = Gen3File(endpoint, auth)

ariba_output_dir = '/home/jovyan/pd/nb_output/tb/ariba/output'
if not os.path.exists(ariba_output_dir):
            os.makedirs(ariba_output_dir)

def query_file(project,chunk,offset,drug_resistant):
    # SDK submission class call peregrine to retrieve data, and convert json to flat json
    params = ""
    for key, value in drug_resistant.items():
        params += ",{}:\"{}\"".format(key, value)
    query = """ 
    {
    subject(project_id:"%s", first:%d, offset:%d, order_by_asc:"submitter_id" %s){
        submitter_id
        samples{
            aliquots{
               read_groups{
                   submitted_unaligned_reads_files{
                       submitter_id
                       file_name
                       object_id
                   }
               } 
            }
        }
    }
    }
    """ %(project,chunk,offset,params)
    subject_res = sub.query(query)
    object_dict = flatten_json.flatten_json(subject_res)
    print(query)
    return object_dict

def parse_json(object_dict,chunk):
    # parse flat json to data frame
    subjects_info = dict()
    for i in range(0,chunk):
        submitter_id = object_dict['data_subject_{0}_submitter_id'.format(i)]
        subjects_info.setdefault(submitter_id,[])
        fastq_submitter_id = object_dict['data_subject_{0}_samples_0_aliquots_0_read_groups_0_submitted_unaligned_reads_files_0_submitter_id'.format(i)]
        fastq_dir = "/home/jovyan/pd/nb_output/tb/fastq_files"
        if not os.path.exists(fastq_dir):
            os.makedirs(fastq_dir)
        filename = fastq_dir + '/' + object_dict['data_subject_{0}_samples_0_aliquots_0_read_groups_0_submitted_unaligned_reads_files_0_file_name'.format(i)]
        object_id = object_dict['data_subject_{0}_samples_0_aliquots_0_read_groups_0_submitted_unaligned_reads_files_0_object_id'.format(i)]
        url = file.get_presigned_url(object_id,protocol="s3")['url']
        if not os.path.isfile(filename):
            print("Downloading %s"%(filename))
            download_file(url,filename)
        else:
            print("%s Exist"%(filename))
        print("************")
        subjects_info[submitter_id].append(filename.split("/")[7])
        subjects_info[submitter_id].append(object_id)
        filename = fastq_dir + '/' + object_dict['data_subject_{0}_samples_0_aliquots_0_read_groups_0_submitted_unaligned_reads_files_1_file_name'.format(i)]
        object_id = object_dict['data_subject_{0}_samples_0_aliquots_0_read_groups_0_submitted_unaligned_reads_files_1_object_id'.format(i)]
        url = file.get_presigned_url(object_id,protocol="s3")['url']
        if not os.path.isfile(filename):
            print("Downloading %s"%(filename))
            download_file(url,filename)
        else:
            print("%s Exist"%(filename))
        print("************")
        sample = filename.split("/")[7]
        sample = sample.split(".")[0]
        sample = sample.split("_")[0]
        subjects_info[submitter_id].append(filename.split("/")[7])
        subjects_info[submitter_id].append(object_id)
        subjects_info[submitter_id].append(sample)
        subjects_info[submitter_id].append(fastq_submitter_id)
    df = pd.DataFrame.from_dict(subjects_info,orient='index')
    # define column names
    df.columns = ['fastq1', 'GUID1','fastq2', 'GUID2', 'SRR',"FASTQ_SI"]
    print(df[["fastq1","GUID1", "fastq2", "GUID2", "SRR"]])
    return df

def download_file(url,filename):
    # call api to download data from url
    req = requests.get(url)
    with open(filename,'wb') as f:
        f.write(req.content)
    
def runAriba(df):
    # loop through each row
    for index, row in df.iterrows():
        file1 = "/home/jovyan/pd/nb_output/tb/fastq_files/" + row['fastq1']
        file2 = "/home/jovyan/pd/nb_output/tb/fastq_files/" + row['fastq2']
        output = ariba_output_dir + "/" + index + ".out"
        print("Processing {}".format(index))
        start = time.time()
        # run ariba command
        subprocess.run(["ariba","run","/home/jovyan/pd/nb_output/tb/ariba/prepareref.out",file1, file2, output])
        end = time.time()
        print("It takes %s sec to complete"%(round(end-start,2)))
        print("****************\n")

def extract_ariba_predict(dir):
    preds = dict()
    if os.path.isdir("/home/jovyan/pd/nb_output/tb/ariba/output/.ipynb_checkpoints"):
        os.rmdir("/home/jovyan/pd/nb_output/tb/ariba/output/.ipynb_checkpoints")
    subfolders = [f.path for f in os.scandir(dir) if f.is_dir()] 
    for p in subfolders:
        subject = p.split("/")[8]
        subject = subject.split(".")[0]
        file = p + "/report.tsv"
        # define file path
        if not os.path.isfile(file):
            file = p + "/" +subject + ".report.tsv"
        else:
            dst = p + "/" +subject + ".report.tsv"
            os.rename(file,dst)
            file = p + "/" +subject + ".report.tsv" 
    cmd = "ariba summary /home/jovyan/pd/nb_output/tb/ariba/out.summary"
    # append file end with report.tsv
    for p in subfolders:
        for file in os.listdir(p):
            if "debug" not in file and file.endswith("report.tsv"):
                cmd += " %s/%s"%(p,file)
    process = subprocess.Popen(cmd, shell=True)
    time.sleep(60)
    for p in subfolders:
        # extract prediction result from each report.tsv and out.summary.csv
        subject = p.split("/")[8]
        subject = subject.split(".")[0]
        file = p + "/" +subject + ".report.tsv" 
        pred = get_prediction_singleSRA(file,"/home/jovyan/pd/nb_output/tb/ariba/out.summary.csv")
        md5sum = md5(file)
        st = os.stat(file)
        pred["md5"] = md5sum
        pred["size"] = st.st_size
        pred["fn"] = subject + ".report.tsv"
        preds[subject] = pred
        df = pd.DataFrame(preds)
    df = df.transpose()
    return df

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def runMykrobe(df):
    # loop through each row
    for index, row in df.iterrows():
        file1 = "/home/jovyan/pd/nb_output/tb/fastq_files/" + row['fastq1']
        file2 = "/home/jovyan/pd/nb_output/tb/fastq_files/" + row['fastq2']
        outfile = "/home/jovyan/pd/nb_output/tb/mykrobeOut" + index + ".out"
        if os.path.isfile(file1) and os.path.isfile(file2):
            print("Processing {}".format(index))
            start = time.time()
            # Run mykrobe for each sample
            subprocess.run(["mykrobe", "predict", index, "tb", "--format", "csv", "-1", file1, file2, "--output", outfile])
            end = time.time()
            print("It takes %s sec to complete"%(round(end-start,2)))
            print("****************\n")

def extract_mykrobe_predict(df):
    preds = dict()
    for index, row in df.iterrows():
        outfile = "mykrobeOut/" + index + ".out"
        pred = get_mykrobe_prediction_singleSRA(outfile)
        pred["md5"] = md5(outfile)
        pred["size"] = os.stat(outfile).st_size
        pred["fn"] = index + ".out"
        preds[index] = pred
    df = pd.DataFrame(preds)
    df = df.transpose()
    return df
        

def get_prediction_singleSRA(ariba_output,ariba_summary):
    #ariba_output='out.run/report.tsv'
    #summary='out_summary.csv'
    results = {}
    df = pd.read_table(ariba_output,  sep="\t")
    df_summary=pd.read_table(ariba_summary,  sep=",")
    n_row=len(df)
    temp={"Ethambutol":'conferring resistance to ethambutol',"Isoniazid":'conferring resistance to isoniazid',"Pyrazinamide":'conferring resistance to pyrazinamide',"Rifampicin":'conferring resistance to rifampicin'}
    for i in range(0,n_row):
        for drug, drugResis in temp.items():
            if drugResis in df["free_text"][i] and (df["cluster"][i]+'.match') in df_summary.columns and len(df_summary.loc[(df_summary['name']==ariba_output)&(df_summary[(df["cluster"][i]+'.match')]=="yes")])==1:
#                 print (drug+":\tR")
                results[drug] = "R"
                del temp[drug]
                break
    
    for drug, drugResis in temp.items():
#         print (drug,':\tS')
        results[drug] = "S"
    return results

    
def get_mykrobe_prediction_singleSRA(mykrobe_output):
     #mykrobe_output='mykrobeOut/result_'+sra+'.csv'
    results = {}
    df = pd.read_csv(mykrobe_output,  sep=",")
    n_row=len(df)
    df=df.sort_values(by=['susceptibility'])
    df= df.reset_index(drop=True)
    for i in range(0,n_row):
#         print (df['drug'][i]+":\t"+df['susceptibility'][i])
        results[df['drug'][i]] = df['susceptibility'][i]
    return results
        
def queryReads(subject):
    # Query for fastq file submitter_id for each subject
    query = '''
    {
        sample(with_path_to:{type:"subject",submitter_id:"%s"}){
            aliquots{
                read_groups{
                    submitted_unaligned_reads_files{
                        submitter_id
                    }
                }
            }
        }
    }
    ''' % subject
    # extract submitted_unaligned_reads submitter_id for that subject
    sample_res = sub.query(query)
    object_dict = flatten_json.flatten_json(sample_res)
    fq1 =object_dict['data_sample_0_aliquots_0_read_groups_0_submitted_unaligned_reads_files_0_submitter_id']
    fq2 =object_dict['data_sample_0_aliquots_0_read_groups_0_submitted_unaligned_reads_files_1_submitter_id']
    return [fq1,fq2]

def submit_results(df,tool):
    # loop through result data frame to submit result to file node sequencing result
    for index, row in df.iteritems():
        pr = queryReads(index)
        submit_id = index + "_seq_res_" + tool.lower()
        d = '''
        [{
          "type": "sequencing_result",
        '''
        # update drug resistant name based on tool that predict the drug resistance
        if tool=="Ariba":
            drugs = ["Ethambutol","Isoniazid","Pyrazinamide","Rifampicin"]
            for drug in drugs:
                if row[drug] == "R":
                    result = "Resistant"
                else:
                    result = "Susceptible"
                aribaDrug = '  "' + str(drug).lower() + '_res_ariba": ' '"' + result + '",' 
                d = d + aribaDrug + '\n'
        if tool == "Mykrobe":
            drugs = ["Amikacin","Capreomycin","Ciprofloxacin","Ethambutol","Isoniazid","Kanamycin","Moxifloxacin","Ofloxacin","Pyrazinamide","Rifampicin","Streptomycin"]
            for drug in drugs:
                if row[drug] == "R":
                    result = "Resistant"
                else:
                    result = "Susceptible"
                mykrobeDrug = '  "' + str(drug).lower() + '_res_mykrobe": ' '"' + result + '",' 
                d = d + mykrobeDrug + '\n'
        # Add other information for submission
        d = d + '''  "experimental_strategy": "Sequencing Experiment",
          "data_type": "cfDNA Sequencing Results",
          "data_format": "TXT",
          "data_category": "Microbial Analysis",
          "file_name": "%s",
          "file_size": %d,
          "md5sum": "%s",
          "submitter_id": "%s",
          "submitted_unaligned_reads_files":[
          {"submitter_id":"%s"},
          {"submitter_id":"%s"}]}]'''%(row["fn"],row["size"],row["md5"],submit_id,pr[0],pr[1])
        sr_json = json.loads(d)
        res = sub.submit_record("TB","PATRIC",sr_json)
        # Print submission success code
        if json.loads(res)['code']==200:
            print("%s is successfully submitted"%(json.loads(res)["entities"][0]['unique_keys'][0]['submitter_id']))
