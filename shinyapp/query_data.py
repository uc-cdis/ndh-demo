#!/usr/bin/env python
from argparse import ArgumentParser
from utils import graphql_api
import json
import shutil
import datetime

study_properties = {
    "study_id": "Study",
    "project_id": "Project"
}

case_properties = {
    "id": "ID",
    "case_id": "Case ID",
    "species": "Species",
    "_follow_ups_count": "Number Visits"
}

demographic_properties = {
    "gender": "Gender",
    "race": "Race",
    "ethnicity": "Ethnicity",
    "vital_status": "Vital Status",
    "year_of_death": "Death Year"
}

follow_up_properties = {
    "_summary_drug_uses_count": "Number Drug Records",
    "_summary_lab_results_count": "Number Lab Records"
}

sample_properties = {
    "_summary_lab_results_count": "Number Lab Records"
}

aliquot_properties = {
    "_mRNA_microarrays_count": "mRNA Array Records"
}

sorted_headers = [
        "id",
        "case_id",
        "project_id",
        "study_id",
        "species",
        "_follow_ups_count",        
        "gender",
        "race",
        "ethnicity",
        "vital_status",
        "year_of_death",
        "_summary_lab_results_count",
        "_summary_drug_uses_count",
        "_mRNA_microarrays_count"
]

excluded_projects = ['ndh-test', 'ndh-dmid-LHV']

# Step for query pagination
step = 20

def parse_cmd_args():
    ''' Read arguments '''

    parser = ArgumentParser()

    parser.add_argument('--keys_file',
                        help='Credentials for the API',
                        default='/srv/shiny-server/ndh/utils/credentials.json')
    parser.add_argument('--api_url',
                        help='URL for the API',
                        default='https://niaid.bionimbus.org/')
    parser.add_argument('--output_file',
                        help='tab-separated text output file',
                        default='data.tsv')

    parser.set_defaults(print_list=False)
    args = parser.parse_args()
    
    return args

def query_data(auth, studies, filename, api_url, keys_file):
    ''' Query data with pagination '''

    data = {}
    for project in projects:
    
        query_txt = """query Cases { _subject_count(project_id: "%s") }""" % project
        counts = graphql_api.query(query_txt, auth, api_url)
        counts = counts["data"]["_subject_count"]
           
        offset = 0  
        data.setdefault(project, {})
        while offset < counts:

          query_txt = """{
              subject(first:%s, offset:%s, order_by_asc: "submitter_id", project_id: "%s"){
                    %s
                    studies(first:0){
                       %s
                    }
                    demographics(first:0){
                       %s
                    }
                    follow_ups(first:0){
                       %s
                    }
                    samples(first:0){
                       %s
                       aliquots(first:0){
                          %s
                       }
                    }                    
                }
              }""" % (step, offset, project,
                 ' '.join(case_properties.keys()).replace('case_id', 'submitter_id'),
                 ' '.join(study_properties.keys()).replace('study_id', 'submitter_id'),  
                 ' '.join(demographic_properties.keys()),
                 ' '.join(follow_up_properties.keys()),
                 ' '.join(sample_properties.keys()),
                 ' '.join(aliquot_properties.keys()))

          itime = datetime.datetime.now()
          data_step = graphql_api.query(query_txt, auth, api_url)
          
          if not data[project]:
              data[project] = data_step["data"]["subject"]
          else:
              data[project] += data_step["data"]["subject"]

          offset += step
          etime = datetime.datetime.now()
          print "Study %s : Submitted (%s / %s) " % (project, str(offset), counts) + str(etime-itime) 

          if offset % 500 == 0:
              print "Saving partial file in %s" % filename
              name_output = get_data_table(data, filename)
              auth = graphql_api.get_api_auth(keys_file, api_url)

    return data

def get_data_table(data, filename):
    ''' Format query output into tab-separated text '''

    dicts = [study_properties, case_properties, demographic_properties, follow_up_properties, sample_properties, aliquot_properties]
    headers = {key:val for d in dicts for key,val in d.items()}

    with open(filename, 'w') as outfile:
 
      # Write headers
      h = []
      for header in sorted_headers:
         h.append(headers[header])
      outfile.write("\t".join(h) + "\n")
 
      # Write data
      for pr in data:              
        for case in data[pr]:
           line = []
           for prop in sorted_headers:
              if prop == 'case_id':
                 line.append(str(case['submitter_id']))
              elif prop == 'project_id':
                 line.append(pr) 
              elif prop == 'study_id':
                 line.append(case['studies'][0]['submitter_id'])              
              elif prop in case:
                 line.append(str(case[prop]))

           if 'demographics' in case and case['demographics']:
               for demo in case['demographics']:  
                   for prop in sorted_headers:
                      if prop in demo:
                         line.append(str(demo[prop]))
           else:
               for value in range(len(demographic_properties)):
                  line.append('None')

           record_counts = {}
           if 'samples' in case and case['samples']:
               for prop in sorted_headers:
                   for sample in case['samples']:
                      if prop in sample:
                         record_counts.setdefault(prop,0)
                         record_counts[prop] += sample[prop]
                      if 'aliquots' in sample and sample['aliquots']:
                         for aliq in sample['aliquots']:
                            if prop in aliq:
                               record_counts.setdefault(prop,0)
                               record_counts[prop] += aliq[prop]

           if 'follow_ups' in case and case['follow_ups']:
               for prop in sorted_headers:
                   for visit in case['follow_ups']:
                      if prop in visit:
                         record_counts.setdefault(prop,0)
                         record_counts[prop] += visit[prop]

           # Add all counts in file
           for prop in sorted_headers:
              counts_headers = follow_up_properties.keys() + sample_properties.keys() + aliquot_properties.keys()
              if prop in record_counts:
                  line.append(str(record_counts[prop]))
              elif prop in counts_headers:
                  line.append('0')

           outfile.write("\t".join(line) + "\n")

    return filename


if __name__ == '__main__':

    args = parse_cmd_args()

    auth = graphql_api.get_api_auth(args.keys_file, args.api_url)  
    projects = graphql_api.get_projects(auth, args.api_url, excluded_projects)  
    data = query_data(auth, projects, args.output_file, args.api_url, args.keys_file)

    filename = get_data_table(data, args.output_file)
