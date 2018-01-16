#!/usr/bin/env python
from argparse import ArgumentParser
from utils import graphql_api
import json
import shutil
import datetime

study_properties = {
    "study_id": "Study",
}

case_properties = {
    "hiv_status": "HIV Status",
    "id": "ID",
    "case_id": "Case ID",
    "_followups_count": "Number Visits"
}

hiv_history_properties = {
    "hiv_history_id": "HIV History ID",
    "status": "Detailed Status",
    "lnegdate": "Last Seronegative Year",
    "fposdate": "First Seropositive Year", 
    "lastcontact": "Last contact",
    "ehaart": "Ever had HAART",
    "lastnohd": "Last HAART free year",
    "frsthaad": "First HAART year"

}

demographic_properties = {
    "gender": "Gender",
    "race": "Race",
    "ethnicity": "Ethnicity",
    "vital_status": "Vital Status",
    "year_of_death": "Death Year"
}

core_visit_properties = {
    "_summary_drug_uses_count": "Number Drug Records",
    "_summary_lab_results_count": "Number Lab Records"
}

sorted_headers = [
        "id",
        "case_id",
        "study_id",
        "hiv_status",
        "_followups_count",   
        "gender",
        "race",
        "ethnicity",
        "vital_status",
        "year_of_death",
        "hiv_history_id",
        "status",
        "lnegdate",
        "fposdate",
        "lastcontact",
        "ehaart",
        "lastnohd",
        "frsthaad",
        "_summary_lab_results_count",
        "_summary_drug_uses_count"
]

excluded_projects = ['ndh-test']

# Step for query pagination
step = 20

def parse_cmd_args():
    ''' Read arguments '''

    parser = ArgumentParser()
    parser.add_argument('--keys_file',
                        help='File for api authorization',
                        default='/home/ubuntu/.secrets')
    parser.add_argument('--copy_file_to_server',
                        help='copies file to object store',
                        action='store_true')
    parser.add_argument('--output_file',
                        help='tab-separated text output file',
                        default='data.tsv')
    #parser.add_argument('--shiny_server',
    #                    help='path to shiny server location',
    #                    default='/srv/shiny-server/sample-apps/bloodpac-exploration/')

    parser.set_defaults(print_list=False)
    args = parser.parse_args()
    
    return args

def query_data(auth, studies, filename):
    ''' Query data with pagination '''

    data = {}
    for study in studies:
    
        query_txt = """query Cases { study(submitter_id: "%s"){ _cases_count }} """ % study
        counts = graphql_api.query(query_txt, auth)
        counts = counts["data"]["study"][0]['_cases_count']

        offset = 0  
        data.setdefault(study, {})
        while offset < counts:

          query_txt = """{
              study(submitter_id: "%s"){
                 cases(first:%s, offset:%s, order_by_asc: "submitter_id"){
                    %s
                    hiv_history_records(first:0){ 
                       %s
                    }
                    demographics(first:0){
                       %s
                    }
                    followups(first:0){
                       %s
                    }
                }
              }
          }""" % (study, step, offset, 
                 ' '.join(case_properties.keys()).replace('case_id', 'submitter_id'), 
                 ' '.join(hiv_history_properties.keys()).replace('hiv_history_id', 'submitter_id'), 
                 ' '.join(demographic_properties.keys()),
                 ' '.join(core_visit_properties.keys()))

          itime = datetime.datetime.now()
          data_step = graphql_api.query(query_txt, auth)
          
          if not data[study]:
              data[study] = data_step["data"]["study"][0]["cases"]
          else:
              data[study] += data_step["data"]["study"][0]["cases"]

          offset += step
          etime = datetime.datetime.now()
          print "Study %s : Submitted (%s / %s) " % (study, str(offset), counts) + str(etime-itime) 

          if offset % 500 == 0:
              print "Saving partial file in %s" % filename
              name_output = get_data_table(data, filename)

    return data

def get_data_table(data, filename):
    ''' Format query output into tab-separated text '''

    dicts = [study_properties, case_properties, hiv_history_properties, demographic_properties, core_visit_properties]
    headers = {key:val for d in dicts for key,val in d.items()}

    with open(filename, 'w') as outfile:
 
      # Write headers
      h = []
      for header in sorted_headers:
         h.append(headers[header])
      outfile.write("\t".join(h) + "\n")
 
      # Write data
      for st in data:              
        for case in data[st]:
           line = []
           for prop in sorted_headers:
              if prop == 'case_id':
                 line.append(str(case['submitter_id']))
              elif prop == 'study_id':
                 line.append(st)
              elif prop == 'hiv_status':
                 if case[prop] == False:
                    line.append("Negative")
                 else:
                    line.append("Positive")
              elif prop in case:
                 line.append(str(case[prop]))

           if 'demographics' in case and case['demographics']:
               for demo in case['demographics']:  
                   for prop in sorted_headers:
                      if prop in demo:
                         line.append(str(demo[prop]))

           if 'hiv_history_records' in case and case['hiv_history_records']:
               for hist in case['hiv_history_records']:  
                   for prop in sorted_headers:
                      if prop == 'hiv_history_id':
                         line.append(str(hist['submitter_id']))
                      elif prop in hist:
                         line.append(str(hist[prop]))

           if 'followups' in case and case['followups']:
               record_counts = {}
               for prop in sorted_headers:
                   for visit in case['followups']:
                      if prop in visit:
                         record_counts.setdefault(prop,0)
                         record_counts[prop] += visit[prop]
               for prop in record_counts:
                   line.append(str(record_counts[prop]))    
 
           outfile.write("\t".join(line) + "\n")

    return filename


if __name__ == '__main__':

    args = parse_cmd_args()

    auth = graphql_api.get_api_auth(args.keys_file)  
    projects = graphql_api.get_projects(auth, excluded_projects)
    studies = []
    for p in projects:
        studies += graphql_api.get_studies(auth, p)
    
    studies = ['WIHS', 'MACSv2']
    data = query_data(auth, studies, args.output_file)

    filename = get_data_table(data, args.output_file)
