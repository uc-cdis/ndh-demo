from cdispyutils.hmac4 import get_auth
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from lifelines.datasets import load_waltons
from lifelines import KaplanMeierFitter
from scipy.stats import ranksums
from operator import add
import numpy as np
import os.path
import datetime
import requests
import glob
import json
import os

import warnings
warnings.filterwarnings('ignore')

summary_order = [
   "_study_count",
   "_subject_count",
   "_demographic_count",
   "_hiv_history_count",
   "_follow_up_count",
   "_summary_socio_demographic_count",    
   "_summary_lab_result_count",
   "_summary_drug_use_count"
]

summary_count_headers = {
    "_subject_count": "Cases",
    "_study_count": "Studies",
    "_demographic_count": "Demographic records",
   "_hiv_history_count": "HIV History records",
   "_follow_up_count": "Visit records",
   "_summary_lab_result_count": "Lab Results records",
   "_summary_drug_use_count": "AIDS Drug records",
   "_summary_socio_demographic_count": "Socio-Demographic records"
}

excluded_studies = ['study-01']

chunk = 50

class MetricsTable(dict):
    ''' Represent metrics tables in HTML format for visualization '''
 
    def _repr_html_(self):
        html = []
        html.append("<table style>")
        html.append("<thead>")  
        html.append("<th>Metric</th>")  
        html.append("<th>Value</th>")
        html.append("</thead>")       
        for key in self:
            html.append("<tr>") 
            html.append("<td>%s</td>" % key)             
            html.append("<td>%s</td>" % self[key])           
            html.append("<tr>") 
        html.append("</table>")        
        
        return ''.join(html)

class SummaryTable(dict):
    ''' Represent result tables in HTML format for visualization '''
 
    def _repr_html_(self):
        html = []
        html.append("<table style>")
        html.append("<thead>")  
        html.append("<th>Category</th>")  
        html.append("<th>Counts</th>")
        html.append("</thead>")       
        for key in summary_order:
            html.append("<tr>") 
            html.append("<td>%s</td>" % summary_count_headers[key])             
            html.append("<td>%s</td>" % self[key])           
            html.append("<tr>") 
        html.append("</table>")        
        
        return ''.join(html)

def add_keys(filename):
    ''' Get auth from our secret keys '''

    global auth
    json_data = open(filename).read()
    keys = json.loads(json_data)
    auth = requests.post('https://niaid.bionimbus.org/user/credentials/cdis/access_token', json=keys)

def query_api(query_txt, variables=None):
    ''' Request results for a specific query '''

    if variables == None:
        query = {'query': query_txt}
    else:
        query = {'query': query_txt, 'variables': variables}

    output = requests.post('https://niaid.bionimbus.org/api/v0/submission/graphql', headers={'Authorization': 'bearer ' + auth.json()['access_token']}, json=query).text
    data = json.loads(output)

    if 'errors' in data:
        print(data)

    return data     

def get_studies(project_id):
    ''' Get list of studies for specific project'''    

    query_txt = """{ study(project_id: "%s"){ submitter_id }}""" % project_id
    data = query_api(query_txt) 
    
    studies = []
    for study in data['data']['study']:
        studies.append(study['submitter_id'])
    
    return studies    
    
def query_summary_counts(project_id):
    ''' Query summary counts for each data type'''
   
    query_txt = "query Counts ($projectID: [String]) { "
    for prop in summary_order:
        query_txt += "%s(project_id: $projectID) " % prop 
    query_txt += "}"
    
    variables = { 'projectID': project_id }

    data = query_api(query_txt, variables) 
    
    table = SummaryTable(data['data'])
    
    return table

def query_summary_field(project, node, field, study_id=None):
    ''' Query summary counts for each data type '''
   
    count_query = """{ _%s_count(project_id:"%s") }""" % (node, project)
    counts = query_api(count_query)['data']['_%s_count' % node]
    offset = 0      
    chunk = 1000
    
    data = {}
    errors = {}
    while offset <= counts:
        query_txt = """query { %s(project_id: "%s", first:%d, offset:%d, order_by_asc: "submitter_id") {%s studies{submitter_id}}} """ % (node, project, chunk, offset, field)  
        output = query_api(query_txt)
        if not data:
            data = output
        else:
            data['data']['subject'] = data['data']['subject'] + output['data']['subject']
        offset += chunk  
      
    summary = {}
    total = []
    for d in data['data'][node]:

        if isinstance(d[field], float):
            d[field] = str(d[field])[:-2]        
        
        if 'studies' in d:
            study = d['studies'][0]['submitter_id']
            if study not in excluded_studies:          
                summary.setdefault(study, {})
                summary[study].setdefault(d[field], 0)
                summary[study][d[field]] += 1
                if d[field] not in total:
                    total.append(d[field])            
        else:
            summary.setdefault(d[field], 0)        
            summary[d[field]] += 1
            
    if study_id != None:
        plot_field_metrics(summary, field)
    else:
        plot_overall_metrics(summary, field, total)       
    
    return summary


def plot_overall_metrics(summary_counts, field, totals):    
    ''' Visualize summary results across projects in a barplot ''' 
    
    results = {}
    projects = {}
    for project in summary_counts:
        
        results[project] = []
        projects.setdefault(project, 0)
            
        for value in totals:
            if value in summary_counts[project]:
                results[project].append(summary_counts[project][value])
                projects[project] += summary_counts[project][value]
            else:
                results[project].append(0)

    N = len(totals)
    positions = np.arange(N) 
    sorted_projects = sorted(projects, key=projects.get, reverse=True)
    bar_size = 0.2
    size_prop = (N/10) + 1
    
    plots = []
    plt.figure(figsize=(8, 4))
    left = [0]*N
    for pr in sorted_projects:
        p = plt.barh(positions, results[pr], bar_size, left, align='center', alpha=0.7)        
        plots.append(p[0])
        left = map(add, left, results[pr])
        
    plt.title('Summary counts by (' + field + ')', fontsize=10*size_prop)
    plt.xlabel('COUNTS', fontsize=10*size_prop)    
    plt.xlim(0, max(left)+5)
    plt.ylabel(field.upper(), fontsize=10*size_prop)  
    plt.yticks(positions, totals, fontsize=10*size_prop)    
    plt.legend(plots, sorted_projects, fontsize=10*size_prop)
           
    plt.show()      
    
    
def compare_lab_results(project_id, variable, tag=None):
    ''' Compare any lab result variable after seropositive conversion'''
        
    filename = '%s.json' % variable
    if os.path.isfile(filename):
        json_data=open(filename).read()
        values = json.loads(json_data)
    else:
        values = {}
        values['Last Seronegative'] = []
        values['First Seropositive'] = []
        studies = get_studies(project_id)
        for study in studies:

            count_query = """{ study(submitter_id:"%s"){ _subjects_count }}""" % study
            counts = query_api(count_query)['data']['study'][0]['_subjects_count']
            offset = 0
            chunk = 50

            data = {}
            errors = {}
            while offset <= counts:
                itime = datetime.datetime.now()
                query_txt = """{ study(submitter_id:"%s"){subjects(first:%d, offset:%d, order_by_asc: "submitter_id"){ submitter_id hiv_history_records{posvis negvis} follow_ups(first:0){
                                           submitter_id 
                                           summary_lab_results{%s}
                                        }}}}""" % (study, chunk, offset, variable)
                output = query_api(query_txt)
                if not data:
                    data = output['data']['study'][0]['subjects']
                else:
                    data = data + output['data']['study'][0]['subjects']
                offset += chunk
                
                etime = datetime.datetime.now()
                print("Query (%s) %s" % (offset, str(etime-itime)))

            if tag == None:
                tag = variable

            for c in data:
                case = c['submitter_id']
                negvis = c['hiv_history_records'][0]['negvis']
                posvis = c['hiv_history_records'][0]['posvis']

                if negvis > 0:
                    visit_id = case + "_" + str(negvis)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    negvalue = visit['summary_lab_results'][0][variable]
                                    if negvalue != None:
                                        values['Last Seronegative'].append(negvalue)              

                if posvis > 0:
                    visit_id = case + "_" + str(posvis)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    posvalue = visit['summary_lab_results'][0][variable]
                                    if posvalue != None:
                                        values['First Seropositive'].append(posvalue)                             

        with open(filename, 'w') as fp:
            json.dump(values, fp)
    
    pvalues = run_statistical_test(values, 'Last Seronegative', 'First Seropositive')
         
    compare_boxplot(values, 'Comparison Seronegative/Seropositive', tag, pvalues)
    
    return values


def compare_survival(project_id, variable, tag=None):
    ''' Compare survival with haart vs not ehaart'''
    
    filename = '%s_survival.json' % variable
    if os.path.isfile(filename):
        json_data=open(filename).read()
        values = json.loads(json_data)
    else:
        values = {}
        studies = get_studies(project_id)
        for study in studies:

            count_query = """{ study(submitter_id:"%s"){ _subjects_count }}""" % study
            counts = query_api(count_query)['data']['study'][0]['_subjects_count']
            offset = 0
            chunk = 100

            data = {}
            errors = {}
            while offset <= counts:
                itime = datetime.datetime.now()
                query_txt = """{ study(submitter_id:"%s"){
                                  subjects(first:%d, offset:%d, order_by_asc: "submitter_id"){ 
                                     submitter_id 
                                     hiv_history_records{%s fposdate} 
                                     demographics{vital_status year_of_death}
                                  }
                                }}""" % (study, chunk, offset, variable)
                output = query_api(query_txt)
                if not data:
                    data = output['data']['study'][0]['subjects']
                else:
                    data = data + output['data']['study'][0]['subjects']
                offset += chunk
                
                etime = datetime.datetime.now()
                print("Query (%s) %s" % (offset, str(etime-itime)))

            if tag == None:
                tag = variable

            for c in data:
                case = c['submitter_id']
                ehaart = c['hiv_history_records'][0][variable]
                vital_status = c['demographics'][0]['vital_status']
                death_year = c['demographics'][0]['year_of_death']
                fposdate = c['hiv_history_records'][0]['fposdate']
                
                if ehaart != None and vital_status != None and fposdate != None:
                    values.setdefault(variable,[])
                    values[variable].append(ehaart)
                    
                    values.setdefault('vital_status',[])
                    values['vital_status'].append(int(vital_status == 'Alive'))
                    
                    values.setdefault('death_year',[])
                    if death_year != None and death_year != 9000:
                        values['death_year'].append(death_year - fposdate)
                    else:
                        values['death_year'].append(2017 - fposdate)
        
        with open(filename, 'w') as fp:
            json.dump(values, fp)
    
    # Prepare time for the two compared groups
    t = np.array(values['death_year'])
    t = t - min(t) + 1
    times = t[t != max(t)]
    fix = np.array(values[variable]) 
    ix = fix[t != max(t)]   
    censors = np.array([1]*len(times))                    
    
    # Plot Kaplan-Meier curve
    kmf = KaplanMeierFitter()    
    kmf.fit(times[ix], censors[ix], label='HAART Treated')
    ax = kmf.plot()
    kmf.fit(times[~ix], censors[~ix], label='Non HAART Treated') 
    kmf.plot(ax=ax)
    ax.set_title(tag)
    ax.set_xlabel("Survival time since seen seropositive (years)")
    
    return times               
    
    
def compare_after_haart(project_id, variable, tag=None):
    ''' Compare any lab result variable after seropositive conversion'''
    
    filename = '%s_haart.json' % variable
    if os.path.isfile(filename):
        json_data=open(filename).read()
        values = json.loads(json_data)
    else:
        values = {'Last HAART Free': [], 'First HAART visit': [], '1 Year Treatment': [], '2 Years Treatment': [], '3 Years Treatment': []}
        studies = get_studies(project_id)
        for study in studies:

            count_query = """{ study(submitter_id:"%s"){ _subjects_count }}""" % study
            counts = query_api(count_query)['data']['study'][0]['_subjects_count']
            offset = 0
            chunk = 50

            data = {}
            errors = {}
            while offset <= counts:
                itime = datetime.datetime.now()
                query_txt = """{ study(submitter_id:"%s"){subjects(first:%d, offset:%d, order_by_asc: "submitter_id"){ submitter_id hiv_history_records{frsthaav lastnohv} follow_ups(first:0){
                                           submitter_id 
                                           summary_lab_results{%s}
                                        }}}}""" % (study, chunk, offset, variable)
                output = query_api(query_txt)
                if not data:
                    data = output['data']['study'][0]['subjects']
                else:
                    data = data + output['data']['study'][0]['subjects']
                offset += chunk
                
                etime = datetime.datetime.now()
                print("Query (%s) %s" % (offset, str(etime-itime)))          

            if tag == None:
                tag = variable

            for c in data:
                case = c['submitter_id']
                lastnohv = c['hiv_history_records'][0]['lastnohv']
                frsthaav = c['hiv_history_records'][0]['frsthaav']

                if lastnohv > 0:
                    visit_id = case + "_" + str(lastnohv)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    value = visit['summary_lab_results'][0][variable]
                                    if value != None:
                                        values['Last HAART Free'].append(value)                        
                
 
                if frsthaav > 0:
                    # First visit with treatment
                    visit_id = case + "_" + str(frsthaav)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    value = visit['summary_lab_results'][0][variable]
                                    if value != None:
                                        values['First HAART visit'].append(value)

                    # After one year of treatment
                    visit_id = case + "_" + str(frsthaav+20)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    value = visit['summary_lab_results'][0][variable]
                                    if value != None:
                                        values['1 Year Treatment'].append(value)
                                        
                    # After two year of treatment
                    visit_id = case + "_" + str(frsthaav+40)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    value = visit['summary_lab_results'][0][variable]
                                    if value != None:
                                        values['2 Years Treatment'].append(value)


                    # After three year of treatment
                    visit_id = case + "_" + str(frsthaav+60)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    value = visit['summary_lab_results'][0][variable]
                                    if value != None:
                                        values['3 Years Treatment'].append(value)                                        

        with open(filename, 'w') as fp:
            json.dump(values, fp)
                                                                              
    pvalues = {}
    pvalues = run_statistical_test(values, 'Last HAART Free', 'First HAART visit')
    pvalues['1 Year Treatment'] = run_statistical_test(values, 'Last HAART Free', '1 Year Treatment')['1 Year Treatment']
    pvalues['2 Years Treatment'] = run_statistical_test(values, 'Last HAART Free', '2 Years Treatment')['2 Years Treatment']
    pvalues['3 Years Treatment'] = run_statistical_test(values, 'Last HAART Free', '3 Years Treatment')['3 Years Treatment']
    
    ordered_columns = ['Last HAART Free', 'First HAART visit', '1 Year Treatment', '2 Years Treatment', '3 Years Treatment']
    compare_boxplot(values, 'Comparison HAART Treatment', tag, pvalues, ordered_columns)
    
    return values
                                                                          

def run_statistical_test(values, base, compare):
    ''' Run statistical test between two conditions'''
    
    pvalues = {}
    test = ranksums(values[base], values[compare])
    pvalues[base] = -1
    pvalues[compare] = test.pvalue
      
    return pvalues


def compare_boxplot(values, title, measure, pvalues, sorted_columns=None):
    ''' Visualize metrics to compare two conditions'''
    
    if sorted_columns != None:
        columns = sorted_columns
    else:
        columns = pvalues.keys() 
      
    scale = 'linear'
    results = []
    medians = []
    for c in columns:
        results.append(values[c])
        medians.append(np.median(values[c]))
        if max(values[c]) - min(values[c]) > 10000:
            scale = 'log'
    
    N = len(results)
    positions = np.arange(N) + 1    
          
    plt.figure(figsize=(3*N, 6))
    plt.boxplot(list(results), patch_artist=True, widths=0.15)
    plt.xticks(positions, columns, fontsize = 14)
    plt.title(title, fontsize = 16)
    plt.ylabel(measure, fontsize = 14) 
    plt.yscale(scale)
    bottom,top = plt.ylim()
    if scale =='linear':
        plt.ylim(0,top*1.1)
    else:
        plt.ylim(0,top*5)       
    
    for p, condition in zip(positions, columns):
        if pvalues[condition] >= 0:
            col = "red"
            if pvalues[condition] < 0.01: col = "green"
            plt.text(p, top*0.95, "p={0:.4f}".format(pvalues[condition]),
                 horizontalalignment='center', color=col, weight="bold")    
    
    plt.plot(positions, medians, 'r--', lw=2)
    
    plt.show()