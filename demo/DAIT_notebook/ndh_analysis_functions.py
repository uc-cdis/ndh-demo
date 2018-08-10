import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from lifelines.datasets import load_waltons
from lifelines import KaplanMeierFitter
from scipy.stats import ranksums
from scipy import interpolate
from operator import add
import numpy as np
import pandas as pd
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
   "_follow_up_count",
   "_sample_count",
   "_aliquot_count", 
   "_summary_lab_result_count",
   "_summary_drug_use_count"
]

summary_count_headers = {
    "_subject_count": "Cases",
    "_study_count": "Studies",
    "_demographic_count": "Demographic records",
    "_follow_up_count": "Visit records",
    "_sample_count": "Samples", 
    "_aliquot_count": "Aliquots", 
    "_summary_lab_result_count": "Lab Results records",
    "_summary_drug_use_count": "Drug records"
}

excluded_studies = ['study-01']
excluded_projects = ['ndh-test', 'ndh-vir-simulation']

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

def add_keys(filename):
    ''' Get auth from our secret keys '''

    global auth
    json_data = open(filename).read()
    keys = json.loads(json_data)
    auth = requests.post('https://niaid.bionimbus.org/user/credentials/cdis/access_token', json=keys)    
    
def get_keys():
    ''' Get auth from internal service '''

    global auth
    auth = requests.get('http://fence-service.default.svc.cluster.local/internal/access_token')
    
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

def get_projects():
    ''' Query list of projects '''
    
    query_txt = """query Project { project(first:0) {project_id}} """
   
    data = query_api(query_txt) 
   
    projects = []
    for pr in data['data']['project']:
        if pr['project_id'] not in excluded_projects:
            projects.append(pr['project_id'])
    projects = sorted(projects)   

    return projects

def query_summary_counts(projects=None):
    ''' Query summary counts for each data type'''
   
    if projects == None:
        projects = get_projects()
    elif not isinstance(projects,list):
        projects = [projects]
       
    dftotal = pd.DataFrame()
    for p in projects:
        query_txt = """query Counts ($projectID: [String]) {"""
        for param in summary_order:
            query_txt += """%s(project_id: $projectID)""" % param
        query_txt += "}" 
        variables = { 'projectID': p}
        data = query_api(query_txt, variables)
        indexes, values = [], []
        for key in summary_order:
            indexes.append(summary_count_headers[key])
            if key in data['data']:
                values.append(data['data'][key])
            else:
                values.append(0)           

        df = pd.DataFrame(values, index=indexes, columns=[p])
        if dftotal.empty:
            dftotal = df
        else:
            dftotal[p] = df[p]

    #dftotal = pd.concat(dftotal)    
 
    return dftotal

def query_summary_field(field, field_node, project_id = None):
    ''' Query summary counts for specific node'''
   
    if project_id != None:
        query_txt = """query { %s(first:0, project_id: "%s") {%s}} """ % (field_node, project_id, field) 
    else:
        query_txt = """query { %s(first:0) {%s project_id}} """ % (field_node, field)    
    
        
    data = query_api(query_txt)
    
    summary = {}
    total = []
    for d in data['data'][field_node]:
        
        if isinstance(d[field], float):
            d[field] = str(d[field])[:-2]        
        
        if 'project_id' in d:
            
            if d['project_id'] in excluded_projects:
                continue
                
            summary.setdefault(d['project_id'], {})
            summary[d['project_id']].setdefault(d[field], 0)
            summary[d['project_id']][d[field]] += 1
            if d[field] not in total:
                total.append(d[field])            
        else:
            summary.setdefault(d[field], 0)        
            summary[d[field]] += 1
    
    if project_id != None:
        plot_field_metrics(summary, field)
    else:
        plot_overall_metrics(summary, field, total)        
    
    return summary

def plot_field_metrics(summary_counts, field):
    ''' Plot summary results in a barplot ''' 
    
    N = len(summary_counts)

    values = []
    types = []

    for n in sorted(summary_counts, key=summary_counts.get, reverse=True):
        values.append(summary_counts[n])
        types.append(n)
           
    plot_bars(values, types, field)


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
    size_prop = 1.2 #(N/10) + 1
    
    plots = []
    plt.figure(figsize=(8, 4))
    left = [0]*N
    for pr in sorted_projects:
        p = plt.barh(positions, results[pr], bar_size, left, align='center', alpha=0.7)        
        plots.append(p[0])
        left = list(map(add, left, results[pr]))
        
    plt.title('Summary counts by (' + field + ')', fontsize=10*size_prop)
    plt.xlabel('COUNTS', fontsize=10*size_prop)    
    plt.xlim(0, max(left)+5)
    plt.ylabel(field.upper(), fontsize=10*size_prop)  
    plt.yticks(positions, totals, fontsize=10*size_prop)    
    plt.legend(plots, sorted_projects, fontsize=10*size_prop, loc='upper right')
           
    plt.show()      
    
def field_distribution(field, field_node, project_id, bins = None, distrib=None, rate=None):
    ''' Plot distribution for one field'''
   
    if project_id != None:
        query_txt = """query { %s(first:0, project_id: "%s") {submitter_id %s}} """ % (field_node, project_id, field) 
    else:
        query_txt = """query { %s(first:0) {submitter_id %s project_id}} """ % (field_node, field)    
         
    data = query_api(query_txt)
       
    summary = {}
    total = []
    for d in data['data'][field_node]:
                
        if isinstance(d[field], float):
            d[field] = str(d[field])[:-2]        
        
        if 'project_id' in d:  
            summary.setdefault(d['project_id'], {})
            summary[d['project_id']].setdefault(d[field], 0)
            summary[d['project_id']][d[field]] += 1
            if d[field] not in total:
                total.append(d[field])            
        else:
            summary.setdefault(d[field], 0)        
            summary[d[field]] += 1    
         
    if len(summary)>10:
        
        accumulated = []
        for d in data['data'][field_node]:
            if d[field] != None:
                accumulated.append(float(d[field]))        
        
        plot_histogram(accumulated, field, bins, distrib, rate)
        
    else:
        
        N = len(summary)

        values = []
        types = []

        for n in sorted(summary, key=summary.get, reverse=True):
            values.append(summary[n])
            types.append(n)
            
        total = sum(values)
        positions = np.arange(N)
        fig, ax = plt.subplots(1, 1, figsize=(3*N, N))

        size_prop = (N/10) + 1
        ax.bar(positions, values, 0.2, align='center', alpha=0.4, color='b')
  
        plt.title('Summary counts by (' + field + ')', fontsize=10*size_prop)
        plt.ylabel('COUNTS', fontsize=10*size_prop)    
        plt.ylim(0, max(values)+5)
        plt.xlabel(field.upper(), fontsize=10*size_prop)  
        plt.xticks(positions, types, fontsize=10*size_prop)         
   
        # fit curve
        if distrib == 'exponential':
            fit_curve = expon.pdf(positions, 0, 1.0/rate)*total
            ax.plot(positions, fit_curve, 'r-', lw=2)
        if distrib == 'uniform':
            fit_curve = [total/float(len(positions))] * len(positions)
            ax.plot(positions, fit_curve, 'r-', lw=2)  

    return data    
    
    
def plot_histogram(accumulated, xlab, bins = None, distrib=None, rate=None):  
        
    # the histogram of the data
    plt.figure(figsize=(8, 4))
    fig, ax = plt.subplots(1, 1)
    n, positions, patches = ax.hist(accumulated, bins, facecolor='b', alpha=0.75)
    total=len(accumulated)

    plt.xlabel(xlab)
    plt.ylabel('Counts')
    plt.title('Histogram of ' + xlab)
    plt.grid(True)

    # fit curve
    if distrib == 'exponential':
        fit_curve = expon.pdf(positions, 0, 1.0/rate)*total
        ax.plot(positions, fit_curve, 'r-', lw=2)
    if distrib == 'uniform':
        fit_curve = [total/float(len(positions))] * len(positions)
        ax.plot(positions, fit_curve, 'r-', lw=2)  

        
def plot_bars(values, types, xlab, errors=[]):         
        
    N = len(values)
    positions = np.arange(N)        
    plt.figure(figsize=(3*N, N))   
    
    size_prop = (N/10) + 1
    if errors:
        plt.bar(positions, values, 0.4, align='center', alpha=0.4, color='cadetblue', yerr=errors, error_kw=dict(elinewidth=3, capsize=20, marker='o', markeredgewidth=3))
    else:
        plt.bar(positions, values, 0.4, align='center', alpha=0.4, color='cadetblue')
        plt.ylim(0, int(max(values)*1.1))
        
    plt.title('Summary counts by (' + xlab + ')', fontsize=10*size_prop)
    plt.ylabel('COUNTS', fontsize=10*size_prop)
    plt.yticks(fontsize=10*size_prop)
    plt.xlabel(xlab.upper(), fontsize=10*size_prop)  
    plt.xticks(positions, types, fontsize=10*size_prop, rotation='vertical')    
    
    if not errors:
        for i, v in enumerate(values):
            plt.text(i-0.05, v, str(v), color='red', fontweight='bold', fontsize=10*size_prop)        
        
    
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

                if negvis != None and negvis > 0:
                    visit_id = case + "_" + str(negvis)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    negvalue = visit['summary_lab_results'][0][variable]
                                    if negvalue != None:
                                        values['Last Seronegative'].append(negvalue)              

                if posvis != None and posvis > 0:
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


def compare_survival(project_id, study_id, variable, tag=None):
    ''' Compare survival with haart vs not ehaart'''
    
    if study_id:
        filename = '%s_%s_survival.json' % (variable,study_id)
    else:
        filename = '%s_survival.json' % (variable)
    
    if os.path.isfile(filename):
        json_data=open(filename).read()
        values = json.loads(json_data)
    else:
        values = {}
        if study_id:
            studies = [study_id]
        else:
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
                print("Query %s (%s) %s" % (study, offset, str(etime-itime)))

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
                    if death_year != None and fposdate != None and death_year < 2100 and fposdate < 2100:
                        values['death_year'].append(death_year - fposdate)                                                 
                    elif fposdate != None and fposdate < 2100:
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

                if lastnohv != None and lastnohv > 0:
                    visit_id = case + "_" + str(lastnohv)
                    if c['follow_ups'] != []:
                        for visit in c['follow_ups']:
                            if visit_id == visit['submitter_id']:
                                if visit['summary_lab_results'] != []:
                                    value = visit['summary_lab_results'][0][variable]
                                    if value != None:
                                        values['Last HAART Free'].append(value)                        
                
 
                if frsthaav != None and frsthaav > 0:
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

def plot_scatterfit(x,y,field,a=None,b=None):  
  
    c = 0
    fig, ax = plt.subplots(1, 1, figsize=(20, 8))  
    colors = ['m', 'y', 'b', 'c',  'r']
    for g in x:
        col = [colors[c]]*len(x[g])      
        ax.scatter(x[g], y[g], c=col, label=g)
        c += 1
        
    ax.legend(fontsize=14)
    ax.grid(True)          
        
    plt.title(field + ' mg/dL by Subject', fontsize=20)
    plt.xlabel('AGE', fontsize=20)    
    plt.ylabel(field.upper(), fontsize=20)  
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.show() 
    
def plot_violin(xvalues, yvalues, field):
    
    plt.figure(figsize=(20, 6))
    
    positions = range(1,len(yvalues)+1)
    data = []
    for g in xvalues:
        data.append(yvalues[g])

    parts = plt.violinplot(data, positions, points=500, widths=0.7, showmeans=True,
                          showextrema=True, showmedians=True, bw_method=0.4)

    for pc in parts['bodies']:
        pc.set_facecolor('red')
    
    plt.ylabel(field.upper(), fontsize=20)  
    plt.xticks(positions, xvalues, fontsize=15)
    
    plt.show()    
    
    
def plot_lab_results(variable, groups, projects=None, timepoints=None):    
    
    if projects == None:
        projects = get_projects()
    elif not isinstance(projects,list):
        projects = [projects]    
        timepoints = [timepoints]
    
    data = {}
    for idx, project in enumerate(projects):
        
        cached_file = project + '_lab_records.json'
        
        if os.path.isfile(cached_file):
            for cfile in glob.glob('./' + project + '_lab_records*'):
                json_data=open(cfile).read()
                output = json.loads(json_data)  
                if not data:
                    data = output
                else:
                    data = data + output  
            
            #print("%s: %s Records" % (project, len(data)))
        else:
            
            # Get counts
            count_query = '{ _subject_count(project_id: "%s")}' % (project)
            counts = query_api(count_query)['data']['_subject_count']       

            # Create query with pagination
            offset = 0
            chunk = 3
            while offset <= counts:

                # Log time
                itime = datetime.datetime.now()            

                query_txt = """{ subject(first:%s, offset:%s, project_id: "%s", order_by_asc: "submitter_id"){
                                    demographics{
                                       %s
                                       human_age_at_index
                                       race
                                    }
                                    hiv_status
                                    project_id
                                    follow_ups(visit_number: %s){
                                       age_at_visit
                                       summary_lab_results{
                                        %s
                                        hdlchol
                                    }
                                }}}""" % (chunk, offset, project, groups, timepoints[idx], variable)

                output = query_api(query_txt)
                if not data:
                    data = output['data']['subject']
                else:
                    data = data + output['data']['subject']
                offset += chunk

                etime = datetime.datetime.now()
                print("Query %s (%s/%s): %s" % (project, offset, counts, str(etime-itime)))

                if offset % 102 == 0:
                    with open(cached_file, 'w') as fp:
                       json.dump(data, fp)
                    add_keys('credentials.json')
    
    neg = 0
    male = 0
    female = 0
    if groups == "race":
        dict_race = {
            "Other": "Other",
            "Unspecified": "Other",
            "Unknown": "Other",
            "Asian": "Asian",
            "Asian/Pacific Islander": "Asian",
            "Black": "Black",
            "White": "White",
            "Multi-racial": "Other",
            "American Indian or Alaskan Native": "American Native"
        }
        # Prepare data for plot
        xvalues = {}
        yvalues = {} 
        for subj in data:
            if 'follow_ups' in subj and len(subj['follow_ups'])>0 and \
               'summary_lab_results' in subj['follow_ups'][0] and \
                len(subj['follow_ups'][0]['summary_lab_results'])>0 \
                and (subj['hiv_status'] == False or subj['hiv_status'] == None):
                    group = dict_race[subj['demographics'][0][groups]]
                    xvalues.setdefault(group,[])
                    yvalues.setdefault(group,[])
                    xval = group
                    yval = subj['follow_ups'][0]['summary_lab_results'][0][variable]
                    if xval != None and yval != None:
                        xvalues[group].append(xval)
                        yvalues[group].append(yval)        
        
        # Show violin plot
        plot_violin(xvalues,yvalues,variable)
    else:
        
        # Prepare data for plot
        xvalues = {}
        yvalues = {} 
        for subj in data:
            if subj['hiv_status'] == False:
                neg += 1
                if subj['demographics'][0][groups] == "male":
                    male += 1
                else:
                    female += 1
            if 'follow_ups' in subj and len(subj['follow_ups'])>0 and \
               'summary_lab_results' in subj['follow_ups'][0] and \
                len(subj['follow_ups'][0]['summary_lab_results'])>0 \
                and (subj['hiv_status'] == False or subj['hiv_status'] == None):
                    group = subj['demographics'][0][groups] + '-' + subj['project_id']
                    xvalues.setdefault(group,[])
                    yvalues.setdefault(group,[])
                    xval = subj['demographics'][0]['human_age_at_index']
                    yval = subj['follow_ups'][0]['summary_lab_results'][0][variable]
                    if xval == None:
                        xval = subj['follow_ups'][0]['age_at_visit']
                    if xval != None and yval != None:
                        xvalues[group].append(xval)
                        yvalues[group].append(yval)        
        
        # Show scatter plot
        plot_scatterfit(xvalues,yvalues,variable)
    
    return(data)