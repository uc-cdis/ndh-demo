import matplotlib.pyplot as plt
import numpy as np
import requests
import json
import os

import warnings
warnings.filterwarnings('ignore')

# Create a list of query variables
summary_order_LHV = [
    "_study_count",
    "_subject_count",
    "_demographic_count",
    "_exposure_count",
    "_follow_up_count",
    "_sample_count",
    "_virus_infection_count",
    "_summary_lab_result_count",
    "_aliquot_count",
    "_mRNA_microarray_count",
    "_mRNA_expression_count",
    "_protein_mass_spectrometry_count",
    "_peptide_expression_count",
    "_protein_expression_count",
    "_lipid_mass_spectrometry_count",
    "_metabolite_mass_spectrometry_count"
]

# Create a dictionary of (query variables:header) pairs
summary_count_headers_LHV = {
    "_subject_count": "Subjects",
    "_study_count": "Studies",
    "_demographic_count": "Demographic records",
    "_exposure_count": "Exposure records",
    "_follow_up_count": "Followup records",
    "_sample_count": "Sample records",
    "_virus_infection_count": "Virus Infection records",
    "_summary_lab_result_count": "Lab Results records",
    "_aliquot_count": "Aliquots",
    "_mRNA_microarray_count": "mRNA Microarray",
    "_mRNA_expression_count": "mRNA Expression",
    "_protein_mass_spectrometry_count": "Protein Mass Spectrometry",
    "_peptide_expression_count": "Peptide Expression",
    "_protein_expression_count": "Protein Expression",
    "_lipid_mass_spectrometry_count": "Lipid Mass Spectrometry",
    "_metabolite_mass_spectrometry_count": "Metabolite Mass Spectrometry"
}

# Construct SummaryTable class with properties: content, summary_order and summary_header. SummaryTable is representated as a table.


class SummaryTable:
    ''' Represent result tables in HTML format for visualization '''

    def __init__(self, content, summary_order, summary_header):
        self.content = content
        self.summary_order = summary_order
        self.summary_header = summary_header

    def _repr_html_(self):
        html = []
        html.append("<table style>")
        html.append("<thead>")
        html.append("<th>Category</th>")
        html.append("<th>Counts</th>")
        html.append("</thead>")
        for item in self.summary_order:
            html.append("<tr>")
            html.append("<td>%s</td>" % self.summary_header[item])
            html.append("<td>%s</td>" % self.content[item])
            html.append("<tr>")
        html.append("</table>")
        return ''.join(html)

# Input "credential.json" downloaded from data common profile page to generate access token for data query


def add_keys(filename):
    ''' Get auth from our secret keys '''

    global auth
    json_data = open(filename).read()
    keys = json.loads(json_data)
    auth = requests.post('https://niaid.bionimbus.org/user/credentials/cdis/access_token', json=keys)

# query data from API using query text, returned data is in json format


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


def query_summary_counts(project_id, summary_order, summary_header):
    ''' Query summary counts for each data type'''

    query_txt = "query Counts ($projectID: [String]) { "
    for prop in summary_order:
        query_txt += "%s(project_id: $projectID) " % prop
    query_txt += "}"

    variables = {'projectID': project_id}

    data = query_api(query_txt, variables)

    table = SummaryTable(data['data'], summary_order, summary_header)

    return table


def query_cell_lab():
    ''' Graphql query text for querying virus, timepoint, virus titers and gRNA for cell line '''
    query_text = '''{
      study(first:0,order_by_asc:"submitter_id",project_id: "ndh-dmid-LHV",with_path_to:{type:"subject",species:"Homo sapiens"}){
          submitter_id
          subjects(first:0){
          samples(first:0){
            composition
            hours_to_collection
            virus_infections{
              strain
              mutation
            }
            summary_lab_results(first:0){
              virus_titer_rep1
              virus_titer_rep2
              virus_titer_rep3
              virus_titer_rep4
              virus_titer_rep5
              virus_NP_gRNA
            }
          }
        }
      }
    }'''
    ''' Save response data in dictionary as structure study{virus{timepoint{virus_titer}}} or study{virus{timepoint{gRNA}}}'''
    data = query_api(query_text)
    studies = dict()
    attributes = ['virus_titer_rep1', 'virus_titer_rep2', 'virus_titer_rep3', 'virus_titer_rep4', 'virus_titer_rep5']
    for entity in data['data']['study']:
        study_id = entity['submitter_id']
        studies.setdefault(study_id, {})
        for subject in entity['subjects']:
            for sample in subject['samples']:
                composition = sample['composition']
                tp = str(sample['hours_to_collection'])
                if composition == "Supernatant":
                    for infection in sample['virus_infections']:
                        virus_id = '_'.join(filter(None, [infection['strain'],infection['mutation']]))
                        if virus_id is not None:
                            virus_id = ("_").join(virus_id.split(" "))
                            studies[study_id].setdefault('Supernatant', {})
                            studies[study_id]['Supernatant'].setdefault(virus_id, {})
                            studies[study_id]['Supernatant'][virus_id].setdefault(tp, [])
                            for lab in sample['summary_lab_results']:
                                for attribute in attributes:
                                    if lab[attribute] is not None:
                                        studies[study_id]['Supernatant'][virus_id][tp].append(np.log10(lab[attribute]))
                elif composition == "Cell":
                    for infection in sample['virus_infections']:
                        virus_id = '_'.join(filter(None, [infection['strain'],infection['mutation']]))
                        if virus_id is not None:
                            virus_id = ("_").join(virus_id.split(" "))
                            studies[study_id].setdefault("Cell", {})
                            studies[study_id]['Cell'].setdefault(virus_id, {})
                            studies[study_id]['Cell'][virus_id].setdefault(tp, [])
                            for lab in sample['summary_lab_results']:
                                if lab['virus_NP_gRNA'] is not None:
                                    studies[study_id]['Cell'][virus_id][tp].append(lab['virus_NP_gRNA'])
    return studies


def query_mouse_titer():
    ''' Graphql query text for querying virus, timepoint and virus titer for mouse '''
    query_txt = '''{
      study(first:0,order_by_asc:"submitter_id",project_id: "ndh-dmid-LHV",with_path_to:{type:"subject",species:"Mus musculus"}){
          submitter_id
          subjects(first:0){
            follow_ups(first:0){
              days_to_follow_up
              submitter_id
              samples(first:0){
                summary_lab_results{
                  titer_PFU_per_gram
                }
                }
            }
            virus_infections(first:0){
              strain
              mutation
            }
        }
      }
    }'''
    ''' Save response data in dictionary as structure study{virus{timepoint{virus_titer}}}'''
    data = query_api(query_txt)
    studies = dict()
    for entity in data['data']['study']:
        study_id = entity['submitter_id']
        studies.setdefault(study_id, {})
        for subject in entity['subjects']:
            for infection in subject['virus_infections']:
                virus_id = '_'.join(filter(None, [infection['strain'],infection['mutation']]))
                studies[study_id].setdefault(virus_id, {})
            for follow_up in subject['follow_ups']:
                day_id = str(follow_up['days_to_follow_up'])
                studies[study_id][virus_id].setdefault(day_id, [])
                for sample in follow_up['samples']:
                    for lab in sample['summary_lab_results']:
                        titer = lab['titer_PFU_per_gram']
                        if titer is not None:
                            studies[study_id][virus_id][day_id].append(np.log10(titer + 1))
    return studies


def plot_cell_titer(studies):
    ''' plot virus titer for cell line. Each line is for a virus. x-axis is timepoint, y-axis is average virus titer with std'''
    data = query_cell_lab()
    all_timepoints = []
    fig = plt.figure(figsize=(6, 4))
    for study in studies:
        virus_list = data[study]["Supernatant"].keys()
        virus_list = filter(None, virus_list)
        for virus in virus_list:
            times = list()
            titer_aves = list()
            titer_stds = list()
            for timepoint in sorted(map(int, data[study]["Supernatant"][virus].keys())):
                tp_key = str(timepoint)
                if data[study]["Supernatant"][virus][tp_key]:
                    titer_aves.append(np.mean(data[study]['Supernatant'][virus][tp_key]))
                    titer_stds.append(np.std(data[study]['Supernatant'][virus][tp_key]))
                    times.append(timepoint)
                    all_timepoints.append(timepoint)

            ax1 = fig.add_subplot(111)
            ax1.errorbar(times, titer_aves, yerr=titer_stds, fmt='-o', label=virus)
            ax1.set_xlabel("Timepoint", fontsize=14)
            ax1.set_ylabel("Virus Titer", fontsize=14)
            ax1.set_xticks(times)
            ax1.legend(loc='best', fancybox=True, framealpha=0.5)
        ax1.set_xticks(all_timepoints)


def plot_mouse_titer(study, virus_list=None):
    ''' plot virus titer for mouse. Each line is for a virus. x-axis is timepoint, y-axis is average virus titer with std'''
    data = query_mouse_titer()
    if virus_list is None:
        virus_list = data[study].keys()
        virus_list = filter(None, virus_list)
    fig = plt.figure(figsize=(6, 4))
    for virus in virus_list:
        times = list()
        titer_aves = list()
        titer_stds = list()
        for timepoint in sorted(map(int, data[study][virus].keys())):
            tp_key = str(timepoint)
            if data[study][virus][tp_key]:
                titer_aves.append(np.mean(data[study][virus][tp_key]))
                titer_stds.append(np.std(data[study][virus][tp_key]))
                times.append(timepoint)

        ax1 = fig.add_subplot(111)
        ax1.errorbar(times, titer_aves, yerr=titer_stds, fmt='-o', label=virus)
        ax1.set_xlabel("Timepoint", fontsize=14)
        ax1.set_ylabel("Virus Titer", fontsize=14)
        ax1.set_xticks(times)
        ax1.legend(loc='best', fancybox=True, framealpha=0.5)


def plot_gRNA(studies):
    ''' plot gRNA for cell line. Each line is for a virus. x-axis is timepoint, y-axis is virus genomic RNA'''
    data = query_cell_lab()
    fig = plt.figure(figsize=(6, 4))
    for study in studies:
        virus_list = data[study]["Cell"].keys()
        virus_list = filter(None, virus_list)
        for virus in virus_list:
            times = list()
            gRNA = list()
            for timepoint in sorted(map(int, data[study]["Cell"][virus].keys())):
                tp_key = str(timepoint)
                if data[study]["Cell"][virus][tp_key]:
                    times.append(timepoint)
                    gRNA.append(data[study]['Cell'][virus][tp_key])
            ax1 = fig.add_subplot(111)
            ax1.plot(times, gRNA, 'o-', label=virus)
            ax1.set_xlabel("Timepoint", fontsize=14)
            ax1.set_ylabel("gRNA", fontsize=14)
            ax1.set_xticks(times)
            ax1.legend(loc='best', fancybox=True, framealpha=0.5)


def query_mouse_weight():
    ''' Graphql query text for querying virus, timepoint and weight_percentage for mouse '''
    query_txt = '''{
      study(first:0,order_by_asc:"submitter_id",project_id: "ndh-dmid-LHV"){
          submitter_id
          subjects(first:0, with_links:"follow_ups"){
            follow_ups(first:0){
              days_to_follow_up
              weight_percentage
            }
            virus_infections(first:0){
              strain
              mutation
            }
        }
      }
    }'''
    ''' Save response data in dictionary as structure study{virus{timepoint{weight_percentage}}}'''
    data = query_api(query_txt)
    included_studies = ['IM101', 'IM102', 'IM103']
    studies = dict()
    for entity in data['data']['study']:
        if entity['submitter_id'] not in included_studies:
            continue
        else:
            study_id = entity['submitter_id']
            studies.setdefault(study_id, {})
            for subject in entity['subjects']:
                for infection in subject['virus_infections']:
                    virus_id = '_'.join(filter(None, [infection['strain'],infection['mutation']]))
                    studies[study_id].setdefault(virus_id, {})
                for follow_up in subject['follow_ups']:
                    time = str(follow_up['days_to_follow_up'])
                    studies[study_id][virus_id].setdefault(time, [])
                    weight_percentage = follow_up['weight_percentage']
                    if weight_percentage:
                        studies[study_id][virus_id][time].append(weight_percentage)
    return studies


def plot_weight_percentage(study):
    ''' plot weight_percentage for mouse. Each line is for a virus. x-axis is timepoint, y-axis is average weight_percentage with std '''
    data = query_mouse_weight()
    fig = plt.figure(figsize=(6, 4))
    virus_list = data[study].keys()
    virus_list = filter(None, virus_list)
    for virus in virus_list:
        times = list()
        weight_aves = list()
        weight_stds = list()
        for timepoint in sorted(map(int, data[study][virus].keys())):
            tp_key = str(timepoint)
            if data[study][virus][tp_key]:
                weight_aves.append(np.mean(list(map(float, data[study][virus][tp_key]))))
                weight_stds.append(np.std(list(map(float, data[study][virus][tp_key]))))
                times.append(timepoint)
        ax1 = fig.add_subplot(111)
        ax1.errorbar(times, weight_aves, yerr=weight_stds, fmt='-o', label=virus)
        ax1.set_xlabel("Timepoint", fontsize=14)
        ax1.set_ylabel("weight_percentage", fontsize=14)
        ax1.set_xticks(times)
        ax1.legend(loc='best', fancybox=True, framealpha=0.5)
