#!/usr/bin/env python2
import requests
import sys
import argparse
import subprocess
import datetime
import pandas as pd
from utils import graphql_api
import time
import json
import csv
import os

# Select arguments:
parser = argparse.ArgumentParser(description="Select input files for Bayesian Graphical Model (BGM) simulation:")
parser.add_argument('-m', '--matrix', required=True, help="UUID of filename for the file (CSV) including the binary matrix with features as columns")
parser.add_argument('--groups', help="Sample classification to do cross-validation", default='./groups.csv')
parser.add_argument('--keys_file', help='Credentials for the API', default='./utils/credentials.json')
parser.add_argument('--api_url', help='URL for the API', default='https://niaid.bionimbus.org/')
parser.add_argument('--output', help='Filename for output', default='bgm.png')
args = parser.parse_args()

# Get file from data model
auth = graphql_api.get_api_auth(args.keys_file, args.api_url) 
if "." in args.matrix:
	query_txt = """query Files { mutation_panel(file_name: "%s"){id file_name}}""" % args.matrix
	data = graphql_api.query(query_txt, auth, args.api_url)
	fileid = data['data']['mutation_panel'][0]['id']
else:
	fileid = args.matrix

matrix = graphql_api.download_data(fileid, auth, args.api_url)

# Read binary matrix
mat = pd.read_csv(matrix)

# Read groups
groups = pd.read_csv(args.groups)

# Prepare subset for cross-validation
for c in list(groups):
   selecidx = [index for index,value in enumerate(groups[c]) if value == 1]
   subset = mat.iloc[selecidx]
   subset.to_csv(c+".csv", index=False, quoting=csv.QUOTE_NONNUMERIC)

# Run BGM with HyPhy tool
os.system("/home/ubuntu/hyphy/HYPHYMP /home/ubuntu/hyphy/virulence/VirBGM.bf")

os.system("cp ./bgm.dot ./data/")
os.system("cp ./edgelist.out ./data/")
os.system("dot -Tpng ./data/bgm.dot > ./data/bgm.png")
