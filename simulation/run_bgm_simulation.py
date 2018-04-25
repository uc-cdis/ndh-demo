#!/usr/bin/env python2
import pandas as pd
import requests
import shutil
import sys
import csv
import os

# Get enviromental variables
inurl = os.environ['INPUT_URL']
outurl = os.environ['OUTPUT_URL']

# Download file from preassigned URL
matrix = 'matrix.csv'
response = requests.get(inurl, stream=True)
with open(matrix, 'wb') as out_file:
    shutil.copyfileobj(response.raw, out_file)

# Read binary matrix
mat = pd.read_csv(matrix)

# Read groups
groups = pd.read_csv('groups.csv')

# Prepare subset for cross-validation
for c in list(groups):
   selecidx = [index for index,value in enumerate(groups[c]) if value == 1]
   subset = mat.iloc[selecidx]
   subset.to_csv(c+".csv", index=False, quoting=csv.QUOTE_NONNUMERIC)

# Run BGM with HyPhy tool
outfile = "img.png"
os.system("/home/ubuntu/hyphy/HYPHYMP /home/ubuntu/hyphy/virulence/VirBGM.bf")
os.system("dot -Tpng bgm.dot > " + outfile)

# Submit file to presigned URL
with open(outfile, 'rb') as data:
    requests.put(outurl, data=data)
