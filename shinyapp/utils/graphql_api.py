from cdispyutils.hmac4 import get_auth
import requests
import time
import json

def get_api_auth(filename, api_url):
  
  access_url = api_url + 'user/credentials/cdis/access_token'
  json_data=open(filename).read()
  keys = json.loads(json_data)
  t = requests.post(access_url, json=keys)

  return t


def query(query_txt, auth, api_url, variables=None):

   api_url = api_url + 'api/v0/submission/graphql'

   if variables == None:
      query = {'query': query_txt}
   else:
      query = {'query': query_txt, 'variables': variables}    

   tries = 0
   while tries < 1:
      output = requests.post(api_url, headers={'Authorization': 'bearer '+ auth.json()['access_token']}, json=query).text

      data = json.loads(output)  

      if 'errors' in data:
         print data

      if not 'data' in data:
         print query_txt
         print data

      tries += 1          

   return data

def get_projects(auth, api_url, excluded=[]):
   
   query_txt = """query Project { project(first:0) {project_id}} """
   
   data = query(query_txt, auth, api_url) 
   
   projects = []
   for pr in data['data']['project']:
      if pr['project_id'] not in excluded:
          projects.append(pr['project_id'])
   projects = sorted(projects)   

   return projects

def get_studies(auth, api_url, project_id):
   
   query_txt = """query Study { study(first:0, project_id: "%s"){ submitter_id}} """ % project_id
   data = query(query_txt, auth, api_url) 
   
   studies = []
   for st in data['data']['study']:
      studies.append(st['submitter_id'])
   
   studies = sorted(studies)   

   return studies


