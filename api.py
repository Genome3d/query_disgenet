#! /usr/bin/env python

#For this example we are going to use the python default http library
import requests
import io
import json
import sys
import os
#from bioc import biocjson
#import pandas as pd
import argparse
import time

import logger


def fetch_dga(genes_df, source, email, password):
    
    #select first col of genes_df
    
    genes_list = genes_df[genes_df.columns[0]].tolist()
    n_genes=str(len(genes_list))
    genes_padded= "%2C%20".join(genes_list)
    #print(genes_padded)
    #Build a dict with the following format, change the value of the two keys your DisGeNET account credentials, if you don't have an account you can create one here https://www.disgenet.org/signup/ 
    auth_params = {"email": email,"password":password}
    api_host = "https://www.disgenet.org/api"
    api_key = None
    s = requests.Session()
    result = "no results returned"
    try:
        r = s.post(api_host+'/auth/', data=auth_params)
        if(r.status_code == 200):
            #Lets store the api key in a new variable and use it again in new requests
            json_response = r.json()
            api_key = json_response.get("token")
      #   print(api_key + "This is your user API key.") #Comment this line if you don't want your API key to show up in the terminal
        # else:
             # print(r.status_code)
             # print(r.text)
    except requests.exceptions.RequestException as req_ex:
        #print(req_ex)
        print("Something went wrong with the request.")

    if api_key:
        #Add the api key to the requests headers of the requests Session object in order to use the restricted endpoints.
        s.headers.update({"Authorization": "Bearer %s" % api_key}) 
        #Lets get all the diseases associated to a gene eg. APP (EntrezID 351) and restricted by a source.
        #gda_response = s.get(api_host+'/gda/gene/', params={'gene':"YAP, P53",'source':source, 'format':format})
        gda_response = s.get(api_host+'/gda/gene/'+genes_padded, params={'source':source, 'format':"tsv"})
        
        # if(gda_response.status_code == 404): 
           # # print(gda_response.status_code)
            # print(str(n_genes) + " genes returned no results") 
            
        # else: 
            # result = gda_response.text
            # #print(gda_response.status_code)
        if(gda_response.status_code != 404):
            result = gda_response.text
            
    if s:
        s.close()
            
    return result

