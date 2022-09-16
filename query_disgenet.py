#! /usr/bin/env python

#For this example we are going to use the python default http library
import requests
import io
import json
import sys
import os
#from bioc import biocjson
import pandas as pd
import argparse
import time
import shutil

import logger
import api

def join_path(*args):
    fp = ''
    for arg in args:
        fp = os.path.join(fp, arg)
    return fp

    
def write_results(df, output_fp, logger): 
        
    out_dir = os.path.dirname(output_fp)
    os.makedirs(out_dir, exist_ok=True)
    logger.write('Writing output...')
    with open(output_fp, 'w') as f:
        print(df, file=f) 
    #df.to_csv(output_fp, sep='\t', index=False)
    
    
def chunk_genes(genes, chunksize, out_dir):
    genes = genes.unique()
    chunks = [genes[i:i+chunksize] for i in range(0, len(genes), chunksize)]
    for i, chunk in enumerate(chunks):
        pd.Series(chunk).to_csv(os.path.join(out_dir, f'{i}.txt'), sep='\t',
                                index=False, header=None)    

def parse_genes(genes_args, logger):
    logger.write('Parsing gene input...')
    df = pd.DataFrame()
    if os.path.isfile(genes_args[0]):
        df = pd.read_csv(genes_args[0], sep='\t')
        if not 'gene' in df.columns:
            sys.exit('Input file has no gene column.')
        df = df[df.gene.str.lower() != 'gene']
        df = df[['gene']].drop_duplicates()
        
    else:
        sys.exit('Input gene file does not exist.')
        #if gene list
        # df = pd.DataFrame({'gene': [i.upper() for i in genes_args]})
        # df = df[['gene']].drop_duplicates()
       
    return df
    
def parse_args():
    parser = argparse.ArgumentParser(
        description='Retrieve DisGeNet annotations of gene-disease associations.')
    parser.add_argument(
        '-e', '--email', required=True, help="DisGeNET user email, if you do not have an account you can create one here https://www.disgenet.org/signup/'.")
    parser.add_argument(
        '-p', '--password', required=True, help="DisGeNET user password, if you do not have an account you can create one here https://www.disgenet.org/signup/'.") 
    parser.add_argument(
        '-g', '--genes', required=True, nargs='+',
        help='Space-separated gene names or a file with gene names in the "gene" column.')
    # parser.add_argument(
        # '-i', '--in-dir', required=True,
        # help='Directory containing files each having PUBMED IDs of interest.')
    # parser.add_argument(
        # '-f', '--format', default='tsv',
        # choices=['tsv', 'json', 'xml'],
        # help='Return type format: json, tsv and xml')
    parser.add_argument(
        '-s', '--source', default='CURATED',
        choices=['CURATED','INFERRED', 'ANIMAL_MODELS', 'ALL', 'BEFREE', 'CGI','CLINGEN','CLINVAR', 
        'CTD_human', 'CTD_mouse', 'CTD_rat','GENOMICS_ENGLAND', 'GWASCAT', 'GWASDB' ,'HPO', 
        'LHGDN','MGD','ORPHANET','PSYGENET','RGD', 'UNIPROT'],
        help='Source of the GDA: CURATED , INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN, CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB ,HPO, LHGDN, MGD ,ORPHANET , PSYGENET , RGD , UNIPROT.')
    # parser.add_argument(
        # '-o', '--output', required=True, help='Filepath to write results.') 
    parser.add_argument(
        '-o', '--output-dir', required=True,
        help='Directory to write results.')
    return parser.parse_args()

if __name__ == "__main__":
    pd.options.mode.chained_assignment = None
    args = parse_args()
    start_time = time.time()
    os.makedirs(args.output_dir, exist_ok=True)
    
    genes_dir = os.path.join(args.output_dir, 'genes')
    os.makedirs(genes_dir, exist_ok=True)
    
    temp_dir = os.path.join(args.output_dir, 'temp')
    os.makedirs(temp_dir, exist_ok=True)
    
    global logger
    logger = logger.Logger(logfile=os.path.join(args.output_dir, 'disgenet.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        if arg == 'password':
            
            logger.write(f'{arg}:\t *********')
        
        else: logger.write(f'{arg}:\t {getattr(args, arg)}')
    logger.write('\n')
    genes = parse_genes(args.genes, logger)
    
    logger.write("number of genes to query: " + str(len(genes)))
    #print("number of genes to query: " + str(len(genes)))
    #logger.write('Number of query genes: ')
    #DisGeNet api accepts a max of 100 genes per query
    if len(genes)>100:
        chunk_genes(genes['gene'], 100, genes_dir)
        logger.write('Retrieving disease-gene associations..')
        
        i=0
        for fp in [fp for fp in os.listdir(genes_dir) if fp.endswith('.txt')]:
            #read genes file as df
            batch_fp = os.path.join(genes_dir, fp)
            df = pd.read_csv(batch_fp, sep='\t')
            #pass df to api.fetch_dga
            dga_results=api.fetch_dga(df, args.source, args.email, args.password)
            
            if (dga_results != "no results returned"):
                write_results(dga_results, join_path(temp_dir, f'temp_dga{i}.txt'), logger)
                i=i+1
    
        dga_results= []
        for fp in [fp for fp in os.listdir(temp_dir) if fp.endswith('.txt')]:
            #concatenate output to one file and save new file
            batch_fp = os.path.join(temp_dir, fp)
            #read first row as header
            pd.set_option('display.max_columns', None)
            df = pd.read_csv(batch_fp, header=0, sep='\t')
            dga_results.append(df)
        dga_results = pd.concat(dga_results)
        dga_results.to_csv(join_path(args.output_dir, 'dga_results.txt'), index=False, sep="\t")
        
        #write_results(dga_results, join_path(args.output_dir, 'dga_results.txt'), logger)
        #delete temp_dir
        shutil.rmtree(temp_dir)
    
    
    #if list has < 100 genes
    else: 
        logger.write('Retrieving disease-gene associations..')
        dga_results=api.fetch_dga(genes, args.source, args.email, args.password)
        #might be issue 
        write_results(dga_results, join_path(args.output_dir, 'dga_results.txt'), logger)
        
        # with open('dga_results.txt', 'w') as f:
            # print(dga_results, file=f) 
        
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
