#!/usr/bin/env python

#import everything that is necessary in order to run the script.
import time
import sqlite3 as sl
from datetime import datetime
import requests
import sys
from tqdm import tqdm

"""
time: used to let to program sleep
sqlite3: used to create and push command to generated database
datetime: used to calculate how much time the retrieval from data and push to database is
requests: used to send requests to the GNOMAD api. 
"""

def fetch(data, url="https://gnomad.broadinstitute.org/api"):
    """
    This part of the code fetches the data from the gnomad database. The online API is a GraphQL query, which is also
    the data given to the fetch-code.The request for the API has the url of the database, the json data and headers.
    These headers are needed, otherwise the API will be overloaded much faster. If the response is 200, the json data
    will be retrieved and returned. If there is a 400 or 502 response, the API refuses the request due to an API
    overload or too much API requests in a short time frame (such as an DDOS attack, which is not preffered). If that
    is the case, there will be a 30 second wait time to clear the API. After this, there will be a new chance to
    retrieve the data. The json format of this data will be returned.

    :param data: json data with GraphQL query with variant or ENSEMBL id
    :param url: API of the gnomad GrapQL database
    :return: json: returns json file
    """
    headers = {'User-Agent': 'Mozilla/5.0'}
    response = requests.post(url, json=data, headers=headers)
    # print(str(response))
    # print("here?")

    #if there is an error with API overload: show it here.
    if "200" in str(response):
        if "json" in response.headers.get('content-type'):
            json = response.json()
            # time.sleep(1)

    else:
        if "400" or "502" in str(response):
            # print("===ERROR===")
            time.sleep(30)
            json = fetch(data)
    return json


def get_variant_list(gene_id, dataset="gnomad_r2_1"):
    """
    In this part of the code, the GraphQL query for the variant list of the gene will be made. This query will be
    formated into json, which will be pushed to the fetch code.

    :param gene_id:ENSEMBL gene id
    :param dataset: dataset in which the API will search.
    :return: variant list of gene in json format
    """
    # Note that this is GraphQL, not JSON.
    fmt_graphql = """
    {
        gene(gene_id: "%s", reference_genome: GRCh37) {
          variants(dataset: %s) {
            consequence
            pos
            rsid
            variant_id: variantId
          }
        }
      }
    """
    # This part will be JSON encoded, but with the GraphQL part left as a
    # glob of text.
    req_variantlist = {
        "query": fmt_graphql % (gene_id, dataset),
        "variables": {}
        }
    response = fetch(req_variantlist)
    return response["data"]["gene"]["variants"]


def get_variant(variant_id, dataset = "gnomad_r2_1"):
    """
       In this part of the code, the GraphQL query for variant will be made. This query will be formated into json,
       which will be pushed to the fetch code.

       :param gene_id: variant id
       :param dataset: dataset in which the API will search.
       :return: variant info in json format
       """
    #graphQL for variant data
    fmt_graphql2 = """
    {
      variant(variantId: "%s", dataset: %s) {
        variantId
        reference_genome
        chrom
        pos
        ref
        alt
        colocatedVariants
        sortedTranscriptConsequences {
          gene_symbol
          major_consequence
          transcript_id
          canonical
          polyphen_prediction
          sift_prediction
          lof
        }
      }
    }
    """

    # This part will be JSON encoded, but with the GraphQL part left as a
    # glob of text.
    req_variant_list = {
        "query": fmt_graphql2 % (variant_id, dataset),
        "variables": {}
    }
    response = fetch(req_variant_list)["data"]["variant"]
    return response


def make_db(con, cursor):
    """
    In this function, the local database is created.
    :param con: connection to local database
    :return:
    """
    #SQL for creating db
    cursor.execute("""
            CREATE TABLE IF NOT EXISTS mutation (
                database_id INTEGER,
                variant_id TEXT,
                rs_id TEXT,
                ref TEXT,
                alt TEXT,
                pos INTEGER,
                transcript_id TEXT,
                ref_genome TEXT,
                chrom VARCHAR(2),
                major_consequence TEXT,
                polyphen_prediction TEXT,
                sift_prediction TEXT,
                loss_of_function TEXT,
                PRIMARY KEY(variant_id, rs_id)
            );
        """)
    con.commit()


if __name__ == "__main__":
    input_info = snakemake.config["ensembl"]
    count = 0
    num_mutations = snakemake.config["number_of_mutations"]

    rs_codes = [i.strip() for i in snakemake.config["rs_codes"].split(",")]

    test = get_variant_list(input_info)  # code for htt gene

    for i in range(len(test)):
        if test[i]['rsid'] in rs_codes:
            test.insert(0, test.pop(i))

    if type(num_mutations) == int:
        test = test[:num_mutations]


    with tqdm(total=len(test)) as pbar:
        with open(snakemake.log[0], "w") as f:
            sys.stderr = sys.stdout = f
            output_info = snakemake.config["db"]+".db"

            print("created logging file")
            print('connecting to database')
            con = sl.connect(output_info)
            print("con done")

            print("Create cursor")
            cursor = con.cursor()
            print("cursor created")

            print('creating database')
            make_db(con, cursor)
            print("created database")

            #get all variants of HTT gene
            print('Get variants for gene')

            #loop through all variants
            print("loop through variant list")
            for i in test:
                print("====New Variant====")
                count+=1
                print("start time:", datetime.utcnow())

                #get variant_id and rs_id, and put in sql_data list.
                sql_data = [count, i['variant_id'], i['rsid']]
                variant_id= i['variant_id']

                #get all of the variant data
                print("Get data from variant")
                variant_data = get_variant(variant_id)

                #get variant data for db in list (what I think is necessary)
                for i in variant_data['ref'], variant_data['alt'], variant_data['pos'], variant_data['sortedTranscriptConsequences'][0]['transcript_id'], variant_data['reference_genome'], variant_data['chrom'], variant_data['sortedTranscriptConsequences'][0]['major_consequence'], variant_data['sortedTranscriptConsequences'][0]['polyphen_prediction'], variant_data['sortedTranscriptConsequences'][0]['sift_prediction'], variant_data['sortedTranscriptConsequences'][0]['lof']:
                    sql_data.append(i)

                print("Insert data into sql database")
                #sql query for inserting in db
                sql = 'INSERT OR IGNORE INTO mutation VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'
                cursor.execute(sql, sql_data)
                con.commit()
                pbar.update(1)

