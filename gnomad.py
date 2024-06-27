"""
This script extract all the variants from GnomAD and creates
a database with these variants.
The variants are extracted using the ensemble ID from the input.
and the rs codes from the user are the first few items in the
database.
The data extracted from each variant are the variant ID, rs number,
consequence, the reference allele, alternative allele, chromosome and
the reference genomen.
Packages/libraries used are:
time: Used to let to program sleep.
sqlite3: Used to create and push command to generated database.
datetime: Used to calculate how much time the retrieval from data
          and push to database is.
requests: Used to send requests to the GNOMAD api.
json: Used to extract data from the config.json file.
"""

#!/usr/bin/env python

# import everything that is necessary in order to run the script.
import time
import requests
import json
import sys
from tqdm import tqdm
import sqlite3 as sl
from datetime import datetime


def fetch(data, url="https://gnomad.broadinstitute.org/api"):
	"""
    This part of the code fetches the data from the gnomad database.
    The online API is a GraphQL query, which is also
    the data given to the fetch-code.The request for the API has the
    url of the database, the json data and headers.
    These headers are needed, otherwise the API will be overloaded
    much faster. If the response is 200, the json data
    will be retrieved and returned. If there is a 400 or 502 response,
    the API refuses the request due to an API
    overload or too much API requests in a short time frame
    (such as an DDOS attack, which is not preferred). If that
    is the case, there will be a 30 second wait time to clear the API.
    After this, there will be a new chance to
    retrieve the data. The json format of this data will be returned.

    :param data: json data with GraphQL query with variant or
    ENSEMBL id
    :param url: API of the gnomad GrapQL database
    :return: json: returns json file
    """
	headers = {'User-Agent': 'Mozilla/5.0'}
	response = requests.post(url, json=data, headers=headers)
	
	# if there is an error with API overload: show it here.
	if "200" in str(response):
		if "json" in response.headers.get('content-type'):
			json = response.json()
	else:
		if "400" or "502" in str(response):
			time.sleep(30)
			json = fetch(data)
	return json


def get_variant_list(gene_id, dataset="gnomad_r2_1"):
	"""
    In this part of the code, the GraphQL query for the variant
    list of the gene will be made. This query will be
    formatted into json, which will be pushed to the fetch code.

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
	# This part will be JSON encoded, but with the
	# GraphQL part left as a glob of text.
	req_variant_list = {
		"query": fmt_graphql % (gene_id, dataset),
		"variables": {}
	}
	response = fetch(req_variant_list)
	return response["data"]["gene"]["variants"]


def get_variant(variant_id, dataset="gnomad_r2_1"):
	"""
       In this part of the code, the GraphQL query for variant will
       be made. This query will be formatted into json,
       which will be pushed to the fetch code.

       :param variant_id: variant id
       :param dataset: dataset in which the API will search.
       :return: variant info in json format
       """
	# graphQL for variant data
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
	
	# This part will be JSON encoded, but with the GraphQL
	# part left as a glob of text.
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
    :param cursor: The cursor for the creation of the db.
    :return:
    """
	# SQL for creating db
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


def create_gnomad_db(c_file):
	"""
	This function creates the database of all the variants of the
	ensemble id from GnomAD.
	"""
	file = open(str(c_file))
	data_config = json.load(file)

	input_info = data_config["ensemble"]
	count = 0

	rs_codes = [rs["rs_code"] for item in data_config if
				"population" in item for haplo in
				data_config[item]["haplotypes"] if "haplotype" in haplo
				for rs in data_config[item]["haplotypes"][haplo][
					"SNPs"]]

	test = get_variant_list(input_info)

	for i in range(len(test)):
		if test[i]['rsid'] in rs_codes:
			test.insert(0, test.pop(i))
			print(i)

	with tqdm(total=len(test)) as pbar:
		with open('logfile_GnomAD.log', "w") as f:
			sys.stderr = sys.stdout = f
			output_info = data_config["db"] + ".db"

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

			# get all variants of HTT gene
			print('Get variants for gene')

			# loop through all variants
			print("loop through variant list")
			for i in test:
				print("====New Variant====")
				count += 1
				print("start time:", datetime.utcnow())

				# get variant_id and rs_id, and put in sql_data list.
				sql_data = [count, i['variant_id'], i['rsid']]
				variant_id = i['variant_id']

				# get all of the variant data
				print("Get data from variant")
				variant_data = get_variant(variant_id)

				# get variant data for db in list (what I think is necessary)
				for i in variant_data['ref'], variant_data['alt'], \
						 variant_data['pos'], \
						 variant_data['sortedTranscriptConsequences'][
							 0]['transcript_id'], variant_data[
							 'reference_genome'], variant_data[
							 'chrom'], \
						 variant_data['sortedTranscriptConsequences'][
							 0]['major_consequence'], \
						 variant_data['sortedTranscriptConsequences'][
							 0]['polyphen_prediction'], \
						 variant_data['sortedTranscriptConsequences'][
							 0]['sift_prediction'], \
						 variant_data['sortedTranscriptConsequences'][
							 0]['lof']:
					sql_data.append(i)

				print("Insert data into sql database")
				# sql query for inserting in db
				sql = 'INSERT OR IGNORE INTO mutation VALUES(?, ?, ?, ?,' \
					  ' ?, ?, ?, ?, ?, ?, ?, ?, ?)'
				cursor.execute(sql, sql_data)
				con.commit()
				pbar.update(1)
	return output_info
