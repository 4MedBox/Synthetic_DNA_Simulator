"""
This script is used to generate the mutated sequences per haplotype
and check if they have the correct CAG repeat interval.
By importing functions from check_CAA_CAG_HD.py.
These functions will be used later to create populations
with different haplotypes.
The packages/libraries used are:
biopython: For extracting the chromosome sequence and snp
information from NCBI.
os: Used to work in different directories and to create one.
xml.etree.elementTree: Used to parse XML strings and create XML files.
Math: Used to round numbers.
random: Used to generate random interchange combinations of
mutations in and between
in and between populations and used to create a seed.
"""
# Used to recombine the SNP mutations.
# import itertools
# import sqlite3 as sl
# import sys
# import Thread
# import datetime
# import tqdm
# import itertools
# Used to read the config json file and extract data from it.
#import json
# from Thread import NewThread #, mutate
# Time is used
# to determine the number of seconds it took to run the scrip
#import time
from Bio import Entrez, SeqIO
import os
import xml.etree.ElementTree as ET
# mutate_ref_sequence.py to create/mutate multiple sequences at ones.
import mutate_ref_sequence
# This script test if the created sequence is huntington and
# has the correct number of CAG repeats
import check_CAA_CAG_HD
import math
import random
import csv


def create_chromosome_seq_file(nc_code, file_name_chrom):
    """
    This function uses the Biopython package to extract the
    sequence of the chromosome using the NC code.
    This DNA sequence is then saved this sequence to a file.
    :param nc_code: The NC code.
    :param file_name_chrom: A string of the file name for the
    chromosome sequence.
    :return: The file/name of the file with the sequence of the
    chromosome.
    """
    handle = Entrez.efetch(db='nuccore', id=nc_code,
                           rettype='fasta')
    record = SeqIO.read(handle, 'fasta')
    seq_chrom = record.seq
    file_chromosome = open(file_name_chrom, 'w')
    file_chromosome.write(str(seq_chrom))
    file_chromosome.close()
    return file_chromosome


def check_chromosome_file(nc_file, nc_code):
    """
    This function checks if the chromosome sequence file
    already exists. If it does not exist, it will create one.
    If it does exist, the user will get a warning.
    :param nc_file: File name of the chromosome sequence.
    :param nc_code: NC code ID.
    :return: The name of the file containing the chromosome sequence.
    (string).
    """
    check_file = f'{os.getcwd()}\{nc_file}'
    isfile = os.path.isfile(check_file)
    if not isfile:
        nc_file = create_chromosome_seq_file(nc_code, nc_file)
        return nc_file
    else:
        print(f"There is an already existing NC file!")
        return nc_file


# Def mutate(data_config, nc_file, amount_seq=2):
# 	"""
# 	NOT IN USE!!!!
# 	This function creates synthetic sequences with the use of
# 	GnomAD database. These sequences are written to the output file.
# 	:Param data_config: The config file.
# 	:Param nc_file: File with the sequence of the chromosome.
# 	:Param amount_seq: Number of sequences of the start population.
# 	"""
# 	With open('seq_gen.log', "w") as f:
# 		f.write('connecting to database and creating cursor\n')
#
# 		#the following code connects to a database and creates a cursor.
# 		con = sl.connect(data_config["db"] + ".db")
# 		cursor = con.cursor()
#
# 		cursor.execute("SELECT * FROM mutation")
# 		# #res = cursor.execute(f"SELECT * FROM mutation WHERE 'rs1742354048' in (rs_id);")
# 		# rows = res.fetchall()
# 		# for row in rows:
# 		# 	print(row)
#
# 		f.write('generate starting population\n')
#
# 		#the following code generates the starting population. In this case: there will be a population with combinations of 20 different mutations.
# #	lis = list(range(1,amount_seq+1)) #snakemake.config['number_of_mutations']
#
# 		f.write('get start and end position\n')
# 		#the following selects the highest and lowest position from N mutations, making sure only a part of the gene will be used for generating a new sequence, instead of the whole gene.
# 		cursor.execute("SELECT MIN(pos), MAX(pos) FROM mutation ")
# 		start, end = [i for i in cursor.fetchall()[0]]
#
# 		f.write('get part of sequence based on start and end position\n')
# 		sequence = open(nc_file, "r").readlines()[0][start:end+1]
#
# 		#a list will be created to keep all the mutations in. This also counts for the file.
# 		f.write('open result-file and population list for mutations\n')
# 		population=0
# 		file = open(data_config["gene"]+"_result.fasta", "a")
#
# 		lijst_threads=list(range(1,data_config["treads"]+1)) #4 threads
#
# 		count=1
# 		f.write('loop through starting population, starting with 1 mutation and adding 1 mutation with each loop\n\n')
#
# 	#for loop: every amount of mutations will be used (1 mutation, 2 mutations, 3 mutations etc.)
# 		while count <= lis[len(lis)-1]:
# 			# print("ja")
# 			f.write('====New Number of Possible Mutations===\n')
# 			f.write('Get all possible combinations with the to be used mutations\n')
# 			#with the N to be used mutations, all the possible combinations will be calculated.
# 			test = list(itertools.combinations(lis, count))
# 			# print(test)
# 			# exit()
#
# 			pbar_len, test_len = len(list(itertools.combinations(lis, count))), len(list(itertools.combinations(lis, count)))
#
# 			#loop through all the possible combinations of mutations.
# 			F.write('loop through all possible combinations of mutations.\n\n')
#
# 			with tqdm(total=pbar_len, desc=str(count)+" mutaties") as pbar:
# 				while test_len > 0:
# 					# print("ja")
# 					f.write('===New Possible Mutated Sequence===\n')
#
# 					#get info from all the mutations in a certain combination.
# 					f.write('get info from all mutations in combination\n')
#
# 					fetch_list = [cursor.execute(
# 						"SELECT database_id, variant_id, ref, alt, pos, rs_id FROM mutation WHERE database_id IN {}".format(i+(0,))).fetchall() for i in test[:len(lijst_threads)] ]
#
# 					f.write('create synthetic person with sequence\n')
# 					threads = [Thread.myThread(number+1, "Thread-" + str(number+1), sequence, fetch_list[number]) for number in
# 							   range(len(test[:len(lijst_threads)]))]
#
# 					result = [[thread.mutations, thread.new_sequence, thread.name] for thread in threads]
#
# 					for gene in result:
#
# 						if data_config["gene"].upper() == "HTT":
# 							survivability.Survival(gene[1])
#
# 						f.write("found 2 repeats\n")
# 						population+=1
#
# 						f.write('Add synthetic person to population.\n')
#
# 						#create (meta-)data for result file, and write to result file.
# 						f.write('Write results to file\n\n')
# 						data = [">seq_",str(population),", (",', '.join(f'{x[::-1][0]}'for x in gene[0]),"), ", str(datetime.datetime.now()),"\n",gene[1],"\n\n"]
# 						file.writelines(data)
# 						f.write("done\n")
#
# 						if type(data_config["number_sequences"]) == int:
# 							if population == data_config["number_sequences"]:
# 								# print(population,"ja")
# 								pbar.update(len(test))
# 								sys.exit()
# 							# exit()
# 					del test[:len(lijst_threads)]
# 					# print(test_len, len(threads))
# 					test_len -= len(threads)
# 					pbar.update(len(lijst_threads))
#
# 				pbar.close()
# 			count += 1
# 			# exit()
#
# 		#print end time (you'll likely walk away while running the program, it takes some time :) )
# 		f.write('Script is done.\n')
# 		f.write(str(datetime.datetime.now())+"\n")
# 	f.close()


def create_directory(name_directory):
    """
    This function creates a directory for the XML files that will
    be downloaded from NCBI.
    :param name_directory: String with the name of the directory.
    :return: String of the path to the directory containing the
    RS XML files.
    """
    global path
    try:
        cwd = os.getcwd()
        path = os.path.join(cwd, name_directory)
        os.mkdir(path)
        return path
    except FileExistsError:
        return path


def extract_mutate_info_ncbi(rs_code, dir):
    """
    This function extracts the SNP information using the rs code and
    efetch Biopython. This is then saved in an XML file
    which is stored in a different directory.
    :param rs_code: The rs_code as ID for the efetch.
    :param dir: Name and path of the directory where the XML
    files for each
    are located for each rs_code.
    :return: String with the XML info of the rs_code.
    """
    # change directory
    Entrez.mail = "A.N.Other@example.com"
    os.chdir(dir)
    handle = Entrez.efetch(db='SNP', id=rs_code, rettype='text')
    xml_str = handle.read()
    myroot = ET.fromstring(xml_str)
    tree = ET.ElementTree(myroot)
    tree.write(f'{rs_code}.xml')
    os.chdir("..")
    return myroot


def check_xml_file(rs_path, rs_xml_file, rs_code):
    """
    This function checks if there is already an excising XML file
    from the fetched XML files (Biopython) from NCBI. If not, the
    extract_mutate_info_ncbi function extracts the SNP info from
    NCBI and creates a file.
    If the file already exists, the file is read and an XML string
    is returned.
    :param rs_path: The pathname to the directory containing all the XML
    files.
    :param rs_xml_file: The name of the XML file containing
    the rs code.
    :param rs_code: String of the RS code.
    :return: The XML string.
    """
    dir_list = os.listdir(rs_path)
    if rs_xml_file in dir_list:
        os.chdir(rs_path)
        file = open(rs_xml_file, "r")
        xml_str = file.read()
        myroot = ET.fromstring(xml_str)
        os.chdir("..")
        return myroot
    else:
        myroot = extract_mutate_info_ncbi(rs_code, rs_path)
        return myroot


def extract_mutate_info_per_haplotype(lijst_snp_allele_hap,
                                      data_config, rs_path):
    """
    This function uses the list of SNPs and alleles (in dict) per
    haplotype to extract information and store it in a 2D list.
    The information is either downloaded (XML string) or read from a
    file (XML). This is done with the check_xml_file function.
    For each rs code, the list is created from the nc_code, the start
    position, the reference allele, the alternative allele, the end
    position, the type of mutation, and the rs_code. The end position for
    SNP/SNV is the same as the start position. For an insertion it is
    the start position + the length of the alternative allele, and
    for a deletion, the end position is the start position + the length
    of the reference sequence.
    :param lijst_snp_allele_hap: List of SNP and allele dictionary
    from a haplotype.
    :param data_config: Name of the configuration file.
    :param rs_path: The path to the directory containing the RS XML files.
    :return: Three lists are returned. The first is the mutation info
    of each RS code in the haplotype (list_mutation).
    The second list is the list of start positions (min_pos_list).
    The last list is the list of the end positions (max_pos_list).
    """
    list_mutations = []
    min_pos_list = []
    max_pos_list = []
    for dict_rs in lijst_snp_allele_hap:
        rs = dict_rs['rs_code']
        allele = dict_rs['allele']
        myroot = check_xml_file(rs_path, f'{rs}.xml', rs)
        spdi_txt = myroot.find('SPDI').text
        type_var = myroot.find('SNP_CLASS').text
        lijst_snp_allele = []
        if data_config["NC_code"] in spdi_txt:
            snp_allele = spdi_txt.strip().split(':')
            lijst_snp_allele.append(snp_allele[0])
            lijst_snp_allele.append(int(snp_allele[1]))
            lijst_snp_allele.append(snp_allele[2])
            lijst_snp_allele.append(allele)
            if len(snp_allele[2]) == len(snp_allele[3]):
                lijst_snp_allele.append(int(snp_allele[1]))
            elif len(snp_allele[2]) < len(snp_allele[3]):
                end = int(snp_allele[1]) + len(snp_allele[3])
                lijst_snp_allele.append(end)
            else:
                end = int(snp_allele[1]) + len(snp_allele[2])
                lijst_snp_allele.append(end)
            lijst_snp_allele.append(type_var)
            lijst_snp_allele.append(rs)
            list_mutations.append(lijst_snp_allele)
            min_pos_list.append(lijst_snp_allele[1])
            max_pos_list.append(lijst_snp_allele[4])
    return list_mutations, min_pos_list, max_pos_list


def get_sequence(nc_file, min_pos, max_pos):
    """
    This function reads the file containing the sequence of the
    chromosome and returns a certain part of the sequence.
    :param nc_file: The name of the chromosome sequence file.
    :param min_pos: The first position (int).
    :param max_pos: The last position (int).
    :return: The sequence and the first position (int) for the
    position correction in mutate_ref_sequence.py
    """
    sequence = open(nc_file, "r").readlines()[0][min_pos:max_pos]
    return sequence, min_pos


def obtain_seq_last_pos_and_CAG_repeat(data_config, max_len, begin):
    """
    This function determines whether the RS mutations are
    before the CAG repeat.
    If so, it will return the sequence between the last position and
    the last position and the CAG repeat.
    :param data_config: The config file.
    :param max_len: The maximum length of the sequence.
    Or end position.
    :param begin: The start position of the sequence.
    """
    sequence = open(f'{data_config["NC_code"]}.fasta', "r").readlines()[0][
               3074681:3211920]
    if max_len < 3074876:
        seq = sequence[max_len - begin:195]
        return seq
    else:
        seq = ''
        return seq


def mutate_sequence(data_config, file, pop_name, hap_name, amount,
                    lis_mutations, sequence, min_pos):
    """
    This function generates the synthetic DNA sequences by mutating the
    original sequence (mutate_ref_sequence). This is done per haplotype.
    If the gene is the HTT gene or a certain variant of the name, then
    it will also add or remove CAG repeats and the CAA codon
    if the user wants that (check_CAA_CAG_HD).
    These reference sequences (chromosome sequence) are then
    mutated and saved to a file.
    :param data_config: The configuration file.
    :param file: Name of the output file.
    :param pop_name: Name of the population from which the haplotype
    comes from.
    :param hap_name: Name of the haplotype.
    :param amount: Integer representing the number of sequences
    that need to be created.
    :param lis_mutations: List of rs mutations per haplotype.
    :param sequence: The DNA sequence of the chromosome
    to be mutated.
    :param min_pos: Integer with the lowest/starting position
    """
    global nc
    thread_num = 0
    lis_mutations = sorted(lis_mutations, key=lambda x: x[1], reverse=True)
    rs = []
    mut = []
    max_len = max([i[4] for i in lis_mutations])
    for item in lis_mutations:
        rs.append(item[5])
        mut.append(item[6])
        nc = item[0]
    for x in range(amount):
        name = f'{hap_name}_{x + 1}_{pop_name}'
        row = [name, pop_name]
        with open('populations.tsv', 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t')
            csvwriter.writerow(row)
        thread_id = f'Thread-{thread_num + 1}'
        t = mutate_ref_sequence.myThread(thread_id, name, sequence, lis_mutations, min_pos)
        if data_config["gene"] in ['HTT', 'htt', 'Huntingtin', 'huntingtin']:
            t_seq = obtain_seq_last_pos_and_CAG_repeat(data_config, max_len, min_pos)
            seq, cag = check_CAA_CAG_HD.sequence(t.new_sequence,
                                                 data_config["CAG repeats"]
                                              ["min_cag"],
                                                 data_config["CAG repeats"]
                                              ["max_cag"],
                                                 data_config["gene_or_subseq"],
                                                 data_config["CAG sub"],
                                                 data_config["CAA"], t_seq)
        else:
            seq = t.new_sequence
            cag = ''
        f = open(file, "a")
        f.write(f">{name}, {pop_name}, {hap_name}, {nc}, "
                f"{mut}, {rs}, sequence length: {len(seq)}, "
                f"{cag}\n{seq}\n")
        f.close()


def gene_or_subsequence(data_config, min_pos_list, max_pos_list):
    """
    This function sets the start and end positions for the reference
     sequence. This is for the HTT gene! This can be changed by
    changing the min_pos_1 and max_pos in the if statement
    to the start and end position of the selected gene on the chromosome.
    In this case, chromosome 4 was used.
    :param data_config: The configuration file.
    :param min_pos_list: List of SNP start positions.
    :param max_pos_list: List of SNP end positions.
    :return: Two integers representing the start (min) position and the
    end position (max).
    """
    if data_config["gene_or_subseq"].lower() == 'g' and \
            data_config["gene"] in ['HTT', 'htt', 'Huntingtin', 'huntingtin']:
        min_pos = 3074681
        max_pos = 3243960
    else:
        if data_config["CAG sub"] == 'y':
            min_pos = 3074681
        else:
            min_pos = min(min_pos_list)
        max_pos = max(max_pos_list)
    return min_pos, max_pos


def create_synthetic_DNA_population(population, data_config, nc_file, rs_path):
    """
    This function calls all the functions necessary to
    initiate the synthetic DNA population.
    :param population: The dictionary of population data from the
    config file.
    :param data_config: The config file.
    :param nc_file: The name of the file containing the chromosome
    sequence.
    :param rs_path: The directory where the XML files of the
    RS codes are saved.
    :return: A 2D list of all mutations.
    """
    haplotype = population.get("haplotypes")
    pop_name = population.get("population name")
    amount = math.ceil(data_config["number_sequences"] * (
            population.get("percentage") / 100))
    list_recombine = []
    for item in haplotype:
        snp_haplotype = haplotype.get(item).get('SNPs')
        factor = math.ceil((haplotype.get(item).get('percentage') / 100)
                           * amount)
        list_mutate_info, min_pos_list, max_pos_list = \
            extract_mutate_info_per_haplotype(snp_haplotype, data_config,
                                              rs_path)
        list_recombine.extend(list_mutate_info)
        min_pos_1, max_pos = gene_or_subsequence(data_config, min_pos_list,
                                                 max_pos_list)
        sequence, min_position = get_sequence(nc_file, min_pos_1, max_pos)
        mutate_sequence(data_config, f'{data_config["gene"]}_result.fasta',
                        pop_name, item,
                        factor, list_mutate_info, sequence, min_position) #1
    return list_recombine


def extract_existing_rs_combs(population):
    """
    This function extracts all RS codes used in this
    population.
    :param population: Dictionary of the populations.
    :return: List of all rs codes.
    """
    rsc = []
    for z in population.get("haplotypes"):
        k = [s.get("rs_code") for s in
             population.get("haplotypes").get(z).get('SNPs')]
        rsc.append(k)
    return rsc

def random_combination(iterable, len_new_com):
    """
    This function randomly generates a new combination of RS
    codes by selecting items from the input list of RS codes.
    :param iterable: The input list of RS codes that can be used
    to create new combinations.
    :param len_new_com: The length of the new combination.
    :return: A 2D list of new combinations.
    """
    pool = tuple(iterable)
    len_pool = len(pool)
    sorting_com = sorted(random.sample(range(len_pool), len_new_com))
    return list(pool[i] for i in sorting_com)


def recombine(list_mutations, rsc, population, file_name, nc_file, data_config,
              hap_name):
    """
    This function recombines all mutations of the populations and
    generates new sequences with these combinations.
    This function also checks if the combination already exists.
    :param list_mutations: A 2D list of all mutations.
    :param rsc: A 2D list of rs codes per haplotype.
    :param population: The name of the population.
    :param file_name: Name of the output file.
    :param nc_file: The Name of the NC file containing the DNA
    sequence of the chromosome.
    :param data_config: The configuration file.
    :param hap_name: Combined name of the haplotypes.
    """
    random.seed(5)
    for x in range(data_config["number_recombine"]):
        item = random_combination(list_mutations,
                                  random.randint(1, len(list_mutations)))
        crs = [i[6] for i in item]
        if crs not in rsc and len(item) > 1:
            min_pos_list = [i[1] for i in item]
            max_pos_list = [i[4] for i in item]
            min_pos_1, max_pos = gene_or_subsequence(data_config,
                                                     min_pos_list,
                                                     max_pos_list)
            sequence, min_position = get_sequence(nc_file, min_pos_1,
                                                  max_pos)
            mutate_sequence(data_config, file_name, population,
                            f'{hap_name}_{x}', 1, item, sequence, min_position)


# if __name__ == '__main__':
# 	startTime = time.time()
# 	# open config file
# 	file = open('config.json')
# 	data_config = json.load(file)
# 	# Here the NC file is created with sequence of the chromosome
# 	nc_file = check_chromosome_file(f'{data_config["NC_code"]}.fasta',
#									data_config["NC_code"])
# 	# Function of creating a directory is called returning a path name
# 	rs_path = create_directory(f'rs_{data_config["gene"]}_xml_files')
# 	population = dict({"population name": "p1", "percentage": 100,
# 					   "haplotypes": {"haplotype_p11": {"percentage": 50,
# 														"SNPs": [{
# 																	 "rs_code": "rs1387211216",
# 																	 "allele": "n"},
# 																 {
# 																	 "rs_code": "rs764562793",
# 																	 "allele": "n"}]},
# 									  "haplotype_p12": {"percentage": 50,
# 														"SNPs": [{
# 																	 "rs_code": "rs1712379289",
# 																	 "allele": "nnnn"},
# 																 {
# 																	 "rs_code": "rs1229057625",
# 																	 "allele": "n"}]}}})
#
# 	# Run the function to create synthetic DNA population
# 	list_recombine = create_synthetic_DNA_population(population, data_config,
# 													 nc_file, rs_path)
# 	# rsc = extract_existing_rs_combs(population)
# 	# combinations = create_all_possible_recombinations(list_recombine)
# 	# combs = check_if_comb_already_exists(combinations, rsc)
# 	# Here the function to recombine the mutations and add them to the already
# 	# created population
# 	# p = population.get("population name")
# 	# recombine(list_recombine, rsc, p,
# 	# 		  f'{data_config["gene"]}_result.fasta', nc_file, data_config,
# 	# 		  f"{'_'.join([i for i in population.get('haplotypes')])}")
# 	executionTime = (time.time() - startTime)
# 	print('Execution time in seconds: ' + str(executionTime))
