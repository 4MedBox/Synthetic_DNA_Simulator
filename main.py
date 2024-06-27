"""
This script executes all the code from the Python scripts necessary
to create a synthetic DNA dataset. This script can create a config
file if the user wishes, or it will use the data from an existing
config file.
The packages/libraries used are
json: Used to read and write JSON configuration files
csv: Used to add metadata to CSV file.
argparse: Used to pass arguments/parameters to a Python script
when running it from the command line.
time: Used to see how long it takes the script to
create all populations
"""

# # Import gnomad function
# from gnomad import create_gnomad_db
# Import functions from user_input script to create config file
from user_input import User_input, input_haplotype
# Import population class to create population
from population import Population
import json
import csv
import time
import argparse
# Import the recombination function from sequence_generator.py
from sequence_generator import recombine


def user_input():
    """
    This function takes the user input and creates a user
    input object.
    The haplotypes are added to the class using the
    input_haplotype function, which is read from the
    user_input.py file.
    Return: Object User_Input named user.
    """
    populations = list(input(f'What populations would you like '
                             f'to add?\nFormat example: '
                             f'population_name_1, population_name2,'
                             f' etc..\n').split(', '))
    input_haplotype_name = []
    input_haplotype_perc = []
    input_snp = []
    percentage_populations = []
    allele_input = []
    for pop in populations:
        while True:
            try:
                perc_pop = int(input(f'What are the percentages of '
                                     f'population {pop} '
                                     f'(0-100 THIS MUST BE AN INTEGER).\n'
                                     f'This percentage is before '
                                     f'recombination meaning the '
                                     f'percentage input is not the '
                                     f'same at the end it can be more!\n'))
                if 0 < perc_pop < 101:
                    percentage_populations.append(perc_pop)
                    break
                else:
                    print('This number is not between 0 and 100!')
                    continue
            except [TypeError, ValueError]:
                print('This is not a number or an integer')
        haplo_name, haplo_perc, snp, allele = input_haplotype(pop)
        input_haplotype_name.append(haplo_name)
        input_haplotype_perc.append(haplo_perc)
        input_snp.append(snp)
        allele_input.append(allele)
    input_user = User_input(populations, percentage_populations,
                            input_haplotype_name, input_haplotype_perc,
                            input_snp, allele_input, "ENSG00000197386",
                            "HTT_db")
    return input_user


def create_config_dict(in_user):
    """
    This function creates the dictionary for the config json file.
    :param in_user: Object of the user input.
    :return: A dictionary with the input from the user for the
    config file.
    """
    config_dict = {"CAG repeats": {"min_cag": in_user.get_cag_min(),
                                   "max_cag": in_user.get_cag_max()},
                   "CAG sub": in_user.get_add_cag_sub(),
                   "CAA": in_user.get_caa(),
                   "recombine": in_user.get_recombine(),
                   "recombine_pop": in_user.get_recombine_pop(),
                   "number_recombine": in_user.get_number_recombine(),
                   "gene_or_subseq": in_user.get_subset_gene(),
                   "ensemble": in_user.get_ensemble(),
                   "gene": in_user.get_gene(),
                   "NC_code": in_user.get_nc_code(),
                   "db": in_user.get_db(),
                   "number_sequences": in_user.get_number_sequences(),
                   "treads": in_user.get_treads()}
    x = 0
    for pop_name, pop_perc, hap_name, hap_perc, \
            rs, allele in zip(in_user.get_populations(),
                              in_user.get_percentage_populations(),
                              in_user.get_haplotypes(),
                              in_user.get_percentage_haplotypes(),
                              in_user.get_snps(), in_user.get_allel()):
        x += 1
        dict_hap = {}
        for hn, hp, frs, fa in zip(hap_name, hap_perc, rs, allele):
            list_snp = []
            for r, a in zip(frs, fa):
                dict_snp = {"rs_code": r, "allele": a}
                list_snp.append(dict_snp)
            dict_hap[hn] = {"percentage": hp, "SNPs": list_snp}
        config_dict[f"population {x}"] = {"population name": pop_name,
                                          "percentage": pop_perc,
                                          "haplotypes": dict_hap}
    return config_dict


def writing_config_file(config_dict, config_file):
    """
    This function writes the config file with the user input.
    :param config_dict: Dictionary with the input from the user.
    :param config_file: Name of the config file (string).
    :return: The name of the config file.
    """
    with open(config_file, "w") as outfile:
        json.dump(config_dict, outfile)
    outfile.close()
    return str(config_file)


def file_from_disk(type_file):
    """
    This function gives the user the option to add a
    consisting config file.
    It also checks if this file exists in the working directory.
    :param type_file: The file name.
    :return: The name of the file.
    """
    while True:
        file = input(f'Give the name of the {type_file} file from disk: ')
        try:
            with open(file, 'r'):
                return file
        except FileNotFoundError:
            print('The file was not found!\n')
            continue


def create_config_file():
    """
    This function allows the user to create a new config file
    or to add an existing one.
    If the user decides to create a file,
    it will ask the input form the
    user (user_input()),
    create a config dictionary (create_config_dic)
    and then create the config file (writing_config_file).
    If the user wants to add an already existing config file,
    then the function config_file_from_disk is used.
    :return: The name of the config file.
    """
    while True:
        options = input("Would you like to create your config file (c)"
                        " or do you have one (h(THIS NEEDS TO BE A JSON "
                        "FILE!!))? ").lower()
        if options in ['c', 'h']:
            if options == 'c':
                in_user = user_input()
                config_dict = create_config_dict(in_user)
                file_c = writing_config_file(config_dict, "config.json")
                return file_c
            else:
                file_h = file_from_disk("config")
                return file_h
        else:
            continue


def config_file():
    """
    This function calls to create the config file or determines if
    the user already has a config file they would like to use.
    The while loop is used to check if the file exits and readable.
    There is also a check within the while loop if the data
    within the config file is correct.
    :return: The name of the config file (string).
    """
    file_config = create_config_file()
    while True:
        with open(file_config, 'r') as fa:
            for line in fa:
                print(line)
        check_file_input = input("Is the config file correct yes (y) "
                                 "or no (n): ").lower()
        if check_file_input in ['y', 'n']:
            if check_file_input == 'y':
                return file_config
            else:
                file_config = create_config_file()
                return file_config
        else:
            continue


# Def gnomad_database(config_file):
# 	"""
# 	This function creates a database by calling on the
# 	create_gnomad_db function in gnomad.py if there is
# 	no database.
# 	If there is a database, The user will get a statement
# 	that there is a database file and a question if they
# 	want to use that.
# 	:Param config_file: The name of the config file (String).
# 	:Return: Name of the created database file (String).
# 	"""
# 	file = open(config_file)
# 	data_config = json.load(file)
# 	check_file = f'{os.getcwd()}\{data_config["db"] + ".db"}'
# 	isfile = os.path.isfile(check_file)
# 	if not isfile:
# 		db_file = create_gnomad_db(config_file)
# 		return db_file
# 	else:
# 		print(
# 			"There is an already existing database for the GnomAD variants")
# 		while True:
# 			create_db = input("Would you like to use the database that "
# 							  "already exists (e) or would you like to "
# 							  "create a new one (n): ").lower()
# 			if create_db in ['e', 'n']:
# 				if create_db == 'n':
# 					db_file = create_gnomad_db(config_file)
# 					return db_file
# 				else:
# 					db_file = file_from_disk("database")
# 					return db_file

def populaties(data_config):
    """
    This function creates the synthetic DNA for each population,
    adding it to the output file. This function also shows the time
    it takes for it to run.
    :param data_config: The config file.
    """
    startTime = time.time()
    pop = [i for i in data_config if 'population' in i]
    all_rs_pop_combs = []
    all_mut = []
    for i in pop:
        x = Population(data_config[i], data_config)
        all_rs_pop_combs.extend(x.rsc())
        all_mut.extend(x.combs)
        if data_config["recombine_pop"] == 'w' \
                and data_config["recombine"] == 'y':
            x.recombine()
    if data_config["recombine_pop"] == 'b' and data_config["recombine"] == 'y':
        recombine(all_mut, all_rs_pop_combs, "com_pop",
                  f"{data_config['gene']}_result.fasta",
                  f"{data_config['NC_code']}.fasta", data_config,
                  "com_hap")
    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))


def create_tsv_file_metatdata():
    """
    This function creates a tsv file with metadata about the samples.
    The first column in this file is the sample name without any
    added symbols.
    The second column is the name of the population
    """
    row_htt_NCBI_header = [['sample', 'pop'],
                           ['NC000004.123074681-3243960',
                            'NC000004.123074681-3243960']]
    with open('populations.tsv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t')
        csvwriter.writerows(row_htt_NCBI_header)


if __name__ == "__main__":
    # Here the parser argument is created
    parser = argparse.ArgumentParser(
        description='Create synthetic DNA dataset')
    # Here are the arguments specified for the config file
    # and database file
    parser.add_argument('--config_file',
                        metavar='cf',
                        type=str,
                        help='This needs to be a json config file\n'
                             'With data about:\n'
                             '- min amount of CAG repeats\n'
                             '- max amount of CAG repeats\n'
                             '- CAG repeat in a subsequence\n'
                             '- removal of CAA\n'
                             '- recombination between or within '
                             'population(s)\n'
                             '- amount of sequences for recombination\n'
                             '- entire gene (g) or subsequence (s)\n'
                             '- ensemble ID\n- gene name\n- db name\n'
                             '- nc code\n- amount of treats\n'
                             '- max amount of sequences\n'
                             '- populations and percentages\n'
                             '- haplotypes and percentages\n'
                             '- SNPs per haplotype')
    # parser.add_argument('--database',
    # 					metavar='db',
    # 					type=str,
    # 					help='This is a created db from GnomAD data')
    # Saving the arguments
    args = parser.parse_args()
    if args.config_file is None:  # or args.database is None
        file_config = config_file()
    # print(file_config)
    # if args.config_file is None:
    # 	# create config file
    # 	file_config = config_file()
    # 	print(file_config)
    # elif args.database is None:
    # 	# name of the config file
    # 	file_config = args.config_file
    # 	print(file_config)
    # 	# create the database
    # 	database = gnomad_database(file_config)
    # 	print(database)
    # elif args.config_file is None and args.database is None:
    # 	file_config = config_file()
    # 	print(file_config)
    # 	database = gnomad_database(file_config)
    # 	print(database)
    else:
        # name of the config file
        file_config = args.config_file
    # # name of the database file
    # database = args.database
    # print(file_config)
    # print(database)

    # Open config file
    f = open(file_config)
    # Load the data from the config file
    data_config = json.load(f)
    # Start creating metadat file
    create_tsv_file_metatdata()
    # Create the populations
    populaties(data_config)
