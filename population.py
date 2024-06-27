"""
This class generates the diverse synthetic DNA sequences of each
haplotype per population.
"""

# The functions used in the Population class
from sequence_generator import check_chromosome_file, create_directory, \
    create_synthetic_DNA_population, recombine, extract_existing_rs_combs


class Population:

    def __init__(self, pop_dict, data_config):
        """
        The init of this class, initializes all the variables
        from this class.
        :param pop_dict: The dictionary of each population.
        :param data_config: The config file.
        """
        self.pop_dict = pop_dict
        self.data_config = data_config
        self.file_nc = self.nc_file(self.data_config)
        self.path_rs = self.rs_path(self.data_config)
        # self.rsc = self.rsc()
        self.combs = self.create_population()

    def get_pop_dict(self):
        """
        This function is able to extract the population dictionary.
        The information in this dictionary is the population name,
        percentage, haplotypes, rs codes per haplotype, allele per rs
        code and the rs code and the percentages of the haplotypes.
        :return: The population dictionary.
        """
        return self.pop_dict

    def set_pop_dict(self, pop_dict):
        """
        This function can be used for the setting/editing
        of the population dictionary.
        :param pop_dict: The new population dict.
        """
        self.pop_dict = pop_dict

    def get_data_config(self):
        """
        This function is able to extract the name of the config file.
        :return: The config file name.
        """
        return self.data_config

    def set_data_config(self, data_config):
        """
        This function can be used for the setting/editing
        of the config file.
        :param data_config: The new name of the config file.
        """
        self.data_config = data_config

    def nc_file(self, data_config):
        """
        This function is able to check if the DNA sequence of the
        chromosome is already stored in a file. This is done by importing
        the check_chromosome_file function from the
        sequence_generator.py.
        :param data_config: Name of config file.
        :return: Name of the file with the DNA sequence of
        the chromosome.
        """
        nc_file = check_chromosome_file(f'{data_config["NC_code"]}.fasta',
                                        data_config["NC_code"])
        return nc_file

    def rs_path(self, data_config):
        """
        This function uses the function imported from
        sequence_generator to create a directory to save
        the XML files of the rs codes from NCBI.
        :param data_config: Name of config file.
        :return: Name of the directory where the RS code XML files
        are saved.
        """
        rs_path = create_directory(f'rs_{data_config["gene"]}_xml_files')
        return rs_path

    def rsc(self):
        """
        This function uses the function imported from
        sequence_generator to create a list of all
        rs codes per haplotype.
        :return: 2D list of rs codes per haplotype.
        """
        rsc = extract_existing_rs_combs(self.pop_dict)
        return rsc

    def create_population(self):
        """
        This function uses several functions imported from
        sequence_generator, which is used to create the initial
        population, to extract the existing RS code combinations,
        create new combinations, and check if these combinations
        already exist.
        :return: Combination list with possible RS code combinations.
        """
        list_recombine = create_synthetic_DNA_population(self.pop_dict,
                                                         self.data_config,
                                                         self.file_nc,
                                                         self.path_rs)
        return list_recombine

    def recombine(self):
        """
        This function uses the imported function to recombine
        the known rs codes between haplotypes and create new
        DNA sequences with the combined haplotypes.
        """
        recombine(self.combs, self.rsc(), self.pop_dict.get("population name"),
                  f'{self.data_config["gene"]}_result.fasta',
                  self.file_nc, self.data_config,
                  f"{'_'.join([i for i in self.pop_dict.get('haplotypes')])}")
