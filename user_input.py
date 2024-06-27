"""
This file asks the user for input and stores it in a class.
The input requested from the user is the populations,
percentages of those populations, haplotypes, percentages
of each haplotype, which RS codes are part of that haplotype
with allele, a minimum number of CAG repeats, and a maximum number
of CAG repeats, whether the user wants to add or remove the CAA codon,
whether the CAG repeat will be included in the subsequence,
the NC code, the ensemble ID, the number of threads,
the number of sequences, the name for the database, whether
the user wants to recombine RS codes within or between populations,
and whether the user wants the whole gene or a subsequence.
The packages/libraries used are
re: is used to extract the SNP or RS number from the
haplotype input file.
"""
import re


class User_input:
    """
    This class collects the user input, partly trough parameters
    and through input questions.
    """

    def __init__(self, populations, percentage_populations, haplotypes,
                 percentage_haplotypes, snps, allele, ensemble, db):
        """
        The init of this class, initializes all the variables
        from this class.
        :param populations: List of population names (strings).
        :param percentage_populations: List of floats, these are the
        percentages of each population.
        :param haplotypes: List op strings of the names of the
        haplotypes (2D list)
        :param percentage_haplotypes: List of floats with percentages
        of each haplotype (2D list).
        :param snps: 3D list of all the rs numbers per haplotype
        per population.
        :param db: Name of a database.
        :param allele: The SNP allele (list).
        :param ensemble: The ensemble ID.
        """
        self.check_list = ['h', 'f', 's', 'g', 'n', 'y', 'b', 'w']
        self.populations = populations
        self.percentage_populations = percentage_populations
        self.haplotypes = haplotypes
        self.percentage_haplotypes = percentage_haplotypes
        self.snps = snps
        self.allel = allele
        self.subset_gene = self.string_check(
            input('Do you want to have a '
                  'entire gene (g) or a '
                  'subsequence (s): '))
        if self.subset_gene == 's':
            self.add_cag_sub = self.string_check(
                input('Do you want to have a '
                      'CAG/CAA repeat in '
                      'the subsequence yes '
                      '(y) or no (n): '))
        else:
            self.add_cag_sub = 'NA'
        self.ensemble = ensemble
        self.gene = input('Which gene would you like to use: ')
        if self.gene in ['HTT', 'htt', 'Huntingtin', 'huntingtin']:
            self.cag_min = self.number_check(
                input('What is the minimum amount '
                      'of CAG repeats: '))
            self.cag_max = self.number_check(
                input('What is the maximum amount '
                      'of CAG repeats: '))
            self.caa = self.string_check(
                input('Would you like to have the '
                      'CAA codon removed yes (y) or '
                      'no (n): '))
        else:
            self.cag_min = 'NA'
            self.cag_max = 'NA'
            self.caa = 'NA'
        self.db = db
        self.nc_code = input('Which NC ID would you like to use: ')
        self.number_sequences = self.number_check(
            input("What number of "
                  "sequences would you "
                  "like to start with: "))
        self.treads = self.number_check(input("What amount of treads "
                                              "would you like to use: "))
        self.recombine = self.string_check(input("Would you like to "
                                                 "interchange the mutations "
                                                 "yes (y) or no (n): "))
        if self.recombine == 'y':
            self.recombine_pop = self.string_check(
                input("Would you like to "
                      "interchange the mutations "
                      "between (b) or within (w) "
                      "populations: "))
            self.number_recombine = self.number_check(
                input("What number "
                      "recombined "
                      "sequences would "
                      "you like (Note "
                      "that there can "
                      "be less sequences "
                      "then this given "
                      "amount in the "
                      "output because "
                      "the created "
                      "randomized "
                      "combination "
                      "already exists)"
                      ": "))
        else:
            self.recombine_pop = "NA"
            self.number_recombine = 0

    def get_populations(self):
        """
        Get population is able to extract the population list.
        :return: List of populations.
        """
        return self.populations

    def set_populations(self, populations):
        """
        This function can set/change the population.
        :param populations: Is a list of populations.
        """
        self.populations = populations

    def get_percentage_populations(self):
        """
        Get percentage population is able to extract the population
        list.
        :return: List of population percentages.
        """
        return self.percentage_populations

    def set_percentage_populations(self, percentage_populations):
        """
        This function can set/change the percentages of the population.
        :param percentage_populations: Is a list op the percentages
        per population.
        """
        self.percentage_populations = percentage_populations

    def get_haplotypes(self):
        """
        Get haplotypes is able to extract the haplotype list.
        :return: List of haplotypes.
        """
        return self.haplotypes

    def set_haplotypes(self, haplotypes):
        """
        This function can set/change the haplotypes of the populations.
        :param haplotypes: Is a list op haplotypes.
        """
        self.haplotypes = haplotypes

    def get_percentage_haplotypes(self):
        """
        Get percentage haplotypes is able to extract the
        percentage haplotypes list.
        :return: List of haplotypes per population.
        """
        return self.percentage_haplotypes

    def set_percentage_haplotypes(self, percentage_haplotypes):
        """
        This function can set/change the list with the percentages of
        the haplotypes.
        :param percentage_haplotypes: Is a list of haplotype
        percentages.
        """
        self.percentage_haplotypes = percentage_haplotypes

    def get_snps(self):
        """
        Get snps is able to extract the SNPs list.
        :return: List of rs number per haplotype per population.
        """
        return self.snps

    def set_snps(self, snps):
        """
        This function can set/change the list of SNPs.
        :param snps: Is a list op SNPs (3D).
        """
        self.snps = snps

    def get_allel(self):
        """
        Get snps is able to extract the SNPs list.
        :return: List of RS code allele per population.
        """
        return self.allel

    def set_allel(self, allel):
        """
        This function can set/change the list of SNPs.
        :param allel: Is a list op allele (3D).
        """
        self.allel = allel

    def get_cag_min(self):
        """
        Extracts the minimum number of CAG repeats.
        :return: Integer of minimum number of CAG repeats.
        """
        return self.cag_min

    def set_cag_min(self, cag_min):
        """
        This function can set/change the minimum number of CAG repeats.
        :param cag_min: The minimum number of CAG
        repeats (integer).
        """
        self.cag_min = cag_min

    def get_cag_max(self):
        """
        Extracts the maximum number of CAG repeats.
        :return: Integer of maximum number of CAG repeats.
        """
        return self.cag_max

    def set_cag_max(self, cag_max):
        """
        This function can set/change the maximum number of CAG repeats.
        :param cag_max: The maximum number of CAG
        repeats (integer).
        """
        self.cag_max = cag_max

    def get_subset_gene(self):
        """
        Extracts if the user wants the entire gene or subsequence of
        the gene.
        :return: String s (subsequence) or g (gene).
        """
        return self.subset_gene

    def set_subset_gene(self, subset_gene):
        """
        This function can set/change the option of the user wants a
        subsequence or the entire gene.
        :param subset_gene: Is string option s (subsequence) or g
        (gene).
        """
        self.subset_gene = subset_gene

    def get_ensemble(self):
        """
        Get ensemble is able to extract the ensemble ID.
        :return: List of populations.
        """
        return self.ensemble

    def set_ensemble(self, ensemble):
        """
        This function can set/change the ensemble ID.
        :param ensemble: Is a string of ensemble ID.
        """
        self.ensemble = ensemble

    def get_gene(self):
        """
        Get gene is able to extract the gene name.
        :return: String of gene name.
        """
        return self.gene

    def set_gene(self, gene):
        """
        This function can set/change the gene name.
        :param gene: Is a string of the gene name.
        """
        self.gene = gene

    def get_db(self):
        """
        Get gene is able to extract the db name.
        :return: String of db name.
        """
        return self.db

    def set_db(self, db):
        """
        This function can set/change the db name.
        :param db: Is a string of the db name.
        """
        self.db = db

    def get_nc_code(self):
        """
        Get NC code is able to extract the NC code.
        :return: String of the NC code.
        """
        return self.nc_code

    def set_nc_code(self, nc_code):
        """
        This function can set/change the NC code.
        :param nc_code: Is a string of the NC code.
        """
        self.nc_code = nc_code

    def get_number_sequences(self):
        """
        Get max number of sequences is able to extract the
        maximum number of sequences.
        :return: Integer of number of sequences.
        """
        return self.number_sequences

    def set_number_sequences(self, number_sequences):
        """
        This function can set/change the number of
        sequences.
        :param number_sequences: The max number
        of sequences (integer).
        """
        self.number_sequences = number_sequences

    def get_treads(self):
        """
        Get cores is able to extract the
        number of cores.
        :return: Integer of number of cores.
        """
        return self.treads

    def set_treads(self, treads):
        """
        This function can set/change the number of cores used.
        :param treads: Is an integer of number of cores.
        """
        self.treads = treads

    def get_add_cag_sub(self):
        """
        Get string if a user wants to add CAG repeats to sequence
        when it is a subsequence.
        :return: String y or n.
        """
        return self.add_cag_sub

    def set_add_cag_sub(self, add_cag_sub):
        """
        This function can set/change the user input for CAG
        repeats is subsequence.
        :param add_cag_sub: String y or n.
        """
        self.add_cag_sub = add_cag_sub

    def get_caa(self):
        """
        Get string value if a user wants to add keep the CAA
        in sequence.
        :return: String y or n.
        """
        return self.caa

    def set_caa(self, caa):
        """
        This function can set/change the user wants to keep the CAA
        ot remove it.
        :param caa: String y or n.
        """
        self.caa = caa

    def get_recombine(self):
        """
        Get string value if the user wants to recombine mutations.
        :return: String y or n.
        """
        return self.recombine

    def set_recombine(self, recombine):
        """
        This function can set/change the user input if the
        user wants recombines sequences.
        :param recombine: String y or n.
        """
        self.recombine = recombine

    def get_recombine_pop(self):
        """
        Get string value if a user wants to recombine mutations
        within a population or between them.
        :return: String b or w.
        """
        return self.recombine_pop

    def set_recombine_pop(self, recombine_pop):
        """
        This function can set/change the user input if the
        user wants to recombine mutations
        within a population or between them.
        :param recombine_pop: String b or w.
        """
        self.recombine_pop = recombine_pop

    def get_number_recombine(self):
        """
        Get the number of sequences the user wants to recombine.
        :return: The number of sequences (integer).
        """
        return self.number_recombine

    def set_number_recombine(self, number_recombine):
        """
        This function can set/change the user input if the
        user wants to alter the number of sequences that
        will be recombined.
        :param number_recombine: The number of sequences (integer).
        """
        self.number_recombine = number_recombine

    @staticmethod
    def number_check(number):
        """
        Checks if the input string is able to become an integer,
        if not, the user will get an error message.
        :param number: A string that will be transformed.
        :return: A float.
        """
        try:
            value = int(number)
            return value
        except ValueError:
            value = 4
            return value

    def string_check(self, string):
        """
        This function is used to check if a string is a string.
        Is a string and checks if it is a specific string or not.
        If it is not, an error is returned and the value s is added.
        Param String: The input is a string.
        :return: A string that is considered correct.
        """
        try:
            value = str(string.lower())
            if value in self.check_list:
                return value
            else:
                print('This is not one of the input options!')
                value = 's'
                return value
        except ValueError:
            print('This is not an string')
            value = 's'
            return value


def haplotype_file(pop):
    """
    This function is used when the user has an input file. It attempts
    to read the file after the user has specified the name as
    input. If the file cannot be opened, it continues the while
    loop. If the file can be opened, the haplotypes, the percentage
    of each haplotype, and the SNPs are placed in a list for later
    use in the user input class later. This function checks if the
    percentage is an interest and if the number is between 0 and 100.
    :param pop: String with the identity/name of the population.
    :return: A list of haplotype names, a list of percentages,
    a 2D list of
    percentages, a 2D list of SNP numbers per haplotype, and a 2D list
    of the allele per SNP.
    """
    while True:
        print('Make sure the format of the file (txt) is as followed\n'
              'This is per population!: \n'
              '- Add one haplotype on a newline:\n'
              '  haplotype name1, percentage (0-100 THIS MUST BE AN INTEGER), '
              'RS code 1 (allele (A, C, G and or T)), '
              'RS code 2 (allele (A, C, G and or T)), etc....\n'
              '  haplotype name2, percentage (0-100 THIS MUST BE AN INTEGER), '
              'RS code 1 (allele (A, C, G and or T)), '
              'RS code 2 (allele (A, C, G and or T)), etc....')
        file = input(f'give the name of the haplotype file for '
                     f'{pop} input from disk: ')
        haplo_name = []
        haplo_perc = []
        snp = []
        allele = []
        try:
            with open(file, 'r') as fa:
                for line in fa:
                    line_list = line.strip().split(', ')
                    haplo_name.append(line_list[0])
                    if 0 < float(line_list[1]) < 101:
                        haplo_perc.append(float(line_list[1]))
                    else:
                        haplo_perc.append('NA')
                    snp_t = [re.sub('[(\[].*?[)\]]', '', i).strip() for
                             i in line_list[2:]]
                    snp.append(snp_t)
                    allele_t = [i.split('(')[1].split(')')[0] for i in
                                line_list[2:]]
                    allele.append(allele_t)
            return haplo_name, haplo_perc, snp, allele
        except FileNotFoundError:
            print('The file was not found!\n')
            continue


def hand_allele_input_check(snp):
    """
    This function checks if the filled in allele are A, C, G or T.
    :param snp: RS code (string).
    :return: The allele string.
    """
    while True:
        allele = input(f'Which allele has {snp}: ').upper()
        count = 0
        for x in allele:
            if x in ['A', 'C', 'G', 'T']:
                count += 1
                if count == len(allele):
                    return allele
                else:
                    continue
            else:
                continue


def input_haplo_hand(pop):
    """
    This function is used if the user decides to fill in the
    information of the haplotypes, percentages and SNPs by hand.
    :param pop: String of the identity/name of the population.
    :return: A list of haplotype names, a list of
    percentages, a 2D list of SNP numbers per haplotype and a 2D list
    of the allele per snp.
    """
    snp_list = []
    haplotype_perc = []
    allele_list = []
    haplotypes = list(
        input(f'Which haplotypes would you like to add in {pop}'
              f'\nFormat example: haplotype_name_1, '
              f'haplotype_name_2, etc..\n').split(', '))
    for haplo in haplotypes:
        while True:
            try:
                perc_hap = float(input(f'What are the percentages of'
                                       f' haplotype {haplo} '
                                       f'(0-100 THIS MUST BE AN INTEGER): '))
                if 0 < perc_hap < 101:
                    haplotype_perc.append(perc_hap)
                    break
                else:
                    print('This number is not between 0 and 100!')
                    continue
            except [TypeError, ValueError]:
                print('This is not a number')
        snp = list(input(f'Which SNP/rs numbers are contained in '
                         f'haplotype {haplo}?\nFormat example: rs_number1,'
                         f' rs_number2, etc..\n').split(', '))
        snp_list.append(snp)
        list_allele = []
        for s in snp:
            allele = hand_allele_input_check(s)
            list_allele.append(allele)
        allele_list.append(list_allele)
    return haplotypes, haplotype_perc, snp_list, allele_list


def input_haplotype(pop):
    """
    This function is used to give the user the option to either fill
    in by hand or by the use of a file. This is in a while for if the
    user gives a wrong option.
    :param pop: String of the identity/name of the population.
    :return: A list of the haplotype names, a list of
    percentages, a 2D list of SNP numbers per haplotype and a 2D list
    of the allele per snp.
    """
    while True:
        input_file_or_hand = input('Do you want to fill in the '
                                   'haplotypes by hand (h) or '
                                   'add a file (f): ').lower()
        if input_file_or_hand == 'f':
            haplo_name, haplo_perc, snp, allele = haplotype_file(pop)
            return haplo_name, haplo_perc, snp, allele
        elif input_file_or_hand == 'h':
            h_haplo_name, h_haplo_perc, h_snp, h_allele = input_haplo_hand(
                pop)
            return h_haplo_name, h_haplo_perc, h_snp, h_allele
        else:
            print('This is not an option!')
            continue
