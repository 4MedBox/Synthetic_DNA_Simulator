# SynthDNASim: Creating diverse synthetic mutated DNA data set
A python tool
for generating synthetic DNA sequences based on the dbSNP database from NCBI.

## Why should I use this project?
There is a lack of genomic data in the biomedical research of rare genetic diseases such as HD.
This data is needed, for example,
to determine the genetic process
that causes these diseases or to identify gene variants.
The lack of data is due to the rarity of the disease and the privacy of genomic data.
Regulatory laws protect these individuals' privacy,
which raises the issue of access to genomic data.
One solution is
to create diverse synthetic DNA data
so that researchers can generate their own DNA data,
which makes the research process much faster and more reproducible
because the DNA data does not belong to anyone.
This is why SynthDNASim was created.

## Setup
SynthDNASim is composed of six different scripts.
These are user_input.py, mutate_ref_sequence.py, check_CAA_CAG_HD.py,
sequence_generator.py, population.py, and main.py.
The main.py script executes the tool.
The tool is based on the GrCh38 reference genome assembly.

#### The main script
The main script is used to execute the tool.
It can be used with the following command,
making sure that all necessary scripts and files are in a working directory:
```commandline
python main.py --config_file name_of_config_file.json
```
The script can also be run in an integrated development environment such as PyCharm.
Make sure that all necessary scripts and files are in one working directory.
This script executes all the code from the files necessary to create a 
synthetic DNA data set.
Main.py can create a config file with user input,
or the user can add an already existing config file.
The config file must be a JSON file!
Packages/libraries used are:
* [json](https://docs.python.org/3/library/json.html): Used to read and write the json config file
* [argparse](https://docs.python.org/3/library/argparse.html): Used to give arguments/parameters to a python script when it is run in a command line
* [time](https://docs.python.org/3/library/time.html): Used to see how long it takes for the script to create all populations
* [csv](https://docs.python.org/3/library/csv.html): Used to write metadata in a tsv file

###### Output (main.py)
This tool creates two files: 
a Fasta file with synthetic DNA sequences 
and a metadata file with sample names and population data. 
An example is HTT_result.fasta.

#### User input (user_input.py)
This script is used to ask and store the user input.
Packages/libraries used are:
* [re](https://docs.python.org/3/library/re.html): Used for the extraction of the snp or rs number from

###### Optional input file (haplotype file)
The text file must contain haplotype names and percentages, 
along with rs codes and alleles.
This is per population!: 
- Add one haplotype on a newline!:

  haplotype name1, percentage (0-100 THIS MUST BE AN INTEGER),
  RS code 1 (allele (A, C, G and or T)),
  RS code 2 (allele (A, C, G and or T)), etc…

  haplotype name2, percentage (0-100 THIS MUST BE AN INTEGER), RS code 1 (allele (A, C, G and or T)), RS code 2 (allele (A, C, G and or T)), etc…
Example: haplotype p21, 20, rs986354601 (T), rs946330590 (T)

###### Output
The user input generates a configuration file called config.json. 
The content of this file is described below. 

###### config file 
`config.json` contains all info:

| Item                  | Value(s)               | Description                                                                               |
|-----------------------|------------------------|-------------------------------------------------------------------------------------------|
| `CAG repeats`         | integers               | The min and maximum number of CAG repeat.                                                 | 
| `CAG sub`             | String y or n          | If the user wants the CAG repeat in the sub sequence.                                     |
| `CAA`                 | String y or n          | If the user wants the CAA codon in the CAG repeat.                                        |
| `recombine`           | String y or n          | If the user wants to recombine mutations.                                                 |
| `recombine_pop`       | String b or w          | If the user wants to recombine mutations from within a population or between populations. |
| `amount_recombine`    | Integer                | The number of sequences the user wants with new recombinations.                           |
| `gene_or_subseq`      | String g or s          | If the user wants an entire gene or just a subsequence.                                   |
| `ensemble`            | String ENSG00000197386 | ENSEMBL code for the HTT-gene on gnomAD.                                                  |
| `gene`                | String HTT             | The name of the gene.                                                                     |
| `NC_code`             | String NC000004.12     | NCBI code for the 4th chromosome of the H. Sapiens (GrChr38)                              |
| `db`                  | HTT                    | Name of the database.                                                                     |
| `max_sequences`       | Integer                | Number of generated sequences needed to be displayed in the result file.                  |
| `treats`              | Integer                | Number of threats used.                                                                   |
| `population`          | String/Integer         | Name of the population, percentage of the population, haplotypes, RS codes and allele.    |


#### Creating a database from variants from gnomad (gnomad.py, this is currently not used!)
The python script asks for the ENSEMBL code according to the GrChr37 version of the 
human chromosome, and the number of mutations. It takes the ENSEMBL code and searches
for the variant list in the ENSEMBL databases. It retrieves this data and searches 
for the mutation data with every retrieved RS-code. With certain parts of this 
mutational information, a new database will be created.

###### Library dependencies
* [time](https://docs.python.org/3/library/time.html): used to wait _x_ seconds, therefore, the gnomAD database will not overload
* [sqlite3](https://docs.python.org/3/library/sqlite3.html): used to create and add to database
* [datetime](https://docs.python.org/3/library/datetime.html): used to calculate how much time from the retrieval from gnomAD data and push to a database is
* [requests](https://pypi.org/project/requests/): used to send requests to the GNOMAD api
* [sys](https://docs.python.org/3/library/sys.html): used to make sure all print statements towards terminal are noted in log file
* [tqdm](https://tqdm.github.io/): progressbar, so the user can see how much time passed and at what percentage the (part of the) script is
* [json](https://docs.python.org/3/library/json.html): Used to read and write json config files

###### Output
The command line above will give the database file ending with `.db`. The name of the
given gene will be in front of the filetype.


#### Generate the synthetic sequences (sequence_generator.py)
This script generates mutated sequences for each haplotype and checks
if they have the correct CAG repeat interval.

###### Library dependencies
* [biopython](https://biopython.org/): To extract the chromosome DNA sequence and SNP info from NCBI
* [os](https://docs.python.org/3/library/os.html): Used for working in different directories and to create a directory for the XML RS code files.
* [xml.etree.ElementTree](https://docs.python.org/3/library/xml.etree.elementtree.html#module-xml.etree.ElementTree): Used to read the XML RS code files
* [math](https://docs.python.org/3/library/math.html): Used to round up numbers
* [random](https://docs.python.org/3/library/random.html):
  Used to create random recombinations of mutations
in and between populations and to create a seed
* [csv](https://docs.python.org/3/library/csv.html): Adding rows to the metadata tsv file

###### Output
The output will be a Fasta file with the gene name and `_result.fasta` behind it. 


#### Creating a synthetic population (population.py)
This class creates the synthetic DNA sequences of each haplotype per population.


#### Mutating the sequence
`mutate_ref_sequence.py` mutates the reference DNA sequence of the chromosome.
It uses the library [threading](https://docs.python.org/3/library/threading.html). 
With this library,
it is possible to run multiple threads alongside each other. 


#### Check the CAG repeats for HD 
`check_CAA_CAG_HD.py` This script checks if a synthetically created HTT gene has the correct CAG repeat interval and can also remove the CAA codon. If the CAG repeats are not in the user interval, it can add or remove them.
###### Packages/libraries used are:
* [random](https://docs.python.org/3/library/random.html): Random is used to randomly add an amount of CAG and to add a seed for reproducibility.

