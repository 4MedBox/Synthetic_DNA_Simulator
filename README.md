# Creating synthetic mutated DNA data
A bash/python script for generating synthetic mutated DNA sequences based on the gnomAD database.

## Why should I use this project?
It is difficult to get personal patient data such as DNA due to privacy concerns or low availability
of data. The gnomAD database uses mutations found in population studies. These 
mutations are publicly available due to the consent each person has given to go along 
with the study. Therefore, it is possible to create a synthetic population based on the data of these mutations. With this tool and the to be generated sequences, researchers can perform tests on 
their tools, or can research certain combinations of mutations, which otherwise 
would not be possible. Please make sure to install all libraries and versions.

## Setup
This pipeline consists of multiple scripts, which can be run independently. All 
scripts are combined into one script, called `Snakefile`. Each script uses files 
and/or databases with specific names or NCBI-/ENSEMBL-codes according to the GrChr37 
version of the genes and their according chromosomes.

### The main script
#### What does it do?
The main script calls the other scripts and retrieves data such as chromosomes and 
mutational data, stores it, and uses it to create synthetic sequences. Itt can be used with the following command:
```commandline
snakemake
```
When the command above is used, the default options from `config.yaml` are used.

If you want to run the complete script with other options than the default options, 
you can use `--config` with the options from the other scripts as described below.

#### Commandline input
```commandline
snakemake get_gene
```

#### Output
The command line above will give the fasta file ending with `.fasta`. The name of 
the given gene will be in front of the filetype.

#### User input
There is also a possibility to give other input than the default input:
```commandline
snakemake get_gene --config NC_code="NC_000023.10" gene="DMD"
```

#### Output
In the example given above, the output will be a fasta file named `DMD.fasta`




### get_gene
#### What does it do?
The bash script asks for the NC code from NCBI according to the GrChr37 version of 
the human chromosome, and the name of the gene. It searches on the NCBI databases for 
the given code, and gives a fastafile back. The bash script is located in the 
file `Snakefile`.

#### Library dependencies
* [esearch](https://www.ncbi.nlm.nih.gov/books/NBK25499/): used to search in the nucleotide NCBI database.
* [efetch](https://www.ncbi.nlm.nih.gov/books/NBK25499/): used to fetch the fasta file from NCBI.
* [egrep](https://www.gnu.org/software/grep/manual/grep.html): deletes the fasta header from the fasta file.
* [tr](https://www.gnu.org/software/coreutils/manual/html_node/tr-invocation.html#tr-invocation): deletes all enters, so the fasta file will be one long string.

#### Commandline input
```commandline
snakemake get_gene
```

#### Output
The command line above will give the fasta file ending with `.fasta`. The name of 
the given gene will be in front of the filetype.

#### User input
There is also a possibility to give other input than the default input:
```commandline
snakemake get_gene --config NC_code="NC_000023.10" gene="DMD"
```

#### Output
In the example given above, the output will be a fasta file named `DMD.fasta`


### create_db
#### What does it do?
The python script asks for the ENSEMBL code according to the GrChr37 version of the 
human chromosome, and the number of mutations. It takes the ENSEMBL code and searches
for the variant list in the ENSEMBL databases. It retrieves this data and searches 
for the mutation data with every retrieved RS-code. With certain parts of this 
mutational information, a new database will be created.  


#### Library dependencies
* [time](https://docs.python.org/3/library/time.html): used to wait _x_ seconds, therefore the gnomAD database will not overload. 
* [sqlite3](https://docs.python.org/3/library/sqlite3.html): used to create and add to database
* [datetime](https://docs.python.org/3/library/datetime.html): used to calculate how much time from  the retrieval from gnomAD data and push to database is
* [requests](https://pypi.org/project/requests/) (version 2.28.2): used to send requests to the GNOMAD api.
* [sys](https://docs.python.org/3/library/sys.html): used to make sure all print statements towards terminal are noted in log file
* [tqdm](https://tqdm.github.io/) (version 4.65.0): progressbar, so the user can see how much time passed and at what percentage the (part of the) script is.

#### Commandline input
```commandline
snakemake create_db
```

#### Output
The command line above will give the database file ending with `.db`. The name of the
given gene will be in front of the filetype.

#### User input
There is also a possibility to give other input than the default input:
```commandline
snakemake create_db --config db="test" rs_codes="rs1445714790"
```

#### Output
In the example given above, the output will be a database file named `test.db`. It 
will have the default number of 20 mutations, including the rs-code `rs1445714790`. It is also possible to give more than one rs-code, as long as it is between quotation marks and separated with comma's.

### create_sequences
#### What does it do?
With this script, the user will create sequences based on the mutations in the given 
database.

#### Library dependencies
* [itertools](https://docs.python.org/3/library/itertools.html): creates combinations of mutations based on the number of mutations needed in that particular loop
* [sqlite3](https://docs.python.org/3/library/sqlite3.html): used to create cursor and connection, to search for mutation info in created database.
* [sys](https://docs.python.org/3/library/sys.html): used to exit the program, is certain threshold in HTT sequence is reached.
* Thread (See "_Other scripts: Thread.py_"): used to thread with 4 cores at the same time, therefore 4 combinations of mutations can be created at the same time.
* [datetime](https://docs.python.org/3/library/datetime.html): used to set a timestamp in the result file to see when a sequence is created. 
* [tqdm](https://tqdm.github.io/) (version 4.65.0):  progressbar, so the user can see how much time passed and at what percentage the (part of the) script is.



#### Commandline input
```commandline
 snakemake create_sequences
```

#### Output
The output will be a fasta file with the name of the used database as name. 

#### User input
There is also a possibility to give other input than the default input:
```commandline
 snakemake create_sequences --config db="HTT_30" number_of_mutations=30
```

#### Output
The output will be a fasta file with the name of the gene, with `_result.fasta` 
behind the gene. With the used options, the top 30 mutations in the given database 
will be used to create sequences.

### Other used files
There are two more files, which are being used to run the pipeline and create the synthetic sequences.

#### config.yaml
`config.yaml` contains all info for the default options. The default options are as followed:

| Syntax                | Default         | Description                                                              |
|-----------------------|-----------------|--------------------------------------------------------------------------|
| `ensembl`             | ENSG00000197386 | ENSEMBL code for the HTT-gene on gnomAD                                  | 
| `gene`                | HTT             | Abbreviation for the Huntington gene                                     |
| `NC_code`             | NC_000004:11    | NCBI code for the 4th chromosome of the H. Sapiens                       |
| `db`                  | HTT             | Abbreviation for the Huntington gene                                     |
| `number_of_mutations` | 20              | Amount of mutations in the database  which are to be used                |
| `rs_codes`            | all             | Which mutations are to be used to generate synthetic sequences           |
| `max_sequences`       | all             | Amount of generated sequences needed to be displayed in the result file. |

When the user uses another variable than the default variable, the default variable will be overriden and the user variable will be used.

 
#### thread.py
`thread.py` speeds up the script with the use of multiple threads. It uses the library [threading](https://docs.python.org/3/library/threading.html) (version ???). 
With this library, it is possible to run multiple threads alongside each other. In this case, the default 
number of to be used cores would be 4, the default setting in `config.yaml`. This can be changed in the command line with `--config threads=2` for example.


#### survivability.py
`survivability.py` checks if a synthetically created HTT gene would survive in the real world. It checks based on the CAG repeat. If there is more than 27 CAG repeats in a sequence, the sequence will be printed on the terminal.

