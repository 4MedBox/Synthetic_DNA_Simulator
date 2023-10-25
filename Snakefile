
configfile: "config.yaml"

rule all:
	input: config["gene"]+"_result.fasta"


rule get_gene:
	output: config["gene"]+".fasta"
	shell:
		"""	
		esearch -db nucleotide -query {config[NC_code]} | efetch -format fasta | egrep -v '>' | tr -d '\n' > {config[gene]}."fasta"

		"""

rule create_db:
	output: config["db"]+".db"
	log: "gnomad.log"
	script:
		"gnomad.py"

rule create_sequences:
	input: config["db"]+".db", config['gene']+".fasta"
	output: config["gene"]+"_result.fasta"
	log: "sequencer.log"
	script:
		"sequence_generator.py"
