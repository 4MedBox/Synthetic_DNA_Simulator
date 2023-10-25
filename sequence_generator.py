import itertools
import sqlite3 as sl
import sys
import Thread
import datetime
from tqdm import tqdm
import survivability


# def mutate(new_sequence, mutations):
#     #create synthetic sequence.
#     f.write('create synthetic sequence')
#     for mutation in mutations:
#         pos = mutation[4] - start
#         ref = mutation[2]
#         alt = mutation[3]
#
#         if pos < len(new_sequence):
#             if new_sequence[pos-1:pos-1 + len(ref)] == ref:
#
#                 if len(ref) == len(alt):
#                     f.write('Mutation: SNP')
#                     new_sequence = new_sequence[:pos-1] + alt + new_sequence[pos:]
#
#
#                 # THIS WORKS, DO NOT TOUCH :)
#                 if len(ref) < len(alt):
#                     f.write("Mutation: Insertion")
#                     new_sequence = new_sequence[:pos-1] + alt + new_sequence[pos:]
#
#
#                 # THIS ALSO WORKS, DO NOT TOUCH :)
#                 if len(ref) > len(alt):
#                     f.write('Mutation: Deletion')
#                     new_sequence = new_sequence[:pos-1] + alt + new_sequence[pos-1 + len(ref):]
#
#     f.write('Sequence created.\n')
#
#
#     return new_sequence


if __name__ == '__main__':
    with open(snakemake.log[0], "w") as f:
        f.write('connecting to database and creating cursor\n')

        #the following code connects to a database and creates a cursor.
        con = sl.connect(snakemake.config['db']+'.db') #gnomad_htt.db for complete HTT db
        cursor = con.cursor()


        f.write('generate starting population\n')
        if snakemake.config['number_of_mutations'] == "all":
            test = cursor.execute("SELECT * FROM mutation;")
            number_of_mutations = 1
            for i in test:
                number_of_mutations+=1
        else:
            number_of_mutations= snakemake.config['number_of_mutations'] + 1

        #the following code generates the starting population. In this case: there will be a population with combinations of 20 different mutations.
        lis = list(range(1, number_of_mutations))

        f.write('get start and end position\n')
        #the following selects the highest and lowest position from N mutations, making sure only a part of the gene will be used of generating a new sequence, instead of the whole gene.
        cursor.execute("SELECT MIN(pos), MAX(pos) FROM mutation ")
        start, end = [i for i in cursor.fetchall()[0]]


        f.write('get part of sequence based on start and end position\n')
        sequence = open(snakemake.config['gene']+".fasta", "r").readlines()[0][start:end+1]

        #a list will be created to keep all the mutations in. This also counts for the file.
        f.write('open result-file and population list for mutations\n')
        population=0
        file = open(snakemake.config["gene"]+"_result.fasta", "a")

        lijst_threads=list(range(1,snakemake.config["cores"]+1)) #4 threads

        count=1
        f.write('loop through starting population, starting with 1 mutation and adding 1 mutation with each loop\n\n')

    #for loop: every amount of mutations will be used (1 mutation, 2 mutations, 3 mutations etc.)
        while count <= lis[len(lis)-1]:

            #while loop below is to get sequence with 5232 mutations.
            # while count <= 5232:
            #     count+=1
            #     print(count)
            # print("ja")
            f.write('====New Number of Possible Mutations===\n')
            f.write('Get all possible combinations with the to be used mutations\n')
            #with the N to be used mutations, all of the possible combinations will be calculated.
            combi_list = list(itertools.combinations(lis, count))
            test = combi_list

            test_len = len(combi_list)

            #loop through all the possible combinations of mutations.
            f.write('loop through all possible combinations of mutations.\n\n')


            with tqdm(total=len(combi_list), desc=str(count)+" mutaties") as pbar:
                while test_len > 0:
                    f.write('===New Possible Mutated Sequence===\n')

                    #get info from all the mutations in a certain combination.
                    f.write('get info from all mutations in combination\n')

                    fetch_list = [cursor.execute(
                        "SELECT database_id, variant_id, ref, alt, pos, rs_id FROM mutation WHERE database_id IN {}".format(i+(0,))).fetchall() for i in test[:len(lijst_threads)] ]

                    f.write('create synthetic person with sequence\n')
                    threads = [Thread.myThread(number+1, "Thread-" + str(number+1), sequence, fetch_list[number]) for number in
                               range(len(test[:len(lijst_threads)]))]

                    result = [[thread.mutations, thread.new_sequence, thread.name] for thread in threads]

                    for gene in result:

                        #if statement below is to go to the survival check if the gene is HTT
                        # if snakemake.config['gene'].upper() == "HTT":
                        #     survivability.Survival(gene[1])

                        f.write("found 2 repeats\n")
                        population+=1

                        f.write('Add synthetic person to population.\n')

                        #create (meta-)data for result file, and write to result file.
                        f.write('Write results to file\n\n')
                        data = [">seq_",str(population),", (",', '.join(f'{x[::-1][0]}'for x in gene[0]),"), ", str(datetime.datetime.now()),"\n",gene[1],"\n\n"]
                        file.writelines(data)
                        f.write("done\n")

                        if type(snakemake.config["max_sequences"]) == int:
                            if population == snakemake.config["max_sequences"]:
                                pbar.update(len(test))
                                sys.exit()
                    del test[:len(lijst_threads)]
                    test_len -= len(threads)
                    pbar.update(len(lijst_threads))

                pbar.close()
            count += 1

        #print end time (you'll likely walk away while running the programm, it takes some time :) )
        f.write('Script is done.\n')
        f.write(str(datetime.datetime.now())+"\n")
    f.close()

