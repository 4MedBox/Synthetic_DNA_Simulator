"""
This script mutates the reference sequences this process is sped up
by the use of multi-threading.
Packages/libraries used are:
threading: To run function multiple times at the same time.
"""
import threading


exitFlag = 0


class myThread (threading.Thread):
    """
    This class is designed to provide an output for mutating multiple
    sequences at once.
    """
    def __init__(self, threadID, name, sequence, mutations, min_pos):
        """
        The init of this class initializes all variables of this
        class.
        It also initializes the threading class.
        :param threadID: Thread number.
        :param name: Thread name.
        :param sequence: The DNA sequence of chromosome 4.
        This sequence is changed depending on the mutations
        the user chooses to add.
        :param mutations: List of all the mutations.
        :param min_pos: The start position of the sequence,
        which will be used to correct the position of the
        mutation.
        """
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.sequence = sequence
        self.mutations = mutations
        self.min_pos = min_pos
        self.new_sequence = mutate(sequence, name, mutations, min_pos)

    def run(self):
        """
        This calls the mutate function, which mutates the sequences.
        """
        mutate(self.sequence, self.name, self.mutations, self.min_pos)

# class NewThread(Thread):
# 	def __init__(self, group=None, target=None, name=None, args=(),
# 				 kwargs=None):
# 		Thread.__init__(self, group, target, name, args, kwargs)
#
# 	def run(self):
# 		self._return = self._target(*self._args, **self._kwargs)
#
# 	def join(self, *args):
# 		Thread.join(self, *args)
# 		return self._return


def mutate(new_sequence, threadName, mutations, min_pos):
    """
    This function mutates the DNA with the mutations desired
    by the user.
    :param new_sequence: The DNA sequence of chromosome 4.
    :param threadName: The name of the thread.
    :param mutations: List of all the mutations.
    :param min_pos: The start position of the sequence
    which will be used to correct the position of the
    mutation.
    :return: The mutated sequence.
    """
    if exitFlag:
        threadName.exit()

    for mutation in mutations:
        ref = mutation[2]
        alt = mutation[3]
        begin_pos = mutation[1] - min_pos

        if len(ref) == len(alt):
            new_sequence = new_sequence[:begin_pos] + alt + new_sequence[
                                                           begin_pos + 1:]
            continue

        if len(ref) < len(alt):
            new_sequence = new_sequence[:begin_pos] + alt + \
                           new_sequence[begin_pos + len(ref):]
            continue

        if len(ref) > len(alt):
            new_sequence = new_sequence[:begin_pos] + alt + \
                           new_sequence[begin_pos + len(ref):]
            continue

    return new_sequence


# def mutate(new_sequence, threadName, mutations):
#    #mutations_del = mutations.copy()
#    if exitFlag:
#       threadName.exit()
#
#    # create synthetic sequence.
#    for mutation in mutations:#while len(mutations_del) != 0:
#       #mutation = mutations_del[0]
#       print(mutation)
#       pos_eind = mutation[4] - 3074681 #this number needs to be gone if the gene is used, number is for HTT. or 0? 3076483 4
#       ref = mutation[2]
#       alt = mutation[3]
#
#       if pos_eind < len(new_sequence):
#          #if new_sequence[pos_eind - 1:pos_eind - 1 + len(ref)] == ref:
#
#             if len(ref) == len(alt):
#                print('snp')
#                # logging.info('Mutation: SNP')
#                #new_sequence = new_sequence[:pos - 1] + alt + new_sequence[pos:]
#                new_sequence = new_sequence[:pos_eind] + alt + new_sequence[pos_eind+1:]
#
#             # THIS WORKS, DO NOT TOUCH :)
#             if len(ref) < len(alt):
#                print('ins')
#                # logging.info("Mutation: Insertion")
#                #new_sequence = new_sequence[:pos - 1] + alt + new_sequence[pos: ]
#                new_sequence = new_sequence[:pos_eind-1] + alt + new_sequence[pos_eind: ]
#
#             # THIS ALSO WORKS, DO NOT TOUCH :)
#             if len(ref) > len(alt):
#                print('del')
#                # logging.info('Mutation: Deletion')
#                #new_sequence = new_sequence[:pos - 1] + alt + new_sequence[pos - 1 + len(ref):] #:pos - 1 en pos -1
#                new_sequence = new_sequence[:pos_eind - 1] + alt + new_sequence[pos_eind - 1 + len(ref):]
#       #del mutations_del[0]
#
#    return new_sequence


#
