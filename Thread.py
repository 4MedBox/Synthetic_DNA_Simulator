import threading

exitFlag = 0


class myThread (threading.Thread):
   def __init__(self, threadID, name, sequence, mutations):
      threading.Thread.__init__(self)
      self.threadID = threadID
      self.name = name
      self.sequence = sequence
      self.mutations = mutations
      self.new_sequence = mutate(sequence, name, mutations)

   def run(self):
      mutate(self.sequence, self.name, self.mutations)




def mutate(new_sequence, threadName, mutations):
   mutations_del = mutations.copy()
   if exitFlag:
      threadName.exit()

   # create synthetic sequence.
   while len(mutations_del) != 0:
      mutation = mutations_del[0]
      pos = mutation[4] - 3076483 #this number needs to be gone if the gene is used, number is for HTT.
      ref = mutation[2]
      alt = mutation[3]

      if pos < len(new_sequence):
         if new_sequence[pos - 1:pos - 1 + len(ref)] == ref:

            if len(ref) == len(alt):
               # logging.info('Mutation: SNP')
               new_sequence = new_sequence[:pos - 1] + alt + new_sequence[pos:]

            # THIS WORKS, DO NOT TOUCH :)
            if len(ref) < len(alt):
               # logging.info("Mutation: Insertion")
               new_sequence = new_sequence[:pos - 1] + alt + new_sequence[pos:]

            # THIS ALSO WORKS, DO NOT TOUCH :)
            if len(ref) > len(alt):
               # logging.info('Mutation: Deletion')
               new_sequence = new_sequence[:pos - 1] + alt + new_sequence[pos - 1 + len(ref):]
      del mutations_del[0]

   return new_sequence


#
