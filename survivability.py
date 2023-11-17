
class Survivability ():
   def __init__(self, seq):
      self.seq = seq


def Survival(seq):
    seq_pos = 0
    cag_list = []
    seq_count = 0

    while seq_pos - 1 != len(seq):
        if seq[seq_pos:seq_pos + 3] == "CAG" or (
                seq[seq_pos - 3:seq_pos] == "CAG" and seq[seq_pos + 3:seq_pos + 6] == "CAG"):
            if seq[seq_pos - 3:seq_pos] != "CAG":
                cag_list.append([seq_pos, 0, 0])
            seq_pos += 3

            cag_list[::-1][0] = [cag_list[::-1][0], 0, cag_list[::-1][0][2] + 3]

            if seq[seq_pos + 3:seq_pos + 6] != "CAG":
                cag_list[::-1][0][1] = seq_pos

        if seq[seq_pos:seq_pos + 3] != "CAG":
            seq_pos += 1
        if seq[seq_pos - 3:seq_pos] == "CAG" and seq[seq_pos:seq_pos + 3] != "CAG":
            seq_count += 1

    for i in cag_list:
        if (i[1] - i[0]) / 3 >= 27:
            return ["mutated HTT gene"]

    return ["non mutated HTT gene"]
