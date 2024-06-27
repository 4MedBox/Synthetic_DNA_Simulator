"""
This script checks if the generated sequence has the correct number
of CAG repeats. If the length of the CAG repeats is not within
the interval desired by the user, the script can add more
CAG repeats or remove CAG repeats.
This script can also remove the CAA codon if the user so desires.
Packages/libraries used:
Random: Used to add a random amount of CAG.
"""

import random


class Survivability:
    def __init__(self, seq):
        self.seq = seq


def remove_CAA(cag_sub, caa, g_or_s, seq):
    """
    This function can remove the CAA codon from the sequence if the
    user wishes.
    :param cag_sub: String y or n if the user wishes to add CAG
    repeats to the subsequence.
    :param caa: String y or n if the user wants to add the CAA codon
    to the sequence.
    :param g_or_s: String g or s if the user wants the entire gene
    or subsequence.
    :param seq: The DNA sequence.
    :return: The modified sequence with or without the CAA
    codon and an additional sequence if the user wants a CAA codon.
    """
    if cag_sub == 'y' or g_or_s == 'g':
        if seq[252:255] == "CAA" and caa == 'n':
            extra_seq = ''
            seq = seq[:252] + seq[258:]
            return seq, extra_seq
        else:
            extra_seq = 'CAACAG'
            return seq, extra_seq
    else:
        extra_seq = ''
        return seq, extra_seq


def add_CAG(g_or_s, cag_sub, cag_repeat, min_cag, max_cag, sequence, caa,
            t_seq):
    """
    This function removes or adds CAG repeats based on user input.
    The CAG repeat in the N-terminus starts at
    and is located in exon 1 (file:///C:/Users/cjens/Downloads/
    A_Compact_Model_of_Huntingtin_Toxicity.pdf).
    Loss of CAA means early onset of HD
    (https://www.thelancet.com/journals/laneur/article/PIIS1474-4422
    (20)30343-4/fulltext).
    :param g_or_s: The entire gene or subsequence.
    :param cag_sub: If the user wants CAG repeats with subsequence.
    :param cag_repeat: Current number of CAG repeats.
    :param min_cag: Minimum number of CAG repeats.
    :param max_cag: Maximum number of CAG repeats.
    :param sequence: The original sequence.
    :param caa: String whether the CAA sequence is added or not.
    :param t_seq: If an extra sequence string is added
    if the RS numbers used only make mutations in the sequence
    before the CAG repeat.
    :return: The sequence and a string with the number of new CAG
    repeats.
    """
    seq, extra_seq = remove_CAA(cag_sub, caa, g_or_s, sequence)
    if g_or_s == 'g' or cag_sub == 'y':
        if cag_repeat < min_cag:
            min_rant = int(max_cag - min_cag)
            max_rant = int(max_cag) - cag_repeat
            # random.seed(5)
            add_cag = random.randint(min_rant, max_rant)
            seq_cag = add_cag*'CAG'
            try:
                seq = seq[:195] + t_seq + seq_cag + extra_seq + seq[195:]
                add_cag += 1
            except TypeError:
                seq = sequence[:195] + t_seq + seq_cag + extra_seq + \
                      sequence[195:]
            cag_string = f"CAG repeat length: {cag_repeat + add_cag}" \
                         f" ({add_cag} CAG repeats were added!)"
        else:
            min_rant = cag_repeat - int(max_cag)
            max_rant = cag_repeat - int(min_cag)
            # random.seed(5)
            random_cag = random.randint(min_rant, max_rant)
            remove_seq = 195 + (random_cag * 3)
            seq = seq[:195] + seq[remove_seq:]
            cag_string = f"CAG repeat length: {cag_repeat - random_cag}" \
                         f" ({random_cag} CAG repeats were removed!)"
        return seq, cag_string
    else:
        cag_string = f"CAG repeat length: 0"
        return seq, cag_string


def sequence(seq, min_cag, max_cag, g_or_s, cag_sub, caa, t_seq):
    """
    This function checks the CAG and CAA repeat length and verifies
    that it is within the user interval.
    :param seq: String of the mutated sequence.
    :param g_or_s: The entire gene or subsequence.
    :param cag_sub: If the user wants CAG repeats with subsequence.
    :param caa: String whether the CAA sequence is added or not.
    :param t_seq: If an extra sequence string is added
    if the RS numbers used only make mutations in the sequence
    before the CAG repeat.
    :param min_cag: The minimum number of CAG repeats.
    :param max_cag: The maximum number of CAG repeats.
    :return: A string with the number of repeats and whether it
    was mutated.
    """
    seq_pos = 0
    cag_list = []
    seq_count = 0

    while seq_pos - 1 != len(seq):
        if seq[seq_pos:seq_pos + 3] == "CAG" or (
                seq[seq_pos - 3:seq_pos] == "CAG" and
                seq[seq_pos + 3:seq_pos + 6] == "CAG"):
            if seq[seq_pos - 3:seq_pos] != "CAG":
                cag_list.append([seq_pos, 0, 0])

            seq_pos += 3

            cag_list[::-1][0] = [cag_list[::-1][0], 0,
                                 cag_list[::-1][0][2] + 3]

            if seq[seq_pos + 3:seq_pos + 6] != "CAG":
                cag_list[::-1][0][1] = seq_pos

        if seq[seq_pos:seq_pos + 3] != "CAG":
            seq_pos += 1

        if seq[seq_pos - 3:seq_pos] == "CAG" and \
                seq[seq_pos:seq_pos + 3] != "CAG":
            seq_count += 1

    list_amount_cag = []
    for i in cag_list:
        num = int((i[1] - i[0]) / 3)
        if min_cag <= (i[1] - i[0]) / 3 <= max_cag:
            cag_length = f'CAG repeat length: {num}' \
                         f' (correct mutated HTT gene)'
            return seq, cag_length
        else:
            if num > 1:
                list_amount_cag.append(num)

    try:
        cag_num = list_amount_cag[0]
    except IndexError:
        cag_num = 0
    seq, cag_length = add_CAG(g_or_s, cag_sub, cag_num, min_cag, max_cag,
                              seq, caa, t_seq)
    return seq, cag_length
