"""
Author: Liz Tucker
Done for CSC314 with Professor Dancik
This script produces the optimal global alignment between
two sequences using the Dynamic Programming method.
"""
from Bio import SeqIO
from pprint import pprint

def get_sequence(filename):
    with open(filename) as handle :
        sequence_iter = SeqIO.parse(handle, "genbank")
        seq_record = next(sequence_iter)
        return seq_record.seq

"""Makes the dynamic programming chart"""
def chart_maker(seq_1, seq_2, len_seq_1,len_seq_2):
    index=0
    #create chart
    chart=[[0 for row in range(len_seq_1 + 2)] 
            for col in range(len_seq_2 + 2)] 
    ######POPULATE SEQUENCES###########
    #populate first horizontal row with sequence 1
    chart[0][2:]=[letter for letter in seq_1]
    #populate first vertical col with second sequence and second vertical col with calculations
    chart[2:]=[[seq_2[index], (-4*(index+1))] + lst[2:] for index,lst in enumerate(chart[2:]) for item in lst[:-len(lst)+1]]
    #populate first horizontal line with calculations
    chart[1][1:]=[(-4*index) for index, entry in enumerate(chart[1][1:])]
    return chart

"""Calculates the dynamic programming chart"""
def calc_global_alignment(chart):
    traceback = []
    index = 1 #vertical index, starting at one to be able to pass the correct value into eval_highest
    v_index=2 #vertical index used to store answers in the chart
    for lst in chart[2:]:
        h_index=2 #horizontal index
        for item in lst[2:]:
            chart[v_index][h_index] = eval_highest(chart[index][h_index-1], chart[v_index][h_index-1], chart[index][h_index], determine_match(chart[v_index][0],chart[0][h_index]))[0]
            traceback.append(eval_highest(chart[index][h_index-1], chart[v_index][h_index-1], chart[index][h_index], determine_match(chart[v_index][0],chart[0][h_index])))
            h_index+=1
        v_index+=1
        index+=1
    return chart, traceback

"""Returns true if two character match, false otherwise"""
def determine_match(char_a, char_b):
    if char_a == char_b: return True
    else: return False

"""Evaluates the highest score for a given spot"""
def eval_highest(p_char, h_gap, v_gap, b_match):
    gap_penalty = -4
    match = 5
    mismatch = -1
    calc = []
    if b_match == True:
        temp = p_char + match
        calc.append(temp)
    else:
        temp = p_char + mismatch
        calc.append(temp)
    cal_1 = h_gap + gap_penalty
    cal_2 = v_gap + gap_penalty
    calc.append(cal_1)
    calc.append(cal_2)
    if max(calc) == calc[0] and b_match == True:
        return calc[0], p_char, 'match'
    elif max(calc) == calc[0] and b_match == False:
        return calc[0], p_char, 'mismatch'
    elif max(calc) == calc[1]:
        return calc[1], h_gap, 'h-gap'
    else: return calc[2], v_gap, 'v-gap'

"""Gets traceback sequence"""
def get_traceback(full_traceback, length, chart):
    traceback = []
    traceback.append(full_traceback[-1])
    while True:
        if len(traceback) > (length-2):
            break
        for item in full_traceback:
            if len(traceback) > (length-2):
                break
            if item[0] == traceback[-1][1]:
                traceback.append(item)
            temp = [item for item in full_traceback if traceback[-1][1] in item]
            for item in temp:
                if item[0]==traceback[-1][1] and item[1]==0:
                    traceback.append(item)
    for item in chart[1][1:]:
        index=1
        if traceback[-1][1] == item and traceback[-1][0] == chart[2][index]:
            traceback.append((chart[1][index], item, 'h-gap'))
    for index, item in enumerate(chart[2:]):
        if traceback[-1][1] == item[1] and traceback[-1][0] == chart[index+3][2]:
            traceback.append((chart[1][index], item[1], 'v-gap'))
    traceback = traceback[::-1]
    return traceback

"""Returns length of longer sequence"""
def get_longer_sequence(seq_1,seq_2):
    lengths = []
    lengths.append(len(seq_1))
    lengths.append(len(seq_2))
    return max(lengths)

"""Returns sequences as lists for printing"""
def determine_alignment(seq_1, seq_2, traceback):
    lines=[]
    seq_1_a = []
    seq_2_a = []
    index = 0
    for item in traceback:
        if item[2]=='match':
            lines.append('|')
            if index >= len(seq_1):
                seq_1_a.append('-')
            else:
                seq_1_a.append(seq_1[index])
            seq_2_a.append(seq_2[index])
        elif item[2]=='v-gap':
            lines.append('')
            if index >= len(seq_1):
                seq_1_a.append('-')
            else:
                seq_1_a.append(seq_1[index])
            seq_2_a.append(seq_2[index])
        elif item[2]=='h-gap':
            lines.append('')
            seq_2_a.append('-')
            seq_1_a.append(seq_1[index])
            seq_2_a.append(seq_2[index])
        else: 
            lines.append('')
            if index >= len(seq_1):
                seq_1_a.append('-')            
            else:
                seq_1_a.append(seq_1[index])
            seq_2_a.append(seq_2[index])
        index+=1
    return seq_1_a, lines, seq_2_a

def format_alignment(alignment):
    seq_1 = str()
    lines= str()
    seq_2 = str()
    for char in alignment[0]:
        seq_1+=char
    lines= ''.join(alignment[1])
    for char in alignment[2]:
        seq_2+=char
    return seq_1, lines, seq_2

def main():
    filename_1 = input('Please enter a sequence filename: ')
    filename_2 = input('Please enter a second sequence filename: ')
    seq_1 = get_sequence(filename_1)
    seq_1=seq_1[0:10]
    seq_2 = get_sequence(filename_2)
    seq_2=seq_2[0:10]
    length = get_longer_sequence(seq_1, seq_2)
    chart = chart_maker(seq_1, seq_2, len(seq_1), len(seq_2))
    chart = calc_global_alignment(chart)[0]
    traceback = calc_global_alignment(chart)[1]
    traceback = get_traceback(traceback, length, chart)
    alignment = determine_alignment(seq_1, seq_2, traceback)
    alignment = format_alignment(alignment)
    print()
    print('------------------Dynamic Programming Matrix------------------')
    pprint(chart)
    print('------Traceback------')
    pprint(traceback)
    print('------Alignment------')
    for item in alignment:
        print(item)

if __name__ == "__main__":
    main()