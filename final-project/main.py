"""
Author: Liz Tucker
Done for CSC314 with Professor Dancik
This script produces the optimal global alignment between
two sequences using the Dynamic Programming method.
"""
from pprint import pprint


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
    chart[2:]=[[seq_2[index], (-4*index)] + lst[2:] for index,lst in enumerate(chart[2:]) for item in lst[:-len(lst)+1]]
    #populate first horizontal line with calculations
    chart[1][1:]=[(-4*index) for index, entry in enumerate(chart[1][1:])]
    return chart

"""Calculates the dynamic programming chart"""
def calc_global_alignment(seq_1, seq_2, chart):
    index=0
    gap_penalty = -4
    match = 5
    mismatch = -1
    c_index=0
    for lst in chart[2:]:
        for item in lst[2:]:
            if chart[0][c_index] == chart[c_index][0]:
                print(chart[c_index][c_index])
        c_index+=1
    return chart

def main():
    seq_1 = 'taca'
    seq_2 = 'ttcag'
    chart = chart_maker(seq_1, seq_2, len(seq_1), len(seq_2))
    chart = calc_global_alignment(seq_1, seq_2, chart)
    pprint(chart)

if __name__ == "__main__":
    main()