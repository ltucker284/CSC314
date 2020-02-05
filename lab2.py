accepted_character = ['T', 'A', 'G', 'C']
character_dict = {
    'T':'A',
    'A':'T',
    'G':'C',
    'C':'G'
}

def check_sequence(sequence: str):
    new_sequence = []
    for character in sequence:
        # print(character)
        if character in accepted_character:
            new_sequence.append(character)
        else:
            new_sequence.append('-')
    new_sequence = ''.join(new_sequence)
    return new_sequence

def sequence_length(sequence: str):
    length = 0
    for character in sequence:
        if character in accepted_character:
            length += 1
    return length

def find_complimentary_sequence(sequence: str):
    temp_list = []
    for character in sequence:
        if character in accepted_character:
            temp_list.append(character_dict.get(character))
    return ''.join(temp_list)

sequence = input('Input your sequence here: ')
sequence = check_sequence(sequence.upper())
formatted_sequence_length = sequence_length(sequence)
complimentary_sequence = find_complimentary_sequence(sequence)
reverse_complimentary_sequence = ''

print('The sequence you entered is: 5\'- '+ sequence +' -3\'. \nThe length of the sequence is: ' + str(formatted_sequence_length))
print('The complimentary sequence is: 5\'- ' + complimentary_sequence +' -3\'')
print('The reverse complimentary sequence is: 5\'- ' + reverse_complimentary_sequence +'-3\'')