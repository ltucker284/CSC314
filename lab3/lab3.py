import os
import re 

def read_file(filename: str) -> str:
    file_path = os.path.join(os.getcwd(), "lab3", filename)
    with open(file_path, 'r') as text_file:
        return text_file.read()

def split_out_genome(file_as_string: str) -> str:
    genome = str()
    temp_list = re.split("([G|C|T|A]*[G|C|T|A]*\n)", file_as_string)
    del temp_list[0]
    for item in temp_list:
        if item != '':
            genome += item.split('\n')[0]
    return genome

def create_genome_dict(key_list, *genomes) -> dict:
    keys_and_genomes = zip(key_list, genomes)
    genome_dict = dict()
    for key, genome in keys_and_genomes:
        genome_dict[key] = genome
    return genome_dict

def file_to_dict(filename: str):
    file_as_string = str()
    samples_dict = dict()
    file_path = os.path.join(os.getcwd(), 'lab3', filename)
    with open(file_path, 'r') as text_file:
        file_as_string = text_file.read()
        key_list = list()
        sample_list = list()
        for item in file_as_string.split('\n'):
            key_list.append(item.split(' ')[0] + ' ' +item.split(' ')[1])
            sample_list.append(item.split(' ')[2])
            keys_and_samples = zip(key_list, sample_list)
            for key, sample in keys_and_samples:
                samples_dict[key] = sample
    return samples_dict

def analyze_samples(sample_dict, genome_dict):
    count_of_human_matches = 0
    count_of_ecoli_matches = 0
    count_of_raoultella_matches = 0
    count_of_zika_matches = 0
    count_of_unknown_matches = 0
    for key in sample_dict:
        for entry in genome_dict:
            if re.search(sample_dict[key], genome_dict[entry]):
                if entry == 'human':
                    count_of_human_matches += 1
                elif entry == 'ecoli':
                    count_of_ecoli_matches += 1
                elif entry == 'raoultella':
                    count_of_raoultella_matches += 1
                elif entry == 'zika':
                    count_of_zika_matches += 1
                else: count_of_unknown_matches += 1
    return count_of_human_matches, count_of_ecoli_matches, count_of_raoultella_matches, count_of_zika_matches, count_of_unknown_matches


samples_dict = file_to_dict('samples.txt')
human = split_out_genome(read_file('human.txt'))
ecoli = split_out_genome(read_file('ecoli.txt'))
raoultella = split_out_genome(read_file('raoultella.txt'))
zika = split_out_genome(read_file('zika.txt'))
key_list = 'human', 'ecoli', 'raoultella', 'zika'

# In order to match the correct genome to the correct key,
# you must pass in the variables in the same order as they are in
# the key list
genome_dict = create_genome_dict(key_list, human, ecoli, raoultella, zika)
print(f"The length of each sequence is as follows:\n"
        + "zika: " + str(len(genome_dict["zika"])) + "\n"
        + "ecoli: " + str(len(genome_dict["ecoli"])) + "\n"
        + "raoultella: " + str(len(genome_dict["raoultella"])) + "\n"
        + "human: " + str(len(genome_dict["human"])))
print(analyze_samples(samples_dict, genome_dict))

