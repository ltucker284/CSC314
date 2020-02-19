import os
import re 

def read_file(filename: str) -> str:
    #join the current working directory with the foldername and filename
    #ran into issues just passing in the string even though the files were 
    #in the same directory
    file_path = os.path.join(os.getcwd(), "lab3", filename)
    with open(file_path, 'r') as text_file:
        return text_file.read()

def split_out_genome(file_as_string: str) -> str:
    genome = str()
    #using regex to remove the header from the file
    temp_list = re.split("([G|C|T|A]*[G|C|T|A]*\n)", file_as_string)
    #delete the header 
    del temp_list[0]
    #splits return lists so we must join all the lines together again
    for item in temp_list:
        #checking for empty strings in the split bc we don't care about them
        if item != '':
            #join all the strings of the genome together again as a string
            genome += item.split('\n')[0]
    return genome

#function to create a dict of the genomes from the files provided by passing in
#a list of keys and an unspecified amount of genomes
def create_genome_dict(key_list, *genomes) -> dict:
    #pair a key with a genome (must pass in the genome at the same index as its key)
    keys_and_genomes = zip(key_list, genomes)
    genome_dict = dict()
    for key, genome in keys_and_genomes:
        #create the dict using the current key and set its value to the corresponding genome
        genome_dict[key] = genome
    return genome_dict

#since the sample file is formatted differently we need to split the lines differently
def file_to_dict(filename: str):
    file_as_string = str()
    samples_dict = dict()
    file_path = os.path.join(os.getcwd(), 'lab3', filename)
    with open(file_path, 'r') as text_file:
        file_as_string = text_file.read()
        key_list = list()
        sample_list = list()
        #split each line by newline character
        for item in file_as_string.split('\n'):
            #split out the sample number and append it to the key list
            key_list.append(item.split(' ')[0] + ' ' +item.split(' ')[1])
            #split out the genome sample and append it to the sample list
            sample_list.append(item.split(' ')[2])
            #pair a key with its genome
            keys_and_samples = zip(key_list, sample_list)
            for key, sample in keys_and_samples:
                #create the dict using the current key and set its value 
                #to the corresponding genome sample
                samples_dict[key] = sample
    return samples_dict

#function to compare the samples to the genomes
def analyze_samples(sample_dict, genome_dict):
    count_of_human_matches = 0
    count_of_ecoli_matches = 0
    count_of_raoultella_matches = 0
    count_of_zika_matches = 0
    for key in sample_dict:
        for entry in genome_dict:
            #use re.search to check the full genome for exact matches of a sample
            if re.search(sample_dict[key], genome_dict[entry]):
                #depending on what entry a match is found in, add +1 to the count of matches
                if entry == 'human':
                    count_of_human_matches += 1
                elif entry == 'ecoli':
                    count_of_ecoli_matches += 1
                elif entry == 'raoultella':
                    count_of_raoultella_matches += 1
                elif entry == 'zika':
                    count_of_zika_matches += 1
    return count_of_human_matches, count_of_ecoli_matches, count_of_raoultella_matches, count_of_zika_matches


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

