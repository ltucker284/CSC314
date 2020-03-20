from Bio import SeqIO # to parse Seq data
#åfrom Bio.Seq import Seq # for Seq

exon_count = 0
first_five = str()

with open('MAOA.gb.txt') as handle :

    # gets an iterator that allows you to go through each entry
    sequence_iter = SeqIO.parse(handle, "genbank")

    # gets next sequence (which is the 1st one here)
    seq_record = next(sequence_iter)

    print('Accession #: ', seq_record.id)
    print('Description: ', seq_record.description)
    print('seq length = ', len(seq_record))

# loop through each of the features
for feature in seq_record.features :

    # print out chromosome
    if feature.type == 'source':
        print('=================== Chromosome ==================')
        print(feature.qualifiers['chromosome'][0])
        chromosome = feature.qualifiers['chromosome'][0]

    # print out gene feature
    if feature.type == 'gene':
        print('=================== Gene ==================')
        print(feature)

    # print out exon 13 (this can be used to help count the exons)
    if feature.type == 'exon': 
        exon_count += 1
        if feature.qualifiers['number'][0] == '13':
            print('=================== Exon 13 ==================')
            print(feature)

    # print out CDS location 
    if feature.type == 'CDS': 
        print('=================== CDS location ==================')
        end_cds_position = feature.location.parts[-1] 
        print("CDS end position:", end_cds_position) 
        first_five_location = feature.location.parts[0]
        first_five = seq_record.seq[first_five_location.start:first_five_location.start+5]
        last_three_location = feature.location.parts[-1]
        last_three = seq_record.seq[last_three_location.end - 9:last_three_location.end]

print('1. The gene is on chromosome:', chromosome)
print('2. The first 5 nucleotides are:', first_five)
print('3. The number of exons contained is:', exon_count)
print('4. The last 3 codons of the protein are:', last_three.transcribe())
print('5. The last 3 codons code for:', last_three.translate())