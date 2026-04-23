import Bio
from Bio import SeqIO
from pathlib import Path
import pandas as pd

gbkFile = Path("/home/abdeali/viralR_test_output/Chlamy_punui/Chlamy_punui_contig.gbk")

gbk = SeqIO.read(open(gbkFile, "r"), "genbank")
annotfile = Path("/home/abdeali/viralR_test_output/Chlamy_punui/Chlamy_punui_contig.tsv")
summ_file = pd.read_table(annotfile, header=0, sep="\t")

print(summ_file.columns)
# print(gbk.features[0].qualifiers['locus_tag'])
#print(f"IDs in file: {gbk.features[0]} ")

feature_keys = [i.qualifiers.get('locus_tag') for i in gbk.features]

# print(summ_file['HMM_hit'].loc[summ_file['query'] == "contig_536_1"].values)

for record in gbk.features:
    query = record.qualifiers.get('locus_tag')[0]
    # print(query)
    # print(summ_file.loc[summ_file['query'] == query , 'HMM_hit'].item() )
    record.qualifiers['product'] = summ_file.loc[summ_file['query'] == query , 'Description'].item()
    
print(gbk)

outputpath = Path("/home/abdeali/viralR_test_output/Chlamy_punui/Chlamy_punui_contig_ncvog_desc.gbk")

with open(outputpath, 'w') as outfile :
    SeqIO.write(gbk, outfile, 'genbank')