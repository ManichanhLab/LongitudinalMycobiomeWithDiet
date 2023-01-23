import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Run fungal functional profiling, Input is the diamond blastx output.\n")

parser.add_argument('-i', '--input',
						dest = "input_path",
						action = "store",
						default = None,
						help = "Input ko table.\n",
						required = True)
parser.add_argument('-p', '--pwy_ko',
						dest = "pwy_path",
						action = "store",
						default = None,
						help = "Input to pwy catalog\n",
						required = True)
parser.add_argument('-o', '--output',
						dest = "output_path",
						action = "store",
						default = "",
						help = "Output path of the pwy table.\n",
						required = True)

option=parser.parse_args()

sample_to_process = option.input_path
pwy_to_ko = option.pwy_path
output=option.output_path



ko_output=pd.read_csv(sample_to_process,delimiter="\t")
ko_output=ko_output.set_axis(['gene', 'counts'], axis=1, inplace=False)
ko_output[['ko','species']] = ko_output['gene'].str.split('|',expand=True)
pwy_table=pd.read_csv(pwy_to_ko,delimiter="\t",names=["ko","pwy","class","type"])
pwy_table=pwy_table[["ko","class"]]
merge=ko_output.sort_values('ko').merge(pwy_table, on='ko', how='left')
final_output=pd.DataFrame()
for p in merge["gene"].unique():
 subset=merge.loc[merge["gene"]==p]
 subset["true_counts"]=subset["counts"]/len(subset)
 subset["pathway"]=subset["class"].astype(str)+"|"+subset["species"].astype(str)
 subset=subset[["pathway","true_counts"]]
 final_output=final_output.append(subset,ignore_index = True,sort=True)
final_output=final_output.replace({"nan": "UNGROUPED"}, regex=True)
final_output=final_output.replace({"\|None": ""}, regex=True)
final_output=final_output.groupby(['pathway']).sum()
final_output.to_csv(output,index=True,sep='\t')