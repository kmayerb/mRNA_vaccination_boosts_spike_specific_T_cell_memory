"""
To compare clonal trajectory consider clonal frequency between adjacent 
samples in participant timeseries. Fisher's Exact Test are efficiently 
computed using fishersapi.

Notes:
  2022-12-15 - This file was originally named (tcr_analyis/new_integrative_fisher_expansion_analysis.py)
  It critically makes the __outer_joined__ files which are used to evaluate 
  clonal expansion between adjacent timepoints.
  
  * FDR value based on Fisher's Exact Testing
  * This script also made a control file which recored all of the outer join.
  * We refactored it here to handle file renamed with PUBID instead of PTIDs

  Original Analysis Done, Nov 28, 2021. 
"""
import os 
import platform
import pandas as pd
from fishersapi import fishers_vec
from statsmodels.stats.multitest import multipletests

# <r> is the resource files with tcrdist_input_files 
# Dec 29, 2021 '/fh/fast/gilbert_p/fg_data/mcelrath_covid/tcrdist3ready_06_04_2021'
# Dec 19, 2022 PUBID refactored, project directory can be repo home folder
project_dir = '/fh/fast/gilbert_p/fg_data/mcelrath_covid/mrna_cd8_tcell_project'
r = os.path.join(project_dir, 'data','tcrb','tcrb_deep_repertoires')
output_dir = os.path.join(project_dir, 'data','tcrb','longitudinal')


assert os.path.isdir(r)

# List all files in <r> that end .ready.tsv
fs = [f for f in os.listdir(r) if f.endswith('ready.tsv')]
# Split the filenames by . then by _ to create a list subject, vsit, TCRB
fs2 = [f.split(".")[0].split("_") for f in fs]

# Assemble a sorted triplet for each subject, such that samples are sorted
# in longtiudinal order by visit number.
d = dict()
for f, info in zip(fs, fs2):
    ptid,visit,chain = info
    ptid = ptid
    visit = int(visit.replace("v",""))
    d.setdefault(info[0], []).append((visit, f))

# <triplets> is a dictionary keyed on subject, with a list of tuple 
# Each tuple (visit, filename), in the ascending order, pre, post1, post2
triplets = {ptid:sorted(item, key = lambda x:x[0]) for ptid, item in d.items() if len(item) == 3}

# Next we define the four potential longitudinal comparisons
all_pre_post1_comparison = list()
all_pre_post2_comparison = list()
all_pre_post1_post2_comparison = list()
all_post1_post2_comparison = list()
for k,v in triplets.items():
    all_pre_post1_comparison.append(['post1_pre', r, v[0][1], v[1][1]])
    all_pre_post2_comparison.append(['post2_pre', r, v[0][1], v[2][1]])
    all_pre_post1_post2_comparison.append(['pre_post1_post2', r, v[0][1], v[1][1], v[2][1] ])
    all_post1_post2_comparison.append(['post2_post1', r, v[1][1], v[2][1]])
# For later analysis we want a key that is based on V Family and CDR3, for matching
# legacy datasets like Emerson.
def add_vfam_key(df):
    df = df[df['v_b_gene'].apply(lambda x: isinstance(x, str))].reset_index(drop = True)
    df['v_family'] = df['v_b_gene'].apply(lambda x : x.split("*")[0].split("-")[0])
    df['key'] = df.apply(lambda x: f"{x['v_family']}+{x['cdr3_b_aa']}", axis = 1)
    return(df)

# This loops through all of the specified comparisons. 
for tag, r,f1,f2 in all_pre_post1_comparison + all_pre_post2_comparison + all_post1_post2_comparison: 
    if not f1.startswith('CO001-897362'):
      continue
    print(f"FROM SOURCE FILES FOUND IN {r}")
    print(f"\tCOMPARING {f2} TO PRIOR SAMPLE {f1}")
    df1 = pd.read_csv(os.path.join(r,f1), sep = "\t") 
    df2 = pd.read_csv(os.path.join(r,f2), sep = "\t") 
    
    df1 = add_vfam_key(df1)
    df2 = add_vfam_key(df2)
    df21 = df2.merge(df1, how = "outer", on = ['subject','cdr3_b_aa','v_b_gene','j_b_gene','cdr3_rearrangement','key'])
    df21['productive_frequency_y'] = df21['productive_frequency_y'].fillna(0)
    df21['templates_y'] = df21['templates_y'].fillna(0)
    df21['productive_frequency_x'] = df21['productive_frequency_x'].fillna(0)
    df21['templates_x'] = df21['templates_x'].fillna(0)
    
    n1 = df1['templates'].sum()
    n2 = df2['templates'].sum()
    df21['a'] = df21['templates_x']
    df21['b'] = n2 - df21['templates_x']
    df21['c'] = df21['templates_y']
    df21['d'] = n2 - df21['templates_y']
    odds, pv = fishers_vec(df21['a'],df21['b'],df21['c'],df21['d'],  alternative="two-sided")
    df21['post_pre_OR'] =  odds
    df21['post_pre_pvalue'] =  pv 
    r,pvc,null1,null2 = multipletests(pvals = pv, method = "fdr_bh")
    df21['post_pre_fdr_bh'] = pvc
    r,pvc,null1,null2 = multipletests(pvals = pv, method = "holm")
    df21['post_pre_holm'] = pvc
    
    # Adds information about 
    #df21m = df21.merge(ref_df, how = "left", on = "key", validate = 'many_to_one' )
    ptid = f1.split("_")[0]
    f1tag = f1.split(".")[0]
    f2tag = f2.split(".")[0]
    fout = os.path.join(output_dir, tag, f"{ptid}__{f2tag}__outer_joined__{f1tag}.tsv")
    #fout = f'/fh/fast/gilbert_p/fg_data/mcelrath_covid/vaccine_samples_vfam_cdr3_vs_emerson_immunerace/{tag}/{ptid}__{f2tag}__outer_joined__{f1tag}.tsv'
    print(fout)
    df21.to_csv(fout, sep = "\t", index = False)

# RECORD ALL OF THE FILES AND OUTER JOINS FOR LATER ANALYSIS STEPS

# <fouts> is a list of all of outer_join files created above
fouts = list()
for tag, r,f1,f2 in all_pre_post1_comparison + all_pre_post2_comparison + all_post1_post2_comparison:
    ptid = f1.split("_")[0]
    f1tag = f1.split(".")[0]
    f2tag = f2.split(".")[0]
    fout = os.path.join(output_dir, tag, f"{ptid}__{f2tag}__outer_joined__{f1tag}.tsv")
    #fout = f'/fh/fast/gilbert_p/fg_data/mcelrath_covid/vaccine_samples_vfam_cdr3_vs_emerson_immunerace/{tag}/{ptid}__{f2tag}__outer_joined__{f1tag}.tsv'
    fouts.append(fout)
    
control_df = pd.DataFrame(all_pre_post1_comparison + all_pre_post2_comparison + all_post1_post2_comparison, 
    columns = ['tag','r','f1','f2'])
control_df['filename'] = fouts
# write the control file
control_df.to_csv(
    os.path.join(project_dir, "control", '2021_12_19_control_df_outerjoins.tsv'),
    #'/fh/fast/gilbert_p/fg_data/mcelrath_covid/vaccine_samples_vfam_cdr3_vs_emerson_immunerace/control/2021_11_28_control_df_outerjoins.tsv', 
    sep="\t",
    index = False)



