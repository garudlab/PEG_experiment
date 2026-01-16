import pandas as pd
import numpy as np
import config 

def process_df_cols(dfx):
    
    X = [s.split("_")[:-1] for s in dfx.columns]
    X2 = []

    for x in X:
        y = [x[0],x[1][0],x[1][1]]
        X2.append(y)

    X2 = np.array(X2)

    dfx.columns = pd.MultiIndex.from_arrays(X2.T)
    dfx.columns.names = ["PEG","Sex","Mouse_num"]

    dfx = dfx.T
    
    for level_to_change in [0,2]:
        dfx.index = dfx.index.set_levels(dfx.index.levels[level_to_change].astype(int), level=level_to_change)

    dfx = dfx.sort_index(level="PEG")
    
    p = dfx.index.levels[0]
    p2 = p + (p == 2)*0.5
    dfx.index = dfx.index.set_levels(p2, level=0)
    
    return(dfx)

def return_metadata():
    
    df_meta = pd.read_csv("PEG_metadata.csv")
    df_meta = df_meta.loc[df_meta["Type"] == "Cecal"]
    df_meta = df_meta.loc[df_meta["Species"].isna()]
    df_meta = df_meta.loc[~df_meta["Mouse_num"].isna()]
    df_meta["PEG"] = df_meta.loc[:,"PEG"].astype(int)
    df_meta["Mouse_num"] = df_meta.loc[:,"Mouse_num"].astype(int)
    X = pd.MultiIndex.from_frame(df_meta[["PEG","Sex","Mouse_num"]])
    df_meta = df_meta.drop(["PEG","Sex","Mouse_num"],axis=1)
    df_meta.index=X

    df_meta = df_meta.sort_index(level="PEG")

    p = df_meta.index.levels[0]
    p2 = p + (p == 2)*0.5
    df_meta.index = df_meta.index.set_levels(p2, level=0)

    return(df_meta)

def split_name(name):
    S = name.split("_")
    
    return("_".join([S[1],S[2]]))


def read_abundance_data(abundance_loc="relative_abundance.txt"):

    df = pd.read_csv(abundance_loc,index_col=0,sep="\t")

    species_list_dic = pd.Series({s:split_name(s) for s in config.good_species})

    df.index = species_list_dic.loc[df.index]

    df.index = df.index.map(species_name_dic)

    df = process_df_cols(df)
    
    return(df)


# species_name_dic = {'Akkermansia_muciniphila':'A_muciniphila',
#  'Bacteroides_ovatus':'B_ovatus',
#  'B_theta':'B_thetaiotaomicron',
#  'Clostridium_sporogenes':'C_sporogenes',
#  'Collinsella_stercoris':'C_stercoris',
#  'Enterococcus_faecalis':'E_faecalis',
#  'Escherichia_coli':'E_coli',
#  'Eubacterium_rectale':'Eubacterium rectale',
#  'Faecalibacterium_prausnitzii': 'F_prausnitzii',
#  'M_intestinale':'G6'}

species_name_dic = {'Akkermansia_muciniphila':'Akkermansia muciniphila',
 'Bacteroides_ovatus':'Bacteroides ovatus',
 'B_theta':'Bacteroides thetaiotaomicron',
 'Clostridium_sporogenes':'Clostridium sporogenes',
 'Collinsella_stercoris':'Collinsella stercoris',
 'Enterococcus_faecalis':'Enterococcus faecalis',
 'Escherichia_coli':'Escherichia coli',
 'Eubacterium_rectale':'Eubacterium rectale',
 'Faecalibacterium_prausnitzii':'Faecalibacterium prausnitzii',
 'M_intestinale':'Muribaculum intestinale'}                
             

species_rev_dic = {item:key for key,item in species_name_dic.items()}