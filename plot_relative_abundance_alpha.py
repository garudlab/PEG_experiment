import pandas as pd
import matplotlib.pyplot as plt
import config
import numpy as np
import matplotlib.patches as mpatches
from scipy import stats
import matplotlib.cm as cm
import seaborn as sns
import matplotlib.lines as mlines
from scipy.stats import linregress
from sklearn.decomposition import PCA
import matplotlib.gridspec as gridspec

def take_triu(df):
    
    N = df.shape[0]
    p=np.triu_indices(N,k=1)
    
    return(df[p])


def split_name(name):
    S = name.split("_")
    
    return("_".join([S[1],S[2]]))

def abbreviate(name):
    
    s = name.split(" ")
    
    s1 = s[0][0].upper()
    s2 = s[1][0].lower()
    
    if name == "Collinsella stercoris":
        s2 = "st"
    elif name == "Clostridium sporogenes":
        s2 = "sp"
        
    return(s1+s2)

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


peg_colors = {s:cm.rainbow(i/5.02) for i,s in enumerate([0,2,5,10,15])}
species_colors = cm.tab20(np.linspace(0, 1, 10))

Z=zip([0, 2.5, 5,10, 15],["#636363", "#31a354", "#ffb404","#e6550d", "#a50f15"])

peg_colors = {z[0]:z[1] for z in Z}

sex_shapes = {"M":"o","F":"D"}
sex_colors = {"M":"k","F":"dodgerblue"}

plt.rcParams['xtick.labelsize']=20
plt.rcParams['ytick.labelsize']=20


plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False

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


######## INSERT PATH TO RELATIVE ABUNDANCE HERE ########
spec_dir = f"{config.base_dir}/midas_output/merged_midas_output/species"
df = pd.read_csv(f"{spec_dir}/relative_abundance.txt",index_col=0,sep="\t")

species_list_dic = pd.Series({s:split_name(s) for s in config.good_species})

df.index = species_list_dic.loc[df.index]

df.index = df.index.map(species_name_dic)

df = process_df_cols(df)

df_mean = df.groupby("PEG").mean().T
df_mean = df_mean/df_mean.sum()


######## FORMATTING METADATA ########

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

######## CALCULATING ALPHA DIVERSITY ########

alpha = (-df.T*np.log(df.T.replace(0,np.nan)).replace(np.nan,0)).sum()

df["osmolality"] = df_meta.loc[df.index,"Osmolality"]
df = df.set_index("osmolality",append=True)

df_meta = df_meta.set_index("Osmolality",append=True)
df_meta["Osmolality"] = df_meta.index.get_level_values("Osmolality")
df_meta.index.names = ["PEG","Sex","Mouse_num","osmolality"]

idx_osm = df_meta.sort_values("Osmolality").index


######## PLOTTING ########

fig = plt.figure(figsize=(18,12),tight_layout=True)
gs = gridspec.GridSpec(10, 2,wspace=0,hspace=0)

ax_alpha = fig.add_subplot(gs[:2, :])
ax_stack = fig.add_subplot(gs[2:, :])

    

ax_alpha.spines["bottom"].set_linewidth(2)
ax_alpha.spines["left"].set_linewidth(2)
ax_stack.spines["left"].set_linewidth(2)
ax_stack.spines["bottom"].set_linewidth(2)

ax_alpha.set_xlim([df_meta["Osmolality"].min()*.99, df_meta["Osmolality"].max()*1.01])
ax_alpha.set_ylabel(r"$\alpha$ diversity",size=30)
ax_alpha.set_ylim([alpha.min()*.8,alpha.max()*1.02])

ax_stack.stackplot(df_meta.loc[idx_osm,"Osmolality"],df.loc[idx_osm].T,
             labels=df.columns,colors = species_colors);

df = df.sort_index(level="osmolality")
df = df.loc[:,df.mean().sort_values().index]

alpha = (-df.T*np.log(df.T.replace(0,np.nan)).replace(np.nan,0)).sum()

taxa = df.columns 
data = df.values
osm = [350] + list(df.index.get_level_values("osmolality")) + [900]


ax_alpha.stackplot(alpha.index.get_level_values("osmolality"),alpha.values,alpha=.3)


for idx in alpha.index:
    peg,sex,mnum,osm = idx
    ax_alpha.scatter(osm,alpha.loc[idx],s=150,
                     color=peg_colors[peg],
                     marker=sex_shapes[sex],
                     edgecolor=sex_colors[sex])
    
        
ax_stack.set_xlabel("Osmolality",size=30)
ax_stack.set_ylabel("Relative abundance",size=30)

ax_stack.set_xlim([df_meta["Osmolality"].min()*.99, df_meta["Osmolality"].max()*1.01])
ax_stack.set_ylim([0,1])

plt.tight_layout()

fig.legend(prop={"size":20},bbox_to_anchor=(1.3,1));
