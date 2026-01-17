import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from scipy.stats import linregress
import os,glob
import config

import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sns
import matplotlib.lines as mlines
from matplotlib.lines import Line2D

import scipy.stats as stats
from sklearn.linear_model import LinearRegression
from OD_curve_analysis import to_seconds
from plotting_utils import orange_to_red_shades,peg_colors,species_colors,sex_shapes,sex_colors
from metadata_utils import split_name,process_df_cols,species_name_dic,species_rev_dic
import metadata_utils
import OD_curve_analysis

fig_dir = config.fig_dir

def read_gr_data(gr_loc="OD_analysis/*/*_growth_rate"):
    
    gr_dic = {}
    normalize = False

    for g in glob.glob(gr_loc):

        spec_name = g.split("/")[1]
        gg = pd.read_csv(g,index_col=0)

        gg = gg.sort_index()

        if normalize:
            gg.index = gg.index - 350
           

        gg = gg.clip(lower=0)

        gr_dic[spec_name] = gg
    
    return(gr_dic)
    
def return_gr_regs(gr_dic):
    
    gr_lr_dic = {}
    for key,growth_rates in gr_dic.items():

        g = "r"

        Lg = OD_curve_analysis.regress_growth_rate(growth_rates["r"])

        rg = Lg.slope

        ag = Lg.intercept

        gr_lr_dic[key] = [rg,ag]

    gr_lr_dic = pd.DataFrame(gr_lr_dic,index=["slope","intercept"]).T
    gr_lr_dic = gr_lr_dic.sort_values('intercept')

    ## Estimated growth rate at PEG = 350
    gr_lr_dic["350"] = gr_lr_dic.slope*350 + gr_lr_dic.intercept
    
    ## Percentage drop and basal growth rate
    for s,k in gr_dic.items():
    
        gr_lr_dic.loc[s,"perc_drop"] = (k.iloc[0]["r"]/(gr_lr_dic.loc[s,"slope"]))
        gr_lr_dic.loc[s,"base_growth"] = k.iloc[0]["r"]

    
    return(gr_lr_dic)

## returns abundance data for either relative abundance or ddPCR
def return_abundance_data_choice(choice):

    if choice == "ddPCR":

        dfA = pd.read_csv("05172024-PEGd_ddPCR_Data_10microbes.tsv",sep="\t")

        L = []
        for idx in dfA.index:
            L.append([float(dfA.loc[idx,"PEG"]),dfA.loc[idx,"Sample"].split("_")[1][0],
                      int(dfA.loc[idx,"Sample"].split("_")[1][1]),dfA.loc[idx,"Species"]])

        L = pd.DataFrame(L,columns=["PEG","Sex","Mouse_num","Species"])

        L = pd.MultiIndex.from_frame(L)    

        df_dd = pd.DataFrame(dfA["Undiluted_Copies/ul"].values,index=L,columns=["ddPCR"])

        df_dd = df_dd.unstack()

        df_dd.columns = df_dd.columns.droplevel(0)

        c = list(df_dd.columns)
        c[2] = 'B theta'
        c[-1] = "M_intestinale"
        df_dd.columns = c

        df_dd["osmolality"] = df_meta["Osmolality"]
        df_dd = df_dd.set_index("osmolality",append=True)

        df = df_dd

        df.columns = [species_name_dic[r.replace(" ","_")] for r in df.columns]
        
        df.columns = df.columns.map(species_rev_dic)
        
    elif choice == "Relative abundance":

        df = metadata_utils.read_abundance_data()

        df.columns = df.columns.map(species_rev_dic)

        df_mean = df.groupby("PEG").mean().T
        df_mean = df_mean/df_mean.sum()

        df["osmolality"] = df_meta["Osmolality"]
        df = df.set_index("osmolality",append=True)

    return(df)
    
def plot_growth_osm_slope_intercept(gr_lr_dic,saveloc):

    fig,axs = plt.subplots(2,1,figsize=(8,16))

    ax = axs[0]
    for idx,row in gr_lr_dic.iterrows():

        xx = np.linspace(300,1200,100)

        if idx == "Eubacterium_rectale":
            continue

        if idx in hatch_species:
            ax.scatter(gr_dic[idx].index,gr_dic[idx]["r"],color=species_colors.loc[idx],edgecolor="red",hatch="xx",
                       zorder=10,s=80)
        else:
            ax.scatter(gr_dic[idx].index,gr_dic[idx]["r"],color=species_colors.loc[idx],edgecolor="k",zorder=10,s=80)


        ax.plot(xx,xx*row[0] + row[1],lw=6,color=species_colors.loc[idx],zorder=1,alpha=.4)

    ax.set_xlabel("PEG Concentration",size=15)
    ax.set_ylabel(r"Growth rate $\lambda$",size=15)
    ax.axhline(0,color="grey")

    ax = axs[1]

    for species in gr_lr_dic.index:
        if species in hatch_species:

            ax.scatter(gr_lr_dic.loc[species,"base_growth"],
                       gr_lr_dic.loc[species,"perc_drop"],
                       color=species_colors.loc[species],s=500,label=species,edgecolor="red",
                       hatch="xx")

        else:

            ax.scatter(gr_lr_dic.loc[species,"base_growth"],
                       gr_lr_dic.loc[species,"perc_drop"],
                       color=species_colors.loc[species],s=500,label=species,edgecolor="k")




    ax.set_xlabel("Base growth rate",size=15)

    ax.set_ylabel("% growth drop per 100mOsm",size=15)


    fig.legend(bbox_to_anchor=(1.1,.9))

    fig.savefig(saveloc)



def plot_growth_abundance_sample(col,df,lr_ab,lines,saveloc):

    fig,ax = plt.subplots(figsize=(8,8))

    ax.set_title(f"Osmolality: {col[-1]}",size=20)

    good_idxs = np.isfinite(df[col]) 

    xx = np.linspace(0,lines[col].max(),100)

    yy = lr_ab.loc[col].slope*xx + lr_ab.loc[col].intercept

    cc = df[col][good_idxs]
    ll = lines[col][good_idxs]


    for species in good_idxs.loc[good_idxs].index:

        ax.scatter(ll.loc[species],cc.loc[species],s=150,edgecolor="k",
                   color=species_colors.loc[species],label=species)

        ax.plot(xx,yy,color="grey",alpha=.5)

    ax.set_ylabel(choice,size=20)

    ax.set_xlabel(r"Predicted growth rate $\lambda$",size=20)

    fig.savefig(saveloc,bbox_inches="tight")
 

def plot_full_community_growth(df,lines,saveloc):
    
    fig,axs = plt.subplots(6,5,figsize=(16,15))
    axs = axs.ravel()

    for i in range(df.shape[1]):

        col = df.columns[i]
        ax=axs[i]
        ax.set_title(int(col[-1]))

        good_idxs = np.isfinite(df[col]) 

        xx = np.linspace(0,lines[col].max(),100)

        yy = lr_ab.loc[col].slope*xx + lr_ab.loc[col].intercept

        cc = df[col][good_idxs]
        ll = lines[col][good_idxs]

        for species in good_idxs.loc[good_idxs].index:

            ax.scatter(ll.loc[species],cc.loc[species],s=150,edgecolor="k",
                       color=species_colors.loc[species],label=species)

            ax.plot(xx,yy,color="grey",alpha=.5)


        plt.subplots_adjust(hspace=0.5)     

        ax.set_xlim([-0.05,.9])


        if choice == "Relative abundance":
            ax.set_ylim([-0.1,0.7])
        else:
            ax.set_ylim([.5*min(df[col]) - .2*1e5,max(df[col])*1.2])

    fig.text(0.42,0.06,r"Predicted growth rate $\lambda$",size=20)

    fig.text(0.06,0.42,choice,size=20,rotation=90)

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label=species,
               markerfacecolor=color, markeredgecolor="k", markersize=20)
        for species, color in species_colors.items()
    ]

    # Add the custom legend
    ax.legend(handles=legend_elements, title='Species',bbox_to_anchor=(2,8.9))
    fig.savefig(saveloc,bbox_inches="tight")

    
def plot_community_growth_lr_slope(lr_ab,choice,saveloc):
    
    fig,ax = plt.subplots(figsize=(8,8))

    for idx in lr_ab.index:
        peg,sex,mnum,osm = idx

        ax.scatter(osm,lr_ab.loc[idx,"slope"],marker=sex_shapes[sex],
                  color=peg_colors[peg],edgecolor=sex_colors[sex],s=200)

    ax.set_ylabel(r"Slope:" + f"\n{choice} x growth rate",size=20)

    ax.set_xlabel("Osmolality",size=20)

    fig.savefig(saveloc,bbox_inches="tight")

    

## predict growth rate at osmolarity given in vitro growth rates at that osmolarity
def return_growth_rate_prediction(gr_lr_dic):
 
    xx = df_meta["Osmolality"].values
    
    lines = {}
    for spec_name,row in gr_lr_dic.iterrows():

        lines[spec_name] = pd.Series(row["slope"]*xx + row["intercept"],xx)

    lines = pd.DataFrame(lines).T
    lines = lines.where(lines > 0,0)
    
    return(lines)

## regress abundance against growth rate
def regress_abundance(df,lines):
    
    lr_ab = {}
    for col in df.columns:

        good_idxs = np.isfinite(df[col])
        lr_ab[col] = linregress(lines[col][good_idxs],df[col][good_idxs])

    lr_ab = pd.DataFrame(lr_ab,index=["slope","intercept","rvalue","pvalue","stderr"])
    lr_ab = lr_ab.T
    lr_ab.index.names = df.columns.names    
    
    return(lr_ab)
 
    
    
fig_dir = config.fig_dir

plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

df_meta = metadata_utils.return_metadata()

species_colors = pd.Series([tuple(s) for s in species_colors],index=species_name_dic.values())
species_colors.index = species_colors.index.map(species_rev_dic)

gr_dic = read_gr_data()

gr_lr_dic = return_gr_regs(gr_dic)

gr_lr_dic.loc["Eubacterium_rectale"] = np.nan



df = metadata_utils.read_abundance_data()

df.columns = df.columns.map(species_rev_dic)

df_mean = df.groupby("PEG").mean().T
df_mean = df_mean/df_mean.sum()


### identify species which increase/decrease in relative abundance btwn PEG 0-15
df_mean_change = ((df_mean[0] - df_mean[15]).sort_values(ascending=False) >= 0)

df = df[df_mean_change.index]

hatch_species = df_mean_change.loc[df_mean_change].index


### plot figure showing linear relationship between growth rate and osmolarity
### and also slope/intercept of that relationship
### NOTE: feel free to modify the slope/intercept code to whatever form you see fit
plot_growth_osm_slope_intercept(gr_lr_dic,saveloc=f"{fig_dir}growth_rate_summary")


choice_list = ["ddPCR","Relative abundance"]

for choice in choice_list:
    
    df = return_abundance_data_choice(choice)
    
    #df.columns = df.columns.map(species_rev_dic)

    ## predict growth rate at each osmolarity
    lines = return_growth_rate_prediction(gr_lr_dic)
    
    #lines.columns = df.columns
    
    #df = df.T
    print(df.head())
    print("\n\n\n")
    print(lines.head())
    
    df = df.sort_index(level="osmolality")
    df = df[lines.index]
    df = df.T
    lines = lines.sort_index(axis=1)

    lines = lines.loc[df.index]
    lines.columns = df.columns


    lr_ab = regress_abundance(df,lines)
    


    
    ### Here, choose which 
    i = -1
    col = df.columns[i]

    colname = '_'.join([str(c) for c in col])
    colname = colname.replace(".0","")

    plot_growth_abundance_sample(col,df,lr_ab,lines,f"{fig_dir}/{colname}_growth_abundance_{choice}")



    plot_full_community_growth(df,lines,f"{fig_dir}growth_rate_community_{choice}")



    plot_community_growth_lr_slope(lr_ab,choice,f"{fig_dir}/community_growth_lr_slope")




