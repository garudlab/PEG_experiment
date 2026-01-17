import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
from scipy.stats import linregress
import os
import plotting_utils
from plotting_utils import orange_to_red_shades

def to_seconds(t):
    h, m, s = map(int, t.split(':'))
    return h * 3600 + m * 60 + s


def return_growth_rates(df_mean,bnd=15.1/60,smoothing=16):

    ## Apply smoothing, and identify timepoint where OD is changing fastest ("lag")
    ## Smoothing takes a rolling mean of "smoothing" number of points
    ## and we compute delta OD on the smoothed data
    ## helps to identify 

    df_mean_roll = df_mean.rolling(smoothing).mean().shift(-smoothing).dropna()
    lag = df_mean_roll.diff().idxmax()

    ## To identify growth rate, perform a linear regression on log transformed OD data in a time window ("bnd")
    ## surrounding the timepoint where OD is changing most rapidly 
    ## no smoothing here

    growth_rates = pd.Series(index=lag.index)
    lr_dic = {}

    for col in growth_rates.index:

        y = np.log(df_mean.loc[lag.loc[col] - bnd:lag.loc[col] + bnd,col])

        lr = linregress(y.index,y.values)
        growth_rates.loc[col] = lr.slope
        lr_dic[col] = lr
        
    return(growth_rates,lr_dic,lag)

## perform regression on growth rates       
def regress_growth_rate(growth_rates):
    
    lr_growth = linregress(growth_rates.index,growth_rates.values)
    
    return(lr_growth)

########### Plotting utilities ###########

def plot_OD_curves(df_mean,lag=None,color_map=None, file_name=None,save=False, out_dir = None,
                   smoothing=5,with_growth=False,lr_dic=None,semilogy=True,bnd=15.1/60):
    
    ## If haven't specified manually, 
    ## make a color map for plotting purposes (orange-low PEG --> red-high PEG)
    if color_map is None:
        color_map = orange_to_red_shades(df_mean.columns.unique())
    
    fig,ax = plt.subplots(figsize=(12,8))

    if not with_growth:
        
        for col in df_mean.columns:

            ax.plot(df_mean.index,df_mean.loc[:,col].rolling(smoothing).mean(),color=color_map[col],label=col)

    else:
        
        if lag is None:
            df_mean_roll = df_mean.rolling(16).mean().shift(-16).dropna()
            lag = df_mean_roll.diff().idxmax()
            
        ## maxx allows us to zoom in on the growth region
        maxx = [] 
        for col in df_mean.columns:

            xy = np.log(df_mean.loc[lag.loc[col] - bnd:lag.loc[col] + bnd,col])

            ax.scatter(xy.index,np.exp(xy.values),color=color_map[col],label=col,zorder=10,s=40)

            ax.plot(df_mean.index,np.exp(np.log(df_mean.loc[:,col].values)),color=color_map[col],zorder=0,alpha=.3)

            m = lr_dic[col].slope
            b = lr_dic[col].intercept

            xx = np.linspace(min(xy.index)-10/60,max(xy.index)+10/60,1000)
            yy = m*xx + b
            ax.plot(xx,np.exp(yy),color="k",zorder=1)            
            maxx.append(xx[-1])
     
        ax.set_xlim([-0.25,max(maxx)*1.2]) 
        
    ax.set_xlabel("Time (hours)",size=20)
    ax.set_ylabel("OD 600",size=20)

    if semilogy:
        ax.semilogy()
        
    fig.legend(bbox_to_anchor=(1,0.875),title="PEG concentration")

    ## saves
    if save:
        if with_growth and file_name is None:
            file_name = "OD_with_growth"
            if out_dir is None:
                fig.savefig(f"{file_name}",bbox_inches="tight")
            else:
                fig.savefig(f"{out_dir}/{file_name}",bbox_inches="tight")
                
        elif with_growth and file_name is not None:
            if out_dir is None:
                fig.savefig(f"{out_dir}/{file_name}",bbox_inches="tight")
            else:
                fig.savefig(f"{out_dir}/{file_name}",bbox_inches="tight")  
                
        elif not with_growth and file_name is None:
            file_name = "OD_raw"
            if out_dir is None:
                fig.savefig(f"{file_name}",bbox_inches="tight")
            else:
                fig.savefig(f"{out_dir}/{file_name}",bbox_inches="tight")                
        
        elif not with_growth and file_name is not None:
            if out_dir is None:
                fig.savefig(f"{out_dir}/{file_name}",bbox_inches="tight")
            else:
                fig.savefig(f"{out_dir}/{file_name}",bbox_inches="tight") 
                
                
def plot_growth_rate_regression(growth_rates,lr_growth,color_map=None,save=False,out_dir=None,file_name=None):
    
    ## If haven't specified manually, 
    ## make a color map for plotting purposes (orange-low PEG --> red-high PEG)
    if color_map is None:
        color_map = orange_to_red_shades(df_mean.columns.unique())
    
    xx = np.linspace(0,max(growth_rates.index)*1.05,1000)

    fig,ax = plt.subplots(figsize=(12,8))

    for osm in growth_rates.index:

        ax.scatter(osm,growth_rates.loc[osm],color=color_map[osm],s=150,edgecolor="k",label=osm)

    lr_growth = linregress(growth_rates.index,growth_rates.values)
    xx = np.linspace(0,max(growth_rates.index)*1.05,1000)

    ag = lr_growth.intercept
    rg = lr_growth.slope
    ax.plot(xx,xx*rg+ag,color="dodgerblue")

    ax.set_xlabel("PEG concentration",size=20)

    ax.set_ylabel(r"Growth rate ($\lambda$)",size=20)

    fig.legend(bbox_to_anchor=(1,.89),title="PEG concentration")    
 
    if save:
        
        if file_name is None:
            file_name = "growth_rate_regression"
            if out_dir is None:
                fig.savefig(f"{file_name}",bbox_inches="tight")
            else:
                fig.savefig(f"{out_dir}/{file_name}",bbox_inches="tight")

plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False


if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Process OD data, plot figures, and output growth rates"
    )

    parser.add_argument(
        "--OD_data",
        help="Path to the OD data file"
    )
    parser.add_argument(
        "--out_dir",
        default=None,
        help="Directory to write output files"
    )
    parser.add_argument(
        "--output_growth",
        action="store_true",
        default=True,
        help="Filename or path for the growth rate output"
    )

    args = parser.parse_args()

    # Access arguments
    OD_data = args.OD_data
    out_dir = args.out_dir
    output_growth = args.output_growth
    
    if out_dir is not None:
        
        from pathlib import Path
        
        out_dir = Path(out_dir)
        # Check / create output directory
        if out_dir.exists():
            if not out_dir.is_dir():
                raise NotADirectoryError(f"{out_dir} exists but is not a directory")
        else:
            out_dir.mkdir(parents=True, exist_ok=True)
        
    ## read in data, format
    df_mean = pd.read_csv(OD_data,index_col=0)

    df_mean.columns,df_mean.index = [float(f) for f in df_mean.columns], [float(f) for f in df_mean.index]
    df_mean = df_mean.dropna()
    
    ## make a colormap for plotting
    color_map = orange_to_red_shades(df_mean.columns.unique())

    growth_rates,lr_dic,lag = return_growth_rates(df_mean)
    lr_growth = regress_growth_rate(growth_rates)

    plot_OD_curves(df_mean,semilogy=False,out_dir=out_dir,save=True)

    plot_OD_curves(df_mean,with_growth=True,lr_dic=lr_dic,out_dir=out_dir,save=True,lag=lag)

    plot_growth_rate_regression(growth_rates,lr_growth,save=True,out_dir=out_dir)
    
    if output_growth:
        growth_rates.to_csv(f"{out_dir}/growth_rates.txt")