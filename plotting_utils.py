import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sns
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
import numpy as np

def red_shades(values):
    # Normalize values between 0 and 1
    norm = (np.array(values) - np.min(values)) / (np.max(values) - np.min(values))
    
    # Map to colors from light red (e.g., rgb(255, 200, 200)) to dark red (e.g., rgb(150, 0, 0))
    colors = {}
    for i in range(len(norm)):
        n = norm[i]
        r = 255
        g = int(200 * (1 - n))
        b = int(200 * (1 - n))
        colors[values[i]] = mcolors.to_hex((r/255, g/255, b/255))
    
    return colors

def orange_to_red_shades(values):
    # Normalize values to [0, 1]
    norm = (np.array(values) - np.min(values)) / (np.max(values) - np.min(values))
    
    # Define RGB for orange and red
    orange = np.array([255, 165, 0])
    red = np.array([255, 0, 0])
    
    colors = {}
    for i in range(len(norm)):
        n = norm[i]
        rgb = (1 - n) * orange + n * red
        
        rgb_normalized = tuple(rgb / 255.1)
                
        colors[values[i]] = mcolors.to_hex(rgb_normalized)
    
    return colors

peg_colors = {s:cm.rainbow(i/5.02) for i,s in enumerate([0,2,5,10,15])}
species_colors = cm.tab20(np.linspace(0, 1, 10))

Z=zip([0, 2.5, 5,10, 15],["#636363", "#31a354", "#ffb404","#e6550d", "#a50f15"])

peg_colors = {z[0]:z[1] for z in Z}

sex_shapes = {"M":"o","F":"D"}
sex_colors = {"M":"k","F":"dodgerblue"}