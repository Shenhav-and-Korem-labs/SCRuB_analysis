import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib
import seaborn as sns
import os
import scipy
global_palette = {'Raw':'#4c72b0', 'Restrictive':'#dd8452',
                  'Decontam':'#55a868', 'SCRUB':'#c44e52', 
                  'Restrictive (Original)':'#dd8452',
                  'Input':'#4c72b0',}
sns.set_theme(font_scale=2)
sns.set(rc={'axes.facecolor':'888888', 'figure.facecolor':'888888'})


def simulate_plate(colors = ['blue', 'red'], 
                  color_negs=True):
    
    N=2#8
    M=2#12
    N=8
    M=12

    nx=(M+1)*100
    ny=(N+1)*100
    data = np.random.random((nx,ny))

    for a in range(M+1):
        for b in range(N+1):
            data[ (100*a - 50):(100*a + 50), (100*b - 50):(100*b + 50)].sort(axis=0)
            data[ (100*a - 50):(100*a + 50), (100*b - 50):(100*b + 50)].sort(axis=1)
            data[ (100*a - 50):(100*a + 50), (100*b - 50):(100*b + 50)] = \
                    np.power(data[ (100*a - 50):(100*a + 50), (100*b - 50):(100*b + 50)], 
                                        np.random.randint(3,5))

    data=data.T


    # top = cm.get_cmap('Oranges', 128)
    # bottom = cm.get_cmap('Blues', 128)
    # newcolors = np.vstack((top(np.linspace(1, .3, 128)),
    #                        bottom(np.linspace(.3, 1, 128))))
    # newcmp = ListedColormap(newcolors, name='OrangeBlue')

    norm=plt.Normalize(-2,2)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    cmap_samp=matplotlib.colors.LinearSegmentedColormap.from_list("", colors+['blue']+['blue'])
    blank_map = matplotlib.colors.LinearSegmentedColormap.from_list("", ['white', 'white'])

    fig,ax = plt.subplots(1,  figsize = (9,6) )
    rad=41
    t_rad=38
    
    # Now, loop through coord arrays, and create a circle at each x,y pair
    for a in range(1, N+1):
        for b in range(1, M+1):
            if (a in [1, N])&(b in [1, M]):
                ta=100*b
                tb=100*a
                
                circ = patches.Polygon(np.vstack([[ta-t_rad, tb-t_rad], 
                                                  [ta+t_rad, tb-t_rad], 
                                                  [ta, tb + ( (N/M)*np.sqrt(2)*t_rad)], [
                                                ta-t_rad, tb-t_rad]]), 
                                       transform=ax.transData, 
                                       linewidth=4, 
                                       edgecolor='black', 
                                       facecolor='none')
                
                ax.add_patch(circ)
                ax.imshow(1-data, clip_path=circ, origin='lower', cmap=[cmap if color_negs else blank_map][0])
            else:
                circ = patches.Circle((100*b,100*a),radius=rad, transform=ax.transData, 
                                       linewidth=4, 
                                      edgecolor='black',
                                      facecolor='none')
                ax.add_patch(circ)
                ax.imshow(data, clip_path=circ, origin='lower', cmap=cmap_samp )
                
            
                
    #         plt.show()
    #     ax.add_patch(circ)

    ax.grid(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlim(50, 1250)
    ax.set_ylim(50, 850)
    plt.axis('off')
    return(ax)