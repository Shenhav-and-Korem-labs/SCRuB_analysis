##############
## Functions to plot results for Figure 1
############


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib_venn import venn3_unweighted
import matplotlib
import seaborn as sns
import os
import scipy
from scipy.stats import wilcoxon

global_palette = {'Raw':'#4c72b0', 'Restrictive':'#dd8452',
                  'Decontam':'#55a868', 'SCRUB':'#c44e52', 'SCRuB':'#c44e52', 
                  'Decontam (Low Biomass)':'darkgreen',
                  'Decontam (LB)':'darkgreen',
                  'Decontam_LB':'darkgreen',
                   'Restrictive (Original)':'#dd8452',
                 'Input':'#4c72b0',
                 'No decontamination':'#4c72b0', 
		'microDecon':'purple'}
sns.set_theme(font_scale=2)

sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'}, 
       font_scale=2)



def plot_f1c(data, out_path='../results/Plots/fig_1/F1_C.pdf', global_palette=global_palette, 
            hide_axes=False):
    plt.subplots(1,
                 figsize=(11.02, 10),
                 dpi=500)
    
    tmp=data.loc[(data.names!='Spatial SCRUB')&(data.well_to_well_level==0)]
    tmp=tmp.groupby(['names', 'contam_level', 'well_to_well_level', 
                                'round'])['jsds'].median().reset_index()
    tmp.names=tmp.names.str.replace('SCRUB', 'SCRuB').str.replace('Input', 'No decontamination')
    
    ax=sns.boxplot( x='contam_level', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order =['No decontamination', 'Restrictive', 'Decontam', 'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                width=.75,
                   palette=global_palette, 
                   fliersize=0
               )
    handles, labels = ax.get_legend_handles_labels()
    
    sns.stripplot( x='contam_level', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam', 'Decontam (Low Biomass)','microDecon', 'SCRuB'],
                  dodge=True,
                  ax=ax, 
                  size=2.5, 
                  color='k'
               )
    
    
    # plot vs scrub significance asterisks
    
    # 'worse than' test
    q=tmp[['contam_level', 'names']].drop_duplicates()
    q['sig_v_raw'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names=='No decontamination' ) ]['jsds'].values, 
        alternative='less').pvalue
                        if row.names!='No decontamination' else 1][0],
           axis=1)
    print(q)

    # 'worse than' test
    q=tmp[['contam_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)

    print(q)

    q['y_val']=1.09
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    
    if sum(q.is_sig > 0):
        sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)','microDecon', 'SCRuB'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon','SCRuB'],
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                   )
    
# #     # 'better than' test
    q=tmp[['contam_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.contam_level==row.contam_level)&
                    (tmp.names=='SCRuB') ]['jsds'].values, 
        alternative='less').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)

    q['y_val']=.003
    q['is_sig'] =  q.sig_v_scrub < 1e-4

    if sum(q.is_sig)>0:
        sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='contam_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                       )
    
    ax.legend_.set_title(None)

    # Put the legend out of the figure
    plt.legend(bbox_to_anchor=(1.01, .7), loc=2, borderaxespad=0.)
    ax.set_title(None)
    ax.set_yscale('log')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.set_yticks([.01, .1, .25, 5, .75, 1])
    plt.ylabel('Median Jensen-Shannon Divergence')
    plt.xlabel('Contamination')
    plt.ylim(2.5e-3, 1.2)
    plt.ylim(0.01, 1.2)
    
    
    ax.legend(handles[:6], labels[:6])
    
    if hide_axes:
        ax.get_legend().remove()
        ax.set_xticks([])
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title(None)
        

    plt.savefig(out_path, 
                dpi=900,
               bbox_inches='tight', 
               format='pdf')
    return(None)

def plot_f1d(data, out_path='../results/Plots/fig_1/F1_D.pdf', global_palette=global_palette, 
            hide_axes=False):
    data=data.loc[data.names!='SCRUB']
    data.loc[data.names=='Spatial SCRUB', 'names'] = 'SCRuB'
    
    plt.subplots(1, #figsize =(2.7/2.54, 4.7/2.54) \
    #                                          if True else (2.5,4.5), 
                 figsize=(11.02, 10),
                 dpi=500)

    tmp=data.loc[(data.contam_level==.05)&(data.well_to_well_level != 0)]
    
    tmp=tmp.groupby(['names', 'well_to_well_level', 
                                'round'])['jsds'].median().reset_index()
    
    tmp.names=tmp.names.str.replace('Input', 'No decontamination')
    
    ax=sns.boxplot( x='well_to_well_level', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam', 'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                width=.75,
                palette=global_palette, 
                fliersize=0
               )
    
    sns.stripplot( x='well_to_well_level', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam', 'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                  dodge=True,
                  ax=ax, 
                  size=2.5, 
                  color='k'
               )
    
    handles, labels = ax.get_legend_handles_labels()
    
    ##  plot vs scrub significance asterisks    
    
    # 'worse than' test
    q=tmp[['well_to_well_level', 'names']].drop_duplicates()
    q['sig_v_raw'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names=='No decontamination' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='No decontamination' else 1][0],
           axis=1)
    print(q)


    
    # 'worse than' significance'
    q=tmp[['well_to_well_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)
    print(q)


    q['y_val']=1.09
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    if sum(q.is_sig > 0):

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)','microDecon', 'SCRuB'],
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                       )
    
        # better than 
    q=tmp[['well_to_well_level', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.well_to_well_level==row.well_to_well_level)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='less').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)

    q['y_val']=.003
    q['is_sig'] =  q.sig_v_scrub < 1e-4

    if q.is_sig.sum() > 0:
        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='well_to_well_level', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                       )
    
    ax.legend_.set_title(None)
    
    # Put the legend out of the figure
    plt.legend(bbox_to_anchor=(1.01, .7), loc=2, borderaxespad=0.)
    ax.set_title(None)
    
    ax.set_yscale('log')
    plt.ylabel('Median Jensen-Shannon Divergence')
    plt.xlabel('Well-to-well leakage')
    ax.set_yticks([.01, .1, .25, 5, .75, 1])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:6], labels[:6])
    
    plt.ylim(2.5e-3, 1.2)
    plt.ylim(0.01, 1.2)
    
    if hide_axes:
        ax.get_legend().remove()
        ax.set_xticks([])
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title(None)
    
    
    plt.savefig(out_path, 
                dpi=900,
                bbox_inches='tight',
                format='pdf')
    return(None)






def make_no_contam_plots(hide_axes=False):
    
    
    ## venn diagram
    taxa_results=pd.read_csv('../results/data/Fig_1/No_contamination_removed_taxa.csv', 
                        index_col=0)
    A=taxa_results.SCRuB.values
    B=taxa_results.microDecon.values
    C=taxa_results['Decontam (LB)'].values

    # Create the numbers for the diagram
    only_a = sum(A & ~B & ~C)
    only_b = sum(B & ~A & ~C)
    only_c = sum(~B & ~A & C)
    a_and_b = sum(A & B & ~C)
    a_and_c = sum(A & ~B & C)
    b_and_c = sum(~A & B & C)
    a_b_and_c  = sum(A & B & C)



    colors=[global_palette[a] for a in \
                                        ['SCRUB', 'microDecon', 'Decontam (LB)']]
    plt.figure(figsize=(10,10))

    venn3_unweighted(subsets = (only_a, only_b,  a_and_b, only_c, a_and_c, b_and_c, a_b_and_c), 
          set_labels = ('SCRuB', 'microDecon', 'Decontam (LB)'), 
                    set_colors=colors)
    
    plt.savefig('../results/Supplementary_figures/Fig_S2_a.pdf', 
                dpi=900,
                bbox_inches='tight',
                format='pdf')
    
    
    ## boxplot
    results=pd.read_csv(
        '../results/data/Fig_1/no_contamination_simulation.csv',
                       index_col=0)
    results=results.groupby(['names', 'rounds', 'smp_tps'])['jsds'].median().reset_index()
    
    plt.figure(figsize=(12,7))
    ax=sns.boxplot(y='jsds', x = 'smp_tps', data=results,
                   hue='names', 
                hue_order=['Decontam (LB)', 'Decontam', 'microDecon', 'SCRUB'],
               linewidth=3.7, 
               palette=global_palette, 
                   fliersize=0
               )

    sns.stripplot(y='jsds', x = 'smp_tps', hue='names', data=results,
                hue_order=['Decontam (LB)', 'Decontam', 'microDecon', 'SCRUB'],
               linewidth=3.7, 
               color='black',
                dodge=True,
                size=2.5
               )
    ax.legend().set_title(None)
    plt.ylabel('Jensen-Shannon Divergence')
    plt.xlabel('Contaminant Biodiversity')
    plt.legend()
    plt.yscale('symlog', linthreshy=0.0001)

    ax.axes.get_xaxis().set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:4], labels[:4])

    if hide_axes:
        ax.get_legend().remove()
        ax.set_xticks([])
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title(None)
    
    plt.savefig('../results/Supplementary_figures/Fig_S2_b.pdf', 
                dpi=900,
                bbox_inches='tight',
                format='pdf')
    return(None)
    

    
    

def noise_level_plot(hide_axes=False):

    data=pd.read_csv('../results/data/Fig_1/Noise_varying_simulations.csv')


    data=data.loc[data.names!='SCRUB']
    data.loc[data.names=='Spatial SCRUB', 'names'] = 'SCRuB'
    data.n_taxa=data.n_taxa.str.replace('med', 'medium').str.capitalize()

    plt.subplots(1, #figsize =(2.7/2.54, 4.7/2.54) \
    #                                          if True else (2.5,4.5), 
                 figsize=(11.5, 10),
                 dpi=500)

    tmp=data.copy()#loc[(data.contam_level==.05)&(data.well_to_well_level != 0)]

    tmp=tmp.groupby(['names', 'well_to_well_level', 'n_taxa',
                                'round'])['jsds'].median().reset_index()

    tmp.names=tmp.names.str.replace('Input', 'No decontamination')

    ax=sns.boxplot( x='n_taxa', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                            'Decontam (Low Biomass)', 'microDecon','SCRuB'],
               order = ['Low', 'Medium', 'High'],
                width=.75,
                palette=global_palette, 
                fliersize=0
               )

    sns.stripplot( x='n_taxa', 
                   y='jsds', hue = 'names', data=tmp, 
               hue_order = ['No decontamination', 'Restrictive', 'Decontam',
                            'Decontam (Low Biomass)', 'microDecon','SCRuB'],
               order = ['Low', 'Medium', 'High'],
                  dodge=True,
                  ax=ax, 
                  size=2.5, 
                  color='k'
               )

    handles, labels = ax.get_legend_handles_labels()

    ##  plot vs scrub significance asterisks    

    # 'worse than' significance'
    q=tmp[['n_taxa', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.n_taxa==row.n_taxa)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.n_taxa==row.n_taxa)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='greater').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)


    q['y_val']=1.09
    q['is_sig'] =  q.sig_v_scrub < 1e-4
    if sum(q.is_sig > 0):

        sns.swarmplot( x='n_taxa', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)','microDecon', 'SCRuB'],
               order = ['Low', 'Medium', 'High'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='n_taxa', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon','SCRuB'],
               order = ['Low', 'Medium', 'High'],
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                       )

        # better than 
    q=tmp[['n_taxa', 'names']].drop_duplicates()
    q['sig_v_scrub'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.n_taxa==row.n_taxa)&
                    (tmp.names==row.names) ]['jsds'].values, 
        tmp.loc[ (tmp.n_taxa==row.n_taxa)&
                    (tmp.names=='SCRuB' ) ]['jsds'].values, 
        alternative='less').pvalue
                        if row.names!='SCRuB' else 1][0],
           axis=1)

    q['y_val']=.003
    q['is_sig'] =  q.sig_v_scrub < 1e-4

    if q.is_sig.sum() > 0:
        sns.swarmplot( x='n_taxa', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
               order = ['Low', 'Medium', 'High'],
                       marker='+',
                      size=25/2, 
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2
                       )

        sns.swarmplot( x='n_taxa', y='y_val', hue='names', 
                        data = q.loc[q.is_sig], 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 
                                   'Decontam (Low Biomass)', 'microDecon', 'SCRuB'],
               order = ['Low', 'Medium', 'High'],
                       marker='x',
                      size=17.5/2,
                        ax=ax,
                        palette=global_palette, 
                      dodge=True, 
                      color='black',
                      linewidth=3.5/2,
                        edgecolor="black"
                       )

    ax.legend_.set_title(None)

    # Put the legend out of the figure
    plt.legend(bbox_to_anchor=(1.01, .7), loc=2, borderaxespad=0.)
    ax.set_title(None)

    ax.set_yscale('log')
    plt.ylabel('Median Jensen-Shannon Divergence')
    plt.xlabel('Noise level')
    ax.set_yticks([.01, .1, .25, 5, .75, 1])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:6], labels[:6], loc='upper center', bbox_to_anchor=(0.5,-0.2),)

    plt.ylim(2.5e-3, 1.2)# .2)# 1.2)
    
    if hide_axes:
        ax.get_legend().remove()
        ax.set_xticks([])
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title(None)
 
    
    plt.savefig('../results/Supplementary_figures/Fig_S2_varying_noise.pdf', 
                dpi=900,
                bbox_inches='tight',
                format='pdf')
    return(None)
