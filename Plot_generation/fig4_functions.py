import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import scipy
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import fdrcorrection

global_palette = {'Raw':'#4c72b0', 'Restrictive':'#dd8452',
                  'Decontam':'#55a868', 'SCRUB':'#c44e52', 
                   'Restrictive (Original)':'#dd8452',
                 'Input':'#4c72b0',}
sns.set_theme(font_scale=2)
sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'})



def n_cont_plot(sim_data, global_palette=global_palette):
    plt.subplots(1, 
                 dpi=500, 
                figsize=(8,8))

    ax=sns.boxplot( x='n_controls', y='jsds', data=sim_data, 
                width=.9,
                   color=global_palette['SCRUB'],
                        boxprops={'lw':1}, 
                        capprops={'color':'k', 'lw':1},
                        whiskerprops={'color':'k', 'lw':1},
                    fliersize=0
               )
    
    sns.swarmplot( x='n_controls', y='jsds', data=sim_data,
                  color='black', 
                  size=3
               )
    
    plt.xlabel(None)
    plt.ylim(0, .01)
    ax.axes.get_xaxis().set_visible(False)
    plt.ylabel('Jensen-Shannon Divergence')
    ax.set(yticklabels=[])
    plt.xlabel(None)
    plt.ylabel(None)
    return(ax)


def plate_format_plot(sim_data, global_palette=global_palette):
    plt.subplots(1, 
                 dpi=500, 
                figsize=(8,8))

    ax=sns.boxplot( x='plate_format', y='jsds', data=sim_data.loc[sim_data.n_controls!=1], 
                width=.9,
                   color=global_palette['SCRUB'],
                        boxprops={'lw':1}, 
                        capprops={'color':'k', 'lw':1},
                        whiskerprops={'color':'k', 'lw':1},
                    fliersize=0
               )
    
    sns.swarmplot( x='plate_format', y='jsds', data=sim_data.loc[sim_data.n_controls!=1], 
                   color='black', 
                 size=3)
    
    # plt.ylabel('Jensen-Shannon Divergence')
    plt.xlabel(None)
    plt.ylim(0, .01)#.08)
    ax.axes.get_xaxis().set_visible(False)
    # ax.axes.get_yaxis().set_visible(False)
    plt.ylabel('Jensen-Shannon Divergence')
    ax.set(yticklabels=[])
    plt.xlabel(None)
    plt.ylabel(None)
    return(ax)


def median_confint_bootstrap(data, confidence=0.95):
    np.random.seed(0)
    n_pts=data.shape[0]
    
    boots= np.array( 
                [np.median(data.values[ np.random.choice(np.arange(n_pts),
                                            int(n_pts*.25), replace=False ) ]) for i in range(2500)]
                    )
    
    return(    np.linspace(np.percentile(boots, 0),
                           np.percentile(boots, 1)
                          )
          )

def one_v_all_test(qq, pf, nc, sig=1e-2):
    return( mannwhitneyu( qq.loc[ (qq.plate_format==pf)&
                          (qq.n_controls==nc)].jsds, 
                  qq.loc[ ((qq.plate_format==pf)&
                          (qq.n_controls==nc)) == False ].jsds,
                         alternative='less'
                ).pvalue 
          )

def one_v_all_test_no_rounding(qq, pf, nc, sig=1e-2):
    return( mannwhitneyu( qq.loc[ (qq.plate_format==pf)&
                          (qq.n_controls==nc)].jsds, 
                  qq.loc[ ((qq.plate_format==pf)&
                          (qq.n_controls==nc)) == False ].jsds,
                         alternative='less'
                ).pvalue
          )

def make_heatmap(qq, out_path = 'tmp.pdf', sig=0.05, hide_axes=False):

    tmp = qq.groupby(['plate_format', 'n_controls'])['jsds'].apply(lambda x: median_confint_bootstrap(x))\
                        .explode().reset_index()

    tmp['indexer'] = tmp.n_controls.astype(str) + '_' + \
                        ( 1 - ( ( tmp.index%150 ) < 10 ).astype(int) ).astype(str) + \
                        ( tmp.index%150 ).astype(str)
    
    print(tmp.jsds.min(), tmp.jsds.max())
    
    plt.figure(figsize=(15,6))
    ax=sns.heatmap( data = tmp.pivot("plate_format", "indexer", "jsds").astype(float),
                cmap='Blues_r',
                cbar=[False if hide_axes else True][0], 
               )
    
    one_r_format=qq[['plate_format', 'n_controls']].drop_duplicates()
    one_r_format['one_v_rest_sig'] = one_r_format.apply(lambda row: 
                                  one_v_all_test_no_rounding(qq, row.plate_format, row.n_controls, 
                                                sig=sig), axis=1)
    ## multiple test correction
    one_r_format.one_v_rest_sig=fdrcorrection(one_r_format.one_v_rest_sig.values)[1] 
    print(one_r_format)
    one_r_format.one_v_rest_sig = ( one_r_format.one_v_rest_sig < sig ) * 100
        

    one_r_format['n_controls'] = 50 * np.log2( one_r_format['n_controls'] ) - 25
    one_r_format['plate_format'] = .5 + one_r_format.plate_format.apply(lambda r: np.where(np.array(['A', 'B', 'C', 'D', 'E', 'F']) == r )[0][0])


    sns.scatterplot( 'n_controls', 'plate_format',  
                    data = one_r_format.loc[one_r_format.one_v_rest_sig > 0 ], 
    #                markers='*'
                    marker='x', s=900/4, color='whitesmoke', 
                    ax=ax,
                    edgecolor="none"
                   )
    
    sns.scatterplot( 'n_controls', 'plate_format',  
                    data = one_r_format.loc[one_r_format.one_v_rest_sig > 0 ], 
    #                markers='*'
                    marker='+', s=1500/4, color='whitesmoke', 
                    ax=ax,
                    edgecolor="none"
                   )
    
    if hide_axes:
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        ax.axes.get_xaxis().set_visible(False)

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.grid(False)
        plt.axis('off')
        
    else:
#         q=['']*49
        plt.xticks([25, 75, 125], labels=['2', '4', '8'], rotation=0)
        plt.yticks(np.linspace(.5, 5.5, 6), labels=['Random',
                           'Cluster in corner',
                           'Column on edge',
                           'Column within plate',
                           'Spread within plate',
                           'Spread across edges'][::-1], 
                  rotation=0)
        plt.xlabel('Number of negative controls')
        plt.ylabel(None)
    
    plt.savefig(out_path, 
                dpi=400,
                bbox_inches='tight', 
                format='pdf'
               )

    
def make_f4_sim_supplemental(sim_results, 
                             out_path='../results/Supplementary_figures/', 
                             hide_axes=False):
    
    
    ## general n control boxplot
    plt.figure(figsize=(12,7))
    ax=sns.boxplot(y='jsds', x = 'n_controls', data=sim_results,
                    color=global_palette['SCRUB'],
                   linewidth=3.7
               )
#     sns.swarmplot(y='jsds', x = 'n_controls', data=sim_results,
#                     color='black'
#                )
    
    ax.legend().set_title(None)
    plt.ylabel('Median Jensen-Shannon divergence')
    plt.xlabel('Number of negative controls')
    
    plt.legend()
    plt.ylim(0, .02)

    if hide_axes:
        ax.set(yticks = [0, .005, .01, .015, .02],
               yticklabels=[])
        ax.get_legend().remove()
        plt.xlabel(None)
        plt.ylabel(None)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
    else:
        ax.set(yticks = [0, .005, .01, .015, .02])
        
        
    plt.savefig(out_path+'F_S6_a.pdf', 
                dpi=400,
                bbox_inches='tight', 
                format='pdf'
               )
    
    
    ## general plate format boxplot
    plt.figure(figsize=(12,7))
    ax=sns.boxplot(y='jsds', x = 'plate_format', data=sim_results,
                    color=global_palette['SCRUB'],
                   linewidth=3.7
               )
    
#     sns.swarmplot(y='jsds', x = 'plate_format', data=sim_results,
#                     color='black'
#                )
    
    ax.legend().set_title(None)
    plt.ylabel('Median Jensen-Shannon divergence')
    plt.xlabel(None)
    
    plt.legend()
    plt.ylim(0, .02)

    if hide_axes:
        ax.set(yticks = [0, .005, .01, .015, .02],
               yticklabels=[])
        ax.get_legend().remove()
        plt.xlabel(None)
        plt.ylabel(None)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
    else:
        ax.set(yticks = [0, .005, .01, .015, .02])
        plt.xticks( np.linspace(0, 5, 6), labels=['Random',
                           'Cluster in corner',
                           'Column on edge',
                           'Column within plate',
                           'Spread within plate',
                           'Spread across edges'][::-1], rotation=90 )
        
        plt.savefig(out_path+'F_S6_b.pdf', 
                dpi=400,
                bbox_inches='tight', 
                format='pdf'
               )
               
    
    
    ## general heatmap
        
    make_heatmap(sim_results, out_path=out_path+'F_S6_c.pdf', hide_axes=hide_axes)
    return(None)




def make_s7_plots():
    results=pd.read_csv('../results/data/Fig_4/Fig_S7_Simulations.csv', index_col=0)
    results=results.loc[results.names=='Spatial SCRUB']
    results=results.groupby(['sample_type', 'round', 'well_to_well_level', 
                             'contam_level', 'n_controls'])['jsds'].median().reset_index()
    results.n_controls=results.n_controls.astype(str)
        
    tmp_idx = -1
    nm_list = ['Gut', 'Skin', 'Vaginal']
    for smp_tp in ['high', 'med', 'low']:    
        tmp_idx += 1 
        for contam_lvl in [.05, .25, .5]:
            
            plt.figure(figsize=(10, 7))
            ax = sns.lineplot('well_to_well_level', 'jsds', hue='n_controls', 
                             data=results.loc[(results.sample_type==smp_tp)&
                                                 (results.contam_level==contam_lvl)], 
                        linewidth=4)

            plt.xlim(0.05, .5)

            plt.legend(title='Number of controls', loc='upper left') 

            plt.ylabel('Median Jensen-Shannon divergence')
            plt.xlabel('Level of well-to-well leakage')
        
            plt.savefig('../results/Supplementary_figures/Fig_S7_{}_{}_contam.pdf'.format(nm_list[tmp_idx], contam_lvl),
                        dpi=500, 
                        bbox_inches='tight', 
                        format='pdf'
                       )








