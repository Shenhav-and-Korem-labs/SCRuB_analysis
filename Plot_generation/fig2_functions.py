import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os
import scipy
import pylab as pl
from matplotlib.colors import SymLogNorm
from sklearn.metrics import roc_curve, auc
from scipy.stats import wilcoxon


global_palette = {'Raw':'#4c72b0', 'Restrictive':'#dd8452',
                  'Decontam':'#55a868', 'SCRUB':'#c44e52', 'SCRuB':'#c44e52', 
                  'Decontam (Low Biomass)':'darkgreen',
                  'Decontam (LB)':'darkgreen',
                  'Decontam_LB':'darkgreen',
                  'Decontam (LB)':'darkgreen',
                   'Restrictive (Original)':'#dd8452',
                  'microDecon':'purple',
                 'Input':'#4c72b0',
                 'No decontamination':'#4c72b0'}

sns.set_theme(font_scale=2)
sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'}, 
       font_scale=2)

def make_f2a_plate_plot(a, datasets, plate_locs, col_idx=6, data_dir = 'results/Fig_2/'):
    plate_df = pd.DataFrame( np.hstack(( plate_locs[a], 
                        datasets[a].index.str.contains('NTC')[:, np.newaxis], 
                         datasets[a].iloc[:, col_idx].values[:, np.newaxis]
                        
                        )),
                        columns=['Row', 'Column', 'Control', 'Abundance'],
            
            )
    
    fig=plt.figure(figsize=(9,6))
    ax=sns.scatterplot(x='Column', 
                       y='Row', 
                       data=plate_df,
                       style = 'Control',
                       hue='Abundance',
                       markers = ['o'] +  ['^']*(a in ['raw', 'truth']),
                       hue_norm=SymLogNorm(vmin=0, vmax=1, linthresh=1e-4),
                       palette='copper', #'viridis',
                       s=1100, 
                      edgecolor='black')
    ax.get_legend().remove()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.grid(False)
    plt.axis('off')
    plt.ylim(-.5, 7.5) # -.75, 7.75)
    with plt.rc_context({'image.composite_image': False}):
        plt.savefig( data_dir + 'F2_A_' + a + '.pdf', 
                    dpi=500,
                    bbox_inches='tight', 
                    format='pdf'
                   )
    return(None)

def plot_f2_a(data_dir='results/Plots/fig_2/'):
    raw=pd.read_csv(data_dir+'Raw.csv', index_col=0)
    scrubbed=pd.read_csv(data_dir+'SCRUB.csv', index_col=0)
    decontam=pd.read_csv(data_dir+'Decontam.csv', index_col=0)
    decontam.index=scrubbed.index

    decontam=decontam/decontam.sum(axis=1)[:, np.newaxis]
    decontam[ np.isnan(decontam) ]=0
    scrubbed=scrubbed/scrubbed.sum(axis=1)[:, np.newaxis]
    raw=raw/raw.sum(axis=1)[:, np.newaxis]


    scrubbed=scrubbed.loc[scrubbed.index.str.contains('NTC')==False]
    decontam=decontam.loc[decontam.index.str.contains('NTC')==False]

    true_elements= raw.loc[raw.index.str.contains('Sourc')].values.argmax(axis=1)
    truth=np.zeros_like(raw.loc[raw.index.str.contains('Sourc')].values)
    truth[ range(true_elements.shape[0]), true_elements ]=1

    truth=pd.concat( [ pd.DataFrame( truth, columns=raw.columns, index=raw.loc[raw.index.str.contains('Sourc')].index ), 
                 raw.loc[ (raw.index.str.contains('Sourc')==False )&
                          (raw.index.str.contains('NTC')==False)]*0 ], axis=0 )
    
    truth=pd.concat([ truth, raw.loc[raw.index.str.contains('NTC')]*0 ], axis=0)    
    
    letters=np.array( ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] )
    tmp_fn = lambda x: np.where(letters == x)[0][0]
    def format_numbered_column(df):
        return( np.vstack(( df.apply( lambda row: 7 - tmp_fn(row.name.split('.')[-1][0]) , axis=1).values, 
                            df.apply( lambda row: int( row.name.split('.')[-1][1:] ) , axis=1).values 
                            )).T 
              )
    datasets = {'raw': raw, 
                'decontam':decontam, 
               'scrubbed':scrubbed, 
               'truth':truth}
    
    plate_locs = { a:format_numbered_column(b) for a,b in datasets.items() }
    
    def make_colorbar_plot(out_path='../results/Plots/fig_2/Colorbar.pdf'):
        a = np.array([[0,1]])
        pl.figure(figsize=(9, 1.5))
        img = pl.imshow(a, cmap='copper', #"viridis", 
                        norm = SymLogNorm(vmin=0, vmax=1, linthresh=1e-4))
        pl.gca().set_visible(False)
        cax = pl.axes([0.1, 0.2, 0.8, 0.6])
        pl.colorbar(orientation="horizontal", cax=cax)
        pl.xticks(fontsize=30)
        pl.savefig(out_path, 
                  dpi=400, 
                  bbox_inches='tight')
        return(None)
    
    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
    make_colorbar_plot(out_path='../results/Plots/fig_2/' + 'Colorbar.pdf')
    sns.set(rc={'axes.facecolor':'888888', 'figure.facecolor':'888888'})
    [make_f2a_plate_plot(a, datasets, plate_locs, col_idx=8, data_dir='../results/Plots/fig_2/') for a in datasets]
    
    return(None)

def plot_f2_b(raw, scrubbed, decontam, path='../results/Plots/fig_2/F2_B.pdf'):
    true_elements= raw.argmax(axis=1)
    values=np.hstack(
            [ raw[:,true_elements].diagonal(), 
              decontam[:, true_elements].diagonal(), 
              scrubbed[:,true_elements].diagonal()
            ]
    )
    n_samps=raw.shape[0]
    names = ['Raw']*n_samps + ['Decontam']*n_samps + ['SCRUB'] * n_samps


    plot_df=pd.DataFrame([names, values]).T
    plot_df.columns = ['Dataset', 'Value']

    plt.subplots(1, #figsize =(2.7/2.54, 4.7/2.54) \
    #                                          if True else (2.5,4.5), 
                 dpi=144, 
                figsize=(2.5,8))
    plt.ylim(0,1)

    ax=sns.boxplot(x='Dataset', 
               y='Value', 

               hue_order = ['Raw', 'Restrictive (Original)', 'Decontam', 'SCRUB'],
               data=plot_df, 
               palette=global_palette)

    sns.swarmplot( x='Dataset', y='Value', data=plot_df, 
               order = ['Raw', 'Decontam', 'SCRUB'],
    #               palette=global_palette, 
                 color='black', 
                 size=5)


    ax.set(yticklabels=[])
    plt.xlabel(None)
    plt.ylabel(None)
    ax.axes.get_xaxis().set_visible(False)



    plt.savefig(path, 
               bbox_inches='tight', 
               dpi=900, 
               format='pdf')
    return(None)


def plot_f2_c(raw, scrubbed, decontam, path='../results/Plots/fig_2/F2_C.pdf'):
    n_samps=raw.shape[0]
    
    true_elements= raw.argmax(axis=1)
    false_indices = np.where( np.eye(raw.shape[0])==0 ) 
    values=np.hstack(
            [ [ raw[i][ np.delete( true_elements, i)].sum() for i in range(decontam.shape[0]) ], 
              [ decontam[i][ np.delete( true_elements, i)].sum() for i in range(decontam.shape[0]) ], 
              [ scrubbed[i][ np.delete( true_elements, i)].sum() for i in range(decontam.shape[0]) ]
            ]
    )

    names = ['Raw']*n_samps + ['Decontam']*n_samps + ['SCRUB'] * n_samps


    plot_df=pd.DataFrame([names, values]).T
    plot_df.columns = ['Dataset', 'Value']

    plt.subplots(1, #figsize =(2.7/2.54, 4.7/2.54) \
    #                                          if True else (2.5,4.5), 
                 dpi=144, 
                figsize=(2.5,8))
    # plt.ylim(0,1)

    ax=sns.boxplot(x='Dataset', 
               y='Value', 
               hue_order = ['Raw', 'Restrictive (Original)', 'Decontam', 'SCRUB'],
               data=plot_df, 
               palette=global_palette)

    plt.ylim(-.0001, .1*10)
    plt.yscale('symlog', linthreshy=0.0001)
    sns.swarmplot( x='Dataset', y='Value', data=plot_df, 
               order = ['Raw', 'Decontam', 'SCRUB'],
    #               palette=global_palette, 
                 color='black', 
                 size=5)


    ax.set(yticklabels=[])
    plt.xlabel(None)
    plt.ylabel(None)
    ax.axes.get_xaxis().set_visible(False)

    plt.savefig(path, 
               bbox_inches='tight', 
               dpi=900, 
               format='pdf')
    return(None)



def plot_f2_d(raw, scrubbed, decontam, path='../results/Plots/fig_2/F2_D.pdf'):
    true_elements= raw.argmax(axis=1)
    truth=np.zeros_like(raw)
    truth[ range(16), true_elements ]=1

    scr_dists=scipy.spatial.distance.jensenshannon(scrubbed.T, truth.T)
    raw_dists=scipy.spatial.distance.jensenshannon(raw.T, truth.T)
    dec_dists=scipy.spatial.distance.jensenshannon(decontam.T, truth.T)

    values=np.hstack(
            [ raw_dists, 
              dec_dists, 
              scr_dists
            ]
    )
    n_fps=raw_dists.shape[0]
    names = ['Raw']*n_fps + ['Decontam']*n_fps + ['SCRUB'] * n_fps


    plot_df=pd.DataFrame([names, values]).T
    plot_df.columns = ['Dataset', 'Value']
    
    plt.subplots(1, #figsize =(2.7/2.54, 4.7/2.54) \
#                                          if True else (2.5,4.5), 
             dpi=144, 
            figsize=(2.5,8))
    plt.ylim(0,1)

    ax=sns.boxplot(x='Dataset', 
               y='Value', 
               hue_order = ['Raw', 'Restrictive (Original)', 'Decontam', 'SCRUB'],
               data=plot_df, 
               palette=global_palette)

    # plt.ylim(-.0001, 1)
    # plt.yscale('symlog', linthreshy=0.0001)
    sns.swarmplot( x='Dataset', y='Value', data=plot_df, 
               order = ['Raw', 'Decontam', 'SCRUB'],
    #               palette=global_palette, 
                 color='black', 
                 size=5)


    ax.set(yticklabels=[])
    plt.xlabel(None)
    plt.ylabel(None)
    ax.axes.get_xaxis().set_visible(False)

    plt.savefig(path, 
               bbox_inches='tight', 
               dpi=900, 
               format='pdf')
    return(None)


def make_cont_truth(data_dir, raw):
    raw_complete=pd.read_csv(data_dir+'Raw.csv', index_col=0).values
    raw_complete=raw_complete/raw_complete.sum(axis=1)[:, np.newaxis]


    consideration=( (raw_complete>.15 ).sum(axis=0) > 0 )

    n_taxa=raw.shape[1]    

    contaminant_truth = ( 1 - (np.arange(0,n_taxa,1) == np.unique( raw.values.argmax(axis=1) ) \
            [:, np.newaxis] ).sum(axis=0) ) 
    
    ecr = ( raw.loc[raw.index.str.contains('Control')] .sum(axis=0) > 0 ).astype(float)
    
    return(contaminant_truth, ecr, consideration)

def make_estimated_cont_plot(estimated_contaminants, 
                             datasets, 
                             path='fig_2/F2_E.pdf',
                             data_dir='results/data/Fig_2/', 
                             split_data=None, 
                                   hide_axes=False):
    
    
    if split_data is None:
        contaminant_truth, estimated_contaminants['Restrictive'], consideration = \
                        make_cont_truth(data_dir, datasets['Raw'] )
        
    else:
        tmp = [ make_cont_truth(data_dir[i], split_data[i][0]['Raw'] ) for i in range(len(data_dir)) ]
        contaminant_truth = np.hstack([a[0] for a in tmp])
        estimated_contaminants['Restrictive'] = pd.Series( np.hstack([a[1] for a in tmp]) )
        consideration = np.hstack([a[2] for a in tmp])

    roc_curves = {a:roc_curve(contaminant_truth[consideration], 
                              estimated_contaminants[a].fillna(0).values.flatten()[consideration] )
                         for a in estimated_contaminants }
    

    cont_name_map= {
        'Restrictive':'Restrictive',
        'decont_cont':'Decontam',
        'decontam_low_biomass_cont':'Decontam (LB)',
        'microdecon_cont':'microDecon',
        'scrub_contaminant':'SCRuB'
        }
    plt.rcParams["font.family"] = "Calibri"

    fig=plt.figure(figsize=(8,8))
    for a in cont_name_map:
        plt.plot(roc_curves[a][0],
                 roc_curves[a][1],
                 c=global_palette[cont_name_map[a]], 
                 linewidth=6, 
                 alpha=.8,
                 label = cont_name_map[a] + ': auROC='+str(round( auc(roc_curves[a][0],
                                                                      roc_curves[a][1]), 2 )), 
                 
                )
    
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.2), prop={'size':23})
    plt.xlim(0,1)
    plt.ylim(0,1)
    
    if hide_axes:
        plt.xticks(np.linspace(0,1,2), labels=[]*2)
        plt.yticks(np.linspace(0,1,2), labels=[]*2)
        plt.legend('',frameon=False)

        
    else:
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.xticks(np.linspace(0,1,2))
        plt.yticks(np.linspace(0,1,2))
        
    plt.savefig(path, 
           bbox_inches='tight', 
           dpi=900, 
           format='pdf')
    return(None)
    
    
def make_f2_swarmplot_data(datasets):
    
    datasets_high_prev = {a:datasets[a].loc[datasets[a].index.str.contains('Sink')].values
                         for a in datasets}

    datasets_low_prev = {a:datasets[a].loc[datasets[a].index.str.contains('Sourc')].values
                         for a in datasets}

    n_high_prev_samps=datasets_high_prev['Raw'].shape[0]
    true_elements_high_prev =  np.array([0]*n_high_prev_samps)
    values=np.hstack([datasets_high_prev[a][:, true_elements_high_prev].diagonal() 
                      for a in datasets_high_prev])

    names = np.repeat( [a for a in datasets_high_prev], n_high_prev_samps) 

    plot_df_prev=pd.DataFrame([names, values]).T
    plot_df_prev.columns = ['Dataset', 'Value']
    plot_df_prev['Group'] = 'High-Prevalence'


    true_elements_low_prev = datasets_low_prev['Raw'].argmax(axis=1)
    

    values=np.hstack([datasets_low_prev[a][:, true_elements_low_prev].diagonal() 
                      for a in datasets_low_prev])
    n_samps= datasets_low_prev['Raw'].shape[0]
    names = np.repeat( [a for a in datasets_low_prev], n_samps) 

    plot_df=pd.DataFrame([names, values]).T
    plot_df.columns = ['Dataset', 'Value']
    plot_df['Group']='Low-Prevalence'
    
    return( pd.concat([plot_df_prev, plot_df]) )
    


def make_f2_swarmplot(datasets, path='results/Plots/fig_2/F2_F.pdf', 
                                   hide_axes=False):
    
    if type(datasets)==dict:
        full_plot_df = make_f2_swarmplot_data(datasets)
    else:
        full_plot_df = pd.concat([make_f2_swarmplot_data(datasets[i][0]) for i in range(len(datasets))])

    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
                'axes.edgecolor':'black', 
               'grid.color': 'black'})
    plt.subplots(1, #figsize =(2.7/2.54, 4.7/2.54) \
    #                                          if True else (2.5,4.5), 
                 dpi=144, 
                figsize=(8, 8))
    full_plot_df.Dataset=full_plot_df.Dataset.str.replace('Raw', 'No decontamination')
    
    print(full_plot_df.groupby(['Dataset', 'Group'])['Value'].apply(lambda x: (x < 1e-4).sum() ) )
    
    plt.ylim(-.05,1.05)
    ax=sns.swarmplot(x='Group',
                      y='Value', 
                      hue='Dataset',
                      data=full_plot_df, 
                      hue_order = ['No decontamination', 'Restrictive', 'Decontam', 'Decontam_LB', 'microDecon', 'SCRUB'],
    #                order = ['Raw', 'Decontam', 'SCRUB'],
                      order = ['Low-Prevalence', 'High-Prevalence'],
                      palette=global_palette, 
    #                   color='black', 
                      dodge=True,
                      size=5)
    

    if hide_axes:
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.get_legend().remove()
        plt.axis('off')
    else:
        plt.ylabel('Relative Abundance')
        plt.xlabel('Taxon Classification')
        ax.legend(title=None)
    
    plt.savefig(path, 
               bbox_inches='tight', 
               dpi=900, 
               format='pdf')
    
    return(None)


def make_supplemental_jsd_swarmplot_data( datasets ):
    
    datasets_high_prev = {a:datasets[a].loc[datasets[a].index.str.contains('Sink')].values
                         for a in datasets}

    datasets_low_prev = {a:datasets[a].loc[datasets[a].index.str.contains('Sourc')].values
                         for a in datasets}
    

    n_high_prev_samps=datasets_high_prev['Raw'].shape[0]
    true_elements_high_prev =  np.array([0]*n_high_prev_samps)
    
    true_elements_high_prev = datasets_high_prev['Raw'].argmax(axis=1)
    truth_high_prev=np.zeros_like(datasets_high_prev['Raw'])
    truth_high_prev[ range(true_elements_high_prev.shape[0]), true_elements_high_prev ]=1
    
    values=np.hstack([scipy.spatial.distance.jensenshannon( datasets_high_prev[a].T, 
                                                            truth_high_prev.T )
                      for a in datasets_high_prev])

    names = np.repeat( [a for a in datasets_high_prev], n_high_prev_samps) 

    plot_df_prev=pd.DataFrame([names, values]).T
    plot_df_prev.columns = ['Dataset', 'Value']
    plot_df_prev['Group'] = 'High-Prevalence'


    true_elements_low_prev = datasets_low_prev['Raw'].argmax(axis=1)
    truth_low_prev=np.zeros_like(datasets_low_prev['Raw'])
    truth_low_prev[ range(true_elements_low_prev.shape[0]), true_elements_low_prev ]=1

    values=np.hstack([ scipy.spatial.distance.jensenshannon( datasets_low_prev[a].T, 
                                                            truth_low_prev.T )
                      for a in datasets_low_prev])
    n_samps= datasets_low_prev['Raw'].shape[0]
    names = np.repeat( [a for a in datasets_low_prev], n_samps) 

    plot_df=pd.DataFrame([names, values]).T
    plot_df.columns = ['Dataset', 'Value']
    plot_df['Group']='Low-Prevalence'
    
    return( pd.concat([plot_df_prev, plot_df]) )



def make_supplemental_jsd_swarmplot(datasets, path='../results/Supplementary_Figures/Well_to_well/F2_JSD.pdf', 
                                   hide_axes=False):
    
    if type(datasets)==dict:
        full_plot_df = make_supplemental_jsd_swarmplot_data(datasets)
    else:
        full_plot_df = pd.concat([make_supplemental_jsd_swarmplot_data(datasets[i][0]) for i in range(len(datasets))])
    
    full_plot_df=full_plot_df.fillna(full_plot_df.Value.max())
    print('low-prevalence')
    
    
    print('SCRuB vs microDecon')
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group=='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='microDecon')&(full_plot_df.Group=='Low-Prevalence')
                                                                                    ].Value.values ) )
    
    

    print('SCRuB vs Decontam')
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group=='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Decontam')&(full_plot_df.Group=='Low-Prevalence')
                                                                                    ].Value.values ) )
    
    print('SCRuB vs Decontam (LB')
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group=='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Decontam_LB')&(full_plot_df.Group=='Low-Prevalence')].Value.values ) )
    
    print('SCRuB vs Restrictive')            
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group=='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Restrictive')&(full_plot_df.Group=='Low-Prevalence')].Value.values ) )
    
    
    print('Decontam (LB) vs Raw')            
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='Decontam_LB')&(full_plot_df.Group=='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Raw')&(full_plot_df.Group=='Low-Prevalence')].Value.values ) )
    
    print('Decontam vs Raw')            
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='Decontam')&(full_plot_df.Group=='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Raw')&(full_plot_df.Group=='Low-Prevalence')].Value.values  ) )
    
    
    
    print('High Prevalence')
    
    print('SCRuB vs microDecon')
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group!='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='microDecon')&(full_plot_df.Group!='Low-Prevalence')
                                                                                    ].Value.values ) )


    
    print('SCRuB vs Decontam')
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group!='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Decontam')&(full_plot_df.Group!='Low-Prevalence')
                                                                                    ].Value.values ) )
    
    print('SCRuB vs Decontam (LB')
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group!='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Decontam_LB')&(full_plot_df.Group!='Low-Prevalence')].Value.values , 
         alternative='less') )
    
    print('SCRuB vs Restrictive')            
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group!='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Restrictive')&(full_plot_df.Group!='Low-Prevalence')].Value.values, alternative='less' ) )
    
    
    print('Decontam vs Decontam (LB')
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='Decontam')&(full_plot_df.Group!='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Decontam_LB')&(full_plot_df.Group!='Low-Prevalence')].Value.values , alternative='less' ))
    
    print('Decontam  vs Restrictive')            
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='Decontam')&(full_plot_df.Group!='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Restrictive')&(full_plot_df.Group!='Low-Prevalence')].Value.values , alternative='less' ))
    
    print('SCRuB vs Raw')            
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='SCRUB')&(full_plot_df.Group!='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Raw')&(full_plot_df.Group!='Low-Prevalence')].Value.values, alternative='less' ) )
    
    print('Decontam vs Raw')            
    print(scipy.stats.wilcoxon(full_plot_df.loc[(full_plot_df.Dataset=='Decontam')&(full_plot_df.Group!='Low-Prevalence')].Value.values, 
                               full_plot_df.loc[(full_plot_df.Dataset=='Raw')&(full_plot_df.Group!='Low-Prevalence')].Value.values, alternative='less' ) )
    
    
          
                               

    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
                'axes.edgecolor':'black', 
               'grid.color': 'black'})
    plt.subplots(1,
                 dpi=144, 
                figsize=(8, 8))
    
    full_plot_df.Dataset=full_plot_df.Dataset.str.replace('Raw', 'No decontamination')
    
    plt.ylim(-.05,1.05)
    ax=sns.swarmplot(x='Group',
                     y='Value', 
                     hue='Dataset',
                     data=full_plot_df, 
                     hue_order = ['No decontamination', 'Restrictive', 'Decontam', 'Decontam_LB', 'microDecon', 'SCRUB'],
    #                order = ['Raw', 'Decontam', 'SCRUB'],
                     order = ['Low-Prevalence', 'High-Prevalence'],
                     palette=global_palette, 
    #                   color='black', 
                     dodge=True,
                     size=5, 
                     )

    if hide_axes:
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        ax.axes.get_xaxis().set_visible(False)
        ax.get_legend().remove()
        plt.axis('off')
        
    else:
        plt.ylabel('Jensen-Shannon divergence')
        plt.xlabel('Taxon classification')
        ax.legend(title=None)

    plt.savefig(path, 
               bbox_inches='tight', 
               dpi=900, 
               format='pdf')
    
    return(None)



