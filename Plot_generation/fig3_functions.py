import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import scipy
#from matplotlib_venn import venn3, venn3_unweighted
from sklearn.metrics import roc_auc_score, auc, roc_curve
from scipy.stats import wilcoxon
# import upsetplot

global_palette = {'Raw':'#4c72b0', 
                  'No decontamination':'#4c72b0', 
                  'Restrictive':'#dd8452',
                  'Decontam':'#55a868',
                          'Decontam (LB)':'darkgreen', 
                          'SCRUB':'#c44e52', 'SCRuB':'#c44e52', 
                  'Restrictive (Original)':'#dd8452',
                  'Input':'#4c72b0',
                  'decontam (Poore et. al)':'#55a868',
                  'decontam (LB)':'darkgreen', 
                  'decontam':'#55a868', 
                   'microDecon':'purple',
                  'Custom':'#FF69B4'}

sns.set_theme(font_scale=2)
sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'})

def plot_cancer_shannon_divs(data, out_path='../results/Plots/fig_3/F3_E.pdf', palette=global_palette, 
                            hide_axes=False):
    plt.subplots(1, #figsize =(2.7/2.54, 4.7/2.54) \
        #                                          if True else (2.5,4.5), 
                 dpi=144, 
                figsize=(3.5,8))

    print(data.Dataset.value_counts())
    
    data.loc[data.Dataset=='Restrictive (Original)', 'Dataset']='Custom'
    
    print(data.Dataset.value_counts())
    
    ax=sns.boxplot( x='Dataset', y='Shannon', data=data, 
                   order = ['Raw', 'Restrictive', 'Decontam',
                            'Decontam (LB)', 'Custom', 'microDecon', 'SCRUB'],
               hue_order = ['Raw', 'Restrictive', 'Decontam', 
                            'Decontam (LB)', 'Custom', 'microDecon', 'SCRUB'],
                width=.9,
                   palette=palette, 
                   fliersize=0
               )
    sns.swarmplot(x='Dataset', y='Shannon', data=data, 
                   order = ['Raw', 'Restrictive', 'Decontam',
                            'Decontam (LB)', 'Custom','microDecon', 'SCRUB'],
                 color='black', size=2.5, 
                  ax=ax
                 )
    
    tmp=data.copy()
    
        # plot vs raw significance asterisks
    q=data[['Dataset']].drop_duplicates()
    q['sig_v_raw'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.Dataset==row.Dataset)]['Shannon'].values, 
        tmp.loc[ (tmp.Dataset=='Raw') ]['Shannon'].values, 
        alternative='less').pvalue
                        if row.Dataset!='Raw' else 1][0],
           axis=1)
    
    print(q)

    q['y_val']=-0.2
    q['is_sig'] =  q.sig_v_raw < 1e-3

    sns.swarmplot( x='Dataset', y='y_val', 
                    data = q.loc[q.is_sig], 
                   order = ['Raw', 'Restrictive', 'Decontam', 'Decontam (LB)', 'Custom', 'microDecon', 'SCRUB'],
               hue_order = ['Raw', 'Restrictive', 'Decontam', 'Decontam (LB)', 'Custom', 'microDecon', 'SCRUB'],
                   marker='+',
                  size=25/2, 
                    ax=ax,
                    palette=global_palette, 
                  dodge=True, 
                  color='black',
                  linewidth=3.5/2
                   )

    sns.swarmplot( x='Dataset', y='y_val',
                    data = q.loc[q.is_sig], 
                   order = ['Raw', 'Restrictive', 'Decontam', 'Decontam (LB)', 'Custom', 'microDecon', 'SCRUB'],
               hue_order = ['Raw', 'Restrictive', 'Decontam', 'Decontam (LB)', 'Custom', 'microDecon', 'SCRUB'],
                   marker='x',
                  size=17.5/2,
                    ax=ax,
                    palette=global_palette, 
                  dodge=True, 
                  color='black',
                  linewidth=3.5/2,
                    edgecolor="black"
                   )
    
    if hide_axes:
        ax.set(yticklabels=[])
        plt.xlabel(None)
        plt.ylabel(None)
        ax.axes.get_xaxis().set_visible(False)
        
    else:
        plt.ylabel('Shannon Diversity')
        plt.xlabel(None)
        plt.xticks(rotation=90)
        
    
    plt.savefig(out_path, 
               bbox_inches='tight', 
               dpi=900, 
               format='pdf')
    return(None)

def make_cancer_venndiagram(out_path='../results/Plots/fig_3/F3_D.pdf',
                            colors=[global_palette[a] for a in \
                                    ['Restrictive (Original)',  'Decontam (LB)', 'SCRUB']]):
    plt.figure(figsize=(10,10))

    # Custom text labels: change the label of group A
    v=venn3_unweighted(subsets=['']*10,
                       set_labels = ('', '', ''), 
                       set_colors=colors,
                       normalize_to=1, 
                       alpha=.8)
   
    plt.savefig(out_path, 
                bbox_inches='tight', 
                dpi=900, 
                format='pdf')
    return(None)





def plot_ICI_rocs(melanoma_rocs, out_path='../results/Plots/fig_3/F3_F.pdf', 
                            hide_axes=False):
    
    melanoma_rocs['Dataset']=melanoma_rocs.Dataset.str.replace('Restictive', 'Restrictive') 
    melanoma_rocs['Dataset']=melanoma_rocs.Dataset.str.replace('SCRUB', 'SCRuB')\
                                    .str.replace('_lb', ' (LB)').str.replace('Microdecon', 'microDecon')
    
    
    melanoma_rocs = melanoma_rocs.loc[melanoma_rocs.Dataset.str.contains('all_at_once')==False]
    
    pal_names=['No decontamination', 'Restrictive', 'Decontam', 'Decontam (LB)',
               'Custom', 'microDecon', 'SCRuB' ]
    
    
    #hue_order=['No decontamination', 'Restrictive', 'Decontam: ', '(LB)',
    #           'Custom', 'microDecon', 'SCRuB' ]
    hue_order=['No decontamination', 'Restrictive', 'Decontam (', '(LB)',
               'Custom', 'microDecon', 'SCRuB' ]

    melanoma_rocs['Dataset']=melanoma_rocs.Dataset.str.replace('Raw', 'No decontamination')
    melanoma_rocs['Dataset'] = melanoma_rocs.Dataset.str.replace(': ', ' (auROC = ') + ')'
    hue_order = [ melanoma_rocs.loc[ melanoma_rocs.Dataset.str.contains(a, regex=False) ].Dataset.values[0]
                    for a in hue_order ]
    
    
    
    melanoma_rocs['round'] = ( (melanoma_rocs.FPR==0) * (melanoma_rocs.TPR==0) ).cumsum()
    
    tmp=melanoma_rocs.copy()
    
    melanoma_rocs['key']=0
    ss=melanoma_rocs.merge( pd.DataFrame( {'full_FPR':np.sort( melanoma_rocs.FPR.unique() ), 
                                 'key':0} ), 
                     on='key')

    melanoma_rocs=ss.loc[ss.FPR<=ss.full_FPR]\
        .groupby(['full_FPR', 'Dataset',  'round'])['TPR'].max().reset_index()
    
    melanoma_rocs.columns=['FPR', 'Dataset', 'round', 'TPR']

    melanoma_rocs=pd.concat( [ melanoma_rocs.loc[melanoma_rocs.Dataset.str.contains(a, 
                                                                               regex=False)]
                              if 'Restrictive' not in a
                                else tmp[list(melanoma_rocs.columns
                                             )].loc[tmp.Dataset.str.contains('Restrictive')]
                              for a in hue_order 
                             ], axis=0).reset_index(drop=True)
    
    melanoma_rocs['Dataset']=melanoma_rocs.Dataset.str.replace('Raw', 'No decontamination')
    pal_names=['No decontamination', 'Restrictive', 'Decontam', 'Decontam (LB)', 'Custom','microDecon', 'SCRuB']
    
    
    plt.rcParams["font.family"] = "Calibri"
    plt.figure(figsize=(8,8))
    ax=sns.lineplot(x='FPR', y='TPR', hue='Dataset', 
                data=melanoma_rocs, 
                palette={a:global_palette[pal_names[i]] for i,a in enumerate(
                                        melanoma_rocs.Dataset.unique()) }, 
                hue_order=hue_order,
                 err_style="band", 
                 ci=95,
                 linewidth=5
                )
    sns.lineplot([0,1], [0,1], color = 'black')
    plt.legend(loc='lower right', prop={'size': 20})

    legend = plt.legend(loc='lower right', prop={'size': 20})
    plt.setp(legend.get_texts(), color='k')

    
    if hide_axes:
        plt.xlabel(None)
        plt.ylabel(None)
        plt.legend(loc='lower right', prop={'size': 20})
        plt.xticks(np.linspace(0,1,2), labels=[]*2)
        plt.yticks(np.linspace(0,1,2), labels=[]*2)
    else:
        plt.xticks(np.linspace(0,1,2))
        plt.yticks(np.linspace(0,1,2))

    plt.ylim(0,1)
    plt.xlim(0,1)


    plt.savefig(out_path, 
                bbox_inches='tight', 
                dpi=900, 
                format='pdf'
               )
    
    return(None)


def plot_plasma_shannons(diversities, out_path='../results/Plots/fig_3/F3_plasma_shannons.pdf', 
                            hide_axes=False):
    plt.subplots(1, 
                 dpi=500, 
                figsize=(2.5,8))
    
    ds=['Raw', 'Restrictive', 'Decontam', 'microDecon', 'SCRuB']

    ax=sns.boxplot( x='Dataset', y='Shannon', data=diversities, 
                   order = ds,
               hue_order = ds,
                width=.9,
                   palette=global_palette, 
                   fliersize=0
               )
    sns.swarmplot(x='Dataset', y='Shannon', data=diversities, 
                   order = ds,
                 color='black', size=2.5
                 )

    
    tmp=diversities.copy()
    
            # print vs restrictive significance
    q=diversities[['Dataset']].drop_duplicates()
    q['sig_v_restrictive'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.Dataset==row.Dataset)]['Shannon'].values, 
        tmp.loc[ (tmp.Dataset=='Restrictive') ]['Shannon'].values, 
        alternative='greater').pvalue
                        if row.Dataset!='Restrictive' else 1][0],
           axis=1)
    
    print(q)
    
        # plot vs raw significance asterisks
    q=diversities[['Dataset']].drop_duplicates()
    q['sig_v_raw'] = q.apply(lambda row: [ wilcoxon( 
                tmp.loc[ (tmp.Dataset==row.Dataset)]['Shannon'].values, 
        tmp.loc[ (tmp.Dataset=='Raw') ]['Shannon'].values, 
        alternative='less').pvalue
                        if row.Dataset!='Raw' else 1][0],
           axis=1)
    
    print(q)

    q['y_val']=-0.25
    q['is_sig'] =  q.sig_v_raw < 1e-3

    sns.swarmplot( x='Dataset', y='y_val', 
                    data = q.loc[q.is_sig], 
                   order = ds,
               hue_order = ds,
                   marker='+',
                  size=25/2, 
                    ax=ax,
                    palette=global_palette, 
                  dodge=True, 
                  color='black',
                  linewidth=3.5/2
                   )

    sns.swarmplot( x='Dataset', y='y_val',
                    data = q.loc[q.is_sig], 
                   order = ds,
               hue_order = ds,
                   marker='x',
                  size=17.5/2,
                    ax=ax,
                    palette=global_palette, 
                  dodge=True, 
                  color='black',
                  linewidth=3.5/2,
                    edgecolor="black"
                   )
    
    plt.xlabel(None)
    if hide_axes:
        ax.set(yticklabels=[])
        plt.ylabel(None)
        ax.axes.get_xaxis().set_visible(False)
    else:
        plt.xticks(rotation=90)
        plt.ylabel('Shannon Diversity')

    plt.savefig(out_path, 
               bbox_inches='tight', 
               dpi=900, 
               format='pdf')
    
    return(None)



def create_roc_dataset(aaa, path, idx=0):

    curve=roc_curve( ( aaa['obs'] == aaa.columns[-1] ).astype(int).values, 
                    aaa[aaa.columns[-1]].values, 
                   drop_intermediate=False)
    
    df = pd.DataFrame({'FPR':curve[0], 
                  'TPR':curve[1], 
                   'Boundary':curve[2],
                   'tp1':path.split('/')[-1].split('_')[0], 
                   'tp2':path.split('/')[-1].split('_')[1], 
                   'AUC':auc(curve[0], curve[1]), 
                   'rnd':idx })
        
    df.FPR= (df.FPR * 2 ).round(1) / 2
    return( df )


def format_rocs(summary):
    summary['key']=0
    ss=summary.merge( pd.DataFrame( {'full_FPR':np.sort( summary.FPR.unique() ), 
                             'key':0} ), 
                 on='key')

    ss=ss.loc[ss.FPR<=ss.full_FPR]\
        .groupby(['full_FPR', 'tp1', 'tp2', 'AUC', 'rnd'])['TPR'].max().reset_index()

    ss=ss[['full_FPR', 'TPR', 'tp1', 'tp2', 'AUC', 'rnd']]
    ss.columns=['FPR', 'TPR', 'tp1', 'tp2', 'AUC', 'rnd']
    return(ss)

def plot_plasma_melanoma_roc(path_to_files, out_path='../results/Plots/fig_3/Plasma_SKCM_v_control.pdf', 
                            hide_axes=False):
    
    subdirs = {'scrubbed_preds/':'SCRuB',
               'restrictive_preds/':"Restrictive",
                'dec_preds/':'Decontam',
                'dec_standard/':'Decontam default', 
               'dec_lb/': 'Decontam (LB)',
               'raw_preds/':'Raw',
               'microdecon_preds/':'microDecon'}
    
    subdirs= { os.path.join(path_to_files + 'trial_0/', a):b for a,b in subdirs.items() }

    full_results={}
    task='SKCM_Control'
    np.random.seed(1)
    full_results[task] ={ b:pd.concat([create_roc_dataset(pd.read_csv(path_to_files + 'trial_{}/'.format(i) + decontam_process.split('/')[-2] + '/' + task + '__Preds.csv'
                               ),
                                                  path_to_files + 'trial_{}/'.format(i) + \
                                                          decontam_process.split('/')[-2] + '/' + task + '__Preds.csv', 
                                                  idx=i)
                                            for i in range(10) ])
                            for decontam_process,b in subdirs.items() }

    for a in full_results['SKCM_Control']:
        if a!='SCRuB':
            print(a)
            print(wilcoxon(full_results['SKCM_Control'][a].groupby('rnd')['AUC'].median().values,
                           full_results['SKCM_Control']['SCRuB'].groupby('rnd')['AUC'].median().values
                          ) )

    for a in full_results:
        for b in full_results[a]:
            full_results[a][b]=format_rocs(full_results[a][b])
            full_results[a][b].columns=['FPR', 'TPR', 'tp1', 'tp2', 'AUC', 'rnd']
            full_results[a][b]['Task']=a
            full_results[a][b]['Dataset']=b
            full_results[a][b]['AUC']=round( full_results[a][b].groupby('rnd')['AUC'].median().median(), 2)
            


    all_roc_plots=pd.concat( [ pd.concat( [full_results[a][b] for b in full_results[a]] ) 
                                 for a in full_results ] ).reset_index(drop=True)

    tmp=all_roc_plots.loc[ (all_roc_plots.Task == 'SKCM_Control')&
                       (all_roc_plots.Dataset.isin(['Raw', 'Decontam', 'SCRuB', 'microDecon', 'Restrictive']))]
    
    tmp['Dataset'] = tmp.Dataset.str.replace('Raw', 'No decontamination')
    
    plt.rcParams["font.family"] = "calibri"
    plt.figure(figsize=(8,8))
    ax=sns.lineplot( "FPR", "TPR", hue='Dataset',
                 data = tmp, 
                 linewidth=4, 
                    palette=global_palette, 
                    hue_order=[ 'No decontamination', 'Restrictive', 'Decontam', 'microDecon', 'SCRuB']
               )

    sns.lineplot([0,1], [0,1], color = 'black')

    ax.legend( 
        labels=[q+')' for q in [' (auROC = '.join( b.astype(str) ) for b in (tmp[['Dataset', 'AUC']].drop_duplicates().values[[3,1,2,4,0], :]
                )] ], 
    prop={'size': 20},
    loc='lower right')
    
    plt.ylim(0,1)
    plt.xlim(0,1)

    if hide_axes:
        plt.xlabel(None)
        plt.ylabel(None)
        plt.xticks(np.linspace(0,1,2), labels=[]*2)
        plt.yticks(np.linspace(0,1,2), labels=[]*2)
        
    else:
        plt.xticks(np.linspace(0,1,2))
        plt.yticks(np.linspace(0,1,2))


        


    
    plt.savefig(out_path, 
                bbox_inches='tight', 
                dpi=900, 
                format='pdf'
               )
        
    return(None)


def plot_plasma_roc_supps(path_to_files, out_path='../results/Supplementary_figures/Plasma_all/', 
                            hide_axes=False):
    subdirs = {'raw_preds/':'Raw',
           'restrictive_preds/':"Restrictive",
           'dec_preds/':'decontam (Poore et. al)',
            'dec_standard/':'decontam', 
           'dec_lb/': 'decontam (LB)',
            'microdecon_preds/':'microDecon',   
           'scrubbed_preds/':'SCRuB'}
    if os.path.exists('../results/Supplementary_figures/Plasma_all')==False:
        os.makedirs('../results/Supplementary_figures/Plasma_all')
    
    tasks= [a.split('/')[-1].split('__')[0] for a in glob.glob(path_to_files + '/trial_0/dec_lb/*__Preds.csv') ]

    
    full_results={}
    for task in tasks:
        np.random.seed(1)
        full_results[task] = { b:pd.concat([create_roc_dataset(pd.read_csv(path_to_files + 'trial_{}/'.format(i) + decontam_process + task + '__Preds.csv'
                               ), #.sample(frac=.8, random_state=i), 
                                                  path_to_files + 'trial_{}/'.format(i) + decontam_process + task + '__Preds.csv', 
                                                  idx=i)
                                            for i in range(10) ])
                                    for decontam_process,b in subdirs.items() }



    all_aucs = {}
    for a in full_results:
        all_aucs[a]={}
        for b in full_results[a]:
            full_results[a][b]=format_rocs(full_results[a][b])
            full_results[a][b].columns=['FPR', 'TPR','tp1', 'tp2', 'AUC', 'rnd']
            full_results[a][b]['Task']=a
            full_results[a][b]['Dataset']=b
            full_results[a][b]['AUC']=round( full_results[a][b].groupby('rnd')['AUC'].median().median(), 2)
            all_aucs[a][b]=b + ': auROC='+str(full_results[a][b].AUC.values[0])


    all_roc_plots=pd.concat( [ pd.concat( [full_results[a][b] for b in full_results[a]] ) 
                                 for a in full_results ] ).reset_index(drop=True)

    plt.rcParams["font.family"] = "calibri"
    for task in all_roc_plots.Task.unique():

        tmp=all_roc_plots.loc[all_roc_plots.Task==task]
        tmp['Dataset'] = tmp.Dataset.str.replace('Raw', 'No decontamination')
        plt.figure(figsize=(8,8))
        ax=sns.lineplot( "FPR", "TPR", hue='Dataset',
                     data = tmp, 
                     linewidth=4, 
                      palette=global_palette
                   )

        sns.lineplot([0,1], [0,1], color = 'black', ax=ax)                            

        ax.legend( 
        labels=[ q+')' for q in [' (auROC = '.join( b.astype(str) ) for b in (tmp[['Dataset', 'AUC']].drop_duplicates().values
                )] ], 
                prop={'size': 17.5}, #20},
                loc='lower right'
                )


        if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.xticks(np.linspace(0,1,2), labels=[]*2)
            plt.yticks(np.linspace(0,1,2), labels=[]*2)
        else:
            plt.xticks(np.linspace(0,1,2))
            plt.yticks(np.linspace(0,1,2))
            
        plt.ylim(0, 1)
        plt.xlim(0, 1)
        
        plt.savefig(out_path + task + '_ROC.pdf', 
                    format='pdf', 
                    dpi=900, 
                    bbox_inches='tight'
                   )
        
        
def plot_tumor_upset():
    df=pd.read_csv('../results/data/Tumor/upset_setup.csv', index_col=0)
    df.columns = df.columns.str.replace('..LB.', ' (LB)')
    upc=upsetplot.from_contents( {a:np.where( df[a] )[0] for a in
        ['Restrictive', 'Decontam',  'Decontam (LB)', 'Custom', 'microDecon', 'SCRuB'] } )

    fig = plt.figure(figsize=(15, 10))
    qq = upsetplot.UpSet(upc , 
             subset_size='count', element_size=None
              )
    qq.plot(fig=fig)
    plt.savefig('../results/Plots/fig_3/upset_plot.pdf', 
               bbox_inches='tight', 
               dpi=900)
    return(None)


    