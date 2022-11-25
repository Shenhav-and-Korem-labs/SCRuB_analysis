import os
from fig1_functions import *
from fig2_functions import *
from fig3_functions import *
from fig4_functions import *
from plate_visualizations import *
import gc
plt.rcParams["font.family"] = "Calibri"

def store_f1_plots(plot_plates=True, hide_axes=False):
    if plot_plates:
        sns.set(rc={'axes.facecolor':'888888', 'figure.facecolor':'888888'})

        simulate_plate(colors=['blue', 'blue'], 
                      color_negs=False)

        with plt.rc_context({'image.composite_image': False}):
            plt.savefig('../results/Plots/fig_1/F1_A_no_contam.pdf', 
                        dpi=500,
                       bbox_inches='tight', 
                       format='pdf')

        simulate_plate(colors=['blue', 'red'], 
                      color_negs=True)

        with plt.rc_context({'image.composite_image': False}):
            plt.savefig('../results/Plots/fig_1/F1_A_with_contam.pdf', 
                        dpi=500,
                       bbox_inches='tight', 
                       format='pdf')
    
    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'}, font_scale=2)
    
    data=pd.read_csv('../results/data/Fig_1/F1_Simulation_Results.csv')
    plot_f1c(data, hide_axes=hide_axes)
    data=pd.read_csv('../results/data/Fig_1/F1_well_leakage_Results.csv')
    plot_f1d(data, hide_axes=hide_axes)
    make_no_contam_plots(hide_axes=hide_axes)
    # noise_level_plot(hide_axes=hide_axes)
    return(None)

def pull_f2_data(data_dir):
    datasets = {}
    estimated_contaminants = {}

    for path in os.listdir(data_dir):
        if path[-4:]=='.csv':
            if '_cont' in path:
                print(path)
                estimated_contaminants[path.split('.')[0]]=pd.read_csv(data_dir + path, index_col=0)
            else:
                a=path.split('.')[0]
                datasets[a]=pd.read_csv(data_dir + path, index_col=0)
                datasets[a]=datasets[a]/datasets[a].sum(axis=1)[:, np.newaxis]
                datasets[a][ np.isnan(datasets[a]) ]=0


    datasets['Raw']=datasets['Raw'].loc[datasets['Raw'].index.str.contains('extractNTC')==False]

    for a in datasets:
        datasets[a].index=datasets['Raw'].index
        datasets[a].columns=datasets['Raw'].columns
        
    return(datasets, estimated_contaminants)

def store_f2_plots(data_dir='results/Fig_2/', hide_axes=False):
    
    if type(data_dir)==list:
        split_data = [ pull_f2_data(a) for a in data_dir ]
        datasets = {a:pd.concat( [ split_data[i][0][a] for i in range(len(split_data) ) ]) for a in split_data[0][0] }
        
        estimated_contaminants = {a:pd.concat( [ split_data[i][1][a] for\
                                                            i in range(len(split_data) ) ]) for a in split_data[0][1] }
    else:
        datasets,estimated_contaminants=pull_f2_data(data_dir)
        split_data=None
        
        
    sns.set(rc={'axes.facecolor':'888888', 'figure.facecolor':'888888'})
    
    plot_f2_a(data_dir= [ data_dir[0] if type(data_dir)==list else data_dir][0])
    
    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'}, font_scale=2)
    
    
    out_path='../results/Plots/fig_2'
    
    make_estimated_cont_plot(estimated_contaminants, datasets, path='../results/Plots/fig_2/F2_E.pdf', data_dir=data_dir, 
                            split_data=split_data, 
                            hide_axes=hide_axes)
    make_f2_swarmplot([ datasets if split_data is None else split_data][0], path=out_path +'/F2_F.pdf', 
                     hide_axes=hide_axes)
    
    
    make_supplemental_jsd_swarmplot([ datasets if split_data is None else split_data][0],
                                      path='../results/Supplementary_figures/Well_to_well/Fig2_JSD.pdf', 
                                   hide_axes=hide_axes)
    

    
    return(None)



def store_f3_plots(run_plates=False, hide_axes=False):
    # plot plates
    if run_plates:
        sns.set(rc={'axes.facecolor':'888888', 'figure.facecolor':'888888'})
        for i,colors in enumerate([['blue','yellow','green', 'red', 'red', 'red'], 
                                   ['blue','yellow','green', 'green'], 
                                   ['blue','yellow'], 
                                   ['blue', 'blue']
                                   ]):
                simulate_plate(colors=colors, 
                                                  color_negs=i!=3)
                with plt.rc_context({'image.composite_image': False}):
                    plt.savefig('../results/Plots/fig_3/f3_a' + str(i+1) + '.png', 
                                dpi=500,
                               )
            
    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'}, font_scale=2)
   
    # plot shannon diversities
    shannon_divs = pd.read_csv('../results/data/Tumor/Shannon_diversities.csv', index_col=0)
    plot_cancer_shannon_divs(shannon_divs,
                            hide_axes=hide_axes)
    
    plasma_shannon_divs = pd.read_csv('../results/data/Fig3_plasma/plasma_shannon_diversities.csv', index_col=0)
    plot_plasma_shannons(plasma_shannon_divs,
                        hide_axes=hide_axes)

    
#     plot SKCM vs Control rocs
    plot_plasma_melanoma_roc('../results/data/Fig3_plasma/',
                            hide_axes=hide_axes)
    plot_plasma_roc_supps('../results/data/Fig3_plasma/',
                         hide_axes=hide_axes)
    
    #plot ICI rocs
    melanoma_rocs=pd.read_csv('../results/data/Tumor/Melanoma_Test_ROCs.csv', index_col=0)
    plot_ICI_rocs(melanoma_rocs, 
                 hide_axes=hide_axes)
    
    # plot_tumor_upset()
    
        
    return(None)


def store_f4_plots(hide_axes=False):
    
    sim_results=pd.read_csv('../results/data/Fig_4/Fig_4_Simulations.csv', index_col=0)
    sim_results=sim_results.loc[sim_results.names=='Spatial SCRUB']
    
    gut_sitch=sim_results.loc[
                        (sim_results.sample_type=='high')
                            ]

    skin_sitch=sim_results.loc[
                        (sim_results.sample_type=='med')
                            ]

    vag_sitch=sim_results.loc[
                        (sim_results.sample_type=='low')
                            ]
    
    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'}, font_scale=2)
    
    print('skin')
    make_heatmap(skin_sitch.groupby(['plate_format', 'n_controls', 'round'])['jsds'].median().reset_index()\
                 .sort_values(['plate_format', 'n_controls', 'jsds'], ascending=False ).reset_index(drop=True),
                 sig=.05,
                 out_path='../results/Plots/fig_4/skin_heatmap.pdf', 
                hide_axes=hide_axes)

    print('vaginal')
    make_heatmap(vag_sitch.groupby(['plate_format', 'n_controls', 'round'])['jsds'].median().reset_index()\
                    .sort_values(['plate_format', 'n_controls', 'jsds'], ascending=False ).reset_index(drop=True), 
                 out_path='../results/Plots/fig_4/vag_heatmap.pdf',
                 sig=.05,
                 hide_axes=hide_axes
                )

    print('gut')
    make_heatmap(gut_sitch.groupby(['plate_format', 'n_controls', 'round'])['jsds'].median().reset_index()\
                    .sort_values(['plate_format', 'n_controls', 'jsds'], ascending=False ).reset_index(drop=True), 
                 out_path='../results/Plots/fig_4/gut_heatmap.pdf',
                 sig=.05,
                 hide_axes=hide_axes
                )
        
        
    make_f4_sim_supplemental( sim_results.groupby(['plate_format', 'n_controls', 'sample_type', 'round'])\
                             ['jsds'].median().reset_index()\
                                .sort_values(['plate_format', 'n_controls', 'jsds'], 
                                            ascending=False ).reset_index(drop=True),
                             hide_axes=hide_axes)
        
        
    make_s7_plots()
    return(None)



def main(hide_axes=False):

    [os.makedirs('../results/Plots/fig_' + str(i+1)) for i in range(4) if os.path.exists('../results/Plots/fig_' + str(i+1))==False]
    store_f1_plots(plot_plates=False, hide_axes=hide_axes)
    gc.collect()
    store_f2_plots(data_dir = ['../results/data/Fig_2/', '../results/data/Fig_2_alternate/'], 
                   hide_axes=hide_axes)
    gc.collect()
    store_f3_plots(hide_axes=hide_axes)
    gc.collect()
    store_f4_plots(hide_axes=hide_axes)
    gc.collect()
    print('Successful run completed!')
    return(None)


if __name__ == "__main__":
    main(hide_axes=False) # set to True to create exact plots from figures
else:
    print('name is not main')





