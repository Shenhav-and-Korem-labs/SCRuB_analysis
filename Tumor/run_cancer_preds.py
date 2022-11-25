import pandas as pd
import numpy as np
import predictions
import pickle
from sklearn.metrics import auc
import warnings
import os
warnings.filterwarnings('ignore')



def pull_cancer_datasets():
    metadata = pd.read_csv('../data/Fig3/full_cancer_metadata.csv')
    metadata['new_SEQ_NAMES'] = metadata.species
    metadata = metadata.drop('species', axis = 1)
    metadata.new_SEQ_NAMES = metadata.new_SEQ_NAMES.str.replace(r'/|-| ', '.')
    metadata.columns = metadata.columns.str.replace(' ', '_')

    scrub = pd.read_csv('../results/data/Tumor/SCRUB_result.csv', index_col = 0)
    scrub = metadata.merge(scrub, on = 'new_SEQ_NAMES')
    
    restrictive = pd.read_csv('../results/data/Tumor/global_fully_restrictive_result.csv', index_col = 0)
    restrictive = metadata.merge(restrictive, on = 'new_SEQ_NAMES')
    
    
    decontam = pd.read_csv('../results/data/Tumor/decontam_result.csv', index_col = 0)
    decontam = metadata.merge(decontam, on = 'new_SEQ_NAMES')
    
    decontam_lb = pd.read_csv('../results/data/Tumor/decontam_low_biomass_result.csv', index_col = 0)
    decontam_lb = metadata.merge(decontam_lb, on = 'new_SEQ_NAMES')
    
    microdecon = pd.read_csv('../results/data/Tumor/microdecon_result.csv', index_col = 0)
    microdecon = metadata.merge(microdecon, on = 'new_SEQ_NAMES')
    

    df_freqs = pd.read_csv('../data/Fig3/df_freqs.csv', index_col = 0).fillna(0)
    df_freq_phylo = pd.read_csv('../data/Fig3/df_filter_tracker.csv', index_col = 0).loc[:,'domain':'species'].merge(
            df_freqs, how = 'inner', left_index=True, right_index = True)

    theirs = df_freq_phylo.loc[df_freq_phylo.species.isna() == False].iloc[:, 7:].T
    theirs.index = theirs.index.str.replace(r'/|-| ', '.')

    theirs_full = scrub.iloc[:, :30].merge(theirs, 
                                          how = 'inner', 
                                          left_on = 'new_SEQ_NAMES', 
                                         right_index = True)

    df_freqs_predecontam = pd.read_csv('../data/Fig3/df_freqs_before_global_contaminants_filter.csv', \
                                       index_col = 0).fillna(0).iloc[:, 7:].T

    df_freqs_predecontam.index = df_freqs_predecontam.index.str.replace(r'/|-| ', '.')

    contaminated = scrub.iloc[:, :30].merge(df_freqs_predecontam, 
                                                  how = 'inner', 
                                                  left_on = 'new_SEQ_NAMES', 
                                                 right_index = True)
    
    
    resp_col = 'Melanoma_(T)_-_Response_to_immune_checkpoint_inhibition_therapy_(NR/R)'
    datasets_dict={'scrub':scrub, 
                   'original':theirs_full, 
                   'restrictive':restrictive,
                   'decontam':decontam, 
                   'decontam_lb':decontam_lb,
                   'microDecon':microdecon,
                   'raw':contaminated
                   }
    
    
    contam_map = pd.read_csv('../data/Fig3/df_freqs_before_global_contaminants_filter.csv', \
                                   index_col = 0).fillna(0).iloc[:, :7]

    original_map = df_freq_phylo.iloc[:,:7]
    original_taxa_info = original_map.reset_index(drop=True)
    raw_taxa_info = contam_map.reset_index(drop=True)


    def group_to_genus(df, 
                       original_taxa_info=original_taxa_info, 
                       raw_taxa_info=raw_taxa_info):
        
        if df.shape[1] < 11200:
            taxa = original_taxa_info
        else:
            taxa = raw_taxa_info

        merged_taxa = pd.concat([taxa, df.iloc[:, 30:].T.reset_index(drop=True)], axis = 1)

        df = pd.concat( [ df.iloc[:, :30], merged_taxa.groupby(['domain', 'phylum', 
                                                                'class', 'order', 
                                                'family', 'genus'])[merged_taxa.columns[7:]].sum().T ], axis=1 )
        return(df)
    
    
    datasets = {a:group_to_genus( b.loc[( b[resp_col].isna()==False )&
                        ( b['Melanoma_(T)_-_body_site'].str.lower().str.contains('lymph|ln') == False )
                       ] ).sort_values('Sample_ID_(WIS)') for a,b in datasets_dict.items()}

    for a,b in datasets.items():
        b['patient_id'] = b['Pat_ID#_(WIS)']
        b.index=b['Sample_ID_(WIS)']
        b['response']=b[resp_col] == 'R'
        b=b.iloc[:, 30:]
        datasets[a]=b.loc[:, b.sum(axis=0)>0]
        
        
    nki_patients=theirs_full.loc[theirs_full.Center_x == 'NKI']['Pat_ID#_(WIS)'].unique()

    return(datasets, nki_patients)

def run_prediction(datasets, 
                   test_pts='infer', 
                   rerun_tuning=False,
                   tuning_results=None,
                   seed=2021):
    
    np.random.seed(seed)
    if test_pts=='infer':
        q=datasets[ [a for a in datasets][0] ]
        test_pts= np.append( q.loc[q.response==False].sample(frac=.20).patient_id.unique(), 
                           q.loc[q.response==True].sample(frac=.20).patient_id.unique()
                         )
    if rerun_tuning:
        tuning_results = { a:predictions.tune_parameters(b, test_pts=test_pts) for a,b in datasets.items() }

        with open('../results/data/Tumor/tuned_params_all_decont_variations.pickle', 'wb') as handle:
            pickle.dump(tuning_results, handle)
        
    
    np.random.seed(seed)
    test_rocs = predictions.run_final_tests(datasets, tuning_results, test_pts)
    return(test_rocs)


def write_melanoma_rocs_nki_test_set(rerun_tuning=False):    

    with open('Tumor/tuned_params_all_decont_variations.pickle', 'rb') as file:
        tuning_results=pickle.load(file)
    file.close()
    
    # load data
    datasets, center_test_pts=pull_cancer_datasets()
    
    # run the predictions
    
    # format results into source data
    test_out=run_prediction(datasets, tuning_results=tuning_results, test_pts=center_test_pts, rerun_tuning=rerun_tuning) # set rerun_tuning if you want to run the k-fold hyperparameter tuning scheme
    test_aucs={ a:np.median( np.array( [ auc(c[0],c[1] )  for c in b] ) ).round(2) for a,b in test_out.items() }
    stacked_rocs = pd.DataFrame(
                        np.vstack( 
                             [np.vstack(( np.hstack([ c[0] for c in b]), 
                                          np.hstack([ c[1] for c in b]), 
                                          np.array( [a+': '+str(test_aucs[a])]*np.hstack([ c[0] for c in b]).shape[0] )
                                            )).T
                               for a,b in test_out.items()]
                                    ), 
                                    columns=['FPR', 'TPR', 'Dataset']
                                )

    stacked_rocs.Dataset=stacked_rocs.Dataset.str.capitalize()
    stacked_rocs.Dataset=stacked_rocs.Dataset.str.replace('Scrub', 'SCRUB').str.replace('Original', 'Custom')
    
    stacked_rocs.to_csv('../results/data/Tumor/Melanoma_Test_ROCs.csv')
    return(None)
    
    
if __name__=='__main__':
    write_melanoma_rocs_nki_test_set(rerun_tuning=False)
    
    
    
    
    
    
    
    
    
    

