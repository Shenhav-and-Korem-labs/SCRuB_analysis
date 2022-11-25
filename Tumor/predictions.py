import numpy as np
import pandas as pd
from itertools import product
import xgboost

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

class PCA_predictor():
    def __init__(self, model, pca_components=20):
        self.model=model
        self.sc = StandardScaler()
        self.pca = PCA(n_components=pca_components)
        
        
    def fit(self, x_, y):
        x = self.sc.fit_transform(x_)
        x = self.pca.fit_transform(x_)
        self.model.fit(x,y)
    
    def predict_proba(self, x_):
        x = self.sc.transform(x_)
        x = self.pca.transform(x_)
        return( self.model.predict_proba(x) )
    
 
    
def format_mat(x_, log_transform, clr_transform):
    x =  (x_.T / x_.sum(axis = 1) ).T
    x[np.isnan(x)]=0
    if log_transform:
        x = np.log10(1+x)
    elif clr_transform:
        x = clr(1+x)  
    return(x)

def preds_one_fold_per_patient(df_,
                               model = RandomForestClassifier(), 
                               test_pts=None,
                               log_transform = False, 
                               clr_transform=False,
                               id_col = 'patient_id',
                               resp_col = 'response',
                               pos_val = 1, 
                               neg_val = 0, 
                               sample = 1, 
                               shuffle=True,
                               n_otu = 100,
                               seed=None):
    df = df_.copy()
    
    df = df.loc[df[id_col].isin(test_pts) == False]
    
    np.random.seed(seed)
    
    if sample!=1:
        df = df.sample(frac=sample)
        
    if shuffle:
        df = df.sample(frac=1)
 
    unique_pts = df[id_col].unique()
    pt_splits = [ (format_mat( df.loc[df[id_col]==a].drop([id_col, resp_col], axis=1).values,
                      log_transform, 
                      clr_transform), 
           df.loc[df[id_col]==a][resp_col].values ) for a in unique_pts ]

    all_trues=[]
    all_preds=[]
    
    for i,pt in enumerate(unique_pts):

        train = np.vstack( [x[0] for j,x in enumerate(pt_splits) if j!=i] )

        feats= list( np.argsort( (train > 0 ).sum(axis=0) )[-n_otu:] )
        
        train=train[:, feats]
        test = pt_splits[i][0][:, feats] 
        train_labels = np.hstack([x[1] for j,x in enumerate(pt_splits) if j!=i])

        test_labels = pt_splits[i][1]

        
        rf = model
        rf.fit(train, train_labels)

        all_trues.append( test_labels )
        all_preds.append( rf.predict_proba(test) )
    
    y = np.hstack(all_trues)
    y_hat = np.vstack(all_preds)[:, 1]
    
    return(roc_curve(y, y_hat), (rf, feats) )


def preds_for_kfold_tuning(df_,
                           test_pts=None,
                           model = RandomForestClassifier(), 
                           log_transform = False, 
                           clr_transform=False,
                           id_col = 'patient_id',
                           resp_col = 'response',
                           pos_val = 1, 
                           neg_val = 0, 
                           sample = 1, 
                           shuffle=True,
                           n_otu = 10,
                           n_folds=5,
                           seed=None):
    df = df_.copy()
    
    df=df.loc[df[id_col].isin(test_pts)==False]
    
    np.random.seed(seed)
    
    if sample!=1:
        df = df.sample(frac=sample)
        
    if shuffle:
        df = df.sample(frac=1)

    t1, t2 = format_mat( df.loc[ (df[resp_col] == pos_val)].drop([id_col, resp_col], axis=1).values, 
                                                                  log_transform, clr_transform ), \
             format_mat( df.loc[ (df[resp_col] == neg_val)].drop([id_col, resp_col], axis=1).values, 
                                                                  log_transform, clr_transform )
        
        
    t1_folds = np.linspace(0, t1.shape[0], n_folds+1, dtype = np.int64)
    t1_folds = [ (t1_folds[i], t1_folds[i+1]) for i in range(n_folds) ]

    t2_folds = np.linspace(0, t2.shape[0], n_folds+1, dtype = np.int64)
    t2_folds = [ (t2_folds[i], t2_folds[i+1]) for i in range(n_folds) ]
    
    
    all_trues = list()
    all_preds = list()
    
    t1_split = [t1[a[0]:a[1]] for a in t1_folds]
    t2_split = [t2[a[0]:a[1]] for a in t2_folds]
    
    for i in range(n_folds):

        t1_train = np.vstack( [x for j,x in enumerate(t1_split) if j!=i] )
        t2_train = np.vstack( [x for j,x in enumerate(t2_split) if j!=i] )

        feats= list( np.argsort( (np.vstack((t1_train, t2_train)) > 0 ).sum(axis=0) )[-n_otu:] )
        t1_train, t2_train = t1_train[:, feats], t2_train[:, feats]
        
        t1_test = t1_split[i][:, feats] 
        t2_test = t2_split[i][:, feats] 

        train_labels = np.hstack((np.zeros(t1_train.shape[0]),
                   np.ones(t2_train.shape[0]) ) )

        test_labels = np.hstack((np.zeros(t1_test.shape[0]),
                   np.ones(t2_test.shape[0]) ) )

        
        full_train = np.vstack((t1_train, t2_train)) 
        full_test = np.vstack((t1_test, t2_test)) 
        
        rf = model
        rf.fit(full_train, train_labels)

        all_trues.append( test_labels )
        all_preds.append( rf.predict_proba(full_test) )
    
    y = np.hstack(all_trues)
    y_hat = np.vstack(all_preds)[:, 1]
    
    return(roc_curve(y, y_hat), (rf, feats) )
    

def preds_for_final_test(df_,
                         test_pts=None,
                       model = RandomForestClassifier(), 
                       log_transform = False, 
                       clr_transform=False,
                       id_col = 'patient_id',
                       resp_col = 'response',
                       pos_val = 1, 
                       neg_val = 0, 
                       sample = 1, 
                       shuffle=True,
                       n_otu = 10,
                       seed=None):
    df = df_.copy()
    
            
    np.random.seed(seed)
    
    if sample!=1:
        df = df.sample(frac=sample)
        
    if shuffle:
        df = df.sample(frac=1)
 
    train=df.loc[df[id_col].isin(test_pts)==False]
    y_train=train[resp_col].values
    test=df.loc[df[id_col].isin(test_pts)].groupby(id_col).sample(n=1)
    y_test = test[resp_col].values
    
    train = format_mat(train.drop([id_col, resp_col], axis=1).values, log_transform, clr_transform)
    test = format_mat(test.drop([id_col, resp_col], axis=1).values, log_transform, clr_transform)
    
    
    feats= list( np.argsort( (train > 0 ).sum(axis=0) )[-n_otu:] )
    
    rf = model
    rf.fit(train[:, feats], y_train)
    
    all_trues=[]
    all_preds=[]
    
    all_trues.append( y_test )
    all_preds.append( rf.predict_proba(test[:, feats]) )
    
    y = np.hstack(all_trues)
    y_hat = np.vstack(all_preds)[:, 1]
    
    return(roc_curve(y_true = y, y_score= y_hat), (rf, feats) )


def tune_parameters(dataset,
                    test_pts, 
                    n_iters=15, 
                    predict_func=preds_for_kfold_tuning, 
                    model=xgboost.XGBClassifier):    
    n_otus = [5, 10, 25, 100]#, 50, 250, 2000]
    pca_components =[None] #[10, 50, None]
    transforms=[None]#, 'log']
    xgb_params = {'max_depth':[3,10], 
                 'n_estimators':[100, 500],#, 250], 
                 'learning_rate':[.1, 1], 
                  'subsample':[.5, .75]
#                   'colsample_bytree':[.5,1],
#                    'colsample_bylevel':[.5, 1]
                 }



    xgb_param_list=[ {'max_depth':a, 
                      'n_estimators':b, 
                      'learning_rate':c, 
                      'subsample':f,
#                       'colsample_bytree':d,
#                     'colsample_bylevel':e,
                  'objective':'binary:logistic', 
                  'eval_metric':'error'} for a,b,c,f in (product(*xgb_params.values())) ]
                    #,d,e in (product(*xgb_params.values())) ]

    best_xgb_auc = 0
    best_xgb_setup = 0
    
    for trans in transforms:
        print('transformation: ', trans)
        for n_o in n_otus:
            print('n_otus: ', n_o)
            for pcas in pca_components:
                print('pcas: ', pcas)
                for xgb_pars in xgb_param_list:
                    print(xgb_pars)
                    if pcas is None:
                        pred_out=[]
                        for seed in range(n_iters):
                            pred_out.append( predict_func(dataset,
                                                          test_pts=test_pts,
                                                          model = model(**xgb_pars),
                                                           shuffle = True, 
                                                           log_transform = trans=='log', 
                                                           resp_col = 'response',
                                                           pos_val = True, 
                                                           neg_val = False, 
                                                           sample = .95,
                                                           n_otu = n_o,
                                                           seed=seed) 
                                           )
                        auc_val = np.array( [auc(a[0][0], a[0][1]) for a in pred_out] ).mean()
                        if auc_val > best_xgb_auc:
                            best_xgb_auc = auc_val
                            best_xgb_setup = {'transformation':trans, 
                                                'n_otus':n_o, 
                                                'pcas':pcas, 
                                                 'xgb_params':xgb_pars}
            
            print(best_xgb_auc)
            print(best_xgb_setup)
    return(best_xgb_auc, best_xgb_setup)



def run_final_tests(datasets, parameters, test_pts):
    test_results={a:[ preds_for_final_test(b, n_otu=parameters[a][1]['n_otus'],
                                           test_pts=test_pts,
                                           seed=i,
                                           sample=1, 
                                           model=xgboost.XGBClassifier(**parameters[a][1]['xgb_params']), 
                                           log_transform=parameters[a][1]['transformation']=='log'
                                          )[0] for i in range(15) ]
                  for a,b in datasets.items() }
    return(test_results)