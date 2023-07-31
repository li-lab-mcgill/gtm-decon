import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error

sns.set(rc={'figure.figsize':(10, 6), 'figure.dpi':350})

def compute_stats(x1, x2):
    '''
    Computes quantitative results

    Args:
        x1 (np.array) : ground truth proportions
        x2 (np.array) : deconvolved proportions (predictions)

    Returns:
        PCC : pearson R
        SCC : spearman R
        CE : cross entropy
        CEW : weighted cross entropy

    '''
    PCC = pearsonr(x1, x2)[0]
    SCC = spearmanr(x1, x2)[0]
    
    t = np.sum(x1)
    t2 = np.sum(x2)

    # normalize props to be same scale
    if abs(np.sum(x2) - 1) > 1e-10 and abs(np.sum(x2) > 1e-10):
        x2 = x2 / np.sum(x2)
    if abs(np.sum(x1) - 1) > 1e-10 and abs(np.sum(x1)) > 1e-10:
        x1 = x1 / np.sum(x1)
        
    
    CE = -np.sum(x1 * np.log(x2.astype('float') + 1e-15))
    CEW = -np.sum(x1 * np.log(x2.astype('float') + 1e-15)) / len(x2)
    

    RMSE = mean_squared_error(x2, x1, squared=False)

        
    MAE = mean_absolute_error(x2, x1)

    return PCC, SCC, CE, CEW, RMSE, MAE

def compile_stats(
    ground_truth,
    methods = [],
    names = [],
):
    '''
    Plots quantitative statistics boxplots. Expect input methods dataframes list to be in same order as names. 

    Args:
        ground_truth (pd.DataFrame) : mixing proportions (samples x celltype)
        methods (list[pd.DataFrame]) : predicted proportions
        names (list[str]) : names of methods (in same order)

    Returns:
        df (pd.DataFrame) : statistics dataframe

    '''


    # populate stats
    stats = []
    for i in ground_truth.index:
        for method_, name_ in zip(methods, names):
            if method_ is not None:
                s = compute_stats(
                    ground_truth.loc[i],
                    method_.loc[i][ground_truth.columns],
                )

                stats.append(
                    tuple(
                        [name_] + list(s)
                    )
                )
    
    df = pd.DataFrame(
        stats,
        columns = [
            "Method", "Pearson R", "Spearman R", "Cross Entropy", "Weighted Cross Entropy", "RMSE", "MAE"
        ]
    )


    return df

def compile_stats_celltypewise(
    ground_truth,
    methods = [],
    names = [],
    spearman = True,
):
    '''
    Computes celltype-specific Pearson R.

    Args:
        ground_truth (pd.DataFrame) : mixing proportions (samples x celltype)
        methods (list[pd.DataFrame]) : predicted proportions
        names (list[str]) : names of methods (in same order)

    Returns:
        df (pd.DataFrame) : statistics dataframe

    '''


    res = []
    for ct in ground_truth:
        ground_truth[ct]

        l = []
        for method in methods:
            if spearman:
                l.append(
                spearmanr(
                    ground_truth[ct],
                    method.loc[ground_truth.index][ct]
                )[0] # * ground_truth.div(ground_truth.sum(axis=1), 0)[ct].mean()
            )
            else:
                l.append(
                    pearsonr(
                        ground_truth[ct],
                        method.loc[ground_truth.index][ct]
                    )[0] # * ground_truth.div(ground_truth.sum(axis=1), 0)[ct].mean()
                )

        res.append(l)
    #     break

    df = pd.DataFrame(res, columns=names, index=ground_truth.columns).T
    return df

def plot_fig2a(path = "../data/fig2_data/"):
    files = sorted(os.listdir(path))

    exp_pairs = (
        ('HBC', 'GSE107011'),
        ('HBC', 'WB'),
        ('Lake', 'ROSMAP'),
        ('PBMC', 'GSE107011'),
        ('PBMC', 'SDY67'),
        ('SegerstolpeSC', 'SegerstolpeBulk')
    )


    COMP_STATS = {f"{ref}_{bulk}" : None for ref, bulk in exp_pairs}
    for ref, bulk in exp_pairs:
        
        methods = []
        names = []
        for file in files:
            if ref in file and bulk in file:
                
                
                if 'proportions' in file:
                    ground_truth = pd.read_csv(os.path.join(path, file), index_col=0)
                else:
                    methods.append(pd.read_csv(os.path.join(path, file), index_col=0))
                    names.append(file.split('_')[2:])
        names = [n[0] if len(n) == 2 else '_'.join(n[:2]) for n in names]
        
        COMP_STATS[f"{ref}_{bulk}"] = compile_stats(
            ground_truth,
            methods,
            names
        )
        
    rename = {'BSEQSC': 'BSEQ-sc',
        'CIBER': 'CIBERSORTx',
        'MUSIC': 'MuSiC',
        'GTM_ALL': 'GTM-ALL',
        'GTM_PP': 'GTM-PP',
        'GTM_HVG': 'GTM-HVG',
        'BISQUE': 'BISQUE',
        'BAYESPRISM': 'BayesPrism',}

    for k in COMP_STATS:
        COMP_STATS[k]['Method'] = COMP_STATS[k]['Method'].replace(rename)

        
    plot_order = [
        'PBMC_SDY67',
        'PBMC_GSE107011',
        'HBC_WB',
        'Lake_ROSMAP',
        'SegerstolpeSC_SegerstolpeBulk'
    ]

    plot_order = [COMP_STATS[p] for p in plot_order]


    titles = ['Ref: PBMC\nBulk: PBMC-1', 
            'Ref: PBMC\nBulk: PBMC-2', 
            'Ref: Blood Cells\nBulk: Whole Blood', 
            'Ref: Frontal Cortex\nBulk: Prefrontal Cortex', 
            'Ref Pancreas\nBulk: Pancreas']


    sns.set(rc={'figure.figsize':(10, 6), 'figure.dpi':350})
    fig, axes = plt.subplots(2, 5, sharex=True)
    order = ['GTM-ALL', 'GTM-PP', 'GTM-HVG', 'BSEQ-sc', 'CIBERSORTx', 'MuSiC', 'BayesPrism', 'BISQUE']


    for i in range(5):
        for j in range(2):
            
            d = plot_order[i]
            d.Method = d.Method.astype('category')
            d.Method = d.Method.cat.set_categories(order)
            d.sort_values('Method')
            
            if j == 0:
                y = 'Spearman R'
            else:
                y = 'RMSE'
                
            sns.boxplot(
                data=d,
                x='Method', y=y, ax = axes[j, i]
            )
            
            axes[j, i].set_xlabel('')
            
            if j == 0:
                axes[j, i].set_title(titles[i])
            if j == 1:
                axes[j, i].tick_params(axis='x', labelrotation=90)
        
            if i != 0:
                axes[j, i].set_ylabel('')
                

    fig.align_ylabels()
    fig.supxlabel('Method')
    fig.suptitle('Quantitative Real Bulk Deconvolution Metrics')
    plt.tight_layout()

    plt.savefig('fig2a.png', bbox_inches='tight')

def plot_fig2c(path = "../data/fig2_data/"):
    def get_acc(df, labels):
        m = np.zeros_like(df.values)
        m[np.arange(len(df)), df.values.argmax(1)] = 1
        df1 = pd.DataFrame(m, columns=df.columns).astype(int)
        
        correct = 0
        incorrect = 0
        
        for i, ct in enumerate(labels):
            if df1.loc[i][ct] == 1:
                correct += 1
            else:
                incorrect += 1
                
        accuracy = correct / (incorrect + correct)

        return accuracy


    path = "/home/mcb/users/zhuang35/projects/GTM/revision_gb/fig2_data/"



    d_acc = {"B_Ex":"B cell",
        "B_NSM":"B cell",
        "B_SM":"B cell",
        "B_naive":"B cell",
        "Basophils":"Basophils",
        "CD4_TE":"T cell",
        "CD4_naive":"T cell",
        "CD8_CM":"T cell",
        "CD8_EM":"T cell",
        "CD8_TE":"T cell",
        "CD8_naive":"T cell",
        "C_mono":"Monocytes",
        "I_mono":"Monocytes",
        "MAIT":"T cell",
        "NC_mono":"Monocytes",
        "NK":"NK cell",
        "Neutrophils":"Neutrophils",
        "Plasmablasts":"Plasmablast",
        "Progenitor":"HSPC",
        "TFH":"T cell",
        "Th1":"T cell",
        "Th1.Th17":"T cell",
        "Th17":"T cell",
        "Th2":"T cell",
        "Treg":"T cell",
        "VD2.":"T cell",
        "VD2..1":"T cell",
        "mDC":"Dendritic cell (Myeloid)",
        "pDC":"Dendritic cell (Plasmacytoid)"}


    music = pd.read_csv(os.path.join(path, 'HBC_GSE107011_MUSIC_decon.csv'), index_col=0)
    bseqsc = pd.read_csv(os.path.join(path, 'HBC_GSE107011_BSEQSC_decon.csv'), index_col=0)
    bisque = pd.read_csv(os.path.join(path, 'HBC_GSE107011_BISQUE_decon.csv'), index_col=0)
    ciber = pd.read_csv(os.path.join(path, 'HBC_GSE107011_CIBERSORTx_decon.csv'), index_col=0)
    GTM_all = pd.read_csv(os.path.join(path, 'HBC_GSE107011_GTM_ALL_decon.csv'), index_col=0)
    GTM_pp = pd.read_csv(os.path.join(path, 'HBC_GSE107011_GTM_PP_decon.csv'), index_col=0)
    GTM_hvg = pd.read_csv(os.path.join(path, 'HBC_GSE107011_GTM_HVG_decon.csv'), index_col=0)
    bayesprism = pd.read_csv(os.path.join(path, 'HBC_GSE107011_BAYESPRISM_decon.csv'), index_col=0)

    keep = [i for i in music.index if 'PBMC' not in i]
        
    music = music.loc[keep]
    ciber = ciber.loc[keep]
    bseqsc = bseqsc.loc[keep]
    bisque = bisque.loc[keep]
    GTM_all = GTM_all.loc[keep]
    GTM_pp = GTM_pp.loc[keep]
    GTM_hvg = GTM_hvg.loc[keep]
    bayesprism = bayesprism.loc[keep]

    labels = bisque.index.str.split('_')
    labels = [i[1:] for i in labels]
    labels = ['_'.join(i) for i in labels]
    bisque['labels'] = labels
    bisque['labels'] = bisque['labels'].replace(d_acc)
    music['labels'] = bisque['labels'].values
    ciber['labels'] = bisque['labels'].values
    bseqsc['labels'] = bisque['labels'].values
    GTM_all['labels'] = bisque['labels'].values
    GTM_pp['labels'] = bisque['labels'].values
    GTM_hvg['labels'] = bisque['labels'].values
    bayesprism['labels'] = bisque['labels'].values

    relabel = {
        'T cell':'Tcells', 
        'B cell':'Bcells', 
        'NK cell':'NKcells', 
    }

    bisque['labels'] = bisque['labels'].replace(relabel)
    music['labels'] = bisque['labels'].values
    ciber['labels'] = bisque['labels'].values
    bseqsc['labels'] = bisque['labels'].values
    GTM_all['labels'] = bisque['labels'].values
    GTM_pp['labels'] = bisque['labels'].values
    GTM_hvg['labels'] = bisque['labels'].values
    bayesprism['labels'] = bisque['labels'].values

    acc_labels = sorted(set(bisque['labels']).intersection(bisque.columns[:-1]))

    GTM_all = GTM_all[GTM_all['labels'].isin(acc_labels)]
    ground_truth = list(GTM_all['labels'])
    GTM_pp = GTM_pp[GTM_pp['labels'].isin(acc_labels)]
    GTM_hvg = GTM_hvg[GTM_hvg['labels'].isin(acc_labels)]
    ciber = ciber[ciber['labels'].isin(acc_labels)]
    bseqsc = bseqsc[bseqsc['labels'].isin(acc_labels)]
    music = music[music['labels'].isin(acc_labels)]
    bisque = bisque[bisque['labels'].isin(acc_labels)]
    bayesprism = bayesprism[bayesprism['labels'].isin(acc_labels)]


    dfs = [GTM_all, GTM_pp, GTM_hvg, ciber, bseqsc, music, bayesprism, bisque]
    names = ['GTM-All', 'GTM-PP', 'GTM-HVG', 'CIBERSORTx', 'BSEQ-sc', 'MuSiC', 'BayesPrism', 'BISQUE']


    l = []
    experiment = "Reference: HBC, Bulk: GSE107011"
    for d, n in zip(dfs, names):
        acc = get_acc(d[acc_labels], ground_truth)
        l.append((n, acc, experiment))

    #####

    d_acc = {"B_Ex":"B_cell",
        "B_NSM":"B_cell",
        "B_SM":"B_cell",
        "B_naive":"B_cell",
        "CD4_TE":"CD4+_T_cell",
        "CD4_naive":"CD4+_T_cell",
        "CD8_CM":"Cytotoxic_T_cell",
        "CD8_EM":"Cytotoxic_T_cell",
        "CD8_TE":"Cytotoxic_T_cell",
        "CD8_naive":"Cytotoxic_T_cell",
        "C_mono":"CD14+_monocyte",
        "NC_mono":"CD16+_monocyte",
        "NK":"Natural_killer_cell",
        "TFH":"CD4+_T_cell",
        "Th1":"CD4+_T_cell",
        "Th1.Th17":"CD4+_T_cell",
        "mDC":"Dendritic_cell",
        "pDC":"Plasmacytoid_dendritic_cell"}


    music = pd.read_csv(os.path.join(path, 'PBMC_GSE107011_MUSIC_decon.csv'), index_col=0)
    bseqsc = pd.read_csv(os.path.join(path, 'PBMC_GSE107011_BSEQSC_decon.csv'), index_col=0)
    ciber = pd.read_csv(os.path.join(path, 'PBMC_GSE107011_CIBERSORTx_decon.csv'), index_col=0)
    GTM_all = pd.read_csv(os.path.join(path, 'PBMC_GSE107011_GTM_ALL_decon.csv'), index_col=0)
    GTM_pp = pd.read_csv(os.path.join(path, 'PBMC_GSE107011_GTM_PP_decon.csv'), index_col=0)
    GTM_hvg = pd.read_csv(os.path.join(path, 'PBMC_GSE107011_GTM_HVG_decon.csv'), index_col=0)
    bayesprism = pd.read_csv(os.path.join(path, 'PBMC_GSE107011_BAYESPRISM_decon.csv'), index_col=0)

    keep = [i for i in music.index if 'PBMC' not in i]


    labels = music.index.str.split('_')
    labels = [i[1:] for i in labels]
    labels = ['_'.join(i) for i in labels]
    music['labels'] = labels
    music['labels'] = music['labels'].replace(d_acc)
    ciber['labels'] = music['labels'].values
    bseqsc['labels'] = music['labels'].values
    GTM_all['labels'] = music['labels'].values
    GTM_pp['labels'] = music['labels'].values
    GTM_hvg['labels'] = music['labels'].values
    bayesprism['labels'] = music['labels'].values



    acc_labels = sorted(set(music['labels']).intersection(music.columns[:-1]))

    GTM_all = GTM_all[GTM_all['labels'].isin(acc_labels)]
    ground_truth = list(GTM_all['labels'])
    GTM_pp = GTM_pp[GTM_pp['labels'].isin(acc_labels)]
    GTM_hvg = GTM_hvg[GTM_hvg['labels'].isin(acc_labels)]
    ciber = ciber[ciber['labels'].isin(acc_labels)]
    bseqsc = bseqsc[bseqsc['labels'].isin(acc_labels)]
    music = music[music['labels'].isin(acc_labels)]
    bisque = bisque[bisque['labels'].isin(acc_labels)]
    bayesprism = bayesprism[bayesprism['labels'].isin(acc_labels)]


    dfs = [GTM_all, GTM_pp, GTM_hvg, ciber, bseqsc, music, bayesprism] #bisque]
    names = ['GTM-All', 'GTM-PP', 'GTM-HVG', 'CIBERSORTx', 'BSEQ-sc', 'MuSiC', 'BayesPrism']#'BISQUE']


    experiment = "Reference: PBMC2, Bulk: GSE107011"
    for d, n in zip(dfs, names):
        acc = get_acc(d[acc_labels], ground_truth)
        l.append((n, acc, experiment))


    #######


    music = pd.read_csv(os.path.join(path, 'HBC_GSE64655_MUSIC_decon.csv'), index_col=0)
    bseqsc = pd.read_csv(os.path.join(path, 'HBC_GSE64655_BSEQSC_decon.csv'), index_col=0)
    bisque = pd.read_csv(os.path.join(path, 'HBC_GSE64655_BISQUE_decon.csv'), index_col=0)
    ciber = pd.read_csv(os.path.join(path, 'HBC_GSE64655_CIBERSORTx_decon.csv'), index_col=0)
    GTM_all = pd.read_csv(os.path.join(path, 'HBC_GSE64655_GTM_ALL_decon.csv'), index_col=0)
    GTM_pp = pd.read_csv(os.path.join(path, 'HBC_GSE64655_GTM_PP_decon.csv'), index_col=0)
    GTM_hvg = pd.read_csv(os.path.join(path, 'HBC_GSE64655_GTM_HVG_decon.csv'), index_col=0)
    bayesprism = pd.read_csv(os.path.join(path, 'HBC_GSE64655_BAYESPRISM_decon.csv'), index_col=0)

    GSE64655 = pd.read_csv(os.path.join(path, 'GSE64655_annot.csv'), header=None)
    GSE64655 = GSE64655[GSE64655[2] != 'PBMC']
    GSE64655 = GSE64655.set_index(0, drop=True)
    GSE64655.columns=['drop', 'label']
    GSE64655 = GSE64655.drop(columns=['drop'])
    d = {
        'Monocyte':'Monocytes',
        'Neutrophil':'Neutrophils',
        'B_cell':'Bcells',
        'T_cell':'Tcells',
        'NK':'NKcells',
    }
    GSE64655['label'] = GSE64655['label'].replace(d)


    acc_labels = sorted(set(GSE64655['label']).intersection(music.columns[:-1]))
    GSE64655 = GSE64655[GSE64655['label'].isin(acc_labels)]

    music = music.loc[GSE64655.index]
    ciber = ciber.loc[GSE64655.index]
    bseqsc = bseqsc.loc[GSE64655.index]
    bisque = bisque.loc[GSE64655.index]
    GTM_hvg = GTM_hvg.loc[GSE64655.index]
    GTM_pp = GTM_pp.loc[GSE64655.index]
    GTM_all = GTM_all.loc[GSE64655.index]
    bayesprism = bayesprism.loc[GSE64655.index]

    dfs = [GTM_all, GTM_pp, GTM_hvg, ciber, bseqsc, music, bayesprism, bisque]
    names = ['GTM-All', 'GTM-PP', 'GTM-HVG', 'CIBERSORTx', 'BSEQ-sc', 'MuSiC', 'BayesPrism', 'BISQUE']


    ground_truth = GSE64655['label'].values
    experiment = "Reference: HBC, Bulk: GSE64655"
    for d, n in zip(dfs, names):
        acc = get_acc(d[acc_labels], ground_truth)
        l.append((n, acc, experiment))

    ########


    music = pd.read_csv(os.path.join(path, 'PBMC_GSE64655_MUSIC_decon.csv'), index_col=0)
    bseqsc = pd.read_csv(os.path.join(path, 'PBMC_GSE64655_BSEQSC_decon.csv'), index_col=0)
    # bisque = pd.read_csv(os.path.join(path, 'PBMC_GSE64655_BISQUE_decon.csv'), index_col=0)
    ciber = pd.read_csv(os.path.join(path, 'PBMC_GSE64655_CIBERSORTx_decon.csv'), index_col=0)
    GTM_all = pd.read_csv(os.path.join(path, 'PBMC_GSE64655_GTM_ALL_decon.csv'), index_col=0)
    GTM_pp = pd.read_csv(os.path.join(path, 'PBMC_GSE64655_GTM_PP_decon.csv'), index_col=0)
    GTM_hvg = pd.read_csv(os.path.join(path, 'PBMC_GSE64655_GTM_HVG_decon.csv'), index_col=0)
    bayesprism = pd.read_csv(os.path.join(path, 'PBMC_GSE64655_BAYESPRISM_decon.csv'), index_col=0)


    GSE64655 = pd.read_csv("/home/mcb/users/zhuang35/projects/GTM/revision_gb/PBMC_Real_Bulk/GSE64655_annot.csv", header=None)
    GSE64655 = GSE64655[GSE64655[2] != 'PBMC']
    GSE64655 = GSE64655.set_index(0, drop=True)
    GSE64655.columns=['drop', 'label']
    GSE64655 = GSE64655.drop(columns=['drop'])
    d = {
        'DC':'Dendritic cell', 
        'Monocyte':'Monocytes',
        'Neutrophil':'Neutrophils',
        'B_cell':'B_cell',
        'T_cell':'T_cell',
        'NK':'Natural_killer_cell',
    }
    GSE64655['label'] = GSE64655['label'].replace(d)

    music['Monocytes'] = music['CD14+_monocyte'] + music['CD16+_monocyte']
    music['T_cell'] = music['CD4+_T_cell'] + music['Cytotoxic_T_cell']
    ciber['Monocytes'] = ciber['CD14+_monocyte'] + ciber['CD16+_monocyte']
    ciber['T_cell'] = ciber['CD4+_T_cell'] + ciber['Cytotoxic_T_cell']
    bseqsc['Monocytes'] = bseqsc['CD14+_monocyte'] + bseqsc['CD16+_monocyte']
    bseqsc['T_cell'] = bseqsc['CD4+_T_cell'] + bseqsc['Cytotoxic_T_cell']
    GTM_all['Monocytes'] = GTM_all['CD14+_monocyte'] + GTM_all['CD16+_monocyte']
    GTM_all['T_cell'] = GTM_all['CD4+_T_cell'] + GTM_all['Cytotoxic_T_cell']
    GTM_pp['Monocytes'] = GTM_pp['CD14+_monocyte'] + GTM_pp['CD16+_monocyte']
    GTM_pp['T_cell'] = GTM_pp['CD4+_T_cell'] + GTM_pp['Cytotoxic_T_cell']
    GTM_hvg['Monocytes'] = GTM_hvg['CD14+_monocyte'] + GTM_hvg['CD16+_monocyte']
    GTM_hvg['T_cell'] = GTM_hvg['CD4+_T_cell'] + GTM_hvg['Cytotoxic_T_cell']
    bayesprism['Monocytes'] = bayesprism['CD14+_monocyte'] + bayesprism['CD16+_monocyte']
    bayesprism['T_cell'] = bayesprism['CD4+_T_cell'] + bayesprism['Cytotoxic_T_cell']

    acc_labels = sorted(set(GSE64655['label']).intersection(music.columns[:-1]))
    GSE64655 = GSE64655[GSE64655['label'].isin(acc_labels)]

    music = music.loc[GSE64655.index]
    ciber = ciber.loc[GSE64655.index]
    bseqsc = bseqsc.loc[GSE64655.index]
    # bisque = bisque.loc[GSE64655.index]
    GTM_hvg = GTM_hvg.loc[GSE64655.index]
    GTM_pp = GTM_pp.loc[GSE64655.index]
    GTM_all = GTM_all.loc[GSE64655.index]
    bayesprism = bayesprism.loc[GSE64655.index]
    

    dfs = [GTM_all, GTM_pp, GTM_hvg, ciber, bseqsc, music, bayesprism]#bisque]
    names = ['GTM-All', 'GTM-PP', 'GTM-HVG', 'CIBERSORTx', 'BSEQ-sc', 'MuSiC', 'BayesPrism']


    ground_truth = GSE64655['label'].values
    experiment = "Reference: PBMC2, Bulk: GSE64655"
    for d, n in zip(dfs, names):
        acc = get_acc(d[acc_labels], ground_truth)
        l.append((n, acc, experiment))



    #####

    df = pd.DataFrame(l, columns=['Method', 'Accuracy', 'Exp'])
    df['Exp1'] = ['\n'.join(i.split(', ')) for i in df['Exp']]
    temp = df['Exp'].str.split(', ')
    df['Ref'] = [i[0] for i in temp]
    df['Bulk'] = [i[1] for i in temp]
    df['Bulk'] = df['Bulk'].str.replace('Bulk: ', '')


    sns.set(rc={'figure.figsize':(7, 3), 'figure.dpi':350})
    fig, axes = plt.subplots(1, 2, sharey=True)

    sns.barplot(
        data=df[df['Exp1'].isin(
            df['Exp1'].unique()[[0, 2]]
        )],
        x='Bulk',
        y='Accuracy',
        hue='Method',
        ax=axes[0],
    )

    axes[0].set_xlabel('Reference: HBC')
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend_.remove()

    handles
    sns.barplot(
        data=df[df['Exp1'].isin(
            df['Exp1'].unique()[[1, 3]]
        )],
        x='Bulk',
        y='Accuracy',
        hue='Method',
        ax=axes[1],
    )
    # axes[1].legend_.bbox_to_anchor([1, 1.05])
    axes[1].legend(handles=handles, labels=labels, title='Method', bbox_to_anchor=(1, .975), fontsize=8)

    # axes[1].legend_.remove()
    axes[1].set_xlabel('Reference: PBMC')
    # plt.legend(bbox_to_anchor=(1, .95))

    plt.suptitle('Purified Bulk PBMC Prediction Accuracy')

    plt.tight_layout()



    plt.savefig('fig2c.png', bbox_inches='tight')

def plot_fig2b_supp(path = "../data/fig2_data/"):
    files = sorted(os.listdir(path))
    gts = {}

    exp_pairs = (
        ('HBC', 'WB'),
        ('Lake', 'ROSMAP'),
        ('PBMC', 'GSE107011'),
        ('PBMC', 'SDY67'),
        ('SegerstolpeSC', 'SegerstolpeBulk')
    )


    COMP_STATS = {f"{ref}_{bulk}" : None for ref, bulk in exp_pairs}
    for ref, bulk in exp_pairs:
        
        methods = []
        names = []
        for file in files:
            if ref in file and bulk in file:
                
                
                if 'proportions' in file:
                    ground_truth = pd.read_csv(os.path.join(path, file), index_col=0)
                    gts[f"{ref}_{bulk}"] = ground_truth
                else:
                    methods.append(pd.read_csv(os.path.join(path, file), index_col=0))
                    names.append(file.split('_')[2:])
        names = [n[0] if len(n) == 2 else '_'.join(n[:2]) for n in names]
        
        COMP_STATS[f"{ref}_{bulk}"] = compile_stats_celltypewise(
            ground_truth,
            methods,
            names
        )

    l = []
    for exp in COMP_STATS.keys():
        
        df = COMP_STATS[exp]
        
        for method in df.index:
            for pcc in df.loc[method]:
                l.append((exp, method, pcc))

    rename = {'BSEQSC': 'BSEQ-sc',
        'CIBER': 'CIBERSORTx',
        'CIBERSORTx' : 'CIBERSORTx',
        'MUSIC': 'MuSiC',
        'GTM_ALL': 'GTM-ALL',
        'GTM_PP': 'GTM-PP',
        'GTM_HVG': 'GTM-HVG',
        'BISQUE': 'BISQUE',
        'BAYESPRISM': 'BayesPrism',}
    order = ['GTM-ALL', 'GTM-PP', 'GTM-HVG', 'BSEQ-sc', 'CIBERSORTx', 'MuSiC', 'BayesPrism', 'BISQUE']

    df = pd.DataFrame(l, columns=['Experiment', 'Method', 'PCC'])
    df['Method'] = df['Method'].map(rename)
    df.Method = df.Method.astype('category')
    df.Method = df.Method.cat.set_categories(order)
    df.sort_values('Method')



    mapping = {
        'PBMC_SDY67': 'Ref: PBMC\nBulk: PBMC-1', 
        'PBMC_GSE107011' : 'Ref: PBMC\nBulk: PBMC-2',
        'HBC_WB' : 'Ref: BC\nBulk: WB', 
        'Lake_ROSMAP':'Ref: FC\nBulk: PFC',
        'SegerstolpeSC_SegerstolpeBulk':  'Ref Pancreas\nBulk: Pancreas'
    }
    df['Experiment'] = df['Experiment'].map(mapping)



    rename = {'BSEQSC': 'BSEQ-sc',
        'CIBER': 'CIBERSORTx',
        'CIBERSORTx' : 'CIBERSORTx',
        'MUSIC': 'MuSiC',
        'GTM_ALL': 'GTM-ALL',
        'GTM_PP': 'GTM-PP',
        'GTM_HVG': 'GTM-HVG',
        'BISQUE': 'BISQUE',
        'BAYESPRISM': 'BayesPrism',}

    for k in COMP_STATS:
        COMP_STATS[k].index = COMP_STATS[k].index.map(rename)

        
    plot_order = [
        'PBMC_SDY67',
        'PBMC_GSE107011',
        'HBC_WB',
        'Lake_ROSMAP',
        'SegerstolpeSC_SegerstolpeBulk'
    ]
    
    gts = [gts[p] for p in plot_order]

    plot_order = [COMP_STATS[p] for p in plot_order]


    titles = ['Ref: PBMC\nBulk: PBMC-1', 
            'Ref: PBMC\nBulk: PBMC-2', 
            'Ref: Blood Cells\nBulk: Whole Blood', 
            'Ref: Frontal Cortex\nBulk: Prefrontal Cortex', 
            'Ref Pancreas\nBulk: Pancreas']



    order = ['GTM-ALL', 'GTM-PP', 'GTM-HVG', 'BSEQ-sc', 'CIBERSORTx', 'MuSiC', 'BayesPrism', 'BISQUE']
    
    for p in plot_order:
        if 'BISQUE' not in p.index:
            p.loc['BISQUE'] = np.nan

    titles = ['Ref: PBMC\nBulk: PBMC-1', 
            'Ref: PBMC\nBulk: PBMC-2', 
            'Ref: BC\nBulk: WB', 
            'Ref: FC\nBulk: PFC', 
            'Ref Pancreas\nBulk: Pancreas']
    

    ground_truth_list = []
    for g in gts:
        g = pd.DataFrame(g.mean(), columns=['val'])
        g = g / g.sum()
        g = g.reset_index()
        ground_truth_list.append(g)

    files = sorted(os.listdir(path))
    gts = {}

    exp_pairs = (
        ('HBC', 'GSE107011'),
        ('HBC', 'WB'),
        ('Lake', 'ROSMAP'),
        ('PBMC', 'GSE107011'),
        ('PBMC', 'SDY67'),
        ('SegerstolpeSC', 'SegerstolpeBulk')
    )


    COMP_STATS = {f"{ref}_{bulk}" : None for ref, bulk in exp_pairs}
    for ref, bulk in exp_pairs:
        
        methods = []
        names = []
        for file in files:
            if ref in file and bulk in file:
                
                
                if 'proportions' in file:
                    ground_truth = pd.read_csv(os.path.join(path, file), index_col=0)
                    gts[f"{ref}_{bulk}"] = ground_truth
                else:
                    methods.append(pd.read_csv(os.path.join(path, file), index_col=0))
                    names.append(file.split('_')[2:])
        names = [n[0] if len(n) == 2 else '_'.join(n[:2]) for n in names]
        
        l = []
        for ct in ground_truth.columns:
            l2 = []
            for method in methods:

                _, _, _, _, RMSE, _ = compute_stats(
                    ground_truth[ct],
                    method.loc[ground_truth.index][ct]
                )
                l2.append(1 / RMSE)

            l.append(l2)

        df = pd.DataFrame(l, columns=names, index=ground_truth.columns).T
        COMP_STATS[f"{ref}_{bulk}"] = df


    rename = {'BSEQSC': 'BSEQ-sc',
        'CIBER': 'CIBERSORTx',
        'CIBERSORTx' : 'CIBERSORTx',
        'MUSIC': 'MuSiC',
        'GTM_ALL': 'GTM-ALL',
        'GTM_PP': 'GTM-PP',
        'GTM_HVG': 'GTM-HVG',
        'BISQUE': 'BISQUE',
        'BAYESPRISM': 'BayesPrism',}

    for k in COMP_STATS:
        COMP_STATS[k].index = COMP_STATS[k].index.map(rename)

        
    plot_order2 = [
        'PBMC_SDY67',
        'PBMC_GSE107011',
        'HBC_WB',
        'Lake_ROSMAP',
        'SegerstolpeSC_SegerstolpeBulk'
    ]
    
    gts = [gts[p] for p in plot_order2]

    plot_order2 = [COMP_STATS[p] for p in plot_order2]


    titles = ['Ref: PBMC\nBulk: PBMC-1',
     'Ref: PBMC\nBulk: PBMC-2',
     'Ref: BC \nBulk: WB',
     'Ref: FC\nBulk: PFC',
     'Ref Pancreas\nBulk: Pancreas']

    
    for p in plot_order2:
        if 'BISQUE' not in p.index:
            p.loc['BISQUE'] = np.nan

    sns.set(rc={'figure.figsize':(11, 7), 'figure.dpi':300})



    po3 = [
        p / p.max().max() for p in plot_order2
    ]

    l = []

    for col in order:
        for v in pd.concat(plot_order, axis=1).T[col]:
            l.append((col, v))
            
    df = pd.DataFrame(l, columns=['Method', 'PCC'])

    l = []

    for col in order:
        for v in pd.concat(po3, axis=1).T[col]:
            l.append((col, v))
            
    df2 = pd.DataFrame(l, columns=['Method', 'Inv-RMSE'])

    sns.set(rc={'figure.figsize':(11, 7), 'figure.dpi':300})

    cell_type_orders = []
    for i in range(5):
        g = ground_truth_list[i]
        g = g.sort_values('val', ascending=False)
        cell_type_orders.append([i for i in g['index'].values if 'unclassified' not in i])

    g1 = ['CD4+ T', 'CD14+ Mono', 'NK Cell', 'Cytotoxic T', 'CD16+ Mono', 'B Cell', 'DC', 'pDC']
    g2 = ['CD4+ T', 'CD14+ Mono', 'Cytotoxic T', 'NK Cell', 'B Cell', 'CD16+ Mono', 'DC', 'pDC']
    g3 = ['Neutrophil', 'T Cell', 'Mono', 'B Cell', 'NK Cell']
    g4 = ['Neuron', 'Endo', 'Astro', 'Oligo', 'Mic']
    g5 = ['Alpha', 'Ductal', 'Beta', 'Gamma', 'Acinar', 'Delta', 'PSC', 'Co-Exp', 'Endothelial', 'Epsilon', 'MHC II', 'Mast']


    width_ratios = [len(i) * 3 for i in cell_type_orders] + [11, 1]

    fig, axes = plt.subplots(
        3, 
        7, 
        gridspec_kw={'height_ratios':[1.5,5,5], 'width_ratios' : width_ratios},
    )
    plt.subplots_adjust(
        wspace=0.1,
        hspace=0.075,
    )
    flat = axes.flatten()



    for i in range(5):
        d = plot_order[i]
        d2 = plot_order2[i]
        d2 = d2 / d2.max().max()
        g = ground_truth_list[i]
        
        sns.barplot(
            data = g.set_index('index').loc[cell_type_orders[i]].reset_index(),
            y='val',
            x='index',
            ax = flat[i],
            linewidth=0,
            color='darkseagreen'

        )

        flat[i].set(ylim=(-0.15, 0.7))
        flat[i].set_xticks([])
        flat[i].set_ylabel('')
        flat[i].set_xlabel('')

        
        sns.heatmap(
            data = d.loc[order][cell_type_orders[i]],
            ax = flat[7 + i],

            cmap='vlag',
            vmin=-1, vmax=1,
            xticklabels=False,
            cbar = i == 4,
            cbar_ax = flat[13],
            cbar_kws={'ticks':[-1, -0.5, 0, 0.5, 1], 'shrink':0.5}
        )

        
        sns.heatmap(
            data = d2.loc[order][cell_type_orders[i]],
            ax = flat[14 + i],
            cmap='Reds',
            vmin=0, vmax=1,
            xticklabels=True, 
            cbar = i == 4,
            cbar_ax = flat[20],
            cbar_kws={'ticks':[-1, -0.5, 0, 0.5, 1]}
            
        )
        
        
        flat[i].set_title(titles[i])
        

        

    for i in range(len(flat)):
        if i < 12:
            flat[i].set_xticks([])
            
        if i % 7 != 0:
            if i in [13, 20]:
                continue
            flat[i].set_yticklabels([])
            

            
    flat[5].axis('off')
    flat[6].axis('off')
    sns.barplot(data=df, x='PCC', y='Method', 
                capsize=0.2, errwidth=1, ax=flat[12], 
                color='indianred', 
            )

    flat[12].xaxis.tick_top()
    flat[12].set_yticklabels([])
    flat[12].set_ylabel('')
    flat[12].set_xlim((0, 0.7))
    flat[12].set_xlabel('')


    sns.barplot(data=df2, x='Inv-RMSE', y='Method', 
                capsize=0.2, errwidth=1, ax=flat[19], 
                color='indianred', 
            )

    flat[19].set_yticklabels([])
    flat[19].set_ylabel('')
    flat[19].set_xlim((0, 0.7))
    flat[19].set_xlabel('')

    flat[0].set_ylabel('Bulk\nProportion')
    flat[7].set_ylabel('PCC')
    flat[14].set_ylabel(r'$RMSE^{-1}$')

    flat[14].set_xticklabels(g1)
    flat[15].set_xticklabels(g2)
    flat[16].set_xticklabels(g3)
    flat[17].set_xticklabels(g4)
    flat[18].set_xticklabels(g5)


    flat[12].tick_params(axis='both', which='major', labelsize=8)
    flat[13].tick_params(axis='both', which='major', labelsize=8)


    flat[19].tick_params(axis='both', which='major', labelsize=8)
    flat[20].tick_params(axis='both', which='major', labelsize=8)


    fig.align_ylabels()
    plt.suptitle('Cell-Type Specific Performance of Bulk Deconvolution')


    plt.savefig('fig2b.png', bbox_inches='tight')

def plot_fig2b(path = "../data/fig2_data/"):
    files = sorted(os.listdir(path))
    gts = {}

    exp_pairs = (
        ('HBC', 'WB'),
        ('Lake', 'ROSMAP'),
        ('PBMC', 'GSE107011'),
        ('PBMC', 'SDY67'),
        ('SegerstolpeSC', 'SegerstolpeBulk')
    )


    COMP_STATS = {f"{ref}_{bulk}" : None for ref, bulk in exp_pairs}
    for ref, bulk in exp_pairs:
        
        methods = []
        names = []
        for file in files:
            if ref in file and bulk in file:
                
                
                if 'proportions' in file:
                    ground_truth = pd.read_csv(os.path.join(path, file), index_col=0)
                    gts[f"{ref}_{bulk}"] = ground_truth
                else:
                    methods.append(pd.read_csv(os.path.join(path, file), index_col=0))
                    names.append(file.split('_')[2:])
        names = [n[0] if len(n) == 2 else '_'.join(n[:2]) for n in names]
        
        COMP_STATS[f"{ref}_{bulk}"] = compile_stats_celltypewise(
            ground_truth,
            methods,
            names,
            spearman=True,
        )

    l = []
    for exp in COMP_STATS.keys():
        
        df = COMP_STATS[exp]
        
        for method in df.index:
            for pcc in df.loc[method]:
                l.append((exp, method, pcc))

    rename = {'BSEQSC': 'BSEQ-sc',
        'CIBER': 'CIBERSORTx',
        'CIBERSORTx' : 'CIBERSORTx',
        'MUSIC': 'MuSiC',
        'GTM_ALL': 'GTM-ALL',
        'GTM_PP': 'GTM-PP',
        'GTM_HVG': 'GTM-HVG',
        'BISQUE': 'BISQUE',
        'BAYESPRISM': 'BayesPrism',}
    order = ['GTM-ALL', 'GTM-PP', 'GTM-HVG', 'BSEQ-sc', 'CIBERSORTx', 'MuSiC', 'BayesPrism', 'BISQUE']

    df = pd.DataFrame(l, columns=['Experiment', 'Method', 'PCC'])
    df['Method'] = df['Method'].map(rename)
    df.Method = df.Method.astype('category')
    df.Method = df.Method.cat.set_categories(order)
    df.sort_values('Method')



    mapping = {
        'PBMC_SDY67': 'Ref: PBMC\nBulk: PBMC-1', 
        'PBMC_GSE107011' : 'Ref: PBMC\nBulk: PBMC-2',
        'HBC_WB' : 'Ref: BC\nBulk: WB', 
        'Lake_ROSMAP':'Ref: FC\nBulk: PFC',
        'SegerstolpeSC_SegerstolpeBulk':  'Ref Pancreas\nBulk: Pancreas'
    }
    df['Experiment'] = df['Experiment'].map(mapping)



    rename = {'BSEQSC': 'BSEQ-sc',
        'CIBER': 'CIBERSORTx',
        'CIBERSORTx' : 'CIBERSORTx',
        'MUSIC': 'MuSiC',
        'GTM_ALL': 'GTM-ALL',
        'GTM_PP': 'GTM-PP',
        'GTM_HVG': 'GTM-HVG',
        'BISQUE': 'BISQUE',
        'BAYESPRISM': 'BayesPrism',}

    for k in COMP_STATS:
        COMP_STATS[k].index = COMP_STATS[k].index.map(rename)

        
    plot_order = [
        'PBMC_SDY67',
        'PBMC_GSE107011',
        'HBC_WB',
        'Lake_ROSMAP',
        'SegerstolpeSC_SegerstolpeBulk'
    ]
    
    gts = [gts[p] for p in plot_order]

    plot_order = [COMP_STATS[p] for p in plot_order]


    titles = ['Ref: PBMC\nBulk: PBMC-1', 
            'Ref: PBMC\nBulk: PBMC-2', 
            'Ref: Blood Cells\nBulk: Whole Blood', 
            'Ref: Frontal Cortex\nBulk: Prefrontal Cortex', 
            'Ref Pancreas\nBulk: Pancreas']

    order = ['GTM-ALL', 'GTM-PP', 'GTM-HVG', 'BSEQ-sc', 'CIBERSORTx', 'MuSiC', 'BayesPrism', 'BISQUE']
    
    for p in plot_order:
        if 'BISQUE' not in p.index:
            p.loc['BISQUE'] = np.nan

    titles = ['Ref: PBMC\nBulk: PBMC-1', 
            'Ref: PBMC\nBulk: PBMC-2', 
            'Ref: BC\nBulk: WB', 
            'Ref: FC\nBulk: PFC', 
            'Ref Pancreas\nBulk: Pancreas']
    

    ground_truth_list = []
    for g in gts:
        g = pd.DataFrame(g.mean(), columns=['val'])
        g = g / g.sum()
        g = g.reset_index()
        ground_truth_list.append(g)

    files = sorted(os.listdir(path))
    gts = {}

    exp_pairs = (
        ('HBC', 'GSE107011'),
        ('HBC', 'WB'),
        ('Lake', 'ROSMAP'),
        ('PBMC', 'GSE107011'),
        ('PBMC', 'SDY67'),
        ('SegerstolpeSC', 'SegerstolpeBulk')
    )


    COMP_STATS = {f"{ref}_{bulk}" : None for ref, bulk in exp_pairs}
    for ref, bulk in exp_pairs:
        
        methods = []
        names = []
        for file in files:
            if ref in file and bulk in file:
                
                
                if 'proportions' in file:
                    ground_truth = pd.read_csv(os.path.join(path, file), index_col=0)
                    gts[f"{ref}_{bulk}"] = ground_truth
                else:
                    methods.append(pd.read_csv(os.path.join(path, file), index_col=0))
                    names.append(file.split('_')[2:])
        names = [n[0] if len(n) == 2 else '_'.join(n[:2]) for n in names]
        
        l = []
        for ct in ground_truth.columns:
            l2 = []
            for method in methods:

                _, _, _, _, RMSE, _ = compute_stats(
                    ground_truth[ct],
                    method.loc[ground_truth.index][ct]
                )
                l2.append(1 / RMSE)

            l.append(l2)

        df = pd.DataFrame(l, columns=names, index=ground_truth.columns).T
        COMP_STATS[f"{ref}_{bulk}"] = df


    rename = {'BSEQSC': 'BSEQ-sc',
        'CIBER': 'CIBERSORTx',
        'CIBERSORTx' : 'CIBERSORTx',
        'MUSIC': 'MuSiC',
        'GTM_ALL': 'GTM-ALL',
        'GTM_PP': 'GTM-PP',
        'GTM_HVG': 'GTM-HVG',
        'BISQUE': 'BISQUE',
        'BAYESPRISM': 'BayesPrism',}

    for k in COMP_STATS:
        COMP_STATS[k].index = COMP_STATS[k].index.map(rename)

        
    plot_order2 = [
        'PBMC_SDY67',
        'PBMC_GSE107011',
        'HBC_WB',
        'Lake_ROSMAP',
        'SegerstolpeSC_SegerstolpeBulk'
    ]
    
    gts = [gts[p] for p in plot_order2]

    plot_order2 = [COMP_STATS[p] for p in plot_order2]


    titles = ['Ref: PBMC\nBulk: PBMC-1',
     'Ref: PBMC\nBulk: PBMC-2',
     'Ref: BC \nBulk: WB',
     'Ref: FC\nBulk: PFC',
     'Ref Pancreas\nBulk: Pancreas']

    
    for p in plot_order2:
        if 'BISQUE' not in p.index:
            p.loc['BISQUE'] = np.nan

    sns.set(rc={'figure.figsize':(11, 7), 'figure.dpi':300})



    po3 = [
        p / p.max().max() for p in plot_order2
    ]

    l = []

    for col in order:
        for v in pd.concat(plot_order, axis=1).T[col]:
            l.append((col, v))
            
    df = pd.DataFrame(l, columns=['Method', 'PCC'])

    l = []

    for col in order:
        for v in pd.concat(po3, axis=1).T[col]:
            l.append((col, v))
            
    df2 = pd.DataFrame(l, columns=['Method', 'Inv-RMSE'])

    sns.set(rc={'figure.figsize':(11, 7), 'figure.dpi':300})

    cell_type_orders = []
    for i in range(5):
        g = ground_truth_list[i]
        g = g.sort_values('val', ascending=False)
        cell_type_orders.append([i for i in g['index'].values if 'unclassified' not in i])

    g1 = ['CD4+ T', 'CD14+ Mono', 'NK Cell', 'Cytotoxic T', 'CD16+ Mono', 'B Cell', 'DC', 'pDC']
    g2 = ['CD4+ T', 'CD14+ Mono', 'Cytotoxic T', 'NK Cell', 'B Cell', 'CD16+ Mono', 'DC', 'pDC']
    g3 = ['Neutrophil', 'T Cell', 'Mono', 'B Cell', 'NK Cell']
    g4 = ['Neuron', 'Endo', 'Astro', 'Oligo', 'Mic']
    g5 = ['Alpha', 'Ductal', 'Beta', 'Gamma', 'Acinar', 'Delta', 'PSC', 'Co-Exp', 'Endothelial', 'Epsilon', 'MHC II', 'Mast']

    sns.set(rc={'figure.figsize':(10, 6), 'figure.dpi':350})
    fig, axes = plt.subplots(2, 5, sharex=True)


    for index, df in enumerate(plot_order):
        temp = []

        for i in df.index:
            for j in df.columns:
                temp.append((i, j, df.loc[i][j]))

        df = pd.DataFrame(temp, columns=['Method', 'Celltype', "Spearman R"])
        
        d = df
        d.Method = d.Method.astype('category')
        d.Method = d.Method.cat.set_categories(order)
        d.sort_values('Method')
        
        sns.boxplot(data=d, x='Method', y="Spearman R", ax = axes.flatten()[index])
        
    for index, df in enumerate(plot_order2):
        temp = []

        for i in df.index:
            for j in df.columns:
                temp.append((i, j, 1/df.loc[i][j]))

        df = pd.DataFrame(temp, columns=['Method', 'Celltype', 'RMSE'])
        
        d = df
        d.Method = d.Method.astype('category')
        d.Method = d.Method.cat.set_categories(order)
        d.sort_values('Method')
        
        sns.boxplot(data=d, x='Method', y='RMSE', ax = axes.flatten()[index+5])
        
    for ax in axes.flatten():
        ax.set_xlabel('')

    for i in range(5):
        axes.flatten()[i].set_title(titles[i])
    for i in range(5, 10):
        axes.flatten()[i].tick_params(axis='x', labelrotation=90)
        
    for i in range(1, 5):
        axes[0, i].set_ylabel('')
        axes[1, i].set_ylabel('')


    fig.align_ylabels()
    fig.supxlabel('Method')
    fig.suptitle('Cell-type Specific Real Bulk Deconvolution Performance')
    plt.tight_layout()

    plt.savefig('fig2b.png', bbox_inches='tight')



if __name__ == "__main__":
    sns.set(rc={'figure.figsize':(10, 6), 'figure.dpi':350})
    plot_fig2a()

    sns.set(rc={'figure.figsize':(10, 6), 'figure.dpi':350})
    plot_fig2b()

    sns.set(rc={'figure.figsize':(11, 7), 'figure.dpi':300})
    plot_fig2b_supp()

    sns.set(rc={'figure.figsize':(10, 6), 'figure.dpi':350})
    plot_fig2c()



    