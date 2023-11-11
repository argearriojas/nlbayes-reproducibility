from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
import statsmodels.formula.api as smf

import pandas as pd
import decoupler as dc
from scipy.stats import false_discovery_control
import matplotlib.pyplot as plt
import upsetplot as usp


def make_figure_3(results):
    hue = 'evd_rnd_p'
    df = results.query('net_rnd_p == 0.00')

    fig, axs = plt.subplots(2, 3, figsize=(12,8))

    plt.sca(axs[0, 0])
    for cat in df[hue].sort_values().unique():
        tmp = df.loc[df[hue] == cat]
        fpr, tpr, thr = metrics.roc_curve(tmp.gt_act, tmp.X)
        auc = metrics.auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f"p={cat:.2f}; auc = {auc:.2f}")
    plt.xlim(-.05,1.05)
    plt.ylim(-.05,1.05)
    plt.plot((0,1), (0,1), c='gray', linestyle='--')
    plt.gca().set_aspect('equal')
    plt.legend()
    plt.title('Data randomization\nROC curve')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.grid(visible=True, which='major', color='lightgray', linestyle='-')
    plt.text(-0.30, 1.13,'A)', {'size':16})

    plt.sca(axs[0, 1])
    for cat in df[hue].sort_values().unique():
        tmp = df.loc[df[hue] == cat]
        ppv, tpr, thr = metrics.precision_recall_curve(tmp.gt_act, tmp.X)
        auc = metrics.auc(tpr, ppv)
        plt.plot(tpr, ppv, label=f"p={cat:.2f}; auc={auc:.2f}")
    plt.xlim(-.05,1.05)
    plt.ylim(-.05,1.05)
    plt.plot((0,1), (1,0), c='gray', linestyle='--')
    plt.gca().set_aspect('equal')
    plt.legend(loc='lower left')
    plt.title('Data randomization\nPrecision vs Recall curve')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.grid(visible=True, which='major', color='lightgray', linestyle='-')
    plt.text(-0.30, 1.13,'B)', {'size':16})

    hue = 'net_rnd_p'
    df = results.query('evd_rnd_p == 0.00')
    plt.sca(axs[1, 0])
    for cat in df[hue].sort_values().unique():
        tmp = df.loc[df[hue] == cat]
        fpr, tpr, thr = metrics.roc_curve(tmp.gt_act, tmp.X)
        auc = metrics.auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f"p={cat:.2f}; auc = {auc:.2f}")
    plt.xlim(-.05,1.05)
    plt.ylim(-.05,1.05)
    plt.plot((0,1), (0,1), c='gray', linestyle='--')
    plt.gca().set_aspect('equal')
    plt.legend()
    plt.title('Network randomization\nROC curve')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.grid(visible=True, which='major', color='lightgray', linestyle='-')
    plt.text(-0.30, 1.13,'C)', {'size':16})

    plt.sca(axs[1, 1])
    for cat in df[hue].sort_values().unique():
        tmp = df.loc[df[hue] == cat]
        ppv, tpr, thr = metrics.precision_recall_curve(tmp.gt_act, tmp.X)
        auc = metrics.auc(tpr, ppv)
        plt.plot(tpr, ppv, label=f"p={cat:.2f}; auc={auc:.2f}")
    plt.xlim(-.05,1.05)
    plt.ylim(-.05,1.05)
    plt.plot((0,1), (1,0), c='gray', linestyle='--')
    plt.gca().set_aspect('equal')
    plt.legend(loc='lower left')
    plt.title('Network randomization\nPrecision vs Recall curve')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.grid(visible=True, which='major', color='lightgray', linestyle='-')
    plt.text(-0.30, 1.13,'D)', {'size':16})

    df = results.query('net_rnd_p == 0.00')
    plt.sca(axs[0, 2])
    sns.violinplot(data=df, x="X", y="evd_rnd_p", hue="gt_act", split=True, cut=0, bw_method="scott", inner=None, density_norm='width', common_norm=False, linewidth=0.75, orient='horiz')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='TF is active\nground truth')
    plt.ylabel('Data randomization level (p)')
    plt.xlabel('Active TF posterior probability')
    plt.title('Posterior distributions\nfor TF activity')
    plt.text(-0.40, -0.9,'E)', {'size':16})

    tmp = df.query('evd_rnd_p == 0.25 and gt_act').sort_values('trg_tot')
    tmp = tmp.groupby('trg_tot').mean(numeric_only=True).reset_index()
    x = tmp.trg_tot
    formula = ('X ~ bs(trg_tot, df=8, degree=3)')
    f = smf.ols(formula=formula, data=tmp).fit().predict
    plt.sca(axs[1, 2])
    sns.scatterplot(data=df.query('evd_rnd_p == 0.25'), x="trg_tot", y="X", hue="gt_act")
    plt.legend(loc='upper left', title='TF is active\nground truth')
    plt.xlabel('TF total number of target genes')
    plt.ylabel('Active TF posterior probability')
    plt.title('Posterior distributions for TF activity\nat 25% data randomization level')
    plt.xscale('log')
    plt.plot(x, f(x), c='C1')
    plt.text(0.85, 1.15,'F)', {'size':16})

    plt.tight_layout()

    filepath = Path('figures/fig3.png')
    filepath.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(filepath, dpi=300)
    plt.show()


def make_figure_5():
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    r = robjects.r
    r.source('r_scripts/nlbayes_utils.R')

    annot = pd.read_csv('data/annot.tsv', sep='\t')
    annot['Gene.ID'] = annot['Gene.ID'].astype(str)
    symbol2id = annot.set_index('Gene.symbol')
    annot = annot.set_index('Gene.ID')

    tT = pd.read_csv('data/myc_experiment_tT.tsv', sep='\t').set_index('Gene.ID')
    tT.index = tT.index.astype(str)
    msk = (tT['adj.P.Val'] < 0.05)
    mat = tT.where(msk, 0).loc[:, ['logFC']].T

    regulon = r['regulonbrca']
    net = r['regulon2long'](regulon)
    net = pandas2ri.rpy2py_dataframe(net).set_index('uid')
    net = net.query('target in @mat.columns')
    mat = mat.loc[:, net.target.unique()]

    e, p = dc.run_ulm(mat, net)
    ulm_df = pd.concat([e, p], axis=0).T
    ulm_df.columns = ['e', 'p']

    e, n, p = dc.run_gsea(mat, net)
    gsea_df = pd.concat([e, n, p], axis=0).T
    gsea_df.columns = ['e', 'n', 'p']

    ulm_df['fdr'] = false_discovery_control(ulm_df.p)
    gsea_df['fdr'] = false_discovery_control(gsea_df.p)

    filename = 'data/oe_myc_on_net_regulonbrca_nlbayes_and_viper.csv'
    nlbayes_df = pd.read_csv(filename)
    nlbayes_df['id'] = symbol2id.loc[nlbayes_df.symbol, 'Gene.ID'].values
    nlbayes_df = nlbayes_df.set_index('id')
    nlbayes_df['viper.pvalue'].fillna(1, inplace=True)

    df = pd.DataFrame({
        'gsea': gsea_df.query('fdr < 0.05').p,
        'ulm': ulm_df.query('fdr < 0.05').p,
        'viper': nlbayes_df.query('`viper.pvalue` < 0.05')['viper.pvalue'],
        'nlbayes': nlbayes_df.query('`posterior.p` > 0.2')['posterior.p']})
    df = df.loc[df.index.isin(annot.index)]
    df.index = annot.loc[df.index, 'Gene.symbol'].values
    df.index.name = 'gene_symbol'
    df = df.reset_index().melt(id_vars=['gene_symbol'], var_name='method', value_name='p').dropna()
    memberships = df.groupby('method')['gene_symbol'].apply(list).to_dict()
    fig = plt.figure()
    upset = usp.UpSet(usp.from_contents(memberships), show_counts=True, sort_by='cardinality')
    upset.style_subsets(present=["nlbayes"], facecolor='blue')
    upset.plot(fig=fig)
    plt.savefig('figures/fig5.png', dpi=300)
    plt.show()
