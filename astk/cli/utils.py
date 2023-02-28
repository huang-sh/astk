"""
astk.cli.utils
~~~~~~~~~~~~~~~~~
This module provide some utility function
"""

from .config import *
import astk.utils._cli_func as ul


@cli_fun.command(help="generate metadata template file")
@click.option('-c1', '--ctrl', "ctrl_file", cls=MultiOption, type=click.Path(exists=True),
                required=True, help="file path for condtion 1(control)")
@click.option('-c2', '--case', "case_file", cls=MultiOption, type=click.Path(exists=True),
                required=True, help="file path for condtion 2(case)")
@click.option('-gn', '--groupName', cls=MultiOption, type=str,
                default=(), help="group name")
@click.option('-repN', '--replicate', cls=MultiOption, type=int,
                help="replicate, number")
@click.option('-o', '--output', required=True, type=click.Path(),
                help='metadata output path')
@click.option('-repN1', '--replicate1', cls=MultiOption, type=int,
                help="replicate1 number(control)")
@click.option('-repN2', '--replicate2', cls=MultiOption, type=int, 
                help="replicate2, number(case)")
@click.option('--condition', cls=MultiOption, type=str, help="condition name, default='ctrl case'")
@click.option('-fn', '--filename', 'is_fn', is_flag=True, help="file name")
@click.option('--split', cls=MultiOption, type=str, help="name split symbol and index")
@click.option('-app', '--app', type=click.Choice(['SUPPA2', 'rMATS']), default="SUPPA2",
                help="application for output; default='SUPPA2'")
def meta(*args, **kwargs):
    if not (pdir:= Path(kwargs["output"]).parent).exists():
        raise UsageError(f"Invalid value for '-o' / '--output': Path '{pdir}' does not exist.")
    replicate = kwargs.pop("replicate")
    replicate1, replicate2 = kwargs.pop("replicate1"), kwargs.pop("replicate2")
    if replicate:
        kwargs["ctrl_rep"] = [int(i) for i in replicate]
        kwargs["case_rep"] = [int(i) for i in replicate]
    elif all([replicate1, replicate2]):
        kwargs["ctrl_rep"] = [int(i) for i in replicate1]
        kwargs["case_rep"] = [int(i) for i in replicate2]
    else:
        raise UsageError("you should set one option: -repN or two options: -repN1, -repN2.")

    if (gn := len(kwargs["groupname"])) > 0:
        case_fn, ctrl_fn = len(kwargs["case_file"]), len(kwargs["ctrl_file"])
        if case_fn > ctrl_fn:
            cfn = case_fn
            crep = kwargs["case_rep"]
        else:
            cfn = ctrl_fn
            crep = kwargs["ctrl_rep"]
        if (len(crep) == 1) and (cfn // crep[0] != gn):
                raise UsageError("Value number for '-gn' / '--groupName' is wrong.")
        elif (len(crep) > 1) and len(crep) != gn:
                raise UsageError("Value number for '-gn' / '--groupName' is wrong.")

    ul.meta(*args, **kwargs)


@cli_fun.command(help="install R packages")
@click.option('-r', '--requirement', is_flag=True, default=False,
                help="install astk requirement R packages")
@click.option('-OrgDb', '--OrgDb', "OrgDb", cls=MultiOption, type=str, 
                default=(), help="install Genome wide annotation package")
@click.option('-cran', '--cran', cls=MultiOption, type=str, 
                default=(), help="install CRAN package")
@click.option('-bioc', '--bioconductor', cls=MultiOption, type=str, 
                default=(),help="install Bioconductor package")
@click.option('-j', '--java',  is_flag=True, help="install java software")
@click.option('-m', '--mirror',  is_flag=True, default=False,
                help="use tsinghua mirrors.")
@click.option('--conda',  is_flag=True, default=False,
                help="install packages via conda")
def install(*args, **kwargs):
    ul.install(*args, **kwargs)


@cli_fun.command(help="get motif from meme file")
@click.argument('motifId', nargs=-1, required=True)   
@click.option('-mm', '--meme', type=click.Path(exists=True), help="meme motif file")
@click.option('-db', '--database', type=click.Choice(['ATtRACT', 'CISBP-RNA']),
                help="RBP motif database")
@click.option('-org', '--organism', help="RBP organism")
@click.option('-o', '--output', required=True, help="output path")
def getmeme(*args, **kwargs):

    ul.getmeme(*args, **kwargs)

@cli_fun.command(name="list", help="list OrgDb")
@click.option('-orgdb', '--OrgDb', "OrgDb", is_flag=True, help="list OrgDb")
@click.option('-rbpsp', '--RBPSp', "RBPSp", is_flag=True, help="RNA binding protein  ")
def list_(*args, **kwargs):
    ul.list_(*args, **kwargs)


# @cli_fun.command(help = "generate ChromHMM anchor file")
# @click.argument('file', type=click.Path(exists=True), required=True)
# @click.option('-o', '--output', required=True, help="file output path")
# @click.option('-idx', '--index', type=int, help="element index")
# @click.option('-si', '--sideIndex', type=(int, int), help="the center of two side index")
# @click.option('-u', '--upstreamOffset', "offset5", type=int, default=0, help="upstream offset")
# @click.option('-d', '--downstreamOffset', "offset3", type=int, default=0, help="downstream offset")
# @click.option('-ss', '--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
# def anchor(*args, **kwargs):
#     ul.anchor(*args, **kwargs)
 

@cli_fun.command(help="retrieve coordinates from splice site region")
@click.option('-i', '--input', 'event_file', type=click.Path(exists=True),
                required=True,  help='AS event file')
@click.option('-o', '--output', required=True, help="output path")
@click.option('-si', '--siteIndex', "site_idx", type=int, default=0,
                help="splice site index. if not set, it will use all splice sites. 1-based")
@click.option('-uw', '--upstreamWidth', "ups_width", type=int, default=150,
                help="flank width of splice site upstream")
@click.option('-dw', '--downstreamWidth', "dws_width", type=int, default=150,
                help="flank width of splice site downstream")
@click.option('-ss', '--spliceSite', "sss", is_flag=True, 
                help="get splice site region window sizes")
@click.option('--interval', type=(int, int), default=(None, None),
                help="interval the between two splice sites, 1-based")
@click.option('-fi', 'fasta', type=click.Path(exists=True),
                help="Input FASTA file. if set, the fasta sequence will be extracted")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file")
def getcoor(*args, **kwargs):
    import astk.utils as ul

    output = kwargs.pop("output")
    interval_idx = kwargs.pop("interval")
    fasta = kwargs.pop("fasta")
    if kwargs["site_idx"]:
        sites = [kwargs["site_idx"] - 1]
    elif all(interval_idx):
        sites = [i-1 for i in sorted(interval_idx)]

    kwargs["site_idx"] = sites
    coord_dic = ul.get_ss_bed(*args, **kwargs)
    if len(coord_dic) == 1:
        df = list(coord_dic.values())[0]
        df.to_csv(output, sep="\t", index=False, header=False)
        print(list(coord_dic.values())[0].head())
    elif len(coord_dic) == 2:
        df = ul.get_ss_range(kwargs["event_file"], *sites, kwargs["app"])
        df.to_csv(output, sep="\t", index=False, header=False)


@cli_fun.command(help="Make the TxDb object")
@click.argument('gtf', type=click.Path(exists=True), required=True)
@click.option('-org', '--organism', required=True, help="organism")
@click.option('-o', '--output', help="file output path")
def mktxdb(*args, **kwargs):
    ul.mkTxDb(*args, **kwargs)


@cli_fun.command(help="get gene ID from AS event file")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), help="file output path")
@click.option('-u', '--unique', is_flag=True, help="only save unique gene ID")
def getgene(*args, **kwargs):
    ul.getgene(*args, **kwargs)


@cli_fun.command(name="merge", help="merge file")
@click.option('-i','--input', 'files' ,cls=MultiOption, type=click.Path(exists=True),
                required=True , help="input files")
@click.option('-o', '--output', type=click.Path(), help="output path")
@click.option('-axis', '--axis', type=click.IntRange(min=0, max=1), default=0,
                help="merge direction, 0 for row merge and 1 for column merge")
@click.option('-rmdup', '--rmdup', type=click.Choice(["index", 'all', 'content']), 
                help="remove duplicate rows")
@click.option('-rmna', '--rmna', is_flag=True, help="remove NA data")                
def sub_merge_files(*args, **kwargs):
    from astk.utils.func import merge_files
    merge_files(*args, **kwargs)


@cli_fun.command(name="regionTest", help="region Statistical tests")
@click.option('-e1', "--event1", 'file1', type=click.Path(exists=True), required=True,
                help="event1 file")
@click.option('-e2', "--event2", 'file2', type=click.Path(exists=True), required=True,
                help="event2 file")
@click.option('-bed', "--bed", type=click.Path(exists=True), required=True,
                help="bed file")
@click.option('-uw', '--upstreamWidth', "ups_width", type=int, default=150,
                help="flank width of splice site upstream")
@click.option('-dw', '--downstreamWidth', "dws_width", type=int, default=150,
                help="flank width of splice site downstream")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file")
@click.option('-o', '--output', type=click.Path(), help="output path")              
def sc_region_test(*args, **kwargs):
    import pyranges as pr
    from pandas import DataFrame
    from astk.utils._getfasta import get_coor_fa
    from astk.utils import get_ss_bed, read_fasta, detect_file_info
    from statsmodels.stats.proportion import proportions_chisquare
    import numpy as np
    import pandas as pd

    file1 = kwargs.pop("file1")
    file2 = kwargs.pop("file2")

    app = kwargs.pop("app")

    if app == "auto":
        app = detect_file_info(file1)["app"]

    # outdir = Path(outdir)
    # Path(outdir).mkdir(exist_ok=True)
    coord_dic1 = get_ss_bed(file1, app=app, **kwargs)
    coord_dic2 = get_ss_bed(file2, app=app, **kwargs)

    bedgr = pr.read_bed(kwargs.pop("bed"))

    event_ls1 = set()
    event_ls2 = set()
    df_score = DataFrame()
    for (ssn1, df1),(ssn2, df2)  in zip(coord_dic1.items(), coord_dic2.items()):
        # ss_dir = outdir / ssn
        # ss_dir.mkdir(exist_ok=True)
        # df.to_csv(ss_dir / f"{ssn}.bed", index=False, header=False, sep="\t")
        df1.columns = ["Chromosome", "Start", "End", "Name",  "Score", "Strand"]
        ss_gr1 = pr.PyRanges(df1)
        df2.columns = ["Chromosome", "Start", "End", "Name",  "Score", "Strand"]
        ss_gr2 = pr.PyRanges(df2)

        a1 = ss_gr1.intersect(bedgr)
        a2 = ss_gr2.intersect(bedgr)

        success_cnts = np.array([a1.df.shape[0], a2.df.shape[0]])
        total_cnts = np.array([df1.shape[0], df2.shape[0]])
        chi2, p_val, cont_table = proportions_chisquare(count=success_cnts, nobs=total_cnts)

        event_ls1.update(a1.df["Name"])
        event_ls2.update(a2.df["Name"])
        print(a1.df.shape[0], df1.shape[0], a2.df.shape[0], df2.shape[0], p_val)

    dpsi_df1 = pd.read_csv(file1, sep="\t", index_col=0)
    dpsi_df2 = pd.read_csv(file2, sep="\t", index_col=0)

    dpsi_df1.loc[event_ls1, ].to_csv(Path(file1).with_suffix(".peak.dpsi"), sep="\t")
    dpsi_df2.loc[event_ls2, ].to_csv(Path(file2).with_suffix(".peak.dpsi"), sep="\t")



@cli_fun.command(name="model", help="model prediction")
@click.option('-i', "--input", 'files', cls=MultiOption, type=click.Path(exists=True), 
                required=True, help="feature files")
@click.option('-cv', "--cv", 'nfold', type=int, help="cross validation fold")                
@click.option('-o', '--output', type=click.Path(), help="output path")              
def sc_model_eval(*args, **kwargs):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import cross_val_predict
    from sklearn.model_selection import  StratifiedKFold,StratifiedShuffleSplit
    from sklearn.metrics import accuracy_score
    from sklearn.metrics import matthews_corrcoef
    from sklearn.metrics import precision_recall_fscore_support
    from sklearn.metrics import confusion_matrix, multilabel_confusion_matrix
    from sklearn.model_selection import LeaveOneOut
    from sklearn.metrics import roc_curve, auc, RocCurveDisplay, plot_roc_curve
    from sklearn.preprocessing import StandardScaler, Normalizer
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.utils import shuffle    
    import pyranges as pr
    from pandas import DataFrame
    from astk.utils._getfasta import get_coor_fa
    from astk.utils import get_ss_bed, read_fasta, detect_file_info
    from statsmodels.stats.proportion import proportions_chisquare
    import numpy as np


    def _metrics(y_true, y_pred):
        cm = confusion_matrix(y_true, y_pred)
        mcm = multilabel_confusion_matrix(y_true, y_pred)
        tn = mcm[:, 0, 0]
        tp = mcm[:, 1, 1]
        fn = mcm[:, 1, 0]
        fp = mcm[:, 0, 1]
        acc = accuracy_score(y_true, y_pred)
        mcc = matthews_corrcoef(y_true, y_pred)
        p, r, f1, s = precision_recall_fscore_support(y_true, y_pred, zero_division=0)
        result = {'tn': tn, 'tp': tp, 'fn': fn, 'fp': fp, 
                    'precision': p, 'recall': r, "f1-score": f1,
                    'acc': acc, 'mcc': mcc, 'cm': cm}
        return result

    # StratifiedShuffleSplit
    # StratifiedKFold
    def cv_metrics(clf, x, y, cv):
        cver = StratifiedKFold(n_splits=cv)
        # cver = StratifiedShuffleSplit(n_splits=cv, test_size=0.2, random_state=0)
        y_pred = cross_val_predict(clf, x, y, cv=cver)
        result = {}
        for idx, (_, test_idx) in enumerate(cver.split(x, y)):
            sub_yt, sub_yp = y[test_idx], y_pred[test_idx]
            result[idx] = _metrics(sub_yt, sub_yp)
        return result

    def cv_roc_curve_plot(clf, x, y, cv):
        tprs = []
        aucs = []        
        mean_fpr = np.linspace(0, 1, len(y))
        fig, ax = plt.subplots()
        cver = StratifiedKFold(n_splits=cv)
        for i, (train_idx, test_idx) in enumerate(cver.split(x, y)):
            clf.fit(x[train_idx], y[train_idx])
            viz = plot_roc_curve(clf, x[test_idx], y[test_idx], ax=ax, name=f"ROC fold {i}", alpha=.3, lw=1)
            interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
            tprs.append(interp_tpr)
            aucs.append(viz.roc_auc)    
        ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
        mean_fpr = np.linspace(0, 1, len(y))
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        ax.plot(mean_fpr, mean_tpr, color='b', lw=2, alpha=.8,
                label=r'Mean ROC (AUC = %0.4f)' % mean_auc)
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                        label=r'$\pm$ %0.4f std. dev.' % std_auc)
        ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
            title="Receiver operating characteristic")
        ax.legend(loc="lower right")
        name =  clf.__class__.__name__
        viz = RocCurveDisplay(fpr=mean_fpr, tpr=mean_tpr,
                                roc_auc=mean_auc, estimator_name=name)
        return viz
    
    def normal_data(data):
        scaler = Normalizer()
        new_data = scaler.fit_transform(data)
        return new_data

    def save_report(metric_dic, report_file):
        with open(report_file, "w", encoding="utf-8") as f:
            head = (f"     {'tp':^5}{'fn':^5}{'fp':^5}{'tn':^5}" +
                    f"{'recall':^8}{'precision':^11}{'f1-score':^11}\n")
            for idx, sm in metric_dic.items():
                f.write(f"{idx:^50}")
                f.write('\n')
                for i, j in enumerate(sm['cm']):
                    f.write(f"{i:<4}")
                    f.write('  '.join(map(str, j.tolist())))
                    f.write('\n')
                f.write('\n')
                f.write(head)
                tp, fn, fp, tn = sm['tp'], sm['fn'], sm['fp'], sm['tn']
                ppr, recall, f1s = sm['precision'], sm['recall'], sm['f1-score']
                cls_i = "{:^5}{:^5}{:^5}{:^5}{:^5}{:^8.2f}{:^11.2f}{:^11.2f}\n"
                acc_i = "acc{:>48.2f}\n"
                mcc_i = "mcc{:>48.2f}"
                for i in range(len(tp)):
                    line = (i, tp[i], fn[i], fp[i], tn[i], recall[i], ppr[i], f1s[i])
                    f.write(cls_i.format(*line))
                f.write(acc_i.format(sm['acc']))
                f.write(mcc_i.format(sm['mcc']))
                f.write("\n")
                f.write("-"*55)
                f.write("\n")
            f.write("\n")
    df_ls = [pd.read_csv(file, index_col=0) for file in kwargs["files"]]
    for idx, sdf in enumerate(df_ls):
        sdf["label"] = idx
    df = pd.concat(df_ls)
    df.index = range(df.shape[0]) 
    X = np.array(df.iloc[:, range(df.shape[1]-1)])
    y = df.iloc[:, -1]
    X = normal_data(X)
    clf = SVC(class_weight="balanced")
    me = cv_metrics(clf, X, y, 3)        
    # me = loo_metrics(clf, X, y)        
    save_report(me, kwargs["output"])
