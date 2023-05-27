"""
astk.cli.ml
~~~~~~~~~~~~~~~~~
This module provide machine learning functions
"""

from click_option_group import optgroup

from .config import *
from astk.ml import *
from astk.lazy_loader import LazyLoader

np = LazyLoader("np", globals(), "numpy")
pd = LazyLoader("pd", globals(), "pandas")


@cli_fun.command(name="hpo", help="Hpyper-parameter optimization")
@click.option('-i', "--input", 'files', cls=MultiOption, type=click.Path(exists=True), 
                required=True, help="feature file for hpyper-parameter optimization")
@click.option('-clf', '--classifier', type=click.Choice(['SVM', 'RF', "KNN"]),
                help="classifier selection, default=SVM")
@click.option('-cv', "--cv", type=int, default=5, help="cross validation fold")
@click.option('-p', "--process", type=int, default=4, help="process number")  
@optgroup.group("SVM classifier parameters")
@optgroup.option('-C', '--C', "C", type=(float, float, int), default=(-5, 5, 11), show_default=True,
                 help="regularization parameter range")
@optgroup.option('-g', '--gamma', type=(float, float, int), default=(-5, 3, 9), show_default=True,
                 help="kernel coefficient range")
@optgroup.option('-k', "--kernel", cls=GroupedMultiOption, default=("linear", "rbf"),
                 type=click.Choice(["linear", "poly", "rbf","sigmoid"]), show_default=True, 
                 help="kernel type to be used in SVM")
@optgroup.group("Random Forest classifier parameters")
@optgroup.option("--n-estimators", type=(int, int, int),
                 help="the number of trees in the forest.")
@optgroup.option('--criterion', cls=GroupedMultiOption, default=("gini", "entropy"), 
                    type=click.Choice(["gini", "entropy", "log_loss"]), show_default=True,
                 help="the criterion to measure the quality of a split")
@optgroup.option('--max-features', cls=GroupedMultiOption, type=click.Choice(["sqrt", "log2"]), 
                 default=("sqrt", "log2"), show_default=True,
                 help="the number of features to consider when looking for the best split")
@optgroup.option('--max-depth', type=(int, int, int), help="the maximum depth of the tree.")
@optgroup.group(" K-nearest neighbors classifier parameters")
@optgroup.option("--n-neighbors", type=(int, int, int), help="number of neighbors to use")
@optgroup.option('--weights', cls=GroupedMultiOption, default=("uniform", "distance"), 
                type=click.Choice(["uniform", "distance"]), show_default=True,
                help="weight function used in prediction.")
@optgroup.option('--algorithm', cls=GroupedMultiOption, 
                 default=("auto", "ball_tree", "kd_tree", "brute"), show_default=True, 
                 help="algorithm used to compute the nearest neighbors")
@optgroup.option('--leaf-size', type=(int, int, int), help="Leaf size passed to BallTree or KDTree.")
@optgroup.option('--metric', cls=GroupedMultiOption, 
                 type=click.Choice(["minkowski", "cityblock", "euclidean", "manhattan"]), 
                 help="metric to use for distance computation. ")
def sc_hpo(*args, **kwargs):
    from sklearn.svm import SVC
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.ensemble import RandomForestClassifier

    X, y = load_feature(kwargs["files"])    
    if kwargs["classifier"] == "SVM":
        param_grid = {
            "C": np.logspace(*kwargs["C"], base=2),
            "gamma": np.logspace(*kwargs["gamma"], base=2),
            "kernel": kwargs["kernel"]
        }
        X = normal_data(X)
        clf = SVC(random_state=42, class_weight="balanced")
    elif kwargs["classifier"] == "RF":
        param_grid = {
            "n_estimators":[int(x) for x in np.linspace(*kwargs["n_estimators"])],
            "max_depth": [int(x) for x in np.linspace(*kwargs["max_depth"])],
            "criterion": kwargs["criterion"],
            "max_features": kwargs["max_features"]
        }
        clf = RandomForestClassifier(random_state=42, class_weight="balanced")
    elif kwargs["classifier"] == "KNN":
        param_grid = {
            "n_neighbors": [int(x) for x in np.linspace(*kwargs["n_neighbors"])],
            "weights": kwargs["weights"],
            "algorithm": kwargs["algorithm"],
            "leaf_size": [int(x) for x in np.linspace(*kwargs["leaf_size"])],
            "metric": kwargs["metric"]
        }
        X = normal_data(X)
        clf = KNeighborsClassifier()
    bps = grid_search(clf, X, y, param_grid, kwargs["cv"],  kwargs["process"])
    for k, v in bps.items():
        print(f"{k}: {v}")


@cli_fun.command(name="eval", help="Model evaluation")
@click.option('-i', "--input", 'files', cls=MultiOption, type=click.Path(exists=True), 
                required=True, help="feature file for hpyper-parameter optimization")
@click.option('-o', "--output", type=click.Path(), required=True, help="output path")
@clf_common_options
def sc_eval(*args, **kwargs):
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.neighbors import KNeighborsClassifier

    X, y = load_feature(kwargs["files"])    
    if kwargs["classifier"] in ("SVM", "KNN"):
        X = normal_data(X)
    clf = choose_clf(**kwargs)
    output = Path(kwargs["output"])
    me = cv_metrics(clf, X, y, kwargs["cv"])
    save_report(me, output)
    viz = cv_roc_curve_plot(clf, X, y, kwargs["cv"])
    outfig = output.with_suffix(".png")
    plt.savefig(outfig, dpi=1000)


@cli_fun.command(name="fs", help="Feature importance score")
@click.option('-i', "--input", 'files', cls=MultiOption, type=click.Path(exists=True), 
                required=True, help="feature files")
@click.option('-o', "--output", type=click.Path(), required=True, help="output path")
@click.option('-m', '--method', type=click.Choice(["permutation", "tree"]), default="permutation",
                show_default=True, help="feature importance compute method")
@clf_common_options
@optgroup.group("Permutation feature importance parameters")
@optgroup.option('-scoring', "--scoring", type=click.Choice(["accuracy", "roc_auc", "f1"]), 
                   default="accuracy", show_default=True, help="scorer to use")
@optgroup.option('-repeatN', '--repeatN', type=int, default=5, show_default=True, 
                   help="number of times to permute a feature.")
@optgroup.option('-gf', "--group-feature", is_flag=True, default=False, show_default=True, 
                help="group features for permutation")
@fig_common_options()
@optgroup.group("Important feature saving")
@optgroup.option('-topN', "--topN", type=int, default=0, show_default=True, 
                   help="the top number feature to save")
def sc_fs(*args, **kwargs):
    from itertools import groupby
    from sklearn.model_selection import StratifiedKFold
    from sklearn.ensemble import RandomForestClassifier

    X, y, colnames = load_feature(kwargs["files"], colname=True)
    if kwargs["group_feature"]:
        col_idx_ls = []
        feature_names = []
        group_list = groupby([i.split("_")[0] for i in colnames])
        idx_num = 0
        for label, group in group_list:
            feature_names.append(label)
            group = list(group)            
            col_idx_ls.append([idx_num+i for i, _ in enumerate(group)])
            idx_num += len(group)
    else:
        feature_names = colnames
        col_idx_ls = range(len(colnames))

    if kwargs["method"] == "tree":
        kwargs["classifier"] == "RF"        
    if kwargs["classifier"] in ("SVM", "KNN"):
        X = normal_data(X)
    clf = choose_clf(**kwargs)
    cver = StratifiedKFold(n_splits=kwargs["cv"])
    importance_ls = []

    for train_idx, test_idx in cver.split(X, y):
        clf.fit(X[train_idx], y[train_idx])
        if kwargs["method"] == "tree":
            importances = clf.feature_importances_
        else:
            result = permutation_importance(
                clf, 
                X[test_idx], 
                y[test_idx], 
                col_idx_ls=col_idx_ls,
                n_repeats=kwargs["repeatn"], 
                scoring=kwargs["scoring"],
                n_jobs=kwargs["process"],
                random_state=42
            )
            importances = result.importances_mean
        importance_ls.append(importances)
    importances = pd.Series(np.mean(np.array(importance_ls), axis=0), index=feature_names)
    fig, ax = plt.subplots(figsize=(kwargs["width"], kwargs["height"]))
    importances.plot.barh(yerr=np.std(importance_ls, axis=0), ax=ax)
    ax.set_title("importance score")
    ax.set_ylabel("feature")
    fig.tight_layout()
    outfig = Path(kwargs["output"])
    outfig.with_suffix(f".{kwargs['figfmt']}")
    plt.savefig(outfig)

    if (topn := kwargs["topn"]) > 0:
        top_features = importances.sort_values(ascending=False)[:topn].index
        for file in kwargs["files"]:
            sdf = pd.read_csv(file, index_col=0) 
            suffix = Path(file).suffix
            out = Path(file).with_suffix(f".top{topn}{suffix}")
            sdf.loc[:, top_features].to_csv(out)
