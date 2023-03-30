from astk.lazy_loader import LazyLoader

pd = LazyLoader("pd", globals(), "pandas")
skms = LazyLoader("model_selection", globals(), "sklearn.model_selection")
skpp = LazyLoader( "preprocessing", globals(), "sklearn.preprocessing")


def grid_search(clf, x, y, param_grid, cv, n_jobs):
    grid = skms.GridSearchCV(clf, cv=cv, n_jobs=n_jobs, param_grid=param_grid)
    clf = grid.fit(x, y)
    print(clf.best_score_)
    return clf.best_params_


def normal_data(data):
    scaler = skpp.Normalizer()
    return scaler.fit_transform(data)


def save_report(metric_dic, report_file):
    acc_ls = []
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
            acc_ls.append(sm['acc'])
            f.write("\n")
            f.write("-"*55)
            f.write("\n")
        f.write("\n")
        mean_acc = sum(acc_ls) / len(acc_ls)
        f.write(f"mean acc{mean_acc:>43.4f}\n")


def load_feature(files, colname=False):
    df_ls = [pd.read_csv(file, index_col=0) for file in files]
    for idx, sdf in enumerate(df_ls):
        sdf["label"] = idx
    df = pd.concat(df_ls)
    df.index = range(df.shape[0])
    x = df.iloc[:, range(df.shape[1]-1)].to_numpy()
    y = df.iloc[:, -1]
    return (x, y, df.columns[:-1]) if colname else (x, y)


def choose_clf(**kwargs):
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.neighbors import KNeighborsClassifier

    if kwargs["classifier"] == "SVM":
        clf = SVC(C=kwargs["C"], 
                  gamma=kwargs["gamma"], 
                  kernel=kwargs["kernel"], 
                  class_weight="balanced", 
                  random_state=42)
    elif kwargs["classifier"] == "RF":
        clf = RandomForestClassifier(
                n_estimators=kwargs["n_estimators"], 
                max_depth=kwargs["max_depth"],
                max_features=kwargs["max_features"],
                criterion=kwargs["criterion"],
                class_weight="balanced_subsample", 
                random_state=42)
    elif kwargs["classifier"] == "KNN":
        clf = KNeighborsClassifier(
                n_neighbors=kwargs["n_neighbors"],
                weights=kwargs["weights"],
                algorithm=kwargs["algorithm"],
                leaf_size=kwargs["leaf_size"],
                metric=kwargs["metric"])
    return clf


def plot_pca(output, data, label):
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Preprocess the data by scaling it to have zero mean and unit variance
    scaler = StandardScaler()
    X = scaler.fit_transform(data)

    # Perform PCA on the preprocessed data
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)    

    # Create a new DataFrame with the principal components and the labels
    pc_df = pd.DataFrame(data=pcs, columns=['PC1', 'PC2'])
    pc_df['label'] = label

    # Plot the results
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize=15)
    ax.set_ylabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize=15)
    ax.set_title('PCA', fontsize=20)

    for target in list(dict.fromkeys(label)):
        indices = pc_df['label'] == target
        ax.scatter(
                pc_df.loc[indices, 'PC1'],
                pc_df.loc[indices, 'PC2'],
                label=target)
    ax.legend()
    plt.savefig(output)
