from astk.lazy_loader import LazyLoader

np = LazyLoader("np", globals(), "numpy")
plt = LazyLoader("plt", globals(), "matplotlib.pyplot")
skmt = LazyLoader("metrics", globals(), "sklearn.metrics")
skms = LazyLoader("model_selection", globals(), "sklearn.model_selection")


def _metrics(y_true, y_pred):
    cm = skmt.confusion_matrix(y_true, y_pred)
    mcm = skmt.multilabel_confusion_matrix(y_true, y_pred)
    tn = mcm[:, 0, 0]
    tp = mcm[:, 1, 1]
    fn = mcm[:, 1, 0]
    fp = mcm[:, 0, 1]
    acc = skmt.accuracy_score(y_true, y_pred)
    mcc = skmt.matthews_corrcoef(y_true, y_pred)
    p, r, f1, s = skmt.precision_recall_fscore_support(y_true, y_pred, zero_division=0)
    return {
        'tn': tn,
        'tp': tp,
        'fn': fn,
        'fp': fp,
        'precision': p,
        'recall': r,
        "f1-score": f1,
        'acc': acc,
        'mcc': mcc,
        'cm': cm,
    }


def cv_metrics(clf, x, y, cv):
    cver = skms.StratifiedKFold(n_splits=cv)
    # cver = StratifiedShuffleSplit(n_splits=cv, test_size=0.2, random_state=0)
    y_pred = skms.cross_val_predict(clf, x, y, cv=cver)
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
    cver = skms.StratifiedKFold(n_splits=cv)
    for i, (train_idx, test_idx) in enumerate(cver.split(x, y)):
        clf.fit(x[train_idx], y[train_idx])

        viz= skmt.RocCurveDisplay.from_estimator(
                clf, x[test_idx], y[test_idx], 
                ax=ax, name=f"ROC fold {i}", alpha=.3, lw=1
        )

        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)    
    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
    mean_fpr = np.linspace(0, 1, len(y))
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = skmt.auc(mean_fpr, mean_tpr)
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
    viz = skmt.RocCurveDisplay(fpr=mean_fpr, tpr=mean_tpr,
                            roc_auc=mean_auc, estimator_name=name)
    return viz
