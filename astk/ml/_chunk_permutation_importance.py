

def permutation_importance(
    estimator,
    X,
    y,
    *,
    col_idx_ls=None,
    scoring=None,
    n_repeats=5,
    n_jobs=None,
    random_state=None,
    sample_weight=None,
    max_samples=1.0,
):
    """Permutation importance for feature evaluation [BRE]_.
    refer to https://github.com/scikit-learn/scikit-learn/blob/9aaed4987/sklearn/inspection/_permutation_importance.py#L102
    
    make little modification for multiple columns index supporting
    """
    import numbers
    import numpy as np
    import pandas as pd
    from sklearn.inspection._permutation_importance import  _weights_scorer
    from sklearn.inspection._permutation_importance import _create_importances_bunch
    from sklearn.inspection._permutation_importance import  _calculate_permutation_scores
    from sklearn.metrics._scorer import _check_multimetric_scoring, _MultimetricScorer
    from sklearn.utils import check_random_state, check_array
    from sklearn.utils.parallel import delayed, Parallel
    from sklearn.metrics import check_scoring
    
    if not hasattr(X, "iloc"):
        X = check_array(X, force_all_finite="allow-nan", dtype=None)
        X = pd.DataFrame(X)
    random_state = check_random_state(random_state)
    random_seed = random_state.randint(np.iinfo(np.int32).max + 1)

    if not isinstance(max_samples, numbers.Integral):
        max_samples = int(max_samples * X.shape[0])
    elif not (0 < max_samples <= X.shape[0]):
        raise ValueError("max_samples must be in (0, n_samples]")

    if callable(scoring):
        scorer = scoring
    elif scoring is None or isinstance(scoring, str):
        scorer = check_scoring(estimator, scoring=scoring)
    else:
        scorers_dict = _check_multimetric_scoring(estimator, scoring)
        scorer = _MultimetricScorer(scorers=scorers_dict)

    baseline_score = _weights_scorer(scorer, estimator, X, y, sample_weight)

    if col_idx_ls is None:
        col_idx_ls = range(X.shape[1])

    scores = Parallel(n_jobs=n_jobs)(
        delayed(_calculate_permutation_scores)(
            estimator,
            X,
            y,
            sample_weight,
            col_idx,
            random_seed,
            n_repeats,
            scorer,
            max_samples,
        )
        for col_idx in col_idx_ls
    )

    if isinstance(baseline_score, dict):
        return {
            name: _create_importances_bunch(
                baseline_score[name],
                # unpack the permuted scores
                np.array([scores[col_idx][name] for col_idx in range(X.shape[1])]),
            )
            for name in baseline_score
        }
    else:
        return _create_importances_bunch(baseline_score, np.array(scores))
