# -*- coding: utf-8 -*-
"""
astk.utils.select
~~~~~~~~~~~~~~~~~
This module provides data filter function.
"""

from pathlib import Path

from pandas import read_csv


def sig_filter(df, **kwargs):
    pval = kwargs.get("pval", 0.05)
    qval = kwargs.get("qval", 1)
    dpsi = kwargs.get("dpsi", 0)
    abs_dpsi = kwargs.get("abs_dpsi", 0)
    app = kwargs.get("app", "auto")

    if app == "SUPPA2":
        pval_col = "pval"
        dpsi_col = "dpsi"
        qval_col = "pval"
    elif app == "rMATS":
        pval_col = "FDR"
        dpsi_col = "IncLevelDifference"
        qval_col = "FDR"

    dpsi_df = df
    pval_keep = dpsi_df.loc[:, pval_col] < pval
    pval_keep = (dpsi_df.loc[:, qval_col] < qval) & pval_keep
    if dpsi > 0:
        keep = (dpsi_df.loc[:, dpsi_col] > dpsi) & pval_keep
    elif dpsi < 0:
        keep = (dpsi_df.loc[:, dpsi_col] < dpsi) & pval_keep
    else:
        keep = pval_keep
    if abs_dpsi > 0:
        keep = (abs(dpsi_df.loc[:, dpsi_col]) > abs_dpsi) & keep
    fdf = dpsi_df.loc[keep, ]
    return fdf


#TODO: generate psi files according dpsi(+,-)
class SigFilter:
    def __init__(self, dpsi_file, out, dpsi, pval, 
                abs_dpsi, psi_file, fmt) -> None:
        self.dpsi_file = dpsi_file
        self.dpsi = dpsi
        self.pval = pval
        self.abs_dpsi = abs_dpsi 
        self.psi_file = psi_file if psi_file else []
        self.fmt = fmt
        self.sep = "," if fmt == "csv" else "\t"
        self.sig_psi = []
        self.set_out(out)
 
    def filter_dpsi(self):
        dpsi_df = read_csv(self.dpsi_file, sep="\t", index_col=0)
        old_col = dpsi_df.columns
        dpsi_df.columns = ["dpsi", "pval"]
        filter_df = sig_filter(dpsi_df, dpsi=self.dpsi, abs_dpsi=self.abs_dpsi, pval=self.pval)

        if self.abs_dpsi > 0: 
            pos_df = filter_df.loc[filter_df["dpsi"] > 0, ]
            neg_df = filter_df.loc[filter_df["dpsi"] < 0, ]
            self.pos_sig_event = pos_df.index
            self.neg_sig_event = neg_df.index
            
            pos_df.columns = old_col
            neg_df.columns = old_col
            pos_df.to_csv(self.pos_sig_dpsi, index=True, sep=self.sep, index_label=False)
            neg_df.to_csv(self.neg_sig_dpsi, index=True, sep=self.sep, index_label=False)

        filter_df.columns = old_col
        filter_df.to_csv(self.sig_dpsi, index=True, sep=self.sep, index_label=False)

        self.sig_event = filter_df.index
                       
    def filter_psi(self):
        for i, pf in enumerate(self.psi_file):
            psi = read_csv(pf, sep="\t")
            psi.index = psi["event_id"]
            if len(set(self.sig_event) & set(psi.index)) > 0:
                sig_psi = psi.loc[self.sig_event, ]
                sig_psi.to_csv(self.sig_psi[i], index=False, sep="\t")
    
    def set_out(self, out):
        sig_psi = []
        if out is not None:
            Path(out).mkdir(exist_ok=True)
            sig_dpsi = Path(out) / Path(self.dpsi_file).name
            for i in self.psi_file:
                sig_psi.append(Path(out) / Path(i).name)
        else:
            self.sig_dpsi = Path(self.dpsi_file)
            for i in self.psi_file:
                sig_psi.append(Path(i))  
        self.sig_dpsi = sig_dpsi.with_suffix(".sig.dpsi")
        if self.abs_dpsi >= 0:
            self.pos_sig_dpsi = sig_dpsi.with_suffix(".sig+.dpsi")
            self.neg_sig_dpsi = sig_dpsi.with_suffix(".sig-.dpsi")
            self.pos_sig_psi = []
            self.neg_sig_psi = []
            for i in sig_psi:
                self.sig_psi.append(i.with_suffix(".sig.psi"))
                self.pos_sig_psi.append(i.with_suffix(".sig+.psi"))
                self.neg_sig_psi.append(i.with_suffix(".sig-.psi"))
        if abs(self.dpsi) >= 0:
            for i in sig_psi:
                self.sig_psi.append(i.with_suffix(".sig.psi"))
                
    def run(self):
        self.filter_dpsi()
        self.filter_psi()


class PsiFilter:

    def __init__(self, file, output, threshold, quantile) -> None:
        self.file = file
        self.output = output
        self.threshold = threshold
        self.quantile = quantile
        self.psi_df = read_csv(file, sep="\t", index_col=0).dropna()
        self.mean_psi = self.psi_df.apply(lambda row: sum(row)/len(row), axis=1)   

    def value_filter(self):
        th = self.threshold
        mp = self.mean_psi
        if th >= 0:
            filter_idx = mp[mp > th].index
        else:
            filter_idx = mp[mp < abs(th)].index
        return self.psi_df.loc[filter_idx, ]
    
    def quantile_filter(self):
        qt = self.quantile
        mp = self.mean_psi
        th = mp.quantile(abs(qt))
        if qt >= 0:
            filter_idx = mp[mp > th].index
        else:
            filter_idx = mp[mp <th].index
        return self.psi_df.loc[filter_idx, ]
    
    def run(self):
        if self.threshold:
            psi_df = self.value_filter()
        elif self.quantile:
            psi_df = self.quantile_filter()
        else:
            psi_df = self.psi_df    
        psi_df.to_csv(self.output, sep="\t", index_label=False)
        