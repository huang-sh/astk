import hashlib
import sys
import json
from pathlib import Path

from astk.types import *
import astk.utils.select as sl
from .AS_event import make_events
from .event_psi import get_ioe_psi
from .gtf_parse import construct_genome
from astk.suppa.lib.diff_tools import multiple_conditions_analysis as mca


def check_gtf_used(gtf):
    gtf = Path(gtf)
    gtf_hash = hashlib.blake2b()
    gtf_size = gtf.stat().st_size
    gtf_hash.update(str(gtf_size).encode('utf-8'))
    f = gtf.open("rb")
    for _ in range(10):
        gtf_hash.update(f.readline().strip())
    hash_val = gtf_hash.hexdigest()
    return hash_val


def gtf_parse_cache(gtf):
    home = Path.home()
    astk_cache_dir = home / f".astk"
    astk_cache_dir.mkdir(exist_ok=True)
    gtf_hash = check_gtf_used(gtf)
    gtf_cache_path = astk_cache_dir / "gtfParse"
    gtf_cache_path.mkdir(exist_ok=True)
    pkl = (gtf_cache_path / gtf_hash).with_suffix(".pkl")
    return pkl


def generate_events(gtf, event_types, output, idtype, event_pos):
    genome = construct_genome(gtf)
    if event_types == "ALL":
        event_types =  ['SE', "A5", "A3", "MX", "RI", 'AF', 'AL']
    else:
        event_types = [event_types]
    make_events(output, genome, event_types, idtype, event_pos)


def diff_splice(
    psi_files: Sequence[FilePath],
    exp_files: Sequence[FilePath],
    reference: FilePath,
    output: FilePath,
    method: str
):
    if str(output).endswith(".dpsi"):
        output = str(output)[:-5]
    mca(method, psi_files, exp_files, reference, 1000, 0, 
        False, True, 0.05, True, False, False, 0, 0, str(output))                 


def read_tpm(file, colname, tpm_col):
    import pandas as pd

    quant_df = pd.read_csv(file, sep = "\t")
    tpm_df = pd.DataFrame({colname: quant_df.iloc[:, tpm_col-1]})
    tpm_df.index = quant_df.iloc[:, 0]
    return tpm_df


def calculate_psi(ioe, tpm_files, tpm_th, tpm_col, tx_col):
    import pandas as pd
    ioe_df = pd.read_csv(ioe, sep="\t")

    psi_dic = {}
    tpm_ls = []
    for tf in tpm_files:
        sf_name = Path(tf).parent.name  # only consider salmon output, now
        tpm_df = read_tpm(tf, sf_name, tpm_col)
        psi_dic[sf_name] = get_ioe_psi(ioe_df, tpm_df, tpm_th=tpm_th)
        tpm_ls.append(tpm_df)
    tpm_df = pd.concat(tpm_ls, axis=1)
    psi_df = pd.DataFrame(psi_dic)
    psi_df.index = ioe_df["event_id"]
    return psi_df, tpm_df


def parse_meta(meta_file: FilePath) -> Dict:

    tpm_dic = {}
    with open(meta_file, "r") as f:
        try:
            meta_dic = json.load(f)
        except json.decoder.JSONDecodeError:
            print("ERROR: JSON file is excepted!")
            sys.exit()
    for gn, gdic in meta_dic.items():
        control_tpms = [(sp["path"], sp["name"]) for sp in gdic["ctrl"]["samples"]]
        treatment_tpms = [(sp["path"], sp["name"]) for sp in gdic["case"]["samples"]]
        tpm_dic[gn] = {"ctrl": control_tpms, "case": treatment_tpms}
    return tpm_dic


def ds_flow(
    meta: FilePath,
    gtf: FilePath,
    etypes: Sequence[str],
    outdir: FilePath, 
    method: str,
    id_type: str,
    pval: float,
    abs_dpsi: float,
    tpm_threshold: float,
    event_pos: str
    ):
    import pandas as pd

    Path(outdir).mkdir(exist_ok=True)
    tpm_dic = parse_meta(meta)
    genome = construct_genome(gtf)
    ref_dir = Path(outdir) / "ref"
    ref_dir.mkdir(exist_ok=True)
    make_events(ref_dir / "annotation", genome, etypes, id_type, event_pos)
    
    tpm_dir =  Path(outdir) / "tpm"
    psi_dir =  Path(outdir) / "psi"
    dpsi_dir =  Path(outdir) / "dpsi"
    sig_dir =  Path(outdir) / ("sig" + str(abs_dpsi).replace(".", ''))
    tpm_dir.mkdir(exist_ok=True)
    psi_dir.mkdir(exist_ok=True)
    dpsi_dir.mkdir(exist_ok=True)
    sig_dir.mkdir(exist_ok=True)
    sig_psi_dir = sig_dir / "psi"
    sig_psi_dir.mkdir(exist_ok=True)

    for gn, gn_dic in tpm_dic.items():
        for et in etypes:
            if event_pos is None:
                es = ""
            else:
                es = f"_{event_pos}"
            ioe = ref_dir / f"annotation{es}_{et}_strict.ioe"
            ioe_df = pd.read_csv(ioe, sep="\t")
            ctrl_tpm_ls = []
            case_tpm_ls = []
            ctrl_psi_ls = []
            case_psi_ls = []
            for (cf, cn), (tf, tn) in zip(gn_dic["ctrl"], gn_dic["case"]):
                cdf = read_tpm(cf, cn, 4)
                tdf = read_tpm(tf, tn, 4)
                ctrl_tpm_ls.append(cdf)
                case_tpm_ls.append(tdf)

                cpsi_df = pd.DataFrame(get_ioe_psi(ioe_df, cdf, tpm_th=tpm_threshold), 
                        columns=cdf.columns, index=ioe_df["event_id"])
                tpsi_df = pd.DataFrame(get_ioe_psi(ioe_df, tdf, tpm_th=tpm_threshold), 
                        columns=tdf.columns, index=ioe_df["event_id"])
                ctrl_psi_ls.append(cpsi_df)
                case_psi_ls.append(tpsi_df)

            ctrl_tpm = pd.concat(ctrl_tpm_ls, axis=1)
            case_tpm = pd.concat(case_tpm_ls, axis=1)
            ctrl_tpm_file = tpm_dir / f"{gn}_c1.tpm"
            case_tpm_file = tpm_dir / f"{gn}_c2.tpm"
            ctrl_tpm.to_csv(ctrl_tpm_file, sep="\t", index_label=False)
            case_tpm.to_csv(case_tpm_file, sep="\t", index_label=False)

            ctrl_psi = pd.concat(ctrl_psi_ls, axis=1)
            case_psi = pd.concat(case_psi_ls, axis=1)
            ctrl_psi_file = psi_dir / f"{gn}_{et}_c1.psi"
            case_psi_file = psi_dir / f"{gn}_{et}_c2.psi"
            ctrl_psi.to_csv(ctrl_psi_file, sep="\t")
            case_psi.to_csv(case_psi_file, sep="\t")
            psi_files = [ctrl_psi_file, case_psi_file]
            exp_files = [ctrl_tpm_file, case_tpm_file]
            dpsi_out = dpsi_dir / f"{gn}_{et}"
            mca(method, psi_files, exp_files, ioe, 1000, 0, 
                False, True, 0.05, True, False, False, tpm_threshold, 0, str(dpsi_out))
            
            sig_dpsi_out = (sig_dir / dpsi_out.name).with_suffix(".sig.dpsi")
            sig_pos_dpsi_out = (sig_dir / dpsi_out.name).with_suffix(".sig+.dpsi")
            sig_neg_dpsi_out = (sig_dir / dpsi_out.name).with_suffix(".sig-.dpsi")
            kwargs = {"abs_dpsi": abs_dpsi, "pval": pval, "app": "SUPPA2"}
            dpsi_df = pd.read_csv(dpsi_out.with_suffix(".dpsi"), sep="\t", index_col=0).dropna()
            old_col = dpsi_df.columns
            dpsi_df.columns = ["dpsi", "pval"]
            df_fil = sl.sig_filter(dpsi_df, **kwargs)
            df_fil.columns = old_col      
            df_fil.to_csv(sig_dpsi_out, sep="\t")

            pos_df = df_fil.loc[df_fil[old_col[0]] > 0, ]
            neg_df = df_fil.loc[df_fil[old_col[0]] < 0, ]
            pos_df.to_csv(sig_pos_dpsi_out, sep="\t")
            neg_df.to_csv(sig_neg_dpsi_out, sep="\t")

            ctrl_psi_file = sig_psi_dir / f"{gn}_{et}_c1.sig.psi"
            case_psi_file = sig_psi_dir / f"{gn}_{et}_c2.sig.psi"
            ctrl_psi.loc[df_fil.index, ].to_csv(ctrl_psi_file, sep="\t")
            case_psi.loc[df_fil.index, ].to_csv(case_psi_file, sep="\t")
