import hashlib
import pickle
from pathlib import Path

from .lib import gtf_parse as gp
from .lib.AS_event import make_events


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
    


def generate_events(gtf, event_types, output, promoterSplit):
    import time
    T1 = time.time()

    gtf_pkl = gtf_parse_cache(gtf)
    
    T21 =time.time()
    print('程序运行时间:%s秒' % ((T21 - T1)))

    # if gtf_pkl.exists():
    #     print("Loading gtf parsing cache data.")
    #     with gtf_pkl.open("rb") as f:
    #         genome = pickle.load(f)
    #     T2 =time.time()
    #     print('程序运行时间:%s秒' % ((T2 - T21)))    
    # else:
    #     genome = gp.construct_genome(gtf)
    #     with gtf_pkl.open("wb") as f:
    #         pickle.dump(genome, f)
    genome = gp.construct_genome(gtf)
    T2 =time.time()
    print('程序运行时间:%s秒' % ((T2 - T21)))
    
    if event_types == "ALL":
        event_types = ["SE", "RI"]
    else:
        event_types = [event_types]

    make_events(event_types, genome, output, promoterSplit)

    T3 =time.time()
    print('程序运行时间:%s秒' % ((T3 - T1)))
    
def cal_psi(alter_tx, total_tx, tpm_df, tpm_th=1):
    alter_tx = [i for i in alter_tx if i in tpm_df.index]
    total_tx = [i for i in total_tx if i in tpm_df.index]  
    if len(total_tx) == 0:
        psi = "nan"
    elif len(alter_tx) == 0 and len(total_tx) > 0:
        psi = 0
    else:
        alter_tx_tpm = sum([tpm_df.at[i, 'TPM'] for i in alter_tx])
        total_tx_tpm = sum([tpm_df.at[i, 'TPM'] for i in total_tx])
        if total_tx_tpm <= tpm_th:
            psi = "nan"
        elif alter_tx_tpm <= 0.0001:
            psi = 0
        else:
            try: 
                psi = alter_tx_tpm / total_tx_tpm
            except ZeroDivisionError:
                psi = 0
    return psi


def ioe_df_row(row, tpm_df, tpm_th):
    alternative_txs = row["alternative_transcripts"].split(",")
    total_txs = row["total_transcripts"].split(",")
    psi = cal_psi(alternative_txs, total_txs, tpm_df, tpm_th=tpm_th)
    return psi


def get_ioe_psi(ioe_df, tpm_df, tpm_th=1):
    psi_ls = ioe_df.apply(ioe_df_row, axis=1, args=(tpm_df, tpm_th))
    return psi_ls


def read_tpm(tpm_file, tpm_col=2, tx_col=0):
    import pandas as pd

    tpm_df = pd.read_csv(tpm_file, sep = "\t", index_col=tx_col)
    cols = list(tpm_df.columns)
    cols[tpm_col] = "TPM"
    tpm_df.columns = cols
    return tpm_df


def calculate_psi(ioe, tpm_files, tpm_th, tpm_col, tx_col):
    import pandas as pd
    ioe_df = pd.read_csv(ioe, sep="\t")

    psi_dic = {}
    for tf in tpm_files:
        sf_name = Path(tf).parent.name  # only consider salmon output, now
        tpm_df = read_tpm(tf, tpm_col=tpm_col, tx_col=tx_col)
        psi_dic[sf_name] = get_ioe_psi(ioe_df, tpm_df, tpm_th=tpm_th)
    
    psi_df = pd.DataFrame(psi_dic)
    psi_df.index = ioe_df["event_id"]
    return psi_df
