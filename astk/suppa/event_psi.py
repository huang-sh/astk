from pathlib import Path


def cal_psi(alter_tx, total_tx, tpm_df, tpm_th=1):
    alter_tx = [i for i in alter_tx if i in tpm_df.index]
    total_tx = [i for i in total_tx if i in tpm_df.index]
    if len(total_tx) == 0:
        psi = "nan"
    elif len(alter_tx) == 0 and len(total_tx) > 0:
        psi = 0
    else:
        alter_tx_tpm = sum([tpm_df.at[i, tpm_df.columns[0]] for i in alter_tx])
        total_tx_tpm = sum([tpm_df.at[i, tpm_df.columns[0]] for i in total_tx])
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
    """Get the PSI value of the event for the row in the ioe df
    """
    alternative_txs = row["alternative_transcripts"].split(",")
    total_txs = row["total_transcripts"].split(",")
    psi = cal_psi(alternative_txs, total_txs, tpm_df, tpm_th=tpm_th)
    return psi


def get_ioe_psi(ioe_df, tpm_df, tpm_th=1):
    """Get the PSI value of the event for each row in the ioe df
    """
    psi_ls = ioe_df.apply(ioe_df_row, axis=1, args=(tpm_df, tpm_th))
    return psi_ls.tolist()


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
