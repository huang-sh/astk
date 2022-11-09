from pathlib import Path


def cal_psi(alter_tx, total_tx, tpm_df, tpm_th=1):
    alter_tx = [i for i in alter_tx if i in tpm_df.index]
    total_tx = [i for i in total_tx if i in tpm_df.index]
    if len(total_tx) == 0 or len(alter_tx) == 0:
        psi = "nan"
    else:
        alter_tx_tpm = sum([tpm_df.at[i, tpm_df.columns[0]] for i in alter_tx])
        total_tx_tpms = [tpm_df.at[i, tpm_df.columns[0]] for i in total_tx]
        if sum(total_tx_tpms) / len(total_tx_tpms) < tpm_th:
            psi = "nan"
        elif alter_tx_tpm <= 0.0001:
            psi = 0
        else:
            try:
                psi = alter_tx_tpm / sum(total_tx_tpms)
            except ZeroDivisionError:
                psi = "nan"
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
