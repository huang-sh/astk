

def get_coor_fa(df, fasta, output, strandedness=False, rna=False):
    from pyfaidx import Fasta

    genes = Fasta(fasta)
    with open(output, "w") as f:
        for row in df.itertuples():
            seq = genes[row[1]][row[2]:row[3]]
            if strandedness and row[6] == "-":
                seq = -seq
            if rna:
                seq.seq = seq.seq.replace("T", "U")
                seq.seq = seq.seq.replace("t", "u")
            f.write(f">{seq.fancy_name}\n")
            f.write(seq.seq + "\n")
    return output
