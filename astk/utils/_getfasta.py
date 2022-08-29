

def get_coor_fa(df, fasta, out, strandedness=False):
    import pybedtools

    if df.shape[1] == 6:
        lines = [f"{v[1]}\t{v[2]}\t{v[3]}\t{v[4]}\t{v[5]}\t{v[6]}" for v in df.itertuples()]
    else:
        lines = [f"{v[1]}\t{v[2]}" for v in df.itertuples()]
    bed = pybedtools.BedTool("\n".join(lines), from_string=True)
    bed.sequence(fi=fasta, fo=out, name=True, s=strandedness)
    return bed.seqfn
