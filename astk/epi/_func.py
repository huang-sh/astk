import os
import sys
import json
import shutil
import subprocess
from pathlib import Path
from functools import partial
from typing import Sequence
from tempfile import NamedTemporaryFile

import astk.utils.func  as ul
from astk.constant import *
from astk.ctypes import *
from astk.lazy_loader import LazyLoader

pd = LazyLoader("pd", globals(), "pandas")


def site_flanking(chrN, site, sam, control_sam=None, window=150, bins=15):
    import pandas as pd
    import pysam

    def signal_func(chrN, start, end, sam, control_sam):
        if control_sam:
            treat = sam.count(chrN, start, end)
            control = control_sam.count(chrN, start, end)
            signal = treat / max(control, 0.9)
        else:
            signal = 1e6*sam.count(chrN, start, end)/sam.mapped
        return signal
    
    psignal_func = partial(signal_func, sam=sam, control_sam=control_sam)

    up_bin_start_idx = sorted(list(range(site+1, site-window, -bins)))[:-1]
    down_bin_end_idx = list(range(site+1, site+window+2, bins))[:-1]

    up_signal = [psignal_func(chrN, i-1, i+bins-1) for i in up_bin_start_idx]
    down_signal = [psignal_func(chrN, i-1, i+bins-1) for i in down_bin_end_idx]
    df = pd.concat([
        pd.DataFrame({"signal": up_signal, "direction": "upstream"}),
        pd.DataFrame({"signal": down_signal, "direction": "downstream"})
    ])
    return df 

#TO-DO make it faster
def epi_signal(out, achor_dic, bam_meta, width, binsize):
    import pandas as pd
    import pysam    

    df = pd.read_csv(bam_meta)
    cds = set(df.condition)
    df_ls = []

    if len(cds) == 1:
        tdf_iter = [sdf for _,sdf in df.iterrows()]
        cdf_iter = [None] * df.shape[0]
    elif cds == {'control', 'treatment'}:
        tdf = df.loc[df.condition == "treatment", ]
        cdf = df.loc[df.condition == "control", ]
        if tdf.shape != cdf.shape:
            raise ValueError("control input number dismatched treatment input")
        tdf_iter = [sdf for _,sdf in tdf.iterrows()]
        cdf_iter = [sdf for _,sdf in cdf.iterrows()]
    for c, t in zip(cdf_iter, tdf_iter):
        csam = None if c is None else pysam.AlignmentFile(c["path"])
        tsam = pysam.AlignmentFile(t["path"])
        for an in achor_dic:
            anchor_df = pd.read_csv(achor_dic[an], sep="\t", header=None)
            for ai, site_row in anchor_df.iterrows():
                chrN, pos = site_row[0], int(site_row[1])
                sdf = site_flanking(chrN, pos, tsam, control_sam=csam, window=width, bins=binsize)
                sdf["event_idx"] = ai
                sdf["mark"] = t["name"]
                sdf["rep"] = t["replicate"]
                sdf["anchor"] = an
                df_ls.append(sdf)
    dfs = pd.concat(df_ls)
    dfs["bin"] = dfs.index
    dfs.to_csv(out, index=False)


def mark(output, celltype, bed, marknum, sep, markindex, markname, stacked):
    from itertools import chain
    import pandas as pd

    if all([sep, markindex]):
        marks = [Path(i).stem.split(sep)[markindex-1] for i in bed]
    else:
        marks = markname

    if len(marks) != len(bed):
        print("mark number dismatchs with bed files ")
        sys.exit(1)
    if len(celltype) == len(marknum):
        cells = [i for i in chain(*[[celltype[i]]*int(marknum[i]) for i in range(len(celltype))])]
    elif len(celltype) == 1:
        cells = [celltype[0] for _ in bed]
    else:
        print("ERROR: -ct/--cellType and -mn/--markNum is not dismatched")
        sys.exit(1)
     
    if len(cells) != len(marks):
        print("-mn/--markNum is wrong")
        sys.exit(1)

    if stacked:
        marks = [f"{cells[idx]}_{marks[idx]}" for idx in range(len(marks))]
        cells = [f"genome" for _ in range(len(marks))]
    df = pd.DataFrame({
        "cell": cells,
        "mark": marks,
        "bed": bed
    })

    df.to_csv(output, sep="\t", index=False, header=False)



CHROMSIZES_dir = Path(__file__).parent / "ChromHMM/CHROMSIZES"
genomes = [i.stem for i in CHROMSIZES_dir.glob("*.txt")]

def LearnState(numstates, markfile, directory, binarydir , outdir, binsize, genome, 
                mx, processor, anchordir, coordir, no_binary, name, stacked, nostrand,defaultcoor):
    ChromHMM_dir = Path(__file__).parent / "ChromHMM"
    ChromHMM_jar = f"java -mx{mx} -jar {ChromHMM_dir/'ChromHMM.jar'}"
    
    uv_params = ""
    if coordir:
        p_coordir = Path(coordir)
        uv_params += f"-u {p_coordir.absolute()} "
        (p_coordir / f"{genome}").mkdir(exist_ok=True)
        if defaultcoor:
            shutil.copytree(ChromHMM_dir/f"COORDS/{genome}", p_coordir/f"{genome}", dirs_exist_ok=True)
        for file in p_coordir.glob("*"):
            if file.is_file(): shutil.copy(file, p_coordir/f"{genome}")
    if anchordir:
        p_anchordir = Path(anchordir)
        uv_params += f"-v {p_anchordir.absolute()}"
        (p_anchordir / f"{genome}").mkdir(exist_ok=True)
        if defaultcoor:
            shutil.copytree(ChromHMM_dir/f"ANCHORFILES/{genome}", p_anchordir/f"{genome}", dirs_exist_ok=True)
        for file in p_anchordir.glob("*"):
            if file.is_file(): shutil.copy(file, p_anchordir/f"{genome}")
    
    stacked_param = "-stacked" if stacked else ""
    chrom_len = ChromHMM_dir / f"CHROMSIZES/{genome}.txt"
    strand_param = "-nostrand" if nostrand else ""
    name_param = f"-i {name}"

    ChromHMM_bin = f"BinarizeBed {stacked_param} -peaks -b {binsize} {chrom_len} {directory} {markfile} {binarydir}"
    ChromHMM_lm = f"LearnModel {uv_params} {strand_param} {name_param} -p {processor} \
                    -b {binsize} {binarydir} {outdir} {numstates} {genome}"
    if not no_binary:
        print(ChromHMM_bin)
        info = subprocess.Popen([*ChromHMM_jar.split(), *ChromHMM_bin.split()])
        info.wait()
    print(ChromHMM_lm)
    info = subprocess.Popen([*ChromHMM_jar.split(), *ChromHMM_lm.split()])
    info.wait()

    rscript = Path(__file__).parent / "R" / "ChromHMM_hm.R"
    for of in Path(outdir).glob("*_overlap.txt"):
        info = subprocess.Popen([ul.Rscript_bin(), str(rscript), of])
    info.wait()


def extract_signal(
    output: FilePath,        
    coor_dic: Dict,
    bigwig: FilePath,
    stype : str,
    nbins: int
):
    """extract bigwig signal from coordinate

    Args:
        output (FilePath): output    
        coor_dic (Dict): AS event splice sites coordinate dict
        bigwig (FilePath): bigwig files
        stype (str): extract bigwig signal type
        nbins (int): bin number

    Returns:
        str: _description_
    """
    import pyBigWig

    def _get_signal(row, bw, stype, nbins):
        values = bw.stats(row[0], row[1], row[2], type=stype, nBins=nbins, exact=True)
        if row[5] == "-":
            values = reversed(values)
        return pd.Series(values)
    
    sdf_ls = []
    for ss, df in coor_dic.items():  
        bw = pyBigWig.open(bigwig)
        _func = partial(_get_signal, bw=bw, stype=stype, nbins=nbins)
        sdf = df.apply(_func, axis=1)
        sdf.columns = [f"{ss}_{i}" for i in sdf.columns]
        sdf_ls.append(sdf)
    dfm = pd.concat(sdf_ls, axis=1)
    dfm.to_csv(output)
