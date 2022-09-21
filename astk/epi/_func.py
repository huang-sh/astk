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
from astk.types import FilePath


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
        if c is None:
            csam = None
        else:
            csam = pysam.AlignmentFile(c["path"])

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



# def epi_sc(output, metadata, anchor, name, width, binsize):

#     names = name if len(name)==len(anchor) else list(range(1, len(anchor)+1))

#     anchor_dic = dict(zip(names, anchor))
#     epi_signal(output, anchor_dic, metadata, width, binsize)


def epihm(output, files, fmt, width, height, resolution):

    rscript = Path(__file__).parent / "R" / "signalHeatmap.R"
    param_dic = {
        "file": files,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def sigcmp(output, file, fmt, width, height, resolution):

    rscript = Path(__file__).parent / "R" / "signalCompare.R"
    param_dic = {
        "file": file,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])



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


def epi_sc(event_file, event_label, bam_file, bam_label, width, bin_size,
            output, normalmethod, paired_end, title, bam_merge, as_type):

    rscript = BASE_DIR / "R" / "epiFeature.R"
    param_dic = {
        "binsize": bin_size,
        "width": width, 
        "bam": bam_file,
        "bamlabel": bam_label,
        "region": event_file, 
        "regionlabel": event_label, 
        "title": title,
        "out": output,
        "markmerge": bam_merge,
        "normalmethod": normalmethod,
        "pairedend": paired_end,
        "ASType": as_type
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def epi_profile(file, output, title, ylim, fmt, width, height, resolution):

    rscript = BASE_DIR / "R" / "epiProfile.R"
    param_dic = {
        "title": title,
        "file": file,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output,
        "ylim": ylim
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    # print(param_ls)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def signal_heatmap(
    output: FilePath,
    event_file: FilePath,
    bw_files: Sequence[FilePath],
    bin_size: int,
    ups_width: int,
    dws_width: int,
    threads: int,
    plot_type: str,
    color_map: str,
    plot_format: str,
    fig_height: float,
    fig_width: float,
    regionsLabel: Sequence[str],
    samplesLabel: Sequence[str]
) -> None:
    from deeptools.computeMatrix import main as cpm
    from deeptools.plotHeatmap import main as plothm
    from pandas import concat

    regionFiles = []

    df_dic = ul.get_evnet_ss_bed(event_file, 150, 150)
    for df in df_dic.values():
        f = NamedTemporaryFile(delete=False)
        f.close()
        df.to_csv(f.name, index=False, header=False, sep="\t")
        regionFiles.append(f.name)
        
    out = Path(output)
    mtx_out = str(out.with_suffix(".mat.gz"))
    fig_out = str(out.with_suffix(f".{plot_format}"))
    mtx_args = [
        "reference-point",
        "-R", *regionFiles,
        "-S", *bw_files,
        "--referencePoint", "center",
        "--beforeRegionStartLength", str(ups_width), 
        "--afterRegionStartLength", str(dws_width), 
        "--binSize", str(bin_size), 
        "--numberOfProcessors", str(threads),
        "--outFileName", mtx_out
    ]
    hm_args = [
        "--matrixFile", mtx_out,
        "--outFileName", fig_out,
        "--plotType", plot_type, 
        "--colorMap", color_map,
        "--heatmapHeight", str(fig_height),
        "--heatmapWidth", str(fig_width),
        "--perGroup",
        "--plotFileFormat", plot_format,
        "--regionsLabel", *regionsLabel,
        #"--samplesLabel", *samplesLabel,
        "--refPointLabel", "SS"
    ]
    cpm(mtx_args)
    plothm(hm_args)
    # [os.unlink(i) for i in regionFiles]
    print(regionFiles)


def signal_metaplot(
    output: FilePath,
    mat_file: FilePath,
    name: str,
    groupname: Sequence[str],
    fig_height: int,
    fig_width: int,
    fig_format: str,
    resolution: int
) -> None:
    from pandas import concat, read_csv
    import gzip

    header_dic = {}
    files = []
    for idx, (fn, file) in enumerate(zip(groupname, mat_file)):

        tempf = NamedTemporaryFile(delete=False)
        tempf.close()

        with gzip.open(file, "rb") as f:
            line = f.readline()
            if line[:3] != b'@{"':
                raise ValueError(f"{file} is wrong!")                
            dic = json.loads(line[1:])
            header_dic[idx] = dic
        
        bds = header_dic[idx]["group_boundaries"]
        sub_df_ls = []
        df = read_csv(file, sep="\t", header=None, comment="@")
        for ai in range(len(bds)-1):    
            sub_df = df.iloc[bds[ai]:bds[ai+1], 6:]
            sub_df.index = range(bds[ai+1]-bds[ai])
            print(sub_df.shape)
            sub_df.columns = [f"{name}_{fn}_a{ai+1}_bin{i+1}" for i in range(sub_df.shape[1])]
            sub_df_ls.append(sub_df)
        score_df = concat(sub_df_ls, axis=1)
        score_df.to_csv(tempf.name, index=False)
        files.append(tempf.name)
        
    rscript = BASE_DIR / "R" / "epiProfile.R"
    param_dic = {
        "title": name,
        "file": files,
        "width": fig_width, 
        "height": fig_height, 
        "resolution": resolution,
        "fmt": fig_format,
        "output": output
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])
