"""
astk.cli.utils
~~~~~~~~~~~~~~~~~
This module provide some utility function
"""

from .config import *
import astk.utils._cli_func as ul


@cli_fun.command(help="generate metadata template file")
@click.option('-c1', '--ctrl', "ctrl_file", cls=MultiOption, type=click.Path(exists=True),
                required=True, help="file path for condtion 1(control)")
@click.option('-c2', '--case', "case_file", cls=MultiOption, type=click.Path(exists=True),
                required=True, help="file path for condtion 2(case)")
@click.option('-gn', '--groupName', cls=MultiOption, type=str,
                default=(), help="group name")
@click.option('-repN', '--replicate', cls=MultiOption, type=int,
                help="replicate, number")
@click.option('-o', '--output', required=True, type=click.Path(),
                help='metadata output path')
@click.option('-repN1', '--replicate1', cls=MultiOption, type=int,
                help="replicate1 number(control)")
@click.option('-repN2', '--replicate2', cls=MultiOption, type=int, 
                help="replicate2, number(case)")
@click.option('--condition', cls=MultiOption, type=str, help="condition name, default='ctrl case'")
@click.option('-fn', '--filename', 'is_fn', is_flag=True, help="file name")
@click.option('--split', cls=MultiOption, type=str, help="name split symbol and index")
@click.option('-app', '--app', type=click.Choice(['SUPPA2', 'rMATS']), default="SUPPA2",
                help="application for output; default='SUPPA2'")
def meta(*args, **kwargs):
    if not (pdir:= Path(kwargs["output"]).parent).exists():
        raise UsageError(f"Invalid value for '-o' / '--output': Path '{pdir}' does not exist.")
    replicate = kwargs.pop("replicate")
    replicate1, replicate2 = kwargs.pop("replicate1"), kwargs.pop("replicate2")
    if replicate:
        kwargs["ctrl_rep"] = [int(i) for i in replicate]
        kwargs["case_rep"] = [int(i) for i in replicate]
    elif all([replicate1, replicate2]):
        kwargs["ctrl_rep"] = [int(i) for i in replicate1]
        kwargs["case_rep"] = [int(i) for i in replicate2]
    else:
        raise UsageError("you should set one option: -repN or two options: -repN1, -repN2.")

    if (gn := len(kwargs["groupname"])) > 0:
        case_fn, ctrl_fn = len(kwargs["case_file"]), len(kwargs["ctrl_file"])
        if case_fn > ctrl_fn:
            cfn = case_fn
            crep = kwargs["case_rep"]
        else:
            cfn = ctrl_fn
            crep = kwargs["ctrl_rep"]
        if (len(crep) == 1) and (cfn // crep[0] != gn):
                raise UsageError("Value number for '-gn' / '--groupName' is wrong.")
        elif (len(crep) > 1) and len(crep) != gn:
                raise UsageError("Value number for '-gn' / '--groupName' is wrong.")

    ul.meta(*args, **kwargs)


@cli_fun.command(help="install R packages")
@click.option('-r', '--requirement', is_flag=True, default=False,
                help="install astk requirement R packages")
@click.option('-OrgDb', '--OrgDb', "OrgDb", cls=MultiOption, type=str, 
                default=(), help="install Genome wide annotation package")
@click.option('-cran', '--cran', cls=MultiOption, type=str, 
                default=(), help="install CRAN package")
@click.option('-bioc', '--bioconductor', cls=MultiOption, type=str, 
                default=(),help="install Bioconductor package")
@click.option('-j', '--java',  is_flag=True, help="install java software")
@click.option('-m', '--mirror',  is_flag=True, default=False,
                help="use tsinghua mirrors.")
@click.option('--conda',  is_flag=True, default=False,
                help="install packages via conda")
def install(*args, **kwargs):
    ul.install(*args, **kwargs)


@cli_fun.command(help="get motif from meme file")
@click.argument('motifId', nargs=-1, required=True)   
@click.option('-mm', '--meme', type=click.Path(exists=True), help="meme motif file")
@click.option('-db', '--database', type=click.Choice(['ATtRACT', 'CISBP-RNA']),
                help="RBP motif database")
@click.option('-org', '--organism', help="RBP organism")
@click.option('-o', '--output', required=True, help="output path")
def getmeme(*args, **kwargs):

    ul.getmeme(*args, **kwargs)

@cli_fun.command(name="list", help="list OrgDb")
@click.option('-orgdb', '--OrgDb', "OrgDb", is_flag=True, help="list OrgDb")
@click.option('-rbpsp', '--RBPSp', "RBPSp", is_flag=True, help="RNA binding protein  ")
def list_(*args, **kwargs):
    ul.list_(*args, **kwargs)


# @cli_fun.command(help = "generate ChromHMM anchor file")
# @click.argument('file', type=click.Path(exists=True), required=True)
# @click.option('-o', '--output', required=True, help="file output path")
# @click.option('-idx', '--index', type=int, help="element index")
# @click.option('-si', '--sideIndex', type=(int, int), help="the center of two side index")
# @click.option('-u', '--upstreamOffset', "offset5", type=int, default=0, help="upstream offset")
# @click.option('-d', '--downstreamOffset', "offset3", type=int, default=0, help="downstream offset")
# @click.option('-ss', '--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
# def anchor(*args, **kwargs):
#     ul.anchor(*args, **kwargs)
 

@cli_fun.command(help="retrieve coordinates from splice site region")
@click.option('-i', '--input', 'event_file', type=click.Path(exists=True),
                required=True,  help='AS event file')
@click.option('-o', '--output', required=True, help="output path")
@click.option('-si', '--siteIndex', "ss_idx", type=int, default=0,
                help="splice site index. if not set, it will use all splice sites. 1-based")
@click.option('-uw', '--upstreamWidth', "ups_width", type=int, default=150,
                help="flank width of splice site upstream")
@click.option('-dw', '--downstreamWidth', "dws_width", type=int, default=150,
                help="flank width of splice site downstream")
@click.option('-ss', '--spliceSite', "sss", is_flag=True, 
                help="get splice site region window sizes")
@click.option('--interval', type=(int, int), default=(None, None),
                help="interval the between two splice sites, 1-based")
@click.option('-fi', 'fasta', type=click.Path(exists=True),
                help="Input FASTA file. if set, the fasta sequence will be extracted")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file")
def getcoor(*args, **kwargs):
    import astk.utils as ul

    output = kwargs.pop("output")
    interval_idx = kwargs.pop("interval")
    fasta = kwargs.pop("fasta")
    if kwargs["site_idx"]:
        sites = [kwargs["site_idx"] - 1]
    elif all(interval_idx):
        sites = [i-1 for i in sorted(interval_idx)]

    kwargs["ss_idx"] = sites
    coord_dic = ul.get_ss_bed(*args, **kwargs)
    if len(coord_dic) == 1:
        df = list(coord_dic.values())[0]
        df.to_csv(output, sep="\t", index=False, header=False)
        print(list(coord_dic.values())[0].head())
    elif len(coord_dic) == 2:
        df = ul.get_ss_range(kwargs["event_file"], *sites, kwargs["app"])
        df.to_csv(output, sep="\t", index=False, header=False)


@cli_fun.command(help="Make the TxDb object")
@click.argument('gtf', type=click.Path(exists=True), required=True)
@click.option('-org', '--organism', required=True, help="organism")
@click.option('-o', '--output', help="file output path")
def mktxdb(*args, **kwargs):
    ul.mkTxDb(*args, **kwargs)


@cli_fun.command(help="get gene ID from AS event file")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), help="file output path")
@click.option('-u', '--unique', is_flag=True, help="only save unique gene ID")
def getgene(*args, **kwargs):
    ul.getgene(*args, **kwargs)


@cli_fun.command(name="merge", help="merge file")
@click.option('-i','--input', 'files' ,cls=MultiOption, type=click.Path(exists=True),
                required=True , help="input files")
@click.option('-o', '--output', type=click.Path(), help="output path")
@click.option('-axis', '--axis', type=click.IntRange(min=0, max=1), default=0,
                help="merge direction, 0 for row merge and 1 for column merge")
@click.option('-rmdup', '--rmdup', type=click.Choice(["index", 'all', 'content']), 
                help="remove duplicate rows")
@click.option('-rmna', '--rmna', is_flag=True, help="remove NA data")
@click.option('-gl', '--groupLabel', cls=MultiOption, help="add label prefix for different files")
def sub_merge_files(*args, **kwargs):
    from astk.utils.func import merge_files

    if grouplabel := kwargs["grouplabel"]:
        if len(grouplabel) != len(kwargs["files"]):
            raise UsageError("-gl/--groupLabel and -i/--input values number muset be same!")
    merge_files(*args, **kwargs)


@cli_fun.command(name="regionTest", help="region Statistical tests")
@click.option('-e1', "--event1", 'file1', type=click.Path(exists=True), required=True,
                help="event1 file")
@click.option('-e2', "--event2", 'file2', type=click.Path(exists=True), required=True,
                help="event2 file")
@click.option('-bed', "--bed", type=click.Path(exists=True), required=True,
                help="bed file")
@click.option('-uw', '--upstreamWidth', "ups_width", type=int, default=150,
                help="flank width of splice site upstream")
@click.option('-dw', '--downstreamWidth', "dws_width", type=int, default=150,
                help="flank width of splice site downstream")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file")
@click.option('-o', '--output', type=click.Path(), help="output path")              
def sc_region_test(*args, **kwargs):
    import pyranges as pr
    from pandas import DataFrame
    from astk.utils._getfasta import get_coor_fa
    from astk.utils import get_ss_bed, read_fasta, detect_file_info
    from statsmodels.stats.proportion import proportions_chisquare
    import numpy as np
    import pandas as pd

    file1 = kwargs.pop("file1")
    file2 = kwargs.pop("file2")

    app = kwargs.pop("app")

    if app == "auto":
        app = detect_file_info(file1)["app"]

    # outdir = Path(outdir)
    # Path(outdir).mkdir(exist_ok=True)
    coord_dic1 = get_ss_bed(file1, app=app, **kwargs)
    coord_dic2 = get_ss_bed(file2, app=app, **kwargs)

    bedgr = pr.read_bed(kwargs.pop("bed"))

    event_ls1 = set()
    event_ls2 = set()
    df_score = DataFrame()
    for (ssn1, df1),(ssn2, df2)  in zip(coord_dic1.items(), coord_dic2.items()):
        # ss_dir = outdir / ssn
        # ss_dir.mkdir(exist_ok=True)
        # df.to_csv(ss_dir / f"{ssn}.bed", index=False, header=False, sep="\t")
        df1.columns = ["Chromosome", "Start", "End", "Name",  "Score", "Strand"]
        ss_gr1 = pr.PyRanges(df1)
        df2.columns = ["Chromosome", "Start", "End", "Name",  "Score", "Strand"]
        ss_gr2 = pr.PyRanges(df2)

        a1 = ss_gr1.intersect(bedgr)
        a2 = ss_gr2.intersect(bedgr)

        success_cnts = np.array([a1.df.shape[0], a2.df.shape[0]])
        total_cnts = np.array([df1.shape[0], df2.shape[0]])
        chi2, p_val, cont_table = proportions_chisquare(count=success_cnts, nobs=total_cnts)

        event_ls1.update(a1.df["Name"])
        event_ls2.update(a2.df["Name"])
        print(a1.df.shape[0], df1.shape[0], a2.df.shape[0], df2.shape[0], p_val)

    dpsi_df1 = pd.read_csv(file1, sep="\t", index_col=0)
    dpsi_df2 = pd.read_csv(file2, sep="\t", index_col=0)

    dpsi_df1.loc[event_ls1, ].to_csv(Path(file1).with_suffix(".peak.dpsi"), sep="\t")
    dpsi_df2.loc[event_ls2, ].to_csv(Path(file2).with_suffix(".peak.dpsi"), sep="\t")

