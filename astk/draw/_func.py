import subprocess
from pathlib import Path

from astk.constant import *
import astk.utils.func as ul



def gseplot(termid, output, rdata, fmt, width, height, resolution):
    
    rscript = BASE_DIR / "R" / "gseaPlot.R"
    param_dic = {
        "termid": termid,
        "fmt": fmt, 
        "RData": rdata,
        "output": output,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def upset(files, output, xlabel, dg, fmt, width, height, resolution):

    rscript = BASE_DIR / "R" / "upset.R"

    param_dic = {
        "file": files,
        "dg": dg,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}"),
        "name": xlabel if xlabel else [str(i) for i  in range(len(files))]
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def volcano(file, output, fmt, width, height, resolution):
    rscript = BASE_DIR / "R" / "volcano.R"
    param_dic = {
        "file": file,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def pca(files, output, fmt, width, height, resolution, groupname):

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    rscript = BASE_DIR / "R" / "pca.R"
    param_dic = {
        "file": files,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}"),
        "groupname": groupname
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def heatmap(files, output, cluster, fmt, width, height, resolution):

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    rscript = BASE_DIR / "R" / "heatmap.R"

    param_dic = {
        "file": files,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "clusterinfo": cluster if cluster else "0",
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])



def barplot(output, files, xlabel, dg, fmt, width, height, resolution):
    rscript = BASE_DIR / "R" / "barplot.R"
    param_dic = {
        "file": files,
        "name": xlabel,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output,
        "dg": dg
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])
