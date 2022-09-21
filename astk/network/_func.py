import sys
import subprocess
from pathlib import Path

import astk.utils.func  as ul
from astk.constant import BASE_DIR, RBP_sp_dic


def coSpliceNet(file, output):
    from collections import Counter
    from pyecharts import options as opts
    from pyecharts.charts import Graph
    import pandas as pd

    link_df = pd.read_csv(file)
    link_sort = link_df.sort_values(by=["cor", "pval"], key=abs, ascending=False)
    link_fil = link_sort.drop_duplicates(subset=["rbp", "event"])

    rbp_symbol = []
    event_id = []
    links = []

    for row in link_fil.itertuples():
        rbp = row[-1].split(":")[0]
        event = row[-2]
        rbp_symbol.append(rbp)
        event_id.append(event)
        links.append({"source": rbp, "target": event})

    rbp_nodes = []
    for k, v in Counter(rbp_symbol).items():
        node = {"name": k, "symbolSize": v, "category": k,
                'label': {'normal': {'show': 'True'}}}
        rbp_nodes.append(node)

    category = [{"name": k} for k in set(rbp_symbol)]
    event_nodes = [{"name": k, "symbolSize": max([v,5])} for k, v in Counter(event_id).items()]
    nodes = rbp_nodes + event_nodes

    c = (
        Graph().add(
            "",
            nodes,
            links,
            category,
            repulsion=50,
            linestyle_opts=opts.LineStyleOpts(curve=0.2),
            label_opts=opts.LabelOpts(is_show=False),
            edge_symbol=[None, 'arrow'],
            edge_symbol_size=4
        )
        .set_global_opts(
            legend_opts=opts.LegendOpts(is_show=False),
            title_opts=opts.TitleOpts(title="Co-splicing network"),
        )
        .render(output)
    )


def co_splice_net(output, fasta, psimeta, tqmeta, organism, database, txdb):
    sp = RBP_sp_dic.get(organism, None)
    if sp is None:
        print("--organism is not valid")
        sys.exit(1)

    rscript = BASE_DIR / "R" / "CoSpliceNet.R"
    csv_out = Path(output).with_suffix(".csv")
    html_out = Path(output).with_suffix(".html")
    param_dic = {
        "output": csv_out,
        "fasta": fasta, 
        "psiMeta": psimeta, 
        "tqMeta": tqmeta,
        "organism": sp,
        "database": database,
        "txdb": txdb
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])
    coSpliceNet(csv_out, html_out)
