import re
import json

from astk.constant import BASE_DIR

elm_out_field = {
    "AS_event": [],
    "elm_id": [],
    "transcript": [],
    "elm": [],
    "match_aa": [],
    
}

def search_elm(file, output, genome):
    import pandas as pd
    from Bio.Seq import Seq

    se_dpsi = pd.read_csv(file, sep="\t", index_col=0)
    elm_cls_file = BASE_DIR / "data/motif/ELM/elm_classes.tsv"
    elm_cls_df = pd.read_csv(elm_cls_file, sep="\t", comment="#")

    exon_seq_js = BASE_DIR / f"data/motif/ELM/{genome}_alter_exon_seq.json"
    with open(exon_seq_js, "r") as f:
        ae_seq_dic = json.load(f)

    for _, row in elm_cls_df.iterrows():
        elm_regex = row["Regex"]
        elm_id = row["ELMIdentifier"]
        p = re.compile(elm_regex)
        for ei in se_dpsi.index:
            ei_dic = ae_seq_dic.get(ei, None)
            if ei_dic is None:
                continue
            strand = ei_dic["strand"]
            for tx, seq in ei_dic["sequence"].items():
                exon_seq = seq
                if strand == "+":
                    aa = Seq(exon_seq).translate()
                elif strand == "-":
                    aa = Seq(exon_seq).reverse_complement().translate()
                if match_elm := p.search(str(aa)):
                    m_str = match_elm.group()
                    elm_out_field["AS_event"].append(ei)
                    elm_out_field["elm"].append(elm_regex)
                    elm_out_field["match_aa"].append(m_str)
                    elm_out_field["elm_id"].append(elm_id)
                    elm_out_field["transcript"].append(tx)

    pd.DataFrame(elm_out_field).to_csv(output)