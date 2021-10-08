# import feature_len as fl
# from . import utils  as ul

def get_anchor_coor(event_id, index, sideindex, width, right_flank, left_flank):
    from . import feature_len as fl

    eid = fl.EventID(event_id)
    if index is None and sideindex is None:
        s, e = eid.alter_element_coor
        anchor =  s + round((e - s) / 2)
    else:
        try:
            if index:
                anchor = eid.coordinates[index-1]
            else:
                s, e = [eid.coordinates[i-1] for i in sideindex]
                anchor =  s + round((e - s) / 2)
        except IndexError:
            print("IndexError")
            exit()
    if width and width > 0:
        start, end = anchor - round(width / 2), anchor + round(width / 2)
        return eid.Chr, start, end
    elif all([right_flank, left_flank]):
        start, end = anchor - left_flank, anchor + right_flank
        return eid.Chr, start, end
    else:
        return eid.Chr, anchor, eid.strand


def gen_anchor_bed(dpsi_file, out, index, sideindex, width, right_flank, left_flank):
    import pandas as pd
    from functools import partial

    wget_anchor_coor = partial(get_anchor_coor, 
            index=index, sideindex=sideindex, width=width, 
            right_flank=right_flank, left_flank=left_flank)
    dpsi_df = pd.read_csv(dpsi_file, sep="\t", index_col=0)
    dpsi_df["event_id"] = dpsi_df.index
    coors = dpsi_df["event_id"].apply(wget_anchor_coor)
    coor_df = pd.DataFrame(coors.tolist())
    coor_df.drop_duplicates(inplace=True)
    coor_df.to_csv(out, index=False, header=False, sep="\t")
                                                                                                                                                                                                      

def install(path):
    import os
    import zipfile
    from urllib import request

    jar = path / "ChromHMM/ChromHMM.jar" 
    if jar.exists():
        print("ChromHMM has installed!")
        exit()

    ChromHMM_url = "http://compbio.mit.edu/ChromHMM/ChromHMM.zip"
    software = path / "ChromHMM.zip"
    with request.urlopen(ChromHMM_url) as f:
        data = f.read()
        with open(software, "wb") as h:
            h.write(data)
    with zipfile.ZipFile(software, 'r') as zip_ref:
        zip_ref.extractall(path)
    os.remove(software)
    
