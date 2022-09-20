# -*- coding: utf-8 -*-
"""
astk.utils.meta_template
~~~~~~~~~~~~~~~~~
This module provides a metadata template table and json output 
inferring from replicates and file paths.
"""

import json
from itertools import chain, repeat
from pathlib import Path

from pandas import DataFrame, concat


class Template:
    def __init__(self, condition):
        self.condition = condition if condition else ["ctrl", "case"]

    def infer_file(self, path, rep):
        rep_num = len(rep)
        path_num = len(path)

        if rep_num == 1:
            if isinstance(fold := path_num // int(rep[0]), int):
                group_num = fold
            else:
                raise ValueError("path file number or replicate number is wrong!")
        else:
            if path_num == sum(rep):
                group_num = rep_num
            else:
                raise ValueError("path file number or replicate number is wrong!")
        return group_num

    def df_generate(self, group_name, rep, path, **kwargs):
        group_num = len(group_name)

        if len(rep) == 1:
            replicates = [int(rep[0])] * group_num
        else:
            replicates = rep
        if len(path) == int(rep[0]):
            path_ls = list(chain(*repeat(path,  group_num)))
        else:
            path_ls = path
        rep_ls = list(chain(*((range(1, i+1)) for i in replicates)))
        group_ls = list(chain(*[repeat(i, r) for r,i in zip(replicates, group_name)]))

        is_fn = kwargs.get("is_fn", None)
        split_param = kwargs.get("split", None)
        
        if is_fn:
            names = [Path(i).stem for i in path_ls] 
        else:
            names = [Path(i).parent.name for i in path_ls]
        if split_param:
            sub_name = lambda idxs, nls: "_".join([nls[int(i)-1] for i in idxs])
            ssym, sidx = split_param[0], split_param[1:]
            names = [sub_name(sidx, n.split(ssym)) for n in names]            

        sdf = DataFrame({
            "group": group_ls, 
            "name": names,
            "path": path_ls,
            "replicate": rep_ls
        })
        return sdf

    def complete_df(self, group, path1, path2, repN1, repN2, **kwargs):
        if len(path1) != 0 and len(path2) != 0:
            group1_num = self.infer_file(path1, repN1)
            group2_num = self.infer_file(path2, repN2)
            group_num = max(group1_num, group2_num)
            if min(group1_num, group2_num) > 1 and group1_num != group2_num:
                raise ValueError("path1 dismatched path2!")
            if len(group) == group_num:
                group_name = group
            elif len(group) == 0:
                group_name = list(range(1, group_num+1))
            else:
                raise ValueError("dismatched path files!")

            df1 = self.df_generate(group_name, repN1, path1, **kwargs)
            df2 = self.df_generate(group_name, repN2, path2, **kwargs)
            df1.insert(1, "condition", self.condition[0])
            df2.insert(1, "condition", self.condition[1])
            df = concat([df1, df2]).sort_values(by=["group"])
  
        else:
            if path1:
                cdname = self.condition[0]
                path = path1
                repN = repN1
            else:
                cdname = self.condition[0]
                path = path2
                repN = repN2            
            group_num = self.infer_file(path, repN)
            if len(group) == group_num:
                group_name = group
            elif len(group) == 0:
                group_name = list(range(1, group_num+1))
            else:
                raise ValueError("group names dismatched path files!")                
            df = self.df_generate(group_name, repN, path, **kwargs)
            
            df.insert(1, "condition", cdname) 
        
        self.df = df
    
    def to_json(self, out):
        df = self.df
        meta_dic = {} 
        for gn, gdf in df.groupby("group"):
            meta_dic[gn] = {cd:{"samples": []} for cd in df["condition"].unique()}
            for _, row in gdf.iterrows():
                name = row["name"]
                rep = row["replicate"]
                path = row["path"]
                sp_dic = {"name": name, "replicate": rep, "path": path}
                meta_dic[gn][row["condition"]]["samples"].append(sp_dic)
        
        outjson = Path(out).with_suffix(".json")

        f = outjson.open(mode="w")
        json.dump(meta_dic, f, indent=4)
        f.close()    

    def to_csv(self, out):
        out = Path(out).with_suffix(".csv")
        self.df.to_csv(out, index=False)
