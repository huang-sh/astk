from astk.event import SuppaEventID
from astk.constant import SS_SCORE_LEN, rMATS_POS_COLS
from astk.types import FilePath
from astk.utils import detect_file_info


class EventCoord:

    def __init__(self, file) -> None:
        self.set_metadata(file)

    def get_ss_flank_bed(
        self,
        idx: int,
        ups_w: int = 150,
        dws_w: int = 150,
        sss: bool = False,
        based0 = True,
        **kwargs
    ):
        """get a splice site flank coordination bed
        """            
        df_ss = self.df_ss.copy()
        slen = SS_SCORE_LEN[self.app][self.etype][idx]
        ps_idx = self.ps_idx 
        ns_idx = self.ns_idx
        exon_w = kwargs.get("exon_width", 0) 
        intron_w = kwargs.get("intron_width", 0)
        ps_idx_pos = df_ss.loc[ps_idx, idx]
        ns_idx_pos = df_ss.loc[ns_idx, idx]
        ups_w, dws_w = self._width_adaptor(slen, ups_w, dws_w, exon_w, intron_w, sss)
        # shift coordiante for keeping consistent
        if abs(slen) == 23:            
            ps_idx_pos -= 1
        elif abs(slen) == 9:
            ns_idx_pos -= 1
        if kwargs.get("split", False):
            df_up_bed = self.df_temp.copy()
            df_dw_bed = self.df_temp.copy()
            ex_ups_w, ex_dws_w = self._exclude_region(slen, **kwargs)
            df_up_bed.loc[ps_idx, "start"] = ps_idx_pos - ups_w
            df_up_bed.loc[ps_idx, "end"] = ps_idx_pos - ex_ups_w
            df_up_bed.loc[ns_idx, "start"] = ns_idx_pos  + ex_ups_w
            df_up_bed.loc[ns_idx, "end"] = ns_idx_pos + ups_w
            df_dw_bed.loc[ps_idx, "start"] = ps_idx_pos + ex_dws_w
            df_dw_bed.loc[ps_idx, "end"] = ps_idx_pos + dws_w
            df_dw_bed.loc[ns_idx, "start"] = ns_idx_pos - dws_w 
            df_dw_bed.loc[ns_idx, "end"] = ns_idx_pos  - ex_dws_w
            if not based0:
                df_up_bed.loc[:, "start"] += 1           
                df_dw_bed.loc[:, "start"] += 1
            return df_up_bed, df_dw_bed
        else:
            df_bed = self.df_temp.copy()
            df_bed.loc[ps_idx, "start"] = ps_idx_pos - ups_w 
            df_bed.loc[ps_idx, "end"] = ps_idx_pos + dws_w
            df_bed.loc[ns_idx, "start"] = ns_idx_pos - dws_w 
            df_bed.loc[ns_idx, "end"] = ns_idx_pos + ups_w
            if not based0:
                df_bed.loc[:, "start"] += 1
            return df_bed

    def get_all_flank_bed(
        self,
        ups_w: int = 150,
        dws_w: int = 150,
        sss: bool = False,
        based0 = True,
        **kwargs
    ):  
        """get all splice site flank coordination bed and return a dict
        """
        dic = {}
        slens = SS_SCORE_LEN[self.app][self.etype]
        site_dic = {23: "3SS", 9: "5SS", -23: "pse_3SS", -9: "pse_5SS"}
        for idx, slen in enumerate(slens):
            sn = site_dic[slen]
            key = f"A{idx+1}_{sn}"
            dic[key] = self.get_ss_flank_bed(idx, ups_w, dws_w, sss, based0, **kwargs)
        return dic

    def _5ss_width(self, ups_w, dws_w, exon_w, intron_w, sss):
        if exon_w and intron_w:
            ups_w, dws_w = exon_w, intron_w
        elif sss:
            ups_w, dws_w = 3, 6
        return ups_w, dws_w

    def _3ss_width(self, ups_w, dws_w, exon_w, intron_w, sss):
        if exon_w and intron_w:
            ups_w, dws_w = intron_w, exon_w
        elif sss:
            ups_w, dws_w = 20, 3
        return ups_w, dws_w

    def _width_adaptor(self, site, ups_w, dws_w, exon_w, intron_w, sss):
        if abs(site) == 23:
            ups_w, dws_w = self._3ss_width(ups_w, dws_w, exon_w, intron_w, sss)
        elif abs(site) == 9:
            ups_w, dws_w = self._5ss_width(ups_w, dws_w, exon_w, intron_w, sss)
        return ups_w, dws_w
    
    def _exclude_region(self, site, **kwargs):
        if kwargs.get("excludeSS", False):
            if site == 23:
                ex_ups_w, ex_dws_w = 20, 3
            elif site == 9:
                ex_ups_w, ex_dws_w = 3, 6
            else:
                ex_ups_w, ex_dws_w = 0, 0
        else:
            ex_ups_w, ex_dws_w = 0, 0
        return ex_ups_w, ex_dws_w
    
    def _ns_adaptor(self):
        pass

    def set_metadata(self, file):
        pass


class SuppaEventCoord(EventCoord):

    def __init__(self, file) -> None:
        super().__init__(file)

    def set_metadata(self, file):
        from pandas import DataFrame, read_csv

        fileinfo = detect_file_info(file)
        self.etype = fileinfo["etype"]
        self.app = fileinfo["app"]

        dpsi_df = read_csv(file, sep="\t", index_col=0)
        dpsi_df.drop_duplicates(inplace=True)
        dpsi_df.dropna(inplace=True)
        dpsi_df["event_id"] = dpsi_df.index

        def _get_coor(event_id):
            ei = SuppaEventID(event_id)
            if ei.strand == "-":
                coords = reversed(ei.coordinates)
            else:
                coords = ei.coordinates
            row = ei.Chr, ei.event_id, ei.strand, *coords
            return row
        event_df = DataFrame(dpsi_df["event_id"].apply(_get_coor).tolist(),
                             index=dpsi_df.index)
        self.df_ss = event_df.iloc[:, 3:]
        self.df_ss.columns = range(self.df_ss.shape[1])
        self.ps_idx = event_df.loc[event_df.iloc[:, 2] == "+", ].index
        self.ns_idx = event_df.loc[event_df.iloc[:, 2] == "-", ].index

        cols = ["seqname", "start", "end", "name", "score", "strand"]
        df_temp = DataFrame(columns=cols, index=dpsi_df.index)
        df_temp["seqname"] = event_df[0]
        df_temp["strand"] = event_df[2]
        df_temp["score"] = 0
        df_temp["name"] = event_df[1]
        self.df_temp = df_temp

        
class rMATSEventCoord(EventCoord):

    def __init__(self, file) -> None:
        super().__init__(file)

    def set_metadata(self, file):
        from pandas import DataFrame, read_csv, concat

        dpsi_df = read_csv(file, sep="\t", index_col=0)
        dpsi_df.drop_duplicates(inplace=True)
        dpsi_df.dropna(inplace=True)
        fileinfo = detect_file_info(file)
        self.etype = fileinfo["etype"]
        self.app = fileinfo["app"]

        self.ps_idx = dpsi_df.loc[dpsi_df["strand"] == "+", ].index
        self.ns_idx = dpsi_df.loc[dpsi_df["strand"] == "-", ].index

        cols = ["seqname", "start", "end", "name", "score", "strand"]
        df_temp = DataFrame(columns=cols, index=dpsi_df.index)
        df_temp["seqname"] = dpsi_df["chr"]
        df_temp["strand"] = dpsi_df["strand"]
        df_temp["score"] = 0
        df_temp["name"] = dpsi_df["geneSymbol"]
        self.df_temp = df_temp

        # tansform 1-based index for unified management
        if self.etype == "SE":
            start_col = ["upstreamES", "exonStart_0base", "downstreamES"]
        elif self.etype in ["A5SS", "A3SS"]:
            start_col = ["longExonStart_0base", "flankingES", "shortES"]
        elif self.etype == "MXE":
            start_col = ["upstreamES", "1stExonStart_0base", "2ndExonStart_0base", "downstreamES"]
        elif self.etype == "RI":
            start_col = ["upstreamES", "downstreamES"]
        dpsi_df.loc[:, start_col] += 1

        ps_df = dpsi_df.loc[self.ps_idx, rMATS_POS_COLS[self.etype]["+"]]
        ns_df = dpsi_df.loc[self.ns_idx, rMATS_POS_COLS[self.etype]["-"]]
        ps_df.columns = range(len(rMATS_POS_COLS[self.etype]["+"]))
        ns_df.columns = range(len(rMATS_POS_COLS[self.etype]["-"]))
        self.df_ss = concat([ps_df, ns_df]).loc[dpsi_df.index, :]


def get_event_coord(event, app):
    if app == "auto":
        app = detect_file_info(event)["app"]
    if app == "SUPPA2":
        coori = SuppaEventCoord(event)
    elif app == "rMATS":
        coori = rMATSEventCoord(event)
    return coori


def get_ss_bed(
    event_file: FilePath, 
    ups_width: int = 150, 
    dws_width: int = 150,
    strand_sp: bool = False,
    sss: bool = False,
    app: str = "auto",
    **kwargs
    ):
    coori = get_event_coord(event_file, app)
    df_dic = coori.get_all_flank_bed(ups_w=ups_width,dws_w=dws_width,sss=sss, **kwargs)
    return df_dic


def get_ss_range(event_file, app):
    from pandas import DataFrame
    
    coori = get_event_coord(event_file, app)
    df_ss = coori.df_ss
    ncol = df_ss.shape[1] - 1
    df_len = DataFrame(index=df_ss.index, columns=range(ncol))
    for idx in range(ncol):
        df_len.iloc[:, idx] = abs(df_ss.iloc[:, idx+1] - df_ss.iloc[:, idx])
    return df_len
