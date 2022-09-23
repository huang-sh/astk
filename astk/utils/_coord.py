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
        df_bed = self.df_temp.copy()
        df_ss = self.df_ss
        slen = SS_SCORE_LEN[self.app][self.etype][idx]
        ps_idx = self.ps_idx
        ns_idx = self.ns_idx
        if sss:
            if slen == 23:
                ups_w, dws_w = 20, 2
            elif slen == 9:
                ups_w, dws_w = 2, 6
        offset = 1 if based0 else 0
        if slen == 23:
            df_bed.loc[ps_idx, "start"] = df_ss.loc[ps_idx, idx] - ups_w - offset
            df_bed.loc[ps_idx, "end"] = df_ss.loc[ps_idx, idx] + dws_w
            df_bed.loc[ns_idx, "start"] = df_ss.loc[ns_idx, idx] - dws_w - offset
            df_bed.loc[ns_idx, "end"] = df_ss.loc[ns_idx, idx] + ups_w
        elif slen == 9:
            df_bed.loc[ps_idx, "start"] = df_ss.loc[ps_idx, idx] - ups_w - offset
            df_bed.loc[ps_idx, "end"] = df_ss.loc[ps_idx, idx] + dws_w

            df_bed.loc[ns_idx, "start"] = df_ss.loc[ns_idx, idx] - dws_w - offset
            df_bed.loc[ns_idx, "end"] = df_ss.loc[ns_idx, idx] + ups_w
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
        for idx, slen in enumerate(slens):
            if slen == 23:
                key = f"A{idx}_3SS"
            elif slen == 9:
                key = f"A{idx}_5SS"
            else:
                continue
            dic[key] = self.get_ss_flank_bed(idx, ups_w, dws_w, sss, based0)
        return dic

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
            row = ei.Chr, ei.gene_id, ei.strand, *coords
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
    app: str = "atuo"
    ):
    coori = get_event_coord(event_file, app)
    df_dic = coori.get_all_flank_bed(ups_w=ups_width,dws_w=dws_width,sss=sss)
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
