from warnings import warn

from astk.event import SuppaEventID
from astk.constant import SSN, SS_SCORE_LEN
from astk.types import FilePath
from astk.utils import sniff_AS_type


class SuppaEventCoord:
    """This class is for suppa2 AS event ID coord extraction
    """
    def __init__(self, event_id):
        self.event = SuppaEventID(event_id)
        self.etype = self.event.AS_type
        self.strand = self.event.strand
        self.coord = self.event.coordinates
    
    def get_alter(self, based0=True):
        s, e = self.event.alter_element_coor
        if based0:
            s -= 1
        return s, e
    
    def get_ss_range(
        self,
        start: int,
        end: int,
        strand_sp: bool = False,
        based0: bool = True,
        **kwargs
    ) -> tuple:
        """start, end index: 1-based
        """
        s_ups = kwargs.get("s_ups", 0)
        s_dws = kwargs.get("s_dws", 0)
        e_ups = kwargs.get("e_ups", 0)
        e_dws = kwargs.get("e_dws", 0)

        if strand_sp and self.strand == "-":
            s = self.coord[::-1][start-1] + s_ups - s_dws
            e = self.coord[::-1][end-1] + e_ups - e_dws
            s, e = e, s
        else:
            s = self.coord[start-1] - s_ups + s_dws
            e = self.coord[end-1] - e_ups + e_dws
        if based0:
            s -= 1
        return s, e

    def get_splicescore_flank(
        self,
        idx: int,
        based0 = True
    ) -> tuple:
        ss_len = SS_SCORE_LEN[self.etype][idx-1]
        if not ss_len:
            msg = f"it doesn't make sense for {self.etype} to choose {idx} when sss=True"
            warn(UserWarning(msg))
        if self.strand == "+":
            gpos = self.coord[idx-1]
            if ss_len == 9:
                s = gpos - 2
                e = gpos + 6
            elif ss_len == 23:
                s = gpos - 20
                e = gpos + 2
            else:
                s = gpos
                e = gpos 
        elif self.strand == "-":
            gpos = self.coord[::-1][idx-1]
            if ss_len == 9:
                e = gpos + 2
                s = gpos - 6
            elif ss_len == 23:
                e = gpos + 20
                s = gpos - 2
            else:
                e = gpos
                s = gpos 
            
        if based0:
            s -= 1
        return s, e

    def get_ss_flank(
        self,
        idx: int,
        ups_w: int = 150,
        dws_w: int = 150,
        strand_sp: bool = False,
        sss: bool = False,
        based0 = True
    ) -> tuple:
        """idx: 1-based
        """ 
        if sss:
            s, e = self.get_splicescore_flank(idx, based0)
            return s, e

        if strand_sp and self.strand == "-":
            gpos = self.coord[::-1][idx-1]
            s = gpos - dws_w            
            e = gpos + ups_w
        else:
            gpos = self.coord[idx-1]
            s = gpos - ups_w   
            e = gpos + dws_w

        if based0:
            s -= 1
        return s, e
    
    def get_ss_flank_bed(
        self,
        idx: int,
        ups_w: int = 150,
        dws_w: int = 150,
        strand_sp: bool = False,
        sss: bool = False,
        based0 = True,
        name_col = True,
        score_col = True,
        strand_col = True
    ) -> tuple:
        s, e = self.get_ss_flank(idx, ups_w, dws_w, strand_sp, sss, based0)
        if all([name_col ,score_col, strand_col]):
            bed = (self.event.Chr, s, e, self.event.event_id, 0, self.strand)
        else:
            bed = (self.event.Chr, s, e)
        return bed


def get_file_ss_bed(
    event_file: FilePath, 
    ups_width: int = 150, 
    dws_width: int = 150,
    strand_sp: bool = False,
    sss: bool = False
    ):
    from pandas import DataFrame, read_csv

    dpsi_df = read_csv(event_file, sep="\t", index_col=0)
    dpsi_df.drop_duplicates(inplace=True)
    dpsi_df["event_id"] = dpsi_df.index
    dic = {}
    etype = sniff_AS_type(event_file)

    def _func(eid, idx, strand_sp, sss, ups_w=ups_width, dws_w=dws_width):
        ci = SuppaEventCoord(eid)
        coords = ci.get_ss_flank_bed(idx, ups_w=ups_w, dws_w=dws_w, strand_sp=strand_sp, sss=sss)
        return coords

    for i in range(1, SSN[etype]+1):
        coors = dpsi_df["event_id"].apply(_func, args=(i, strand_sp, sss))
        dic[f"A{i}"] = DataFrame(coors.tolist())
    return dic
