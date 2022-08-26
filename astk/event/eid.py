

class SuppaEventID:
    """This class is for suppa2 AS event ID parsing
    """

    def __init__(self, event_id):
        self.event_id = event_id
        self.parse_id(event_id)
        self.parse_AS_coor()

    def parse_id(self, event_id):
        gene_id, tmp = event_id.split(";")
        AS_type, Chr, *coordinates, strand = tmp.split(":")
        self.gene_id = gene_id
        self.AS_type = AS_type
        self.Chr = Chr
        self.strand = strand
        self.coordinate_str = "-".join(coordinates)

    @property
    def alter_element_len(self):
        alter_len = self.alter_element_coor[1] - self.alter_element_coor[0]
        return alter_len

    def parse_AS_coor(self): 

        AS_func_dic = {
            "SE": self._parse_SE_coor,
            "MX": self._parse_MX_coor,
            "A5": self._parse_SS_coor,
            "A3": self._parse_SS_coor,
            "RI": self._parse_RI_coor,
            "AF": self._parse_FL_coor,
            "AL": self._parse_FL_coor
        }
        self.coordinates = AS_func_dic[self.AS_type]()

    def _parse_SE_coor(self): 
        nums = self.coordinate_str.split("-")
        e1, s2, e2, s3 = map(int, nums)
        self.alter_element_coor = (s2, e2)
        return (e1, s2, e2, s3)
    
    def _parse_MX_coor(self):
        nums = self.coordinate_str.split("-")
        e1, s2, e2, s4, _, s3, e3, _ = map(int, nums)
        self.alter_element_coor = (s2, e2)
        return (e1, s2, e2, s3, e3, s4)

    def _parse_SS_coor(self):
        nums = self.coordinate_str.split("-")
        A5_fs = all([self.strand == "+", self.AS_type == "A5"])
        A3_rs = all([self.strand == "-", self.AS_type == "A3"])

        A5_rs = all([self.strand == "-", self.AS_type == "A5"])
        A3_fs = all([self.strand == "+", self.AS_type == "A3"])

        if A5_fs or A3_rs:
            e2, s3, e1, _ = map(int, nums)
            coors = (e1, e2, s3)
            self.alter_element_coor = (e1, e2)
        elif A5_rs or A3_fs:
            e1, s2, _, s3 = map(int, nums)
            coors = (e1, s2, s3)
            self.alter_element_coor = (s2, s3)
        else:
            coors = ()
        return coors

    def _parse_RI_coor(self):
        nums = self.coordinate_str.split("-")
        s1, e1, s2, e2 = map(int, nums)
        self.alter_element_coor = (e1, s2)
        return (s1, e1, s2, e2)

    def _parse_FL_coor(self):
        nums = self.coordinate_str.split("-")
        AF_fs = all([self.strand == "+", self.AS_type == "AF"])
        AL_rs = all([self.strand == "-", self.AS_type == "AL"])

        AF_rs = all([self.strand == "-", self.AS_type == "AF"])
        AL_fs = all([self.strand == "+", self.AS_type == "AL"])

        if AF_fs or AL_rs:
            s1, e1, s3, s2, e2, _ = map(int, nums)
            coors = (s1, e1, s2, e2, s3)
            self.alter_element_coor = (s1, e1)
        elif AF_rs or AL_fs:
            e1, s2, e2, _, s3, e3 = map(int, nums)
            coors = (e1, s2, e2, s3, e3)
            self.alter_element_coor = (s3, e3)
        else:
            coors = ()
        return coors

