



class EventID:

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


##### they are old function for length parsing

def SE_exon_len(event_id):
    start, end = event_id.split("-")[1].split(":")
    return int(end) - int(start)

def MX_exon_len(event_id):
    _, exon1, _, exon2 = event_id.split("-")[:4]
    s1, e1 = exon1.split(":")
    s2, e2 = exon2.split(":")
    len1 = int(e1) - int(s1)
    len2 = int(e2) - int(s2)
    return len1, len2

def A5_exon_len(event_id):
    splice1, splice2 = event_id.split(":")[2: 4]
    if event_id.endswith("+"):
        end = int(splice1.split("-")[0])
        start = int(splice2.split("-")[0])
    elif event_id.endswith("-"):
         end = int(splice2.split("-")[1])
         start = int(splice1.split("-")[1])
    return end - start

def A3_exon_len(event_id):
    splice1, splice2 = event_id.split(":")[2: 4]
    if event_id.endswith("+"):
        end = int(splice2.split("-")[1])
        start = int(splice1.split("-")[1])
    elif event_id.endswith("-"):
         end = int(splice1.split("-")[0])
         start = int(splice2.split("-")[0])
    return end - start

def RI_intron_len(event_id):
    start, end = event_id.split(":")[3].split("-")
    return int(end) - int(start)

def AF_exon_len(event_id):
    if event_id.endswith("+"):
        s1, e1 = event_id.split("-")[0].split(":")[-2:]
        s2, e2 = event_id.split("-")[1].split(":")[1:]
    elif event_id.endswith("-"):
        s1, e1 = event_id.split("-")[-2].split(":")[:2]
        s2, e2 = event_id.split("-")[1].split(":")[:2]
    len1 = int(e1) - int(s1)
    len2 = int(e2) - int(s2)
    return len1, len2

def AL_exon_len(event_id):
    if event_id.endswith("+"):
        s1, e1 = event_id.split("-")[2].split(":")[:2]
        s2, e2 = event_id.split("-")[1].split(":")[:2]
    elif event_id.endswith("-"):
        s1, e1 = event_id.split("-")[0].split(":")[-2:]
        s2, e2 = event_id.split("-")[1].split(":")[1:]
    len1 = int(e1) - int(s1)
    len2 = int(e2) - int(s2)
    return len1, len2

def AS_len(event_id):
    as_type = event_id.split(";")[1][:2]
    len_func_dic = {
        "SE": SE_exon_len,
        "MX": MX_exon_len,
        "A5": A5_exon_len,
        "A3": A3_exon_len,
        "RI": RI_intron_len,
        "AF": AF_exon_len,
        "AL": AL_exon_len
    }
    func = len_func_dic[as_type]
    f_len = func(event_id)
    if isinstance(f_len, tuple) and len(f_len) == 2:
        res = {"len": f_len[0], "len2":  f_len[1]}
    elif isinstance(f_len, int):
        res = {"len": f_len}
    return res


if "__main__" == "__name__":
    se_test = "ENSMUSG00000024169.17;SE:chr17:25244920-25247778:25247953-25251059:+"
    mx_test = "ENSMUSG00000000037.18;MX:chrX:160004745-160006153:160006236-160008873:160004745-160007457:160007540-160008873:+"

    a5_test_f = "ENSMUSG00000000134.18;A5:chrX:7629452-7631746:7629299-7631746:+"
    a5_test_r = "ENSMUSG00000000131.16;A5:chr7:125708072-125709086:125708072-125709126:-"

    a3_test_f = "ENSMUSG00000000127.16;A3:chr17:64280152-64288475:64280152-64288478:+"
    a3_test_r = "ENSMUSG00000000131.16;A3:chr7:125704106-125704273:125703961-125704273:-"

    ri_test = "ENSMUSG00000000031.17;RI:chr7:142129268:142129851-142129940:142130267:-"

    af_test_f = "ENSMUSG00000000049.12;AF:chr11:108234180:108234350-108286128:108271990:108272222-108286128:+"
    af_test_r = "ENSMUSG00000000103.13;AF:chrY:2150121-2150216:2150346:2150121-2170149:2170409:-"

    al_test_f = "ENSMUSG00000000037.18;AL:chrX:160029363-160038026:160038655:160029363-160039562:160039882:+"
    al_test_r = "ENSMUSG00000000384.16;AL:chr11:6570002:6570143-6570742:6570180:6570235-6570742:-"

    a = EventID(al_test_r)
    print(a.Chr, a.gene_id, a.coordinates, a.strand, a.alter_element_len)

    # a = SE_exon_len(se_test)
    # b = MX_exon_len(mx_test)
    # c_f = A5_exon_len(a5_test_f)
    # c_r = A5_exon_len(a5_test_r)

    # d_f = A3_exon_len(a3_test_f)
    # d_r = A3_exon_len(a3_test_r)

    # e = RI_intron_len(ri_test)

    # f_f = AF_exon_len(af_test_f)
    # f_r = AF_exon_len(af_test_r)

    # g_f = AL_exon_len(al_test_f)
    # g_r = AL_exon_len(al_test_r)

    # print(a)
    # print(b)
    # print(c_f, c_r)
    # print(d_f, d_r)
    # print(e)

    # print(f_f)
    # print(f_r)

    # print(g_f, g_r)