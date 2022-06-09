"""
Modified version of suppa.lib.gtf_store module.
"""

import sys



class Genome:
    def __init__(self):
        self.genes = {}

    def add_gene(self, gene):
        self.genes.setdefault(gene.id, gene)
    
    def add_transcript(self, tx, gene_id):
        self.genes[gene_id].add_transcript(tx)

    def add_exon(self, exon, gene_id, transcript_id):
        self.genes[gene_id].transcripts[transcript_id].add_exon(exon)

    def sort_transcripts(self):
        """
        For each gene in genome sorts transcripts (for further
        precessing).
        """
        for gene, _, _ in self.__iter__():
            gene.sort_transcripts()


class Gene:
    __slots__ = ("seqname", "chr", "start", "end", "strand", "id", "attr",
                "transcripts")

    def __init__(self, id, seqname, start, end, strand, **kwargs):
        self.seqname = seqname
        self.chr = seqname
        self.start = start
        self.end = end
        self.strand = strand
        self.id = id
        self.attr = kwargs
        self.transcripts = {}

    def add_transcript(self, transcript):
        self.transcripts.setdefault(transcript.id, transcript)

    def sort_transcripts(self):
        txs = self.transcripts
        tx_ids = sorted(txs, key=lambda tx: tx.start)
        self.transcripts = {txs[tx_id] for tx_id in tx_ids}
    
    def __str__(self):
        return self.id


class Transcript:
    __slots__ = ("seqname", "chr", "start", "end", "strand", "id", "attr",
                "transcripts", "exons", "TSS", "cds")

    def __init__(self, id, seqname, start, end, strand, **kwargs):
        self.chr = seqname
        self.seqname = seqname
        self.start = start
        self.end = end
        self.strand = strand
        self.id = id
        self.attr = kwargs
        self.transcripts = {}
        self.exons = {}
        self.TSS = start if strand == "+" else end
        self.cds = (float('inf'), float('-inf'))

    def add_exon(self, exon):
        """
        Adds exon to transcript
        """
        self.exons.setdefault(exon.id, exon)
        self.cds = (min([self.cds[0], exon.start]), max([self.cds[1], exon.end]))

    def sort_exons(self):
        exons = self.exons
        exon_ids = sorted(exons, key=lambda exon: exons[exon].start)
        self.exons = {exon_id: exons[exon_id] for exon_id in exon_ids}
    
    def get_exon_coordinates(self):
        self.sort_exons()
        coordinates = ((eid, (ec.start, ec.end)) for eid, ec in self.exons.items())
        return coordinates
    

from collections import namedtuple

Exon = namedtuple('Exon',('id, chr, start, end, strand, attr'))

# class Exon:
#     # __slots__ = ("seqname", "chr", "start", "end", "strand", "id", "attr")

#     def __init__(self, id, seqname, start, end, strand, **kwargs):
#         self.seqname = seqname
#         self.chr = seqname
#         self.start = start
#         self.end = end
#         self.strand = strand
#         self.id = id
#         self.attr = kwargs


GTF_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]


# def parse_gtf_line(line, idx):
#     line = line.strip().split('\t')
#     if len(line) != 9:
#         print(f'Missmatch in number of Fields. skipping line {idx}')
#         return    
#     info_dic = {k:v for k, v in zip(GTF_COLUMNSv[:8], line[:8])}
#     attributes = [att.split() for att in line[-1][:-1].split('; ')]
#     for att in attributes:
#         info_dic.update({att[0]: att[1].strip('"')})
#     return info_dic

def parse_gtf_line(line, feature=["gene", "exon", "transcript"]):
    line = line.strip().split('\t')
    if len(line) != 9:
        print('Missmatch in number of Fields. skipping line')
        return 
    if line[2] not in feature:
        return   
    info_dic = {k:v for k, v in zip(GTF_COLUMNS[:8], line[:8])}
    attributes_str = line[-1][:-1].replace('"', "")
    attributes_ls = [att.split() for att in attributes_str.split('; ')]
    info_dic.update(dict(attributes_ls))
    return info_dic

import mmap

def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def construct_genome(gtf):
    from tqdm import tqdm

    genome = Genome()
    num_lines = get_num_lines(gtf)
    with open(gtf, 'r') as handle:
        for line in tqdm(handle, desc="Parsing GTF annotation", total=num_lines):
            if line.startswith('#'):
                continue
            dic = parse_gtf_line(line)
            if not dic:
                continue
            seqname = dic.pop("seqname")
            start = int(dic.pop("start"))
            end = int(dic.pop("end"))
            strand = dic.pop("strand")        
            if dic["feature"] == "gene":
                gene_id = dic.pop("gene_id")
                # c += 1
                gene = Gene(gene_id, seqname, start, end, strand, **dic)
                genome.add_gene(gene)
                # genome.genes.setdefault(gene.id, gene)
                
            elif dic["feature"] == "transcript":
                transcript_id = dic.pop("transcript_id")
                tx = Transcript(transcript_id, seqname, start, end, strand, **dic)
                genome.add_transcript(tx, tx.attr["gene_id"])
                # genome.genes[tx.attr["gene_id"]].add_transcript(tx)

            if dic["feature"] == "exon":
                exon_id = dic.pop("exon_id")
                exon = Exon(exon_id, seqname, start, end, strand, dic)
                # genome.add_exon(exon, gene_id, exon.attr["transcript_id"])
                genome.genes[gene_id].transcripts[exon.attr["transcript_id"]].add_exon(exon)
    return genome