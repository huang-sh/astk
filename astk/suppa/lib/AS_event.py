"""
Modified version of suppa.lib.event module.
"""
from pathlib import Path
from typing import TypeVar
from typing import Sequence
from typing import Optional
# Gene = TypeVar("Gene")

from .gtf_parse import Exon, Gene, Transcript
from .tools import EWriter


from typing import Sequence

class ASEvent:

    def __init__(self,
        splicein_exon: Sequence["Exon"] = [],
        spliceout_exon: Sequence["Exon"] = [],
        splicein_tx: Sequence["Transcript"] = [],
        spliceout_tx: Sequence["Transcript"] = [],       
    ):
        self.splicein_exon = splicein_exon
        self.spliceout_exon = spliceout_exon
        self.splicein_tx = set(splicein_tx)
        self.spliceout_tx = set(spliceout_tx) 
        self.transcript = {}
        if splicein_exon and spliceout_exon:
            self.validate()

    def set_splicein_exon(self, exons: Sequence["Exon"]):
        self.splicein_exon = exons
        if self.splicein_exon and self.spliceout_exon:
            return self.validate()

    def set_spliceout_exon(self, exons: Sequence["Exon"]):
        self.spliceout_exon = exons
        if self.splicein_exon and self.spliceout_exon:
            return self.validate()

    def add_splicein_tx(self, transcript: "Transcript"):
        self.splicein_tx.add(transcript)

    def add_spliceout_tx(self, transcript: "Transcript"):
        self.spliceout_tx.add(transcript)
    
    def validate(self):
        pass

    def set_basic_info(self, exons: Sequence["Exon"]):
        pass


class SkippingExonEvent(ASEvent):

    def __init__(
        self, 
        splicein_exon: Sequence["Exon"] = [],
        spliceout_exon: Sequence["Exon"] = [],
        splicein_tx: Sequence["Transcript"] = [],
        spliceout_tx: Sequence["Transcript"] = [],        
    ):
        super().__init__(
            splicein_exon,
            spliceout_exon,
            splicein_tx,
            spliceout_tx
        )
        self.etype = 'SE'
        self.status = "infer"
        if splicein_exon:
            self.set_basic_info(splicein_exon)

    def validate(self):
        si1, si2, si3 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        fbool = [
                si1.end == so1.end, 
                si3.start == so2.start,
                so1.end < si2.start < si2.end < so2.start
        ]
        if all(fbool):        
            self.status = "review"
            return self

    def set_basic_info(self, exons: Sequence["Exon"]):
        e1, e2, e3 = exons
        self.strand =  e1.strand
        self.id = f"{e1.chr}:{e1.end}-{e2.start}:{e2.end}-{e3.start}:{e1.strand}"


class RetainedIntronEvent(ASEvent):

    def __init__(
        self, 
        splicein_exon: Sequence["Exon"] = [],
        spliceout_exon: Sequence["Exon"] = [],
        splicein_tx: Sequence["Transcript"] = [],
        spliceout_tx: Sequence["Transcript"] = [],        
    ):
        super().__init__(
            splicein_exon,
            spliceout_exon,
            splicein_tx,
            spliceout_tx
        )
        self.etype = 'RI'
        self.status = "infer"
        if spliceout_exon:
            self.set_basic_info(spliceout_exon)
       
    def validate(self):
        e1, e2 = self.spliceout_exon
        e = self.splicein_exon[0] 
        fbool = [e.start == e1.start, e.end == e2.end, e1.end < e2.start]
        if all(fbool):
            self.status = "review"
            return self
    
    def set_basic_info(self, exons: Sequence["Exon"]):
        e1, e2 = exons
        self.strand =  e1.strand
        self.id = f"{e1.chr}:{e1.start}:{e1.end}-{e2.start}:{e2.end}:{e2.strand}"


class AlternativeSplicing:

    def __init__(self, gene: "Gene"):
        self.gene = gene
        self.gene_id = gene.id
        self.chr = gene.chr
        self.AS_events = {}
        self.start_codon_events = {}
        self.stop_codon_events = {}

    def construct_events(self):
        for transcript in self.gene.transcripts.values():
            self.add_putative_events(transcript)
        for transcript in self.gene.transcripts.values():
            self.add_real_events(transcript)

    def clear_invlaid(self):
        events_ids = list(self.AS_events.keys())
        for eid in events_ids:
            event = self.AS_events[eid]
            if event.status == "infer":
                del self.AS_events[eid]
        self.inner_AS_events = self.AS_events.copy()

    def filter_promoter_event(self):
        """filter AS events that exon overlap with promoter
        """
        self.promoter_event = {}
        events_id = list(self.inner_AS_events.keys())
        for eid in events_id:
            event = self.inner_AS_events[eid]
            for tx in event.splicein_tx:
                if event.strand == "+":
                    if int(tx.TSS) + 2000 > int(event.splicein_exon[0].start):
                        self.promoter_event.setdefault(eid, self.inner_AS_events.pop(eid))
                        break
                else:
                    if int(event.splicein_exon[-1].end) + 2000 > int(tx.TSS):
                        self.promoter_event.setdefault(eid, self.inner_AS_events.pop(eid))
                        break

    def filter_start_codon(self):
        """filter AS events that exon overlap with promoter
        """
        events_id = list(self.inner_AS_events.keys())
        for eid in events_id:
            event = self.inner_AS_events[eid]
            for tx in event.splicein_tx:
                if event.strand == "+":
                    if event.splicein_exon[0].start == tx.cds[0]:
                        if eid not in self.start_codon_events:
                            self.start_codon_events.setdefault(eid, self.inner_AS_events.pop(eid))
                else:
                    if event.splicein_exon[-1].end == tx.cds[1]:
                        if eid not in self.start_codon_events:
                            self.start_codon_events.setdefault(eid, self.inner_AS_events.pop(eid))
                
    def filter_stop_codon(self):
        """filter AS events that exon overlap with promoter
        """
        events_id = list(self.inner_AS_events.keys())
        for eid in events_id:
            event = self.inner_AS_events[eid]
            for tx in event.splicein_tx:
                if event.strand == "+":
                    if event.splicein_exon[-1].end == tx.cds[1]:
                        if eid not in self.stop_codon_events:
                            self.stop_codon_events.setdefault(eid, self.inner_AS_events.pop(eid))
                else:
                    if event.splicein_exon[0].start == tx.cds[0]:
                        if eid not in self.stop_codon_events:
                            self.stop_codon_events.setdefault(eid, self.inner_AS_events.pop(eid))

    def to_ioe_list(self, event: "ASEvent"):
        # print(event)
        splicein_tx = [i.id for i in event.splicein_tx]
        spliceout_tx = [i.id for i in event.spliceout_tx]
        all_tx = splicein_tx + spliceout_tx
        full_event_id = f"{self.gene_id};{self.etype}:{event.id}"
        line = [self.chr, self.gene_id, full_event_id, ",".join(splicein_tx), ",".join(all_tx)]
        return line

    def to_ioe(self, type: Optional[str] = None):
        if type is None:
            events = self.AS_events.copy()
            events.update(self.start_codon_event)
            events.update(self.stop_codon_event)
        elif type == "start_codon":
            events = self.start_codon_events
        elif type == "stop_codon":
            events = self.stop_codon_events
        elif type == "inner":
            events = self.inner_AS_events
        else:
            events = {}
        lines = [self.to_ioe_list(e) for e in events.values()]
        return lines

        # data = {"seqname":[], "gene_id":[], "event_id":[], "alternative_transcripts":[], "total_transcripts": []}
        # if self.AS_events.values():
        #     ioe_df = pd.DataFrame([self.to_ioe_list(event) for event in self.AS_events.values()])
        #     ioe_df.columns = list(data.keys())
        # else:
        #     ioe_df = pd.DataFrame(data)
        # return ioe_df

    def add_putative_events(self, transcript):
        pass

    def add_real_events(self, transcript):
        pass


class SkippingExon(AlternativeSplicing):
    """
    Skipped exon event
    """
    def __init__(self, gene):
        super().__init__(gene)
        self.etype = 'SE'  
        self.construct_events()
        self.clear_invlaid()

    def add_putative_events(self, transcript: "Transcript"):
        exons = sorted(transcript.exons.values(), key=lambda e: e.start)
        if len(exons) < 3:
            return
        for i in range(len(exons) - 2):
            event = SkippingExonEvent(splicein_exon=exons[i:i+3])
            self.AS_events.setdefault(event.id, event)
            self.AS_events[event.id].add_splicein_tx(transcript)

    def add_real_events(self, transcript: "Transcript"):
        exons = sorted(transcript.exons.values(), key=lambda e: e.start)
        if len(exons) < 2:
            return
        for eid, event in self.AS_events.items():
            for i in range((len(exons)-1)):
                if event.set_spliceout_exon(exons[i:i+2]):
                    event.add_spliceout_tx(transcript)
                    break


class RetainedIntron(AlternativeSplicing):
    """
    Retained Intron event
    """
    def __init__(self, gene: "Gene"):
        super().__init__(gene)
        self.etype = 'RI'  
        self.construct_events()
        self.clear_invlaid()

    def add_putative_events(self, transcript: "Transcript"):
        exons = sorted(transcript.exons.values(), key=lambda e: e.start)
        for i in range(len(exons) - 1):
            event = RetainedIntronEvent(spliceout_exon=exons[i:i+2])
            self.AS_events.setdefault(event.id, event)
            self.AS_events[event.id].add_spliceout_tx(transcript)

    def add_real_events(self, transcript: "Transcript"):
        for event in self.AS_events.values():
            for exon in transcript.exons.values():
                if event.set_splicein_exon([exon]):
                    event.add_splicein_tx(transcript)
                    break


def process_events(my_gene, event_cls, output, promoterSplit):
    """
    Generates specified event occurrences in gene and writes the GTF/IOE output
    """
    gene_event = event_cls(my_gene)
    print(gene_event.AS_events)
    # if promoterSplit:
    #     gene_event.filter_promoter_event()
    #     pout = Path(output).with_suffix(".tss.ioe")
    #     gene_event.export_promoter_events_ioe(pout)

    #     output = Path(output).with_suffix(".gb.ioe")
    # gene_event.export_events_ioe(output)


def make_events(events, genome, output, AS_split, b_type="S", th=10):
    """
    Controls event creation and writing for each gene
    """
    boundary = 'variable_{}'.format(th) if b_type == 'V' else 'strict'

    # my_ioe_writer = EWriter(events, output_name, 'ioe', boundary)

    event_cls_dic = {"SE": SkippingExon, "RI": RetainedIntron}
    event_cls_dic = {k: v for k, v in event_cls_dic.items() if k in events}
    import time
    T1 = time.time()
    handle_ls = []
    if "startCodon" in AS_split:
        sc_ioe_writer = EWriter(event_cls_dic.keys(), f"{output}_start", "ioe", boundary)
        handle_ls.append(sc_ioe_writer)
    if "stopCodon" in AS_split:
        se_ioe_writer = EWriter(event_cls_dic.keys(), f"{output}_stop", "ioe", boundary)
        handle_ls.append(se_ioe_writer)
    if "inner" in AS_split:
        mid_ioe_writer = EWriter(event_cls_dic.keys(), f"{output}_inner", "ioe", boundary)
        handle_ls.append(mid_ioe_writer)
    if not AS_split:
        ioe_writer = EWriter(event_cls_dic.keys(), output, "ioe", boundary)
        handle_ls.append(ioe_writer)

    from tqdm import tqdm
    for t, event_cls in event_cls_dic.items():
        for gene in tqdm(genome.genes.values(), desc=f"Calculating {t} events"):
            gene_event = event_cls(gene)
            gene_event.filter_start_codon()
            gene_event.filter_stop_codon()

            if "startCodon" in AS_split:    
                for line in gene_event.to_ioe(type="start_codon"):
                    sc_ioe_writer.write("\t".join(line), t)
            if "stopCodon" in AS_split:
                for line in gene_event.to_ioe(type="stop_codon"):
                    se_ioe_writer.write("\t".join(line), t)
            if "inner" in AS_split:
                for line in gene_event.to_ioe(type="inner"):
                    mid_ioe_writer.write("\t".join(line), t)
            if not AS_split:
                for line in gene_event.to_ioe():
                    ioe_writer.write("\t".join(line), t)
    
    for i in handle_ls:
        i.close()
    

        # ioe_dfs = [event_cls(gene).to_ioe() for gene in genome.genes.values()]
        # df = pd.concat(ioe_dfs)
        # output = Path(output).with_suffix(".ioe")
        # df.to_csv(output, sep="\t")

        # for gene in tqdm(genome.genes.values(), desc=f"Calculating {t} events"):
        #     # print(gene)
        #     process_events(gene, event_cls, output, promoterSplit)
    # my_ioe_writer.close()
    # T2 =time.time()
    # print('程序运行时间:%s秒' % ((T2 - T1)))


# class ARSS(Event):
#     """
#     Alternative 'Right' Splice Site
#     """
#     def __init__(self, gene, *_):
#         super().__init__(self, gene)
#         if self.gene.strand == '-':
#             self.etype = 'A3'
#         else:
#             self.etype = 'A5'
#         self.positive_coords = {}
#         self.negative_coords = {}

#         self.construct_events(gene)

#     def add_putative_events(self, exons, transcript):
#         if len(exons) < 2:
#             return
#         for i in range(len(exons) - 1):
#             firstexon = exons[i]
#             secondexon = exons[i + 1]
#             search_coords = (firstexon[1], secondexon[0])
#             self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

#     def add_real_events(self, exons, transcript):
#         if len(exons) < 2:
#             return
#         for i in range(len(exons) - 1):
#             firstexon = exons[i]
#             secondexon = exons[i + 1]
#             self.create_arss(firstexon[0], firstexon[1], secondexon[0], transcript)

#     def create_arss(self, coord1, coord2, coord3, transcript):
#         gene = self.gene
#         for neg1, neg2 in self.negative_coords:
#             if neg2 == coord3 and (coord1 < neg1 < coord2):
#                 eventid = '{}:{}-{}:{}-{}:{}'.format(gene.chr, coord2, coord3, neg1, neg2, gene.strand)
#                 self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
#                 self.negative_ids[eventid] = self.negative_coords[(neg1, neg2)]

#     def export_events_gtf(self, edge):
#         """
#         Generator of GTF lines
#         """
#         strand = self.gene.strand
#         for event in self.positive_ids:
#             full_event = '{}:{}'.format(self.etype, event)
#             e_vals = full_event.replace('-', ':').split(':')

#             line1 = self.gtf_string.format(int(e_vals[4]) - edge, e_vals[2], strand, full_event, full_event,
#                                            'alternative1')
#             yield line1, self.etype

#             line2 = self.gtf_string.format(e_vals[3], int(e_vals[3]) + edge, strand, full_event, full_event,
#                                            'alternative1')
#             yield line2, self.etype

#             line3 = self.gtf_string.format(int(e_vals[4]) - edge, e_vals[4], strand, full_event, full_event,
#                                            'alternative2')
#             yield line3, self.etype

#             line4 = self.gtf_string.format(e_vals[5], int(e_vals[5]) + edge, strand, full_event, full_event,
#                                            'alternative2')
#             yield line4, self.etype


# class ALSS(Event):
#     """
#     Alternative 'Left' Splice Site
#     """
#     def __init__(self, gene, *_):
#         super().__init__(self, gene)
#         if self.gene.strand == '-':
#             self.etype = 'A5'
#         else:
#             self.etype = 'A3'
#         self.positive_coords = {}
#         self.negative_coords = {}

#         self.construct_events(gene)

#     def add_putative_events(self, exons, transcript):
#         if len(exons) < 2:
#             return
#         for i in range(len(exons) - 1):
#             firstexon = exons[i]
#             secondexon = exons[i + 1]
#             search_coords = (firstexon[1], secondexon[0])
#             self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

#     def add_real_events(self, exons, transcript):
#         if len(exons) < 2:
#             return
#         for i in range(len(exons) - 1):
#             firstexon = exons[i]
#             secondexon = exons[i + 1]
#             self.create_alss(firstexon[1], secondexon[0], secondexon[1], transcript)

#     def create_alss(self, coord1, coord2, coord3, transcript):
#         gene = self.gene
#         for neg1, neg2 in self.negative_coords:
#             if neg1 == coord1 and (coord2 < neg2 < coord3):
#                 eventid = '{}:{}-{}:{}-{}:{}'.format(gene.chr, coord1, coord2, neg1, neg2, gene.strand)
#                 self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
#                 self.negative_ids[eventid] = self.negative_coords[(neg1, neg2)]

#     def export_events_gtf(self, edge):
#         """
#         Generator of GTF lines
#         """
#         strand = self.gene.strand
#         for event in self.positive_ids:
#             full_event = '{}:{}'.format(self.etype, event)
#             e_vals = full_event.replace('-', ':').split(':')

#             line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
#                                            'alternative1')
#             yield line1, self.etype

#             line2 = self.gtf_string.format(e_vals[3], int(e_vals[5]) + edge, strand, full_event, full_event,
#                                            'alternative1')
#             yield line2, self.etype

#             line3 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
#                                            'alternative2')
#             yield line3, self.etype

#             line4 = self.gtf_string.format(e_vals[5], int(e_vals[5]) + edge, strand, full_event, full_event,
#                                            'alternative2')
#             yield line4, self.etype


# class ALTL(Event):
#     """
#     Alternative on 'left' side of the gene
#     """
#     def __init__(self, gene, *_):
#         super().__init__(self, gene)
#         if self.gene.strand == '-':
#             self.etype = 'AL'
#         else:
#             self.etype = 'AF'
#         self.positive_coords = {}
#         self.negative_coords = {}

#         self.construct_events(gene)

#     def add_putative_events(self, exons, transcript):
#         if len(exons) < 2:
#             return
#         firstexon = exons[0]
#         secondexon = exons[1]
#         search_coords = (firstexon[0], firstexon[1], secondexon[0])
#         self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

#     def add_real_events(self, exons, transcript):
#         if len(exons) < 2:
#             return
#         firstexon = exons[0]
#         secondexon = exons[1]
#         self.create_altl(firstexon[0], firstexon[1], secondexon[0], transcript)

#     def create_altl(self, coord1, coord2, coord3, transcript):
#         """
#         Find and create alternative events
#         """
#         gene = self.gene
#         for neg1, neg2, neg3 in self.negative_coords:
#             if neg3 == coord3 and coord2 < neg1:
#                 eventid = '{}:{}:{}-{}:{}:{}-{}:{}'.format(gene.chr, coord1, coord2, coord3,
#                                                            neg1, neg2, neg3, gene.strand)
#                 self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
#                 self.negative_ids[eventid] = self.negative_coords[(neg1, neg2, neg3)]

#     def export_events_gtf(self, edge):
#         """
#         Generator of GTF lines
#         """
#         strand = self.gene.strand
#         for event in self.positive_ids:
#             full_event = '{}:{}'.format(self.etype, event)
#             e_vals = full_event.replace('-', ':').split(':')

#             line1 = self.gtf_string.format(e_vals[2], e_vals[3], strand, full_event, full_event,
#                                            'alternative1')
#             yield line1, self.etype

#             line2 = self.gtf_string.format(e_vals[4], int(e_vals[4]) + edge, strand, full_event, full_event,
#                                            'alternative1')
#             yield line2, self.etype

#             line3 = self.gtf_string.format(e_vals[5], e_vals[6], strand, full_event, full_event,
#                                            'alternative2')
#             yield line3, self.etype

#             line4 = self.gtf_string.format(int(e_vals[7]), int(e_vals[7]) + edge, strand, full_event, full_event,
#                                            'alternative2')
#             yield line4, self.etype


# class ALTR(Event):
#     """
#     Alternative on 'Right' side of the gene
#     """
#     def __init__(self, gene, *_):
#         super().__init__(self, gene)
#         if self.gene.strand == '-':
#             self.etype = 'AF'
#         else:
#             self.etype = 'AL'
#         self.positive_coords = {}
#         self.negative_coords = {}

#         self.construct_events(gene)

#     def add_putative_events(self, exons, transcript):
#         if len(exons) < 2:
#             return
#         firstexon = exons[-2]
#         secondexon = exons[-1]
#         search_coords = (firstexon[1], secondexon[0], secondexon[1])
#         self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

#     def add_real_events(self, exons, transcript):
#         if len(exons) < 2:
#             return
#         firstexon = exons[-2]
#         secondexon = exons[-1]
#         self.create_altr(firstexon[1], secondexon[0], secondexon[1], transcript)

#     def create_altr(self, coord1, coord2, coord3, transcript):
#         """
#         Find and create alternative events
#         """
#         gene = self.gene
#         for neg1, neg2, neg3 in self.negative_coords:
#             if neg1 == coord1 and neg2 > coord3:
#                 eventid = '{}:{}-{}:{}:{}-{}:{}:{}'.format(gene.chr, coord1, coord2, coord3,
#                                                            neg1, neg2, neg3, gene.strand)
#                 self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
#                 self.negative_ids[eventid] = self.negative_coords[(neg1, neg2, neg3)]

#     def export_events_gtf(self, edge):
#         """
#         Generator of GTF lines
#         """
#         strand = self.gene.strand
#         for event in self.positive_ids:
#             full_event = '{}:{}'.format(self.etype, event)
#             e_vals = full_event.replace('-', ':').split(':')

#             line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
#                                            'alternative2')
#             yield line1, self.etype

#             line2 = self.gtf_string.format(e_vals[3], e_vals[4], strand, full_event, full_event,
#                                            'alternative2')
#             yield line2, self.etype

#             line3 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
#                                            'alternative1')
#             yield line3, self.etype

#             line4 = self.gtf_string.format(e_vals[6], e_vals[7], strand, full_event, full_event,
#                                            'alternative1')
#             yield line4, self.etype


# class MXE(Event):
#     """
#     Mutually exclusive exons event
#     """
#     def __init__(self, gene, *_):
#         super().__init__(self, gene)
#         self.etype = 'MX'
#         self.positive_coords = {}
#         self.negative_coords = {}

#         self.construct_events(gene)

#     def add_putative_events(self, exons, transcript):
#         if len(exons) < 3:
#             return
#         for i in range(len(exons) - 2):
#             firstexon = exons[i]
#             midexon = exons[i + 1]
#             lastexon = exons[i + 2]
#             search_coords = firstexon[1], midexon[0], midexon[1], lastexon[0]
#             self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

#     def add_real_events(self, exons, transcript):
#         if len(exons) < 3:
#             return
#         for i in range((len(exons) - 2)):
#             firstexon = exons[i]
#             midexon = exons[i + 1]
#             lastexon = exons[i + 2]
#             search_coords = firstexon[1], midexon[0], midexon[1], lastexon[0]

#             for eventid, neg_coord in self.coord_match(search_coords):
#                 self.negative_ids[eventid] = self.negative_ids.get(eventid, []) + [transcript]
#                 self.positive_ids[eventid] = self.negative_coords[neg_coord]

#     def coord_match(self, new_coords):
#         """
#         Creates and returns event IDs that match tupled with
#         alternative IDs
#         """
#         all_events = []
#         gene = self.gene

#         for orig_cord in self.negative_coords:
#             compare_f = orig_cord[0] + orig_cord[-1]
#             compare_t = new_coords[0] + new_coords[-1]
#             mm_f = orig_cord[2]
#             mm_t = new_coords[1]

#             if compare_f == compare_t and mm_f < mm_t:
#                 eventid = ('{}:{}-{}:{}-{}:{}-{}:{}-{}:{}'
#                            .format(gene.chr, orig_cord[0], orig_cord[1], orig_cord[2], orig_cord[3],
#                                    new_coords[0], new_coords[1], new_coords[2], new_coords[3], gene.strand))
#                 all_events.append([eventid, orig_cord])
#         return all_events

#     def export_events_gtf(self, edge):
#         """
#         Generator of GTF lines
#         """
#         strand = self.gene.strand
#         for event in self.positive_ids:
#             full_event = '{}:{}'.format(self.etype, event)
#             e_vals = full_event.replace('-', ':').split(':')

#             line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
#                                            'alternative1')
#             yield line1, self.etype

#             line2 = self.gtf_string.format(e_vals[3], e_vals[4], strand, full_event, full_event,
#                                            'alternative1')
#             yield line2, self.etype

#             line3 = self.gtf_string.format(e_vals[9], int(e_vals[9]) + edge, strand, full_event, full_event,
#                                            'alternative1')
#             yield line3, self.etype

#             line4 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
#                                            'alternative2')
#             yield line4, self.etype

#             line5 = self.gtf_string.format(e_vals[7], e_vals[8], strand, full_event, full_event,
#                                            'alternative2')
#             yield line5, self.etype

#             line6 = self.gtf_string.format(e_vals[9], int(e_vals[9]) + edge, strand, full_event, full_event,
#                                            'alternative2')
#             yield line6, self.etype


# class RI(Event):
#     """
#     Retained intro event
#     """
#     def __init__(self, gene, *_):
#         super().__init__(self, gene)
#         self.etype = 'RI'

#         self.construct_events(gene)

#     def add_putative_events(self, exons, transcript):
#         gene = self.gene
#         for i in range(len(exons) - 1):
#             firstexon = exons[i]
#             secondexon = exons[i + 1]
#             eventid = '{}:{}:{}-{}:{}:{}'.format(gene.chr, firstexon[0], firstexon[1], secondexon[0],
#                                                  secondexon[1], gene.strand)

#             self.negative_ids[eventid] = self.negative_ids.get(eventid, []) + [transcript]
#             search_coord = (firstexon[0], secondexon[1])
#             search = self.alt_search.get(search_coord, set())
#             search.add(eventid)
#             self.alt_search[search_coord] = search

#     def add_real_events(self, exons, transcript):
#         for i in range(len(exons)):
#             search_coord = exons[i]

#             if search_coord in self.alt_search:
#                 eventids = self.alt_search[search_coord]
#                 for eventid in eventids:
#                     self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]

#     def export_events_gtf(self, *_):
#         """
#         Generator of GTF lines
#         """
#         strand = self.gene.strand
#         for event in self.positive_ids:
#             full_event = '{}:{}'.format(self.etype, event)
#             e_vals = full_event.replace('-', ':').split(':')

#             line1 = self.gtf_string.format(e_vals[2], e_vals[3], strand, full_event, full_event, 'alternative2')
#             yield line1, self.etype

#             line2 = self.gtf_string.format(e_vals[4], e_vals[5], strand, full_event, full_event, 'alternative2')
#             yield line2, self.etype

#             line3 = self.gtf_string.format(e_vals[2], e_vals[5], strand, full_event, full_event, 'alternative1')
#             yield line3, self.etype

