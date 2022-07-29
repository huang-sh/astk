"""
astk.suppa.AS_event
~~~~~~~~~~~~~~~~~~~~
Modified version of suppa.lib.event module.
"""
from itertools import product
from typing import Sequence, Optional

from .gtf_parse import Exon, Gene, Transcript, Genome
from .lib.tools import EWriter


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
        self.status = "infer"
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


class SEEvent(ASEvent):

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
            self.set_basic_info()
            return self

    def set_basic_info(self):
        e1, e2, e3 = self.splicein_exon
        self.strand =  e1.strand
        self.id = f"{e1.chr}:{e1.end}-{e2.start}:{e2.end}-{e3.start}:{e1.strand}"


class RIEvent(ASEvent):

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
       
    def validate(self):
        e1, e2 = self.spliceout_exon
        e = self.splicein_exon[0] 
        fbool = [e.start == e1.start, e.end == e2.end, e1.end < e2.start]
        if all(fbool):
            self.status = "review"
            self.set_basic_info()
            return self
    
    def set_basic_info(self):
        e1, e2 = self.spliceout_exon
        self.strand =  e1.strand
        self.id = f"{e1.chr}:{e1.start}:{e1.end}-{e2.start}:{e2.end}:{e2.strand}"


class MXEvent(ASEvent):

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
        self.etype = 'MX'

    def validate(self):
        si1, si2, si4 = self.splicein_exon
        so1, so3, so4 = self.spliceout_exon

        # use <= for get same result to suppa2, 
        # but it should to consider to classify a new AS type
        fbool = [
                si1.end == so1.end, 
                si4.start == so4.start,
                si2.start < si2.end < so3.start <= so3.end  
        ]
        if all(fbool):        
            self.status = "review"
            self.set_basic_info()
            return self

    def set_basic_info(self):
        e1, e2, e4 = self.splicein_exon
        e1, e3, e4 = self.spliceout_exon
        self.strand =  e1.strand
        self.id = (f"{e1.chr}:{e1.end}-{e2.start}:{e2.end}-"
                   f"{e4.start}:{e1.end}-{e3.start}:{e3.end}-{e4.start}:{self.strand}")

    def set_tid(self, exons: Sequence["Exon"]):
        """just for AS event inferring
        """
        e1, e2, e3 = exons
        # coor_str = ":".join([f"{e.start}-{e.end}" for e in exons])
        coor_str = f"S-{e1.end}:{e2.start}-{e2.end}:{e3.start}-E"
        self.tid = f"{self.etype}:{coor_str}"
    
    def set_boundary():
        pass


class A5Event(ASEvent):

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
        self.etype = 'A5'

    def validate(self):
        """It is not exact when exon is not adjacent.
        """
        si1, si2 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        if so1.strand == "+":
            fbool = [
                si2.start == so2.start, 
                si1.start < so1.end < si1.end < so2.start,

            ]
        elif so1.strand == "-":
            fbool = [
                si1.end == so1.end, 
                si1.end < si2.start < so2.start < si2.end
            ]
        else:
            fbool = []
        if all(fbool):        
            self.status = "review"
            self.set_basic_info()
            return self

    def set_basic_info(self):
        si1, si2 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        self.strand = si1.strand
        self.id = f"{si1.chr}:{si1.end}-{si2.start}:{so1.end}-{so2.start}:{self.strand}"


class A3Event(ASEvent):

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
        self.etype = 'A3'

    def validate(self):
        """It is not exact when exon is not adjacent.
        """
        si1, si2 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        if so1.strand == "-":
            fbool = [
                si2.start == so2.start, 
                si1.start < so1.end < si1.end < si2.start
            ]
        elif so1.strand == "+":
            fbool = [
                si1.end == so1.end, 
                si1.end < si2.start < so2.start < si2.end
            ]
        else:
            fbool = []
        if all(fbool):        
            self.status = "review"
            self.set_basic_info()
            return self

    def set_basic_info(self):
        si1, si2 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        self.strand = si1.strand
        self.id = f"{si1.chr}:{si1.end}-{si2.start}:{so1.end}-{so2.start}:{self.strand}"


class AFEvent(ASEvent):

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
        self.etype = 'AF'

    def validate(self):
        si1, si2 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        if so1.strand == "+":
            fbool = [
                si2.start == so2.start, 
                si1.end < so1.start < so1.end < si2.start
            ]
        elif so1.strand == "-":
            fbool = [
                si1.end == so1.end, 
                si1.end < so2.start < so2.end < si2.start
            ]
        else:
            fbool = []
        if all(fbool):        
            self.status = "review"
            self.set_basic_info()
            return self

    def set_basic_info(self):
        si1, si2 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        self.strand = si1.strand
        if self.strand == "+":
            self.id = (f"{si1.chr}:{si1.start}:{si1.end}-{si2.start}:"
                        f"{so1.start}:{so1.end}-{so2.start}:{self.strand}")
        elif self.strand == "-":
            self.id = (f"{si1.chr}:{so1.end}-{so2.start}:{so2.end}:"
                        f"{si1.end}-{si2.start}:{si2.end}:{self.strand}")
 
 
class ALEvent(ASEvent):

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
        self.etype = 'AL'

    def validate(self):
        si1, si2 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        if si1.strand == "-":
            fbool = [
                si2.start == so2.start, 
                si1.end < so1.start < so1.end < si2.start
            ]
        elif si1.strand == "+":
            fbool = [
                si1.end == so1.end, 
                si1.end < so2.start < so2.end < si2.start
            ]
        else:
            fbool = []
        if all(fbool):        
            self.status = "review"
            self.set_basic_info()
            return self

    def set_basic_info(self):
        si1, si2 = self.splicein_exon
        so1, so2 = self.spliceout_exon
        self.strand = si1.strand
        if self.strand == "-":
            self.id = (f"{si1.chr}:{si1.start}:{si1.end}-{si2.start}:"
                        f"{so1.start}:{so1.end}-{so2.start}:{self.strand}")
        elif self.strand == "+":
            self.id = (f"{si1.chr}:{so1.end}-{so2.start}:{so2.end}:"
                        f"{si1.end}-{si2.start}:{si2.end}:{self.strand}")


class AlternativeSplicing:

    def __init__(self, gene: "Gene"):
        self.gene = gene
        self.gene_id = gene.id
        self.chr = gene.chr
        self.AS_events = {}
        self.FT_events = {}
        self.LT_events = {}

    def filter_promoter_event(self):
        """filter AS events that exon overlap with promoter
        """
        self.promoter_event = {}
        events_id = list(self.inner_AS_events.keys())
        for etid in events_id:
            event = self.inner_AS_events[etid]
            for tx in event.splicein_tx:
                if event.strand == "+":
                    if int(tx.TSS) + 2000 > int(event.splicein_exon[0].start):
                        self.promoter_event.setdefault(etid, self.inner_AS_events.pop(etid))
                        break
                else:
                    if int(event.splicein_exon[-1].end) + 2000 > int(tx.TSS):
                        self.promoter_event.setdefault(etid, self.inner_AS_events.pop(etid))
                        break

    def filter_FT(self):
        """filter AS events that exon overlap with promoter
        """
        events_tid = list(self.inner_AS_events.keys())
        for etid in events_tid:
            event = self.inner_AS_events[etid]
            for tx in event.splicein_tx:
                if event.strand == "+":
                    if event.splicein_exon[0].start == tx.cds[0]:
                        if etid not in self.FT_events:
                            self.FT_events.setdefault(etid, self.inner_AS_events.pop(etid))
                else:
                    if event.splicein_exon[-1].end == tx.cds[1]:
                        if etid not in self.FT_events:
                            self.FT_events.setdefault(etid, self.inner_AS_events.pop(etid))
                
    def filter_LT(self):
        """filter AS events that exon overlap with promoter
        """
        events_id = list(self.inner_AS_events.keys())
        for etid in events_id:
            event = self.inner_AS_events[etid]
            for tx in event.splicein_tx:
                if event.strand == "+":
                    if event.splicein_exon[-1].end == tx.cds[1]:
                        if etid not in self.LT_events:
                            self.LT_events.setdefault(etid, self.inner_AS_events.pop(etid))
                else:
                    if event.splicein_exon[0].start == tx.cds[0]:
                        if etid not in self.LT_events:
                            self.LT_events.setdefault(etid, self.inner_AS_events.pop(etid))

    def to_ioe_list(self, event: "ASEvent"):

        splicein_tx = [i.id for i in event.splicein_tx]
        spliceout_tx = [i.id for i in event.spliceout_tx]
        all_tx = splicein_tx + spliceout_tx
        full_event_id = f"{self.gene_id};{self.etype}:{event.id}"
        line = [self.chr, self.gene_id, full_event_id, ",".join(splicein_tx), ",".join(all_tx)]
        return line
 
    def to_ioe(self, type: Optional[str] = None):
        if type is None:
            events = self.AS_events.copy()
            events.update(self.FT_events)
            events.update(self.LT_events)
        elif type == "FTE":
            events = self.FT_events
        elif type == "LTE":
            events = self.LT_events
        elif type == "inner":
            events = self.inner_AS_events
        else:
            events = {}
        lines = [self.to_ioe_list(e) for e in events.values()]
        return lines


class SkippingExon(AlternativeSplicing):
    """
    Skipped exon event
    """
    def __init__(self, gene):
        super().__init__(gene)
        self.etype = 'SE'  
        self.construct_events()

    def _infer_splicein_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = [tx for tx in self.gene.transcripts.values() if len(tx.exons) > 2]
        for tx in transcripts:
            exons = sorted(tx.exons.values(), key=lambda e: e.start)
            for i in range((len(exons)-2)):
                tid = f"{exons[i].end}-{exons[i+2].start}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[i:i+3], "tx": tx})
        return pe_dic

    def _infer_spliceout_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = [tx for tx in self.gene.transcripts.values() if len(tx.exons) > 1]
        for tx in transcripts:
            exons = sorted(tx.exons.values(), key=lambda e: e.start)
            for i in range((len(exons)-1)):
                tid = f"{exons[i].end}-{exons[i+1].start}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[i:i+2], "tx": tx})
        return pe_dic

    def construct_events(self):
        psi_dic = self._infer_splicein_events()
        pso_dic = self._infer_spliceout_events()
        shared_keys = set(psi_dic.keys()) & set(pso_dic.keys())
        for tid in shared_keys:
            for se, so in product(psi_dic[tid], pso_dic[tid]):
                if se["tx"] == so["tx"]:
                    continue
                event = SEEvent(splicein_exon=se["exons"],
                                spliceout_exon=so["exons"])
                if event.status == "review":
                    self.AS_events.setdefault(event.id, event)
                    self.AS_events[event.id].splicein_tx.add(se["tx"])
                    self.AS_events[event.id].spliceout_tx.add(so["tx"])

        self.inner_AS_events = self.AS_events.copy()


class RetainedIntron(AlternativeSplicing):
    """
    Retained Intron event
    """
    def __init__(self, gene: "Gene"):
        super().__init__(gene)
        self.etype = 'RI'  
        self.construct_events()

    def _infer_splicein_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = self.gene.transcripts.values()
        for tx in transcripts:
            exons = sorted(tx.exons.values(), key=lambda e: e.start)
            for i in range(len(exons)):
                tid = f"{exons[i].start}-{exons[i].end}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[i:i+1], "tx": tx})
        return pe_dic

    def _infer_spliceout_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = [tx for tx in self.gene.transcripts.values() if len(tx.exons) > 1]
        for tx in transcripts:
            exons = sorted(tx.exons.values(), key=lambda e: e.start)
            for i in range((len(exons)-1)):
                tid = f"{exons[i].start}-{exons[i+1].end}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[i:i+2], "tx": tx})
        return pe_dic

    def construct_events(self):
        psi_dic = self._infer_splicein_events()
        pso_dic = self._infer_spliceout_events()
        shared_keys = set(psi_dic.keys()) & set(pso_dic.keys())
        for tid in shared_keys:
            for se, so in product(psi_dic[tid], pso_dic[tid]):
                if se["tx"] == so["tx"]:
                    continue
                event = RIEvent(splicein_exon=se["exons"],
                                spliceout_exon=so["exons"])
                if event.status == "review":
                    self.AS_events.setdefault(event.id, event)
                    self.AS_events[event.id].splicein_tx.add(se["tx"])
                    self.AS_events[event.id].spliceout_tx.add(so["tx"])

        self.inner_AS_events = self.AS_events.copy()


class MutuallyExclusiveExon(AlternativeSplicing):
    """
    Mutually Exclusive Exon
    """
    def __init__(self, gene: "Gene"):
        super().__init__(gene)
        self.etype = 'MX'  
        self.construct_events()

    def _infer_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = [tx for tx in self.gene.transcripts.values() if len(tx.exons) > 2]
        for tx in transcripts:
            exons = sorted(tx.exons.values(), key=lambda e: e.start)
            for i in range((len(exons)-2)):
                tid = f"{exons[i].end}-{exons[i+2].start}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[i:i+3], "tx": tx})
        return pe_dic

    def construct_events(self):
        events = self._infer_events()
        for item in events.values():
            for si, so in product(item, item):
                if si["tx"] == so["tx"]:
                    continue
                event = MXEvent(splicein_exon=si["exons"],
                                spliceout_exon=so["exons"])
                if event.status == "review":
                    self.AS_events.setdefault(event.id, event)
                    self.AS_events[event.id].splicein_tx.add(si["tx"])
                    self.AS_events[event.id].spliceout_tx.add(so["tx"])
        self.inner_AS_events = self.AS_events.copy()


class Alternative5SS(AlternativeSplicing):
    """
    Alternative 5' Splice-site
    """
    def __init__(self, gene: "Gene"):
        super().__init__(gene)
        self.etype = 'A5'  
        self.construct_events()

    def _infer_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = [tx for tx in self.gene.transcripts.values() if len(tx.exons) > 1]
        if strand == "+":
            for tx in transcripts:
                exons = sorted(tx.exons.values(), key=lambda e: e.start)
                for i in range((len(exons)-1)):
                    tid = f"{exons[i+1].start}:{strand}"
                    pe_dic.setdefault(tid, [])
                    pe_dic[tid].append({"exons": exons[i:i+2], "tx": tx})
        else:
            for tx in transcripts:
                exons = sorted(tx.exons.values(), key=lambda e: e.start)
                for i in range((len(exons)-1)):
                    tid = f"{exons[i].end}:{strand}" 
                    pe_dic.setdefault(tid, [])
                    pe_dic[tid].append({"exons": exons[i:i+2], "tx": tx})
        return pe_dic

    def construct_events(self):
        events = self._infer_events()
        for item in events.values():
            for si, so in product(item, item):
                if si["tx"] == so["tx"]:
                    continue
                event = A5Event(splicein_exon=si["exons"],
                                spliceout_exon=so["exons"])
                if event.status == "review":
                    self.AS_events.setdefault(event.id, event)
                    self.AS_events[event.id].splicein_tx.add(si["tx"])
                    self.AS_events[event.id].spliceout_tx.add(so["tx"])
        self.inner_AS_events = self.AS_events.copy()


class Alternative3SS(AlternativeSplicing):
    """
    Alternative 3' Splice-site
    """
    def __init__(self, gene: "Gene"):
        super().__init__(gene)
        self.etype = 'A3'
        self.construct_events()

    def _infer_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = [tx for tx in self.gene.transcripts.values() if len(tx.exons) > 1]
        if strand == "-":
            for tx in transcripts:
                exons = sorted(tx.exons.values(), key=lambda e: e.start)
                for i in range((len(exons)-1)):
                    tid = f"{exons[i+1].start}:{strand}"
                    pe_dic.setdefault(tid, [])
                    pe_dic[tid].append({"exons": exons[i:i+2], "tx": tx})
        else:
            for tx in transcripts:
                exons = sorted(tx.exons.values(), key=lambda e: e.start)
                for i in range((len(exons)-1)):
                    tid = f"{exons[i].end}:{strand}" 
                    pe_dic.setdefault(tid, [])
                    pe_dic[tid].append({"exons": exons[i:i+2], "tx": tx})
        return pe_dic

    def construct_events(self):
        events = self._infer_events()
        for item in events.values():
            for si, so in product(item, item):
                if si["tx"] == so["tx"]:
                    continue
                event = A3Event(splicein_exon=si["exons"],
                                spliceout_exon=so["exons"])
                if event.status == "review":
                    self.AS_events.setdefault(event.id, event)
                    self.AS_events[event.id].splicein_tx.add(si["tx"])
                    self.AS_events[event.id].spliceout_tx.add(so["tx"])
        self.inner_AS_events = self.AS_events.copy()


class AlternativeFirstExon(AlternativeSplicing):
    """
    Alternative First Exon
    """
    def __init__(self, gene: "Gene"):
        super().__init__(gene)
        self.etype = 'AF'  
        self.construct_events()

    def _infer_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = [tx for tx in self.gene.transcripts.values() if len(tx.exons) > 1]
        if strand == "+":
            for tx in transcripts:
                exons = sorted(tx.exons.values(), key=lambda e: e.start)
                tid = f"{exons[1].start}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[:2], "tx": tx})
        else:
            for tx in transcripts:
                exons = sorted(tx.exons.values(), key=lambda e: e.start)
                tid = f"{exons[-2].end}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[-2:], "tx": tx})
        return pe_dic

    def construct_events(self):
        events = self._infer_events()
        for item in events.values():
            for si, so in product(item, item):
                if si["tx"] == so["tx"]:
                    continue
                event = AFEvent(splicein_exon=si["exons"],
                                spliceout_exon=so["exons"])
                if event.status == "review":
                    self.AS_events.setdefault(event.id, event)
                    self.AS_events[event.id].splicein_tx.add(si["tx"])
                    self.AS_events[event.id].spliceout_tx.add(so["tx"])
        self.inner_AS_events = self.AS_events.copy()


class AlternativeLastExon(AlternativeSplicing):
    """
    Alternative Last Exon
    """
    def __init__(self, gene: "Gene"):
        super().__init__(gene)
        self.etype = 'AL'
        self.construct_events()

    def _infer_events(self):
        pe_dic = {}
        strand = self.gene.strand
        transcripts = [tx for tx in self.gene.transcripts.values() if len(tx.exons) > 1]
        if strand == "-":
            for tx in transcripts:
                exons = sorted(tx.exons.values(), key=lambda e: e.start)
                tid = f"{exons[1].start}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[:2], "tx": tx})
        else:
            for tx in transcripts:
                exons = sorted(tx.exons.values(), key=lambda e: e.start)
                tid = f"{exons[-2].end}:{strand}"
                pe_dic.setdefault(tid, [])
                pe_dic[tid].append({"exons": exons[-2:], "tx": tx})
        return pe_dic

    def construct_events(self):
        events = self._infer_events()
        for item in events.values():
            for si, so in product(item, item):
                if si["tx"] == so["tx"]:
                    continue
                event = ALEvent(splicein_exon=si["exons"],
                                spliceout_exon=so["exons"])
                if event.status == "review":
                    self.AS_events.setdefault(event.id, event)
                    self.AS_events[event.id].splicein_tx.add(si["tx"])
                    self.AS_events[event.id].spliceout_tx.add(so["tx"])
        self.inner_AS_events = self.AS_events.copy()


def make_events(output: str,
                genome: "Genome",            
                events: Sequence[str],                
                AS_split: Sequence[str], 
                b_type: str = "S", 
                th: int = 10
    ):
    """
    Controls event creation and writing for each gene
    """
    boundary = 'variable_{}'.format(th) if b_type == 'V' else 'strict'

    event_cls_dic = {
                    "MX": MutuallyExclusiveExon,
                    "SE": SkippingExon, "RI": RetainedIntron, 
                    "A5": Alternative5SS, "A3": Alternative3SS,
                    "AF": AlternativeFirstExon, "AL": AlternativeLastExon
                }
    event_cls_dic = {k: v for k, v in event_cls_dic.items() if k in events}

    handle_ls = []
    if "FTE" in AS_split:
        sc_ioe_writer = EWriter(event_cls_dic.keys(), f"{output}_FT", "ioe", boundary)
        handle_ls.append(sc_ioe_writer)
    if "LTE" in AS_split:
        se_ioe_writer = EWriter(event_cls_dic.keys(), f"{output}_LT", "ioe", boundary)
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
            gene_event.filter_FT()
            gene_event.filter_LT()

            if "FTE" in AS_split:
                for line in gene_event.to_ioe(type="FTE"):
                    sc_ioe_writer.write("\t".join(line), t)
            if "LTE" in AS_split:
                for line in gene_event.to_ioe(type="LTE"):
                    se_ioe_writer.write("\t".join(line), t)
            if "inner" in AS_split:
                for line in gene_event.to_ioe(type="inner"):
                    mid_ioe_writer.write("\t".join(line), t)
            if not AS_split:
                for line in gene_event.to_ioe():
                    ioe_writer.write("\t".join(line), t)
    
    for i in handle_ls:
        i.close()