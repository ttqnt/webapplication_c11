from app import *
from Bio.Blast import NCBIWWW, NCBIXML
#from Bio import Entrez, Medline

# class voor alignmentpresentatiemethoden, om ze overzichtelijker op de webpagina te kunnen presenteren.
class AlignmentPresentation:
    def __init__(self, query, sbjct):
        self.query = query
        self.sbjct = sbjct
        self.match = self.setMatch()
        self.full_seq = self.full_seq()

    def setMatch(self):
        m = ""
        i= 0
        while len(self.query)-1 >= i:
            if self.query[i] == self.sbjct[i]:
                m += "|"
                i += 1
            else:
                m += " "
                i += 1
        return m

    def full_seq(self):
        alignment = ""
        i = 0
        while len(self.query)-76 >= i:
            alignment += self.query[i:i+75] + "\n" + self.match[i:i+75] + "\n" + self.sbjct[i:i+75] + "\n\n"
            i+=75
        alignment += self.query[i:len(self.query)-1] + "\n" + self.match[i:len(self.query)-1] + "\n" + self.sbjct[i:len(self.query)-1]
        return alignment

#blast een inputsequentie met qblast, en returns een met NCBIXML parsed bestand van de alignments.
def blast_seq(input_seq, seq_type):
    if seq_type == "nuc":
        result_handle = NCBIXML.read(NCBIWWW.qblast("blastn", "nt", input_seq))
    else:
        result_handle = NCBIXML.read(NCBIWWW.qblast("blastp", "nr", input_seq))
    return result_handle

#returns nested dictionary van de ingevoerde sequenties op volgorde van invoeren,
#de accessienummers van de alignments met bijbehorende annotatie.
#(handled dus meerdere sequenties uit een .fasta bestand.)
def fastaBlastRecords(blast_record):
    blast_records = {}
    i = 1
    for record in blast_record:
        alignments = {}
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if ">" in alignment.hit_def:
                    title = alignment.hit_def[:alignment.hit_def.index(">")]
                else:
                    title = alignment.hit_def
                alignments[alignment.accession] = {
                    "accessionnr" : alignment.accession,
                    "hit_id" : alignment.hit_id,
                    "hit_def" : title,
                    "length" : alignment.length,
                    "expect" : hsp.expect,
                    "score" : hsp.score,
                    "identities" : hsp.identities,
                    "positives": hsp.positives,
                    "gaps" : hsp.gaps,
                    "query" : hsp.query,
                    "sbjct" : hsp.sbjct
                }
                alignObject = AlignmentPresentation(hsp.query, hsp.sbjct)
                alignments[alignment.accession].update({"match" : alignObject.match, "full_seq" : alignObject.full_seq})
        blast_records[i] = alignments
        i+=1
    return blast_records

#returns dictionary van alignment.accessienummers met annotatiegegevens van de alignment.
def seqBlastRecords(blast_record):
    alignments = {}
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if ">" in alignment.hit_def:
                title = alignment.hit_def[:alignment.hit_def.index(">")]
            else:
                title = alignment.hit_def
            alignments[alignment.accession] = {
                "accessionnr" : alignment.accession,
                "hit_id" : alignment.hit_id,
                "hit_def" : title,
                "length" : alignment.length,
                "expect" : hsp.expect,
                "score" : hsp.score,
                "identities" : hsp.identities,
                "positives": hsp.positives,
                "gaps" : hsp.gaps,
                "query" : hsp.query,
                "sbjct" : hsp.sbjct
                }
            alignObject = AlignmentPresentation(hsp.query, hsp.sbjct)
            alignments[alignment.accession].update({"match" : alignObject.match, "full_seq" : alignObject.full_seq})
    return alignments
