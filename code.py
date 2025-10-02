# Install Biopython  final
!pip install biopython

from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW, NCBIXML

# Setting up email for NCBI Entrez
Entrez.email = "validemail@mail.com"

# Fetching nucleotide sequence by accession number
accession = "NM_000321"
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

print("Accession:", record.id)
print("Organism:", record.annotations['organism'])


dna_seq = record.seq
print("\n DNA sequence:\n", dna_seq)

# Transcribe mRNA sequence
mrna_seq = dna_seq.transcribe()
print("\n mRNA sequence:\n", mrna_seq)

# Complementary strand length
complement_seq = dna_seq.complement()
print("\nComplementary strand:\n", complement_seq)

# Translate the sequence
seq_len = len(dna_seq)
trimmed_len = seq_len - (seq_len % 3)
dna_seq_trimmed = dna_seq[:trimmed_len]
protein_seq = dna_seq_trimmed.translate()
print("\n Translated protein sequence:\n", protein_seq)

# GC Content calculation
gc_percent = round(gc_fraction(dna_seq) * 100, 2)
print("\nGC Content (%):", gc_percent)

# BLASTn
print("\nRunning BLASTn on sequence ")
result_handle = NCBIWWW.qblast("blastn", "nt", str(dna_seq))

with open("blastn_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
print("BLASTn results saved to blastn_result.xml")

#  Displaying top 3 alignments
def print_top_alignments(blast_xml_file, top_n=3):
    with open(blast_xml_file) as result_handle:
        blast_record = NCBIXML.read(result_handle)
        print(f"\nTop {top_n} Alignments for BLASTn:\n")

        count = 0
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print("Query:", hsp.query)
                print("Match:", hsp.match)
                print("Subject:", hsp.sbjct)
                print("Sequence:", alignment.title)
                print("Length:", alignment.length)
                print(f"Score: {hsp.score}, E-value: {hsp.expect}")
                print()
                count += 1
                if count == top_n:
                    return

print_top_alignments("blastn_result.xml", top_n=3)