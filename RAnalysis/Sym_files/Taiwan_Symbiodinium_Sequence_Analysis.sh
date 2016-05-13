# Manuscript Title: Integrating physiological and symbiotic signatures distinguishes innate functionality and environmental responsiveness in reef corals
# Authors: HM Putnam, PJ Edmunds, TY Fan, RW Lee, EA Lenz, J Lemus, LH Wang, DM Yost, and RD Gates
# Last Updated: January 3, 2016
#-----------------------------------------------------------------------------
# This script processes ITS2 sequences generated from 454 pyrosequencing. Sequences provided by Research and Testing Laboratory (Lubbock, TX) are trimmed using cutadapt and processed using SymTyper, a custom bioinformatic pipeline for Symbiodinium ITS2 sequence analysis. "New" sequences identified by SymTyper are subsequently clustered using QIIME and representative sequences are compared to the custom BLAST database. Sequences phylogenetically placed at internal tree nodes by SymTyper are contextualized by identifying their closest descendant leaf taxa and corresponding distances.
# Input files to run this script are 1) Taiwan.fasta, and 2) Taiwan.ids.
# Final output files of this script are 1) breakdown.tsv 2) new_otu_table.tsv 3) new_otu_clades.txt and 4) TableS1.csv
#-----------------------------------------------------------------------------
# Trim primers
# Count number of starting sequences
count_seqs.py -i Taiwan.fasta  # 312974  : Taiwan.fasta (Sequence lengths (mean +/- std): 316.9799 +/- 45.3669)
# Trim forward primers using cutadapt
#   Allow error rate of 15% (~3 indels/mismatches), minimum overlap of 18
#   Remove sequences that do not have forward primer sequence
cutadapt -g GTGAATTGCAGAACTCCGTG -e 0.15 -O 18 --discard-untrimmed Taiwan.fasta -o Taiwan_trimF.fasta
# Trim forward primers again (may be multiple internal primer sequences), but do not discard sequences that do not contain primer
cutadapt -g GTGAATTGCAGAACTCCGTG -e 0.15 -O 18 Taiwan_trimF.fasta -o Taiwan_trimF2.fasta
# Trim reverse primers using cutadapt
#   Allow error rate of 15% (~3 indels/mismatches), minimum overlap of 18
#   Do not remove sequences that do not have reverse primer sequence
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 -O 18 Taiwan_trimF2.fasta -o Taiwan_trimF2_trimR.fasta
# Trim reverse primers again for any remaining internally
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 -O 18 Taiwan_trimF2_trimR.fasta -o Taiwan_trimmed.fasta
# Count sequences remaining after primer trimming
count_seqs.py -i Taiwan_trimmed.fasta  # 312537  : Taiwan_trimmed.fasta (Sequence lengths (mean +/- std): 280.7927 +/- 40.5028)

#-----------------------------------------------------------------------------
# SymTyper bioinformatic pipeline
symTyper.py -t 24 clade -i Taiwan_trimmed.fasta -s Taiwan.ids -e 1e-20 --hmmdb ~/symTyper/database/HMMER_ITS2_DB/All_Clades.hmm
# Changed e-value cutoff from 1e-05 to 1e-20 in blastall command in ProgramRunner.py script before running subtype command
symTyper.py -t 24 subtype -s Taiwan.ids -H hmmer_hits/ -b blast_output/ -r blastResults/ -f fasta --blastdb ~/symTyper/database/blast_DB/ITS2_Database_04_23_13.fas
symTyper.py -t 20 resolveMultipleHits -s Taiwan.ids -m blastResults/MULTIPLE/ -c resolveMultiples/
symTyper.py -t 20 buildPlacementTree -c resolveMultiples/correctedMultiplesHits/corrected -n  ~/symTyper/database/clades_phylogenies -o placementInfo
symTyper.py -t 20 stats --outputs_dir . -i Taiwan_trimmed.fasta --out_file ./outputfile
symTyper.py -t 20 makeTSV --outputs_dir .

#-----------------------------------------------------------------------------
# Process SymTyper "new" sequences de novo
# Extract "new" sequences from SymTyper output
cat blastResults/NEW/*.out > new_names.txt
gawk 'NR==FNR { a[">"$1]; next } ($1 in a) { print; getline; print }' new_names.txt Taiwan_trimmed.fasta > Taiwan_new_seqs.fasta
count_seqs.py -i Taiwan_new_seqs.fasta  # 2858
# Reformat for QIIME
sed 's/::/_/g' Taiwan_new_seqs.fasta > Taiwan_new_seqs_for_QIIME.fasta

# Cluster "new" sequences de novo at 97% similarity using QIIME
pick_otus.py -i Taiwan_new_seqs_for_QIIME.fasta -s 0.97 -o new_otus  
# Create "new" OTU table
make_otu_table.py -i new_otus/Taiwan_new_seqs_for_QIIME_otus.txt -o new_otus/new_otu_table.biom
biom convert -i new_otus/new_otu_table.biom -o new_otus/new_otu_table.tsv --to-tsv
# Pick representative sequences (most abundant) for new OTUs
pick_rep_set.py -i new_otus/Taiwan_new_seqs_for_QIIME_otus.txt -f Taiwan_new_seqs_for_QIIME.fasta -m most_abundant -o new_otus/rep_set.fasta
# Get clade assignments for denovo representative set
awk '/>/{sub(/_/, "::"); print $2}' new_otus/rep_set.fasta > new_otus/rep_set_names.fasta
for i in `seq 1 $(grep -c ">" new_otus/rep_set.fasta)`
do
seqname=$(sed -n "$i"p new_otus/rep_set_names.fasta)
samplename=$(echo $seqname | cut -d ":" -f 1)
clade=$(grep -w "$seqname" hmmer_parsedOutput/$samplename/HIT | cut -f 4 | cut -d "_" -f 1)
otu=$(grep -w `echo $seqname | sed 's/::/_/'` new_otus/rep_set.fasta | cut -d " " -f 1 | sed 's/>//')
echo $otu $clade >> new_otu_clades.txt
done	
# Blast de novo representative set to get closely related type information
assign_taxonomy.py -i new_otus/rep_set.fasta -r /symTyper/database/blast_DB/ITS2_Database_04_23_13.fas -m blast -o new_otus/blast_taxonomy
# Extract list of blast hits and e-values
paste -d'_' <(sort new_otu_clades.txt | cut -d' ' -f2 | cut -c6) <(sort new_otus/blast_taxonomy/rep_set_tax_assignments.txt | cut -f1,3,4) | awk '{print $1","$3","$2}' | sed 's/enovo/:/g' > denovo_blastid.txt




#-----------------------------------------------------------------------------
# Get closest leaf taxa for SymTyper internal nodes
python
import sys
import os
from ete2 import Tree

orig_stdout = sys.stdout
f = file('node_closest_leaves.txt', 'w')
sys.stdout = f

At = Tree('placementInfo/A/placed_clade_A.nwk')
Ct = Tree('placementInfo/C/placed_clade_C.nwk')
Dt = Tree('placementInfo/D/placed_clade_D.nwk')
Ft = Tree('placementInfo/F/placed_clade_F.nwk')
Gt = Tree('placementInfo/G/placed_clade_G.nwk')

for node in At.traverse("postorder"):
	print node.name
	node.get_closest_leaf()

for node in Ct.traverse("postorder"):
	print node.name
	node.get_closest_leaf()

for node in Dt.traverse("postorder"):
	print node.name
	node.get_closest_leaf()

for node in Ft.traverse("postorder"):
	print node.name
	node.get_closest_leaf()

for node in Gt.traverse("postorder"):
	print node.name
	node.get_closest_leaf()

sys.stdout = orig_stdout
f.close()
exit()

# Extract list of closest leaves and their distances from internal nodes
paste -d'_' <(awk '/^I_/{print; getline; print}' node_closest_leaves.txt | grep -oe "'[^']*'" -oe '^I_.*$' -oe '\s[0-9][^)]*' | sed 's/'\''//g;s/ //g' | sed ':a;N;$!ba;s/\n/,/g;s/,I_/\nI_/g' | cut -d, -f2 | cut -c1) <(awk '/^I_/{print; getline; print}' node_closest_leaves.txt | grep -oe "'[^']*'" -oe '^I_.*$' -oe '\s[0-9][^)]*' | sed 's/'\''//g;s/ //g' | sed ':a;N;$!ba;s/\n/,/g;s/,I_/\nI_/g') | sed 's/_I_/_i:/g' > internal_leafid.txt

#-----------------------------------------------------------------------------
# Concatenate internal node and denovo cluster identifications for Supplementary Table S1
cat denovo_blastid.txt internal_leafid.txt | sort > TableS1.csv
