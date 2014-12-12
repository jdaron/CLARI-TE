This readme should help you get started with "clariTE.pl", which is for 
TE annotation

--------------------------------------------------------------------------------
Setting up -------------------------------------------------------------------
--------------------------------------------------------------------------------

i) Using clariTE.pl
====================

"clariTE.pl" is a program that post process the xm output file of RepeatMasker.


ii) Usage
=========

clariTE.pl -fasta <fasta_file> -LTR <position of LTR> -classi <classification> -gene <embl> <xm file>

Command line example :
clariTE.pl -LTR clariTeRep_LTRposition.txt -gene v443_0001.annoGene.embl -classi clariTeRep_classification.txt -fasta v443_0001.fas -dir [your directory] v443_0001.fas.out.xm

iii) REQUIREMENTS
=================

To recap, "clariTE.pl" requires:

  - the output file ".xm" from RepeatMasker
  - a tabular file containing the classification of the library use in RepeatMasker :
  - a tabular file containing the position of the LTR retrotransposons present in the library
  - embl file containing the gene annotation of the sequence (in triAnnot format).
  - a FASTA-format file containing your query sequence (multi-FASTA file is not accepted)

iv) Input FILE
==============

1) xm output file from RepeatMasker

example :
---------
1725 20.7  2.6  0.8 v443_2791 5 393 (40687) C ctgB_rep_0246#Unknown (7481) 396 1 *
4791 21.8  0.8  3.3 v443_2791 378 1704 (39376) + ctgF_rep_0077#Unknown 322 1615 (7557) *
5649 11.4  0.4  0.4 v443_2791 799 1704 (39376) + ctgF_rep_0055#Unknown 722 1627 (11086) 
15719  9.1  0.5  0.7 v443_2791 6858 9147 (31933) + ctgF_rep_0055#Unknown 5829 8115 (4598) 
9263 22.2  0.9  2.8 v443_2791 7600 9932 (31148) + ctgF_rep_0077#Unknown 2241 4528 (4644) *

2) Tabular file of the classification 

example :
---------
#sequenceID	family
TREP100	RLC_famc1.1
TREP1000	DTT_famn14
TREP1001	DTT_famn1

3) Tabular file of the position of LTR in the LTR retrotransposons 
Note : the total number of LTR retrotransposons is not obligated to be present in this file.

example :
---------
#sequenceID	family	LTRtype	start	end
ctgH_rep_0307	RLG_famc3	LTR5	1	3945
ctgH_rep_0307	RLG_famc3	LTR3	7931	12091
sr2_rep_0226	RLC_famc9	LTR5	1	214
sr2_rep_0226	RLC_famc9	LTR3	4972	5185

4) Gene annotation in EMBL format TriAnnot
Note : only the locus_tag and the id is require to run clariTE.pl, other tag (such as blastp_file...) are not necessary.

example :
---------
ID   unknown; SV 1; linear; unassigned DNA; STD; UNC; 1411106 BP.
XX
AC   unknown;
XX
XX
FT   CDS             join(141960..142006,142121..142147,142248..142370,142493..142739,142850..142873)
FT                   /locus_tag="v443_0002_EXONERATE_BLASTX_protOSA_6"
FT                   /blastp_file="v443_0002_EXONERATE_BLASTX_protOSA_141960_142873_Q5JM42_Match_0005_mRNA_CDS.bltp"
FT                   /id="v443_0002_EXONERATE_BLASTX_protOSA_141960_142873_Q5JM42_Match_0005_mRNA_joinedCDS"
FT                   /note="Similar_to: hypothetical_protein"
FT                   /note="BestBlastHit: B9EZI3_ORYSJ TrEMBL databank Putative uncharacterized protein - %25id: 91.67 - hcov: 13.78 - qcov: 100.00"
FT                   /note="Status: High Confidence"
FT   CDS             complement(join(143435..144154,144239..144363,145030..145267))
FT                   /locus_tag="v443_0002_EXONERATE_BLASTX_validated_9"
FT                   /expressed
FT                   /blastp_file="v443_0002_EXONERATE_BLASTX_validated_143435_145267_AFR_02_CAT01_3_Match_0001_mRNA_CDS.bltp"
FT                   /id="v443_0002_EXONERATE_BLASTX_validated_143435_145267_AFR_02_CAT01_3_Match_0001_mRNA_joinedCDS"
FT                   /note="Similar_to: putative_function - F2CSA4_HORVD TrEMBL databank Predicted protein OS Hordeum vulgare var distichum PE 2 SV 1"
FT                   /note="BestBlastHit: F2CSA4_HORVD TrEMBL databank Predicted protein - %25id: 96.12 - hcov: 100.56 - qcov: 100.00"
FT                   /note="Function_coverage: 94.71"
FT                   /note="Function_identity: 97.94"
FT                   /note="Function_target: F2CSA4 22 361"
FT                   /note="Status: High Confidence"


5) FASTA file
Note : only one sequence could be contained it the FASTA file, multi fasta file are accepted

v) Help
=======

#### clariTE.pl ####
#
# CREATED:    2012-05-31
# LAST MODIF: 01/02/2014
# AUTHOR:     Josquin Daron (INRA Clermont-Ferrand)
# VERSION:    1
#
# This script is used to post process a embl output of RepeatMasker
# 

USAGE:
       clariTE.pl -fasta <fasta_file> -LTR <position of LTR> -classi <classification> -gene <embl> <xm file>

       ### OPTIONS ###

       -h|--help:             print this help
       -LTR:                  tabulation file of the LTR annotations
       -classi:               corresponding file sequence id to family id
       -fasta:                fasta file of the annotated sequence
       -gene:                 embl file of the gene annotation
       -dir:                  directory name
       -v:                    verbosity (3,4)
       

vi) Output format
=================

	1) In STDERR

--> read all input File
--> step overlaping feature
--> step merge collinear feature
--> step Join
--> /home/v443_0001.fas_annoTE.embl embl file created
--> failure is not an option

	2) In embl file

Example embl output :
---------------------
ID   unknown; SV 1; linear; unassigned DNA; STD; UNC; 1411106 BP.
XX
AC   unknown;
XX
XX
FH   Key             Location/Qualifiers
FH
FT   repeat_region   complement(join(10..373,1824..2754))
FT                   /parent="mp1" # parent id 
FT                   /status="fragmented" # status of the predicton (complete or fragmented)
FT                   /range="3854..4785,4881..5242," # position on the reference sequence used to annotated the prediction
FT                   /post="RLG_famc36 8671bp 3854..5242" # tag
FT                   /compo="RLG_famc36 100.00 " # composition of the prediction
FT                   /id="3_v443_0001" # id of the prediction 
FT                   /copie="ctgD_rep_0264" # sequenceID in the library used to annotated the prediction

vii) Licence
============

Copyright or © or Copr. Josquin Daron INRA-GDEC 01/09/2014
 
email: josquin.daron@clermont.inra.fr

This software is a computer program whose purpose is to predict Transposable 
Elements (TEs) in complexe genome such as wheat. The program correct raw similarity
search result from RepeatMasker.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 
 
As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.


