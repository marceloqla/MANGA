# MANGA
### Mutation ANalysis Guided by Annotations - Version: Ubá 1.0 BETA

-----
## Requirements:
1. Python 3.6 or higher
2. [MAFFT](https://mafft.cbrc.jp/alignment/software/) and/or [PRANK](https://www.ebi.ac.uk/research/goldman/software/prank)
3. [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/)
4. [CD-HIT](http://cd-hit.org/)
5. [MAESTRO](https://pbwww.services.came.sbg.ac.at/?page_id=416)
6. Libraries:
   * numpy==1.19.4
   * networkx==2.5
   * xmltramp2==3.1.1
   * ProDy==1.11
   * requests==2.18.4
   * scipy==1.5.4
   * shap==0.40.0
   * matplotlib==3.1.1
   * beautifulsoup4==4.10.0
   * Bio==1.3.3
   * scikit_learn==1.0.1

## Installation:

1. Make sure Python 3.6 or higher, and programs from Requirements (1-4) are installed.
2. Download or Clone this repository.
4. Run *pip install requirements.txt*.
5. Edit scripts:
   * In script *2-get_homologs_align.py*

     Edit lines 326, 327 and 328, inserting MAFFT, Prank and CD-HIT running paths or aliases

   > MAFFT_PATH = "EDIT ME"
   > 
   > PRANK_PATH = "EDIT ME"
   > 
   > CD_HIT_PATH = "EDIT ME"
   
   * In script *4-get_annotations.py*

     Edit line 1515 inserting DSSP running path or alias
   > DSSP_COMMAND = 'EDIT ME'
   
   * In script * 5-gen_predict_reports.py*

     Edit lines 278 and 279, inserting MAESTRO and MAESTRO's config.xml running paths or aliases
    > MAESTRO_COMMAND_LINE = "EDIT ME"
    > MAESTRO_CONFIG_XML = "EDIT ME"

## Running:

1. Run the code below replacing:
   * *UniProtCode* with the UniProt code of your protein of interest
   
   > 
   > #### *python 1-get_homologs_align.py UniProtCode*
   > 
 
   This will generate an *UniProtCode_seqstr.pickle* file, which will be used later.
   
   And the *seqfasta* folder with a *UniProtCode.fasta* file containing the fasta sequence of your protein of interest.
   
   **The _./seqfasta/UniProtCode.fasta_ file will be used in the next script.**

---

2. Run the code below replacing:
   * *UniProtCode* with the UniProt code of your protein of interest
   
   * and *YourEmailHere* with your email (required for running jobs using the EBI API)
   
   > 
   > #### *python 2-get_homologs_align.py YourEmailHere UniProtCode ./seqfasta/UniProtCode.fasta*
   > 

 
   This will generate several folders and files under the *PHMMER_JOBS* created folder.
   
   The file we are interested in is named:
   
   *./PHMMER_JOBS/UniProtCode/msa/seqs_aln_mafftauto_cln_maxid_refseq.fasta* 
   
   You can edit the *2-get_homologs_align.py* for yielding MAFFT-LINSI or PRANK alignments.
   
   If so, the file of interest could also be named:
   
   *./PHMMER_JOBS/UniProtCode/msa/seqs_aln_mafftlinsi_cln_maxid_refseq.fasta* 
   
   or:
   
   *./PHMMER_JOBS/UniProtCode/msa/seqs_aln_prank_cln_maxid_refseq.fasta* 

---


3. Run the code below replacing:
   * *UniProtCode* with the UniProt code of your protein of interest
   
   * and *MSA_file* with the MSA file name described in the previous section)
   
   > 
   > #### *python 3-get_cons_coevo.py MSA_file UniProtCode**
   > 
   
   This will generate an *UniProtCode_conscoevo.pickle* file, which will be used later.
   
   Several folders and files will also be created under the *CONSCOEVO_JOBS* folder.

---

4. Run the code below replacing:
   * *UniProtCode* with the UniProt code of your protein of interest
   
   * *PDB_code* with the selected PDB entry which your protein of interest is mapped to.
   
   * *Chain* with the selected PDB entry chain of your protein of interest (see above)
   
   All *PDB_code* and *Chain* listed in the UniProt Structure mappings are allowed (see: https://www.uniprot.org/uniprot/YourUniProtCode#structure)
   
   * *UniProtCode_seqstr.pickle* with the file generated in the first section of this tutorial
   
   * and *UniProtCode_conscoevo.pickle* with the file generated in the previous section of this tutorial
   
   > 
   > #### *python 4-get_annotations.py UniProtCode, PDB_code, Chain, UniProtCode_seqstr.pickle, UniProtCode_conscoevo.pickle*
   > 

   This will generate an *UniProtCode_PDB_code:Chain_annomap.pickle* file, which will be used in the last step.
   
   Several folders and files will also be created under the *MAPPING_JOBS* folder.

---

5. Run the code below replacing:
   * *UniProtCode_seqstr.pickle* with the file generated in the first section of this tutorial
   * *UniProtCode_conscoevo.pickle* with the file generated in the third section of this tutorial
   * and *UniProtCode_PDB_code:Chain_annomap.pickle* with the file generated in the fourth section of this tutorial
   
   > 
   > #### *python 5-gen_predict_reports.py UniProtCode_seqstr.pickle, UniProtCode_conscoevo.pickle, UniProtCode_PDB_code:Chain_annomap.pickle**
   > 
   
   Several folders and files will also be created under the *PREDICT_JOBS* folder.
   
   **The final results are located under the _RESULTS/UniProtCode_PDB_code:Chain_ folder**
   
   By opening the *prediction.html* file, quantitative and qualitative reports of each possible mutation for the protein and structure of interest can be explored interactively.
   


