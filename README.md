# PacBio_pbIsoCollapse

This workflow is designed to further correct and reduce redundant isoforms generated from [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake). Although PacBio reads have 3´-polyA signals, they may not be true signals due to internal poly-A primings or other errors. Besides, since the 5’ends of these isoforms are not accurate, they may not be collapsed by cDNA_Cupcake.

## 1. Dependencies:

This tool requires the installation of **bedtools**. You can install it via conda:

`conda install -c bioconda bedtools`

## 2. Installation:

You can clone this GitHub repository, and then add its path to your $PATH variable.

`git clone https://github.com/lhuang3s/PacBio_pbIsoCollapse.git`

`export PATH=<path_to_PacBio_pbIsoCollapse>/bin:$PATH`

## 3. Usages:

We will use the input files in the **test_data** folder as an example. To run the workflow, you need the following inputs:

- **Input bed12 file**. Generated by **cDNA_Cupcake**. Isoform names must be in **PB.id.n** or **Name.n** format.
- **A list of unique gene/cluster names**, the PB.id or Name must match the name prefix in the bed12 file.
- **Isoform abundance file** with at least two space-/tab- delimited columns: 
  PB.id.n/Name.n and an integer count, no header. It can also contain multiple count columns, representing each data
  If not provided, the workflow will generate this file and treat all counts as 1.
- **Prefix of the output files**

The command will be like:

`Usage: pbIsoCollapse.sh [Annotation.bed12] [List|PB.id] [Abundance] [OutputPrefix]`

`pbIsoCollapse.sh Test.realPB13.bed12 Test.realPB13.bed12.list Test.realPB13.bed12.abd2.txt testrunPB13`

## 4. Output files:

`OutputPrefix.bed12`  Final collapsed and corrected isoforms, sorted by genomic locations

`OutputPrefix.Original2Final.matching`  Matching relationships between the original isoform names and the collapsed isoform names in OutputPrefix.bed12

`OutputPrefix.abundance.txt`  Provide the merged abundances for each collapsed and corrected isoform
