# BIOL 7180-EA1: Scripting for Biologist

# Group Members
| Name                | Username          |Email|
| ------------------- | ----------------- |--------------|
| Melika Shiran       | Melikashiran      |mzg0146@auburn.edu|
| Sean Onileowo      | Seun-Onileowo     |sao0027@auburn.edu|
| Surma Mohiudden Meem        | surma-meem        |sum0005@auburn.edu|
| Prince Mensah Ansah | princemensahansah |pma0020@auburn.edu|

# Project Title
Germline PRMD9 expression across selected taxa

## Project Description

Proper segregation of homologous chromosomes requires homologous recombination, which occurs during meiosis. The recombination process occurs at specific binding sites determined by PRMD9 (**DNA-binding histone methyltransferase**). PRDM9 has four functional domains conserved evolutionarily, serving as regulators of transcription. Importantly, organisms with functional gene copies of PRDM9 are known to have recombination hotspots and specific binding sequences. Whereas organisms without it initiate recombination at the promoter region. Evolutionary studies suggest that some non-avian reptile species might have lost this protein through evolution or evolved to regulate non-canonical recombination patterns. For instance, Baker et al. 2017 found orthologs of PRDM9 across distantly related vertebrates. Most taxa had partially conserved domains of the proteins. Moreover, taxa with complete conservation of all domains were those with rapid evolution of the protein’s binding affinity. The protein,PRDM9,is well-studied in mammals and other primates, and its role has been properly documented. However, studies on the presence and functions of PRDM9 (**DNA-binding histone methyltransferase**) in non-avian reptiles are still lacking. Given what is known about the protein and its functions in mammals and other primates, and its role in meiosis, it is likely to be differentially expressed across various tissues. Homologous recombination occurs during meiosis, thus, PRDM9 must be differentially upregulated in germline tissues. PRDM9 will be expressed in the germline in taxa with functional copies of the gene, whereas taxa that have lost the gene through evolution will not have it expressed. By analyzing the transcriptome of germline tissues (gonads, etc.) of non-avian reptiles, PRDM9’s expression can be ascertained.

## Method
# Methods
To test the hypothesis, we will use gonad and a control tissue (liver/heart/skin) transcriptome data from the following species: Pogona vitticeps (Georges et al. 2015), Mauremys mutica (Yuan et al. 2021), Trachemys scripta elegans (Hatkevich et al. 2025), Crotalus tigris, Hemicordylus capensis, Sceloporus undulatus, Pantherophis guttatus, Erythrolamprus reginae, Rhineura floridana, Podarcis muralis, Naja naja, Anniella stebbinsi, and Pelodiscus sinensis (Zhu et al. 2022). This RNAseq data is publicly available through NCBI. For each dataset, the paired-end .fastq R1 and R2 files will be downloaded from NCBI and brought onto the Alabama Supercomputer for analysis. After quality assessment and trimming, the reads will be aligned to the reference genome (.fasta) and annotation (.gff) files. Along with mapping, a gene count matrix file will be produced (.csv), to assess gene expression. 
The project will employ `Bash scripting` for batch processings and `Python` for data visualization


## Data availability





