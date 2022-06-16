# miRNApaper
## description

Pipelines and analyses for Shiekhattarlab's paper "The Integrator complex regulates microRNA abundance through RISC loading",
currently available at https://www.biorxiv.org/content/10.1101/2021.09.21.461113v1

- folder: conda_environments: share the used environments (yml files) to import in miniconda
- sh files:
\t subsample_fastq.sh: To account for varying read numbers and simplify data normalization, we merged multiple sets (if available) and randomly subsampled all samples to 30 million reads. 
\t smrna_mapping.sh: describes treatement of fastq-files from adapter removal, removal of repetitive regions, mapping against the human genome, to bigwig file generation
\t mirdeep2_mapping.sh: fastq-files were specifically mappend to human miRNAs. 
- folder: smRNA_Deseq2: analyses of miRNA abundance 
- folder: pri-miRNA_mapping: files and analyses to assess primary transcript abundance
- folder: GOterms_miRNAmachinery: definition of miRNA-related GOterms and analyses of expression changes
- folder: s4U_figures: analyses of s4U metabolic labeling experiment
- folder: s4U_AGO_RIP_figures: analyses of s4U-labeling followed by AGO2 RIP
- folder: ints11_ago2_bound_miRNA_characteristics: quantification of miRNAs bound by INTS11 eCLIP or AGO2 RIP