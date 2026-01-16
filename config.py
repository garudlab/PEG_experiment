base_dir="/u/project/ngarud/Garud_lab/metagenomic_fastq_files/PEG_experiment"
ann_dir=f"{base_dir}/Reference_genomes"
midas_db_dir=f"{base_dir}/midas_db"
raw_data_dir = f"{base_dir}/raw_data"
base_snps_dir = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/PEG_experiment/midas_output/merged_midas_output/snps"
base_genes_dir = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/PEG_experiment/midas_output/merged_midas_output/genes"
fig_dir = "figures/"

good_species = ["2022July7_B_theta_VPI5482_PROKKA_01062023",
                "2023Jan2_Akkermansia_muciniphila_DSM_22959_PROKKA_01062023",
                "2023Jan2_Bacteroides_ovatus_ATCC_8483_PROKKA_01062023",
                "2023Jan2_Clostridium_sporogenes_ATCC_15579_PROKKA_01062023",
                "2023Jan2_Collinsella_stercoris_DSM_13279_PROKKA_01062023",
                "2023Jan2_Enterococcus_faecalis_TX1322_PROKKA_01062023",
                "2023Jan2_Escherichia_coli_BW25113_strain_K_12_PROKKA_01062023",
                "2023Jan2_Eubacterium_rectale_ATCC_22656_PROKKA_01062023",
                "2023Jan2_Faecalibacterium_prausnitzii_A2_165_PROKKA_01062023",
                "2023Jan2_M_intestinale_G6_hybridAssembly_PROKKA_01062023"]


species_name_dic = {'Akkermansia_muciniphila':'Akkermansia muciniphila',
 'Bacteroides_ovatus':'Bacteroides ovatus',
 'B_theta':'Bacteroides thetaiotaomicron',
 'Clostridium_sporogenes':'Clostridium sporogenes',
 'Collinsella_stercoris':'Collinsella stercoris',
 'Enterococcus_faecalis':'Enterococcus faecalis',
 'Escherichia_coli':'Escherichia coli',
 'Eubacterium_rectale':'Eubacterium rectale',
 'Faecalibacterium_prausnitzii':'Faecalibacterium prausnitzii',
 'M_intestinale':'Muribaculum intestinale'}                
             