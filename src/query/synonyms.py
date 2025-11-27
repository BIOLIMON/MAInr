"""
Synonym dictionaries for query expansion in MAInr.
"""

CONDITION_SYNONYMS = {
    # Abiotic Stress
    "drought": ["drought", "water deficit", "water stress", "dehydration", "osmotic stress"],
    "salt": ["salt stress", "salinity", "NaCl", "sodium chloride", "saline", "ionic stress"],
    "cold": ["cold stress", "chilling", "freezing", "low temperature", "frost"],
    "heat": ["heat stress", "high temperature", "thermal stress", "elevated temperature"],
    "uv": ["UV", "UV-B", "UV-A", "ultraviolet", "light stress"],
    "oxidative": ["oxidative stress", "H2O2", "hydrogen peroxide", "ROS", "reactive oxygen"],
    "heavy_metal": ["heavy metal", "cadmium", "lead", "mercury", "aluminum", "metal toxicity"],
    "flooding": ["flooding", "waterlogging", "submergence", "hypoxia", "anoxia"],
    "nutrient_deficiency": ["nutrient deficiency", "nitrogen", "phosphorus", "potassium", "iron", "starvation"],
    
    # Biotic Stress
    "pathogen": ["pathogen", "infection", "disease", "bacterial", "fungal", "viral"],
    "herbivory": ["herbivory", "insect", "pest", "aphid", "caterpillar"],
    
    # Hormones
    "auxin": ["auxin", "IAA", "indole acetic acid"],
    "cytokinin": ["cytokinin", "CK", "6-BA", "kinetin"],
    "abscisic_acid": ["abscisic acid", "ABA"],
    "ethylene": ["ethylene", "ACC"],
    "jasmonic_acid": ["jasmonic acid", "JA", "jasmonate", "MeJA"],
    "salicylic_acid": ["salicylic acid", "SA", "salicylate"],
}

EXPERIMENT_SYNONYMS = {
    "rna_seq": ["RNA-seq", "RNAseq", "transcriptome", "transcriptomic", "RNA sequencing"],
    "single_cell": ["single-cell", "scRNA-seq", "single cell RNA-seq", "10x Genomics"],
    "microarray": ["microarray", "gene chip", "expression array", "Affymetrix"],
    "small_rna": ["small RNA", "miRNA", "sRNA", "microRNA"],
    "wgs": ["whole genome sequencing", "WGS", "genome sequencing"],
    "chip_seq": ["ChIP-seq", "ChIPseq", "chromatin immunoprecipitation"],
    "atac_seq": ["ATAC-seq", "ATACseq", "chromatin accessibility"],
    "bisulfite": ["bisulfite", "WGBS", "methylation", "BS-seq", "methylome"],
}

ORGANISM_SYNONYMS = {
    "tomato": ["Solanum lycopersicum", "tomato", "Lycopersicon esculentum"],
    "potato": ["Solanum tuberosum", "potato"],
    "arabidopsis": ["Arabidopsis thaliana", "Arabidopsis"],
    "rice": ["Oryza sativa", "rice"],
    "maize": ["Zea mays", "maize", "corn"],
    "wheat": ["Triticum aestivum", "wheat"],
    "soybean": ["Glycine max", "soybean"],
}
