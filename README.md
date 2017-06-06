# Overview

This software package predicts a human isoform interactome from a human reference interactome using the domain-based method DIIP: Domain-based Isoform Interactome Prediction. Given a reference interactome consisting of experimentally observed reference interactions between reference proteins, we first annotate reference protein with structural domains by scanning their sequences using HMMER hmmscan with an E-value cutoff of 10e-5. Next, we annotate interactions in the reference interactome with domain-domain interactions (DDIs), where an interaction is annotated with a DDI if one of the interacting proteins is annotated with one of the domains of the DDI and the other interacting protein is annotated with the other domain of the same DDI. We then predict an isoform interactome from the DDI-annotated reference interactome by predicting isoform-specific interactions for the alternative isoforms of the reference proteins in the DDI-annotated interactome. Given a reference interaction between two reference proteins annotated with a DDI, we predict that an alternative isoform of one of the interacting proteins loses interaction with the partner of that protein if it loses the DDI annotation, otherwise the interaction is retained by that alternative isoform. Alternative isoforms are annotated with structural domains also using HMMER hmmscan with an E-value cutoff of 10e-5.

# Code descriptions

### Code: isoform_partners_go_coexpr.m

This code first predicts the isoform interactome from the reference interactome, either HI-II-14 or IntAct interactome. Then it calculates GO association similarity and co-expression for pairs of reference proteins interacting with the same subset of isoforms of the same gene, pairs of reference proteins interacting with different subsets of isoforms of the same gene, and pairs of reference proteins interacting with protein products of different genes. GO similarity is calculated as the fraction (Jaccard similarity index) of GO terms shared by the two proteins. Co-expression is calculated using Pearson's correlation coefficient.

Data files used by this code:
	
	If using HI-II-14 reference interactome:
	- HI-II-14.tsv (HI-II-14 reference interactome)
	- HI-II-14_spEntrezMap.tab (SwissProt ID to Entrez ID mapping)
		* column 1 (name: From): SwissProt ID
		* column 2 (name: To):   Entrez ID
	
	If using IntAct reference interactome:
	- intact.txt: IntAct interactions
	- IntAct_spGeneMap.txt (SwissProt-to-GeneName mapping)
		* column 1 (name: From): SwissProt ID
		* column 2 (name: To):   Gene Name
	
	Other files:
	- 3did_flat.txt (3did domain-domain interactions)
	- domine_interactions.xlsx (DOMINE domain-domain interactions)
	- hmmer_canonical_full.domtab.txt (hmmscan output table for human reference proteins)
	- hmmer_isoforms.domtab.txt (hmmscan output table for human alternative isoforms)
	- human_protein_sequences.tab (human reference protein sequences)
	- human_protein_can_iso_sequences.fasta (human alternative isoform sequences)
	- gene_association.goa_ref_human.xlsx (Gene Ontology associations)
		May split into two files if MATLAB cannot load at once:
		- gene_association.goa_ref_human_a.xlsx
		- gene_association.goa_ref_human_b.xlsx	
	- E-MTAB-513.tsv.txt (gene tissue expression)

### Code: HI-II-14_isoform_partners_disSubNetSim.m

This code first predicts the isoform interactome from the HI-II-14 reference interactome. Then it calculates disease subnetwork similarity for pairs of reference proteins interacting with the same subset of isoforms of the same gene, pairs of reference proteins interacting with different subsets of isoforms of the same gene, and pairs of reference proteins interacting with protein products of different genes. Disease subnetwork similarity is calculated as the fraction (Jaccard similarity index) of disease subnetworks shared by two proteins, where two proteins share a disease subnetwork if each protein or its interaction partner in the HI-II-14 reference interactome is associated with the disease.

Data files used by this code:
	
	- HI-II-14.tsv (HI-II-14 reference interactome)
	- HI-II-14_spEntrezMap.tab (SwissProt ID to Entrez ID mapping)
		* column 1 (name: From): SwissProt ID
		* column 2 (name: To):   Entrez ID
	- 3did_flat.txt (3did domain-domain interactions)
	- domine_interactions.xlsx (DOMINE domain-domain interactions)
	- hmmer_canonical_full.domtab.txt (hmmscan output table for human reference proteins)
	- hmmer_isoforms.domtab.txt (hmmscan output table for human alternative isoforms)
	- human_protein_sequences.tab (human reference protein sequences)
	- human_protein_can_iso_sequences.fasta (human alternative isoform sequences)
	- curated_gene_disease_associations.tsv	(gene disease associations from DisGeNET)
	- mapa_4_uniprot_crossref.tsv (SwissProt-to-GeneName mapping)
		* column 1 (name: UniProtKB):	SwissProt ID
		* column 2 (name: GENE_SYMBOL): Gene Name

### Code: IntAct_isoform_partners_disSubNetSim.m

This code first predicts the isoform interactome from the IntAct reference interactome. Then it calculates disease subnetwork similarity for pairs of reference proteins interacting with the same subset of isoforms of the same gene, pairs of reference proteins interacting with different subsets of isoforms of the same gene, and pairs of reference proteins interacting with protein products of different genes. Disease subnetwork similarity is calculated as the fraction (Jaccard similarity index) of disease subnetworks shared by the two proteins, where two proteins share a disease subnetwork if each protein or its interaction partner in the high-quality HI-II-14 reference interactome is associated with the disease.

Data files used by this code:
	
	- HI-II-14.tsv (HI-II-14 reference interactome)
	- HI-II-14_spEntrezMap.tab (SwissProt ID to Entrez ID mapping)
		* column 1 (name: From): SwissProt ID
		* column 2 (name: To):   Entrez ID
	- intact.txt: IntAct interactions
	- IntAct_spGeneMap.txt (SwissProt-to-GeneName mapping)
		* column 1 (name: From): SwissProt ID
		* column 2 (name: To):   Gene Name
	- 3did_flat.txt (3did domain-domain interactions)
	- domine_interactions.xlsx (DOMINE domain-domain interactions)
	- hmmer_canonical_full.domtab.txt (hmmscan output table for human reference proteins)
	- hmmer_isoforms.domtab.txt (hmmscan output table for human alternative isoforms)
	- human_protein_sequences.tab (human reference protein sequences)
	- human_protein_can_iso_sequences.fasta (human alternative isoform sequences)
	- curated_gene_disease_associations.tsv	(gene disease associations from DisGeNET)
	- mapa_4_uniprot_crossref.tsv (SwissProt-to-GeneName mapping)
		* column 1 (name: UniProtKB):	SwissProt ID
		* column 2 (name: GENE_SYMBOL): Gene Name

### Code: validation.m

This code validates the performance of the domain-based isoform interaction prediction method on the experimental dataset of Yang et al. (2016). Interactions between Orfeome proteins and newly cloned reference proteins are used as the reference interactome and domain-domain interactions are mapped onto these reference interactions. Interactions for the newly cloned alternative isoforms are then predicted from the DDI-annotated reference interactions using two methods:

- First method uses reference interactions having full DDI annotations only. A reference interaction has a full DDI annotation if the newly cloned protein contains an interacting domain of a DDI and its orfeome interaction partner contains the other interacting domain of the same DDI. An alternative isoform of the newly cloned reference protein is then predicted to lose the interaction if it loses all DDI annotations with the orfeome partner, otherwise the interaction is retained.

- Second method uses reference interactions that have either a full DDI annotation or a partial DDI annotation. A reference interaction has a partial DDI annotation if the newly cloned protein contains an interacting domain of a DDI whereas its orfeome interaction partner does not contain the other interacting domain of the DDI. An alternative isoform of the newly cloned reference protein is then predicted to lose the interaction if it loses all of its interacting domains, otherwise the interaction is retained.

# Data file links

- HI-II-14.tsv
	- http://interactome.dfci.harvard.edu/H_sapiens/download/HI-II-14.tsv

- intact.txt
	- ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip

- 3did_flat.txt
	- http://3did.irbbarcelona.org/download/current/3did_flat.gz

- domine_interactions.xlsx
	- http://domine.utdallas.edu/Domine-2.0/domine-tables-2.0.zip
	- rename interaction.xlsx as domine_interactions.xlsx

- gene_association.goa_ref_human.xlsx
	- ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN///goa_human.gaf.gz
	- remove descriptions in beginning and end of file
	- save as excel file: gene_association.goa_ref_human.xlsx

- E-MTAB-513.tsv.txt
	- https://www.ebi.ac.uk/gxa/experiments/E-MTAB-513/

- curated_gene_disease_associations.tsv
	- http://www.disgenet.org/web/DisGeNET/menu/downloads#gdascurated
