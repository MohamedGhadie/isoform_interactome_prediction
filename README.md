# Overview

This software package predicts a human isoform interactome from a human reference interactome using the domain-based method DIIP: Domain-based Isoform Interactome Prediction. Given a reference interactome consisting of experimentally observed reference interactions between reference proteins, we first annotate reference protein with structural domains by scanning their sequences using HMMER hmmscan with an E-value cutoff of 10e-5. Next, we annotate interactions in the reference interactome with domain-domain interactions (DDIs), where an interaction is annotated with a DDI if one of the interacting proteins is annotated with one of the domains of the DDI and the other interacting protein is annotated with the other domain of the same DDI. We then predict an isoform interactome from the DDI-annotated reference interactome by predicting isoform-specific interactions for the alternative isoforms of the reference proteins in the DDI-annotated interactome. Given a reference interaction between two reference proteins annotated with a DDI, we predict that an alternative isoform of one of the interacting proteins loses interaction with the partner of that protein if it loses the DDI annotation, otherwise the interaction is retained by that alternative isoform. Alternative isoforms are annotated with structural domains also using HMMER hmmscan with an E-value cutoff of 10e-5.

# File descriptions:

File: isoform_partners_go_coexpr.m

This script first predicts the isoform interactome from the reference interactome, either HI-II-14 or IntAct interactome. Then it calculates GO association similarity and co-expression for pairs of reference proteins interacting with the same subset of isoforms of the same gene, pairs of reference proteins interacting with different subsets of isoforms of the same gene, and pairs of reference proteins interacting with protein products of different genes. GO similarity is calculated as the fraction (Jaccard similarity index) of GO terms shared by the two proteins. Co-expression is calculated using Pearson's correlation coefficient.

File: HI-II-14_isoform_partners_disSubNetSim.m

This script first predicts the isoform interactome from the HI-II-14 reference interactome. Then it calculates disease subnetwork similarity for pairs of reference proteins interacting with the same subset of isoforms of the same gene, pairs of reference proteins interacting with different subsets of isoforms of the same gene, and pairs of reference proteins interacting with protein products of different genes. Disease subnetwork similarity is calculated as the fraction (Jaccard similarity index) of disease subnetworks shared by two proteins, where two proteins share a disease subnetwork if each protein or its interaction partner in the HI-II-14 reference interactome is associated with the disease.

File: IntAct_isoform_partners_disSubNetSim.m

This script first predicts the isoform interactome from the IntAct reference interactome. Then it calculates disease subnetwork similarity for pairs of reference proteins interacting with the same subset of isoforms of the same gene, pairs of reference proteins interacting with different subsets of isoforms of the same gene, and pairs of reference proteins interacting with protein products of different genes. Disease subnetwork similarity is calculated as the fraction (Jaccard similarity index) of disease subnetworks shared by the two proteins, where two proteins share a disease subnetwork if each protein or its interaction partner in the high-quality HI-II-14 reference interactome is associated with the disease.

File: validation.m

This script validates the performance of the domain-based isoform interaction prediction method on the experimental dataset of Yang et al. (2016). Interactions between Orfeome proteins and newly cloned reference proteins are used as the reference interactome and domain-domain interactions are mapped onto these reference interactions. Interactions for the newly cloned alternative isoforms are then predicted from the DDI-annotated reference interactions using two methods:

- First method uses reference interactions having full DDI annotations only. A reference interaction has a full DDI annotation if the newly cloned protein contains an interacting domain of a DDI and its orfeome interaction partner contains the other interacting domain of the same DDI. An alternative isoform of the newly cloned reference protein is then predicted to lose the interaction if it loses all DDI annotations with the orfeome partner, otherwise the interaction is retained.

- Second method uses reference interactions that have either a full DDI annotation or a partial DDI annotation. A reference interaction has a partial DDI annotation if the newly cloned protein contains an interacting domain of a DDI whereas its orfeome interaction partner does not contain the other interacting domain of the DDI. An alternative isoform of the newly cloned reference protein is then predicted to lose the interaction if it loses all of its interacting domains, otherwise the interaction is retained.

