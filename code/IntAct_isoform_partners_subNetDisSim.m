% This script first predicts the isoform interactome from the IntAct reference
% interactome. Then it calculates disease subnetwork similarity for pairs
% of reference proteins interacting with the same subset of isoforms of the same
% gene, pairs of reference proteins interacting with different subsets of isoforms of
% the same gene, and pairs of reference proteins interacting with protein products of
% different genes. Disease subnetwork similarity is calculated as the
% fraction (Jaccard similarity index) of disease subnetworks shared by two
% proteins, where two proteins share a disease subnetwork if each protein or
% its interaction partner in the HI-II-14 reference interactome is associated
% with the disease.

load_processed_data = 1;
save_processed_data = 1;
calculate_pvalues = 1;
processed_data_dir = 'interactome_processed/';
processed_data_exists = 0;

disp('Loading and processing HI-II-14 interactome data');
interactome = 'HI-II-14';
interactomeFile = 'HI-II-14.xlsx';
isoformInteractomeFile = [processed_data_dir 'HI-II-14_isoform_interactome'];
spEntrezMapFile = 'HI-II-14_spEntrezMap.tab';
numTimesReported = 1;
removeSelfInteractions = 1;
[rolland_I, rolland_PPIs, rolland_spID, genes] = load_interactome(interactomeFile, [], spEntrezMapFile, numTimesReported, removeSelfInteractions, interactome);
rolland_numGenes = size(rolland_I,1);

if (load_processed_data == 1) && (exist([processed_data_dir 'HI-II-14_data.mat'], 'file') == 2)
    disp('Loading number of common partners for each pair of reference proteins in HI-II-14 interactome');
    load([processed_data_dir 'HI-II-14_data.mat'], 'numCommonPartners');
    rolland_numCommonPartners = numCommonPartners;
    clear numCommonPartners
else
    disp('Calculating number of common partners for each pair of reference proteins in HI-II-14 interactome');
    %rolland_numCommonPartners = get_num_common_partners(rolland_I);
end

% load gene_disease associations
tdfread('curated_gene_disease_associations.tsv');
clear geneId sourceId score description diseaseName NofPmids NofSnps

geneName = strtrim(num2cell(geneName,2));
diseaseId = strtrim(num2cell(diseaseId,2));
disNames = unique(diseaseId);
numDis = length(disNames);

% load SwissProt-to-geneName mapping table
tdfread('mapa_4_uniprot_crossref.tsv');
disSpGeneMap = [strtrim(num2cell(UniProtKB,2)) strtrim(num2cell(GENE_SYMBOL,2))];
clear UniProtKB GENE_SYMBOL

fprintf('\nCreating gene-disease association matrix for HI-II-14 interactome\n');
rolland_geneDisI = zeros(rolland_numGenes,numDis);
for i = 1:rolland_numGenes
    spid = rolland_spID{i};
    if ~isempty(spid)
        gn = unique(disSpGeneMap(strcmpi(disSpGeneMap(:,1),spid),2));
        dis = unique(diseaseId(ismember(geneName,gn)));
        if ~isempty(dis)
            rolland_geneDisI(i,ismember(disNames,dis)) = 1;
        end
    end
end
clear spid gn dis i

rolland_netDisProfile = zeros(rolland_numGenes,size(rolland_geneDisI,2));
for i = 1:rolland_numGenes
    neighbors = find(rolland_I(i,:));
    rolland_netDisProfile(i,:) = sum(rolland_geneDisI(sort([i neighbors]),:),1)>0;
end
rolland_geneNumDis = sum(rolland_netDisProfile,2);

% Now calculate subnetwork disease similarity for IntAct proteins using
% their subnetwork disease profiles in HI-II-14

if (load_processed_data == 1) && (exist([processed_data_dir 'IntAct_data.mat'], 'file') == 2)
    disp(['Loading IntAct interactome processed data from file ' processed_data_dir 'IntAct_data.mat']);
    processed_data_exists = 1;
    load([processed_data_dir 'IntAct_data.mat']);
else
    disp('Loading and processing IntAct interactome data');
    interactome = 'IntAct';
    interactomeFile = 'intact.txt';
    isoformInteractomeFile = [processed_data_dir 'IntAct_isoform_interactome'];
    spGeneMapFile = 'IntAct_spGeneMap.tab.txt';
    numTimesReported = 2;
    removeSelfInteractions = 1;
    [I, PPIs, spID, genes] = load_interactome(interactomeFile, spGeneMapFile, [], numTimesReported, removeSelfInteractions, interactome);
    outputfile = [processed_data_dir 'IntAct_data.mat'];
end

if processed_data_exists == 0
    numPPI = sum(sum(triu(I)));
    numGenes = size(I,1);
    did3File = '3did_flat.txt';
    domineFile = 'interaction.xlsx';
    processedDDIfile = [processed_data_dir 'DDIs.txt'];
    disp('Loading domain-domain interactions');
    if exist(processedDDIfile, 'file') == 2
        tdfread(processedDDIfile);
        DDIs = cell(size(dom1,1),2);
        for i = 1:size(dom1,1)
            DDIs(i,:) = {strtrim(dom1(i,:)), strtrim(dom2(i,:))};
        end
        clear dom1 dom2
    elseif (exist(did3File, 'file') == 2) && (exist(domineFile, 'file') == 2)
        DDIs = load_3did_and_domine_DDIs(did3File, domineFile);
        fid = fopen(processedDDIfile,'w');
        fprintf(fid,'dom1\tdom2\n');
        for i = 1:size(DDIs,1)
            fprintf(fid,[DDIs{i,1} '\t' DDIs{i,2} '\n']);
        end
        fclose(fid);
    else
        disp('Domain-domain interaction file not found. Exiting script');
        return
    end
    domains = unique([DDIs(:,1); DDIs(:,2)]);
    numDom = length(domains);
    
    if exist('hmmer_canonical_full.domtab.txt', 'file') == 2
        disp('Loading structural domains for reference proteins');
        [refDomMap, domRefPos] = load_hmmscan_ref('hmmer_canonical_full.domtab.txt');
    else
        disp('Reference-isoform Domain-mapping file (hmmer_isoforms.domtab.txt) not found. Exiting script');
        return
    end
    
    if exist('hmmer_isoforms.domtab.txt', 'file') == 2
        disp('Loading structural domains for alternative isoforms');
        [domAltMap, domAltPos] = load_hmmscan_alt('hmmer_isoforms.domtab.txt');
    else
        disp('Alternative-isoform Domain-mapping file (hmmer_isoforms.domtab.txt) not found. Exiting script');
        return
    end
    
    if exist('human_protein_sequences.tab', 'file') == 2
        disp('Loading reference protein sequences');
        [headers, prSeq] = load_ref_protein_seq('human_protein_sequences.tab');
    else
        disp('Human reference protein sequence file (human_protein_sequences.tab) not found. Exiting script');
        return
    end
    
    if exist('human_protein_can_iso_sequences.fasta', 'file') == 2
        disp('Loading alternative isoform sequences');
        isoSeq = load_alt_protein_seq('human_protein_can_iso_sequences.fasta');
    else
        disp('Human alternative protein sequence file (human_protein_can_iso_sequences.fasta) not found. Exiting script');
        return
    end
    isoNames = unique(isoSeq(:,1));
    
    disp('Counting number of alternative isoforms per gene');
    [altIsoforms, numAltIsoforms, maxIsoform] = getIsoforms(spID, isoNames);
    
    disp('Identifying isoform interacting domains');
    isoInterDomains = get_isoform_interacting_domains(spID, maxIsoform, domains, refDomMap, domAltMap);
    
    disp('Creating domain interaction matrix');
    domI = create_DDI_matrix(domains, DDIs);
    
    disp('Creating domain-protein mapping matrix');
    domPrI = create_domain_protein_matrix(domains, refDomMap, spID);
    
    disp('Counting number of domain-domain interactions mapping to each protein-protein interaction');
    numDDImap = create_numDDI_matrix(I, domPrI, domI);
    
    disp('Calculating number of common partners for each pair of reference proteins');
    numCommonPartners = get_num_common_partners(I);
    
    disp('Predicting isoform interactome');
    [inIsoInteractome, ref_ref_interactions, alt_ref_interactions, ~, ~] = predict_isoform_interactome(spID,domains,PPIs,I,domI,domPrI,numDDImap,altIsoforms,isoInterDomains,interactome,isoformInteractomeFile);
    
    disp('Comparing isoform interaction profiles');
    [numLosingIsoforms, fracLosingIsoforms, numIsoPairs, fracDiffIsoProfilePerGene, avgIsoProfileDistPerGene] = compare_isoform_interactions(I,domI,domPrI,numDDImap,altIsoforms,numAltIsoforms,isoInterDomains);
    
    if save_processed_data == 1
        if exist(processed_data_dir, 'dir') ~= 7
            mkdir(processed_data_dir);
        end
        save(outputfile,'interactome','interactomeFile','isoformInteractomeFile', 'numTimesReported', 'removeSelfInteractions',...
            'PPIs', 'I', 'spID', 'genes', 'numPPI', 'numGenes', 'DDIs', 'domains', 'numDom', 'refDomMap',...
            'domRefPos', 'domAltMap', 'domAltPos', 'headers', 'prSeq', 'isoSeq', 'isoNames', 'altIsoforms', 'numAltIsoforms',...
            'maxIsoform', 'isoInterDomains', 'domI', 'domPrI', 'numDDImap', 'numCommonPartners', 'inIsoInteractome',...
            'ref_ref_interactions', 'alt_ref_interactions', 'numLosingIsoforms', 'fracLosingIsoforms', 'numIsoPairs',...
            'fracDiffIsoProfilePerGene', 'avgIsoProfileDistPerGene');
        disp(['Interactome processed data saved to file ' outputfile]);
    end
end

fprintf('\nIntAct Predicted isoform interactome statistics:\n');
disp([num2str(length(unique(ref_ref_interactions(:,1:2)))) ' reference proteins and ' num2str(length(unique(alt_ref_interactions(:,1)))) ' alternative isoforms']);
disp([num2str(size(ref_ref_interactions,1)) ' known reference interactions']);
disp([num2str(sum(strcmpi(alt_ref_interactions(:,3),'retained'))) ' predicted retained isoform interactions']);
disp([num2str(sum(strcmpi(alt_ref_interactions(:,3),'lost'))) ' predicted lost isoform interactions']);
disp([num2str(sum(numAltIsoforms(sum(numDDImap)>0)>0)) ' genes have at least two isoforms (including reference)']);
disp([num2str(sum(numLosingIsoforms>0)) ' genes have at least one isoform losing an interaction']);
fprintf('\n');

[pairs,selected,nodiffTargets,subsetTargets,changeoverTargets] = label_iso_partner_pairs(spID,I,domI,domPrI,isoInterDomains,numDDImap,maxIsoform);

nodiff_pairs = pairs(selected(:,1) & ~selected(:,2) & ~selected(:,3),:);
diff_pairs = pairs(~selected(:,1) & (selected(:,2)|selected(:,3)),:);

[numNodiffPairsRepeated, nodiff_pairs] = remove_gene_pair_duplicates(nodiff_pairs,genes);
disp([num2str(numNodiffPairsRepeated) ' protein pairs with identical isoform interaction profiles occur more than once with different uniProt IDs: duplicates removed']);

[numDiffPairsRepeated, diff_pairs] = remove_gene_pair_duplicates(diff_pairs,genes);
disp([num2str(numDiffPairsRepeated) ' protein pairs with different isoform interaction profiles occur more than once with different uniProt IDs: duplicates removed']);

disp([num2str(size(nodiff_pairs,1)) ' protein pairs interacting with the same subset of isoforms of the same gene']);
disp([num2str(size(diff_pairs,1)) ' protein pairs interacting with different subsets of isoforms of the same gene']);

fprintf('\nCreating gene-disease association matrix for genes in IntAct interactome\n');
geneDisI = zeros(numGenes,numDis);
for i = 1:numGenes
    spid = spID{i};
    if ~isempty(spid)
        gn = unique(disSpGeneMap(strcmpi(disSpGeneMap(:,1),spid),2));
        dis = unique(diseaseId(ismember(geneName,gn)));
        if ~isempty(dis)
            geneDisI(i,ismember(disNames,dis)) = 1;
        end
    end
end
clear spid gn dis i
geneNumDis = sum(geneDisI,2);

rolland_numNeighbors = sum(rolland_I,2);
intact_numNeighbors = sum(I,2);
intact_disSubNetProfiles = zeros(size(geneDisI));
geneNum = zeros(numGenes,1);
for i = 1:numGenes
    sp = spID{i};
    if ~isempty(sp)
        num = find(strcmpi(rolland_spID,sp),1);
        if ~isempty(num)
            geneNum(i) = num;
        end
        if geneNum(i)>0
            intact_disSubNetProfiles(i,:) = (geneDisI(i,:)>0) | (rolland_netDisProfile(geneNum(i),:)>0);
        else
            intact_disSubNetProfiles(i,:) = geneDisI(i,:)>0;
        end
    end
end
geneNumDisSubNet = sum(intact_disSubNetProfiles,2);

diff_disSim = [];
diff_numDis = [];
for i = 1:size(diff_pairs,1)
    sp1 = spID{diff_pairs(i,1)};
    sp2 = spID{diff_pairs(i,2)};
    if ~isempty(sp1) && ~isempty(sp2)
        geneNum1 = geneNum(diff_pairs(i,1));
        geneNum2 = geneNum(diff_pairs(i,2));
        if geneNum1>0 && geneNum2>0
            nb1 = rolland_I(geneNum1,:)>0;
            nb2 = rolland_I(geneNum2,:)>0;
            nb1(geneNum1) = 1;
            nb2(geneNum2) = 1;
            p1_disProfile = (geneDisI(diff_pairs(i,1),:) + sum(rolland_geneDisI(nb1,:),1)) > 0;
            p2_disProfile = (geneDisI(diff_pairs(i,2),:) + sum(rolland_geneDisI(nb2,:),1)) > 0;
        elseif geneNum1>0
            nb1 = rolland_I(geneNum1,:)>0;
            nb1(geneNum1) = 1;
            p1_disProfile = (geneDisI(diff_pairs(i,1),:) + sum(rolland_geneDisI(nb1,:),1)) > 0;
            p2_disProfile = geneDisI(diff_pairs(i,2),:) > 0;
        elseif geneNum2>0
            nb2 = rolland_I(geneNum2,:)>0;
            nb2(geneNum2) = 1;
            p1_disProfile = geneDisI(diff_pairs(i,1),:) > 0;
            p2_disProfile = (geneDisI(diff_pairs(i,2),:) + sum(rolland_geneDisI(nb2,:),1)) > 0;
        else
            p1_disProfile = geneDisI(diff_pairs(i,1),:) > 0;
            p2_disProfile = geneDisI(diff_pairs(i,2),:) > 0;
        end
        p1num = sum(p1_disProfile);
        p2num = sum(p2_disProfile);
        if p1num>0 && p2num>0
            diff_disSim = [diff_disSim 1-sum(p1_disProfile~=p2_disProfile)/sum(p1_disProfile|p2_disProfile)];
            diff_numDis = [diff_numDis; p1num p2num];
        end
    end
end

nodiff_disSim = [];
nodiff_numDis = [];
for i = 1:size(nodiff_pairs,1)
    sp1 = spID{nodiff_pairs(i,1)};
    sp2 = spID{nodiff_pairs(i,2)};
    if ~isempty(sp1) && ~isempty(sp2)
        geneNum1 = geneNum(nodiff_pairs(i,1));
        geneNum2 = geneNum(nodiff_pairs(i,2));
        if geneNum1>0 && geneNum2>0
            nb1 = rolland_I(geneNum1,:)>0;
            nb2 = rolland_I(geneNum2,:)>0;
            nb1(geneNum1) = 1;
            nb2(geneNum2) = 1;
            p1_disProfile = (geneDisI(nodiff_pairs(i,1),:) + sum(rolland_geneDisI(nb1,:),1)) > 0;
            p2_disProfile = (geneDisI(nodiff_pairs(i,2),:) + sum(rolland_geneDisI(nb2,:),1)) > 0;
        elseif geneNum1>0
            nb1 = rolland_I(geneNum1,:)>0;
            nb1(geneNum1) = 1;
            p1_disProfile = (geneDisI(nodiff_pairs(i,1),:) + sum(rolland_geneDisI(nb1,:),1)) > 0;
            p2_disProfile = geneDisI(nodiff_pairs(i,2),:) > 0;
        elseif geneNum2>0
            nb2 = rolland_I(geneNum2,:)>0;
            nb2(geneNum2) = 1;
            p1_disProfile = geneDisI(nodiff_pairs(i,1),:) > 0;
            p2_disProfile = (geneDisI(nodiff_pairs(i,2),:) + sum(rolland_geneDisI(nb2,:),1)) > 0;
        else
            p1_disProfile = geneDisI(nodiff_pairs(i,1),:) > 0;
            p2_disProfile = geneDisI(nodiff_pairs(i,2),:) > 0;
        end
        p1num = sum(p1_disProfile);
        p2num = sum(p2_disProfile);
        if p1num>0 && p2num>0
            nodiff_disSim = [nodiff_disSim 1-sum(p1_disProfile~=p2_disProfile)/sum(p1_disProfile|p2_disProfile)];
            nodiff_numDis = [nodiff_numDis; p1num p2num];
        end
    end
end
clear sp1 sp2 geneNum1 geneNum2 p1num p2num p1_disProfile p2_disProfile

diffGeneDisSim_file = [processed_data_dir 'IntAct_diffGene_disSubNetSim.mat'];

if (load_processed_data == 1) && (exist(diffGeneDisSim_file, 'file') == 2)
    disp('Loading disease subnetwork similarity results for pairs of proteins -');
    disp(['interacting with products of different genes from file ' diffGeneDisSim_file]);
    load(diffGeneDisSim_file);
else
    disp('Calculating disease subnetwork similarity for pairs of proteins -');
    disp('interacting with products of different genes');
    maxsum = cumsum(1:numGenes-1);
    maxnum = maxsum(numGenes-1);
    diffGene_disSim = zeros(1,maxnum);
    diffGene_pairs = zeros(maxnum,2);
    c = 0;
    for i = 1:numGenes-1
        if ~mod(i,1000)
            disp(['- Calculated disease subnetwork similarity for ' num2str(i) ' genes out of ' num2str(numGenes) ' genes']);
        end
        sp1 = spID{i};
        if ~isempty(sp1)
            if geneNumDisSubNet(i) > 0
                for j = i+1:numGenes
                    if geneNumDisSubNet(j) > 0
                        sp2 = spID{j};
                        if ~isempty(sp2)
                            if numCommonPartners(i,j) == 0
                                if geneNum(i)>0 && geneNum(j)>0
                                    if rolland_numCommonPartners(geneNum(i),geneNum(j)) > 0
                                        continue
                                    end
                                end
                                c = c + 1;
                                diffGene_disSim(c) = 1-sum(intact_disSubNetProfiles(i,:)~=intact_disSubNetProfiles(j,:))/sum(intact_disSubNetProfiles(i,:)|intact_disSubNetProfiles(j,:));
                                diffGene_pairs(c,:) = [i j];
                            end
                        end
                    end
                end
            end
        end
    end
    diffGene_disSim = diffGene_disSim(1:c);
    diffGene_pairs = diffGene_pairs(1:c,:);
    clear i j c maxsum maxnum
end

if save_processed_data == 1
    save(diffGeneDisSim_file, 'diffGene_disSim', 'diffGene_pairs');
    disp('Disease similarity results for pairs of proteins interacting with -');
    disp(['protein products of different genes were saved in file ' diffGeneDisSim_file]);
end

disp('Mean fraction of diseases shared by protein pairs interacting with:');
disp(['1. same isoforms of the same gene (identical profiles): ' num2str(mean(nodiff_disSim))]);
disp(['2. different isoforms of the same gene (different profiles): ' num2str(mean(diff_disSim))]);
disp(['3. protein products of different genes: ' num2str(mean(diffGene_disSim))]);

% bootstrap test
if calculate_pvalues == 1
    itr = 100000;
    disp('Significance in difference between protein pairs with identical and different isoform interaction profiles:');
    bootMeanSim = bootstrapSim(nodiff_disSim,diff_disSim,itr);
    p = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_disSim)-mean(diff_disSim)))/size(bootMeanSim,1);
    disp(['bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p)]);
    
    itr = 100;
    disp('Significance in difference between protein pairs with different isoform interaction profiles');
    disp('and protein pairs interacting with protein products of different genes:');
    bootMeanSim = bootstrapSim(diff_disSim,diffGene_disSim,itr);
    p = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_disSim)-mean(diffGene_disSim)))/size(bootMeanSim,1);
    disp(['bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p)]);
end

plotdata = [mean(nodiff_disSim); mean(diff_disSim); mean(diffGene_disSim)];
sderror = [std(nodiff_disSim)/sqrt(length(nodiff_disSim)); std(diff_disSim)/sqrt(length(diff_disSim));...
    std(diffGene_disSim)/sqrt(length(diffGene_disSim))];
figure
hold on
hb = bar(plotdata,0.6,'FaceColor','c','LineWidth',1.5,'EdgeColor','k');
ylim([0 0.05]);
set(gca,'XTick',1:3,'XTickLabel',{' Interacting with\newlinethe same subset\newline  of isoforms of\newline the same gene',...
    ' Interacting with\newlinedifferent subsets\newline  of isoforms of\newline the same gene',...
    '   Interacting with\newline  protein products\newline of different genes'},'FontSize',16);
ylabel('Fraction of disease subnetworks\newlineshared by pairs of proteins','FontSize',22);
set(gca,'YGrid','on');
set(gca,'tickDir','out');
box off
pause(0.1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    for j = 1:length(xData)
        line([xData(j) xData(j)], [plotdata(j,ib)-sderror(j,ib) plotdata(j,ib)+sderror(j,ib)],'Color','k','LineWidth',1.5);
    end
end
clear hb ib j xData plotdata sderror
