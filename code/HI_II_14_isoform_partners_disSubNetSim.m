% This script first predicts the isoform interactome from the HI-II-14 reference
% interactome. Then it calculates disease subnetwork similarity for pairs
% of reference proteins interacting with the same subset of isoforms of the same
% gene, pairs of reference proteins interacting with different subsets of isoforms of
% the same gene, and pairs of reference proteins interacting with protein products of
% different genes. Disease subnetwork similarity is calculated as the
% fraction (Jaccard similarity index) of disease subnetworks shared by two
% proteins, where two proteins share a disease subnetwork if each protein or
% its interaction partner in the HI-II-14 reference interactome is associated
% with the disease.

% select interactome: HI-II-14 or IntAct
interactome = 'HI-II-14';

% keep only interactions reported this many times or more
numTimesReported = 1;

% remove self interactions,: 1 for yes, 0 otherwise
removeSelfInteractions = 1;     

% laod processed data (if exists) from .mat file:
% 1 for yes, 0 to process interactome data from scratch
load_processed_data = 1;

% save interactome processed data to .mat file: 
% 1 for yes, 0 otherwise
save_processed_data = 1;

% Calculate p-values for results: 1 for yes, 0 otherwise
calculate_pvalues = 1;

% processed data directory where interactome processed data will be saved
processed_data_dir = 'interactome_processed/';

% load an process interactome data
[I, PPIs, spID, genes, domains, DDIs, refDomMap, domRefPos, domAltMap, domAltPos, prSeq, ...
isoSeq, isoNames, altIsoforms, numAltIsoforms, maxIsoform, isoInterDomains, domI, domPrI, ...
numDDImap, numCommonPartners, inIsoInteractome, ref_ref_interactions, alt_ref_interactions, ...
numLosingIsoforms, fracLosingIsoforms, numIsoPairs, fracDiffIsoProfilePerGene, ...
avgIsoProfileDistPerGene, allDist] ...
= process_interactome(interactome, load_processed_data, save_processed_data, ...
                      numTimesReported, removeSelfInteractions, processed_data_dir);

% display results from domain-resolved interactome and predicted isoform interactome 
display_isoform_interactome_results(I,numDDImap,numAltIsoforms,numIsoPairs,inIsoInteractome,...
                                    ref_ref_interactions,alt_ref_interactions,numLosingIsoforms,...
                                    fracLosingIsoforms,fracDiffIsoProfilePerGene,avgIsoProfileDistPerGene, allDist);

% number of genes
numGenes = size(I,1);

% label pairs of reference proteins based on isoform interaction profiles
[pairs,selected,nodiffTargets,subsetTargets,changeoverTargets] = label_iso_partner_pairs(spID,I,domI,domPrI,isoInterDomains,numDDImap,maxIsoform);

% select pairs with identical isoform interaction profiles
nodiff_pairs = pairs(selected(:,1) & ~selected(:,2) & ~selected(:,3),:);
% select pairs with different isoform interaction profiles
diff_pairs = pairs(~selected(:,1) & (selected(:,2)|selected(:,3)),:);

disp([num2str(size(nodiff_pairs,1)) ' protein pairs interacting with the same subset of isoforms of the same gene']);
disp([num2str(size(diff_pairs,1)) ' protein pairs interacting with different subsets of isoforms of the same gene']);

%-------------------------------------------------------------------------
% load gene_disease associations
tdfread('curated_gene_disease_associations.tsv');
clear geneId sourceId score description diseaseName NofPmids NofSnps

geneName = strtrim(num2cell(geneName,2));
diseaseId = strtrim(num2cell(diseaseId,2));
disNames = unique(diseaseId);
numDis = length(disNames);

% load SwissProt-to-geneName mapping table provided by DisGeNET
tdfread('mapa_4_uniprot_crossref.tsv');
disSpGeneMap = [strtrim(num2cell(UniProtKB,2)) strtrim(num2cell(GENE_SYMBOL,2))];
clear UniProtKB GENE_SYMBOL

% create gene-disease association matrix
fprintf('\nCreating gene-disease association matrix\n');
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

% number of diseases per gene
numDis = sum(geneDisI,2);
disp([num2str(sum(numDis>0)) ' genes associated with at least one disease']);
disp([num2str(sum(numDis==0)) ' genes associated with no disease']);

% plot disease distribution
figure
hist(numDis(numDis>0),20);
xlabel('Number of diseases');
ylabel('Number of genes');
title('Number of diseases per gene');

%-------------------------------------------------------------------------
% calculate disease subnetwork similarity for pairs of proteins interacting
% with different subsets of isoforms of the same gene
diff_disSim = [];
for i = 1:size(diff_pairs,1)
    nb1 = I(diff_pairs(i,1),:)>0;
    nb2 = I(diff_pairs(i,2),:)>0;
    nb1(diff_pairs(i,1)) = 1;
    nb2(diff_pairs(i,2)) = 1;
    p1_disProfile = sum(geneDisI(nb1>0,:),1)>0;
    p2_disProfile = sum(geneDisI(nb2>0,:),1)>0;
    p1num = sum(p1_disProfile);
    p2num = sum(p2_disProfile);
    if p1num>0 && p2num>0
        diff_disSim = [diff_disSim 1-sum(p1_disProfile~=p2_disProfile)/sum(p1_disProfile|p2_disProfile)];
    end
end

%-------------------------------------------------------------------------
% calculate disease subnetwork similarity for pairs of proteins interacting
% with the same subset of isoforms of the same gene
nodiff_disSim = [];
for i = 1:size(nodiff_pairs,1)
    nb1 = I(nodiff_pairs(i,1),:)>0;
    nb2 = I(nodiff_pairs(i,2),:)>0;
    nb1(nodiff_pairs(i,1)) = 1;
    nb2(nodiff_pairs(i,2)) = 1;
    p1_disProfile = sum(geneDisI(nb1>0,:),1)>0;
    p2_disProfile = sum(geneDisI(nb2>0,:),1)>0;
    p1num = sum(p1_disProfile);
    p2num = sum(p2_disProfile);
    if p1num>0 && p2num>0
        nodiff_disSim = [nodiff_disSim 1-sum(p1_disProfile~=p2_disProfile)/sum(p1_disProfile|p2_disProfile)];
    end
end
clear p1num p2num p1_neighbors p2_neighbors p1_disProfile p2_disProfile

%-------------------------------------------------------------------------
% file directory to load disease subnetwork similarity for protein pairs
% interacting with protein products of different genes
diffGeneDisSim_file = [processed_data_dir 'HI-II-14_diffGene_disSubNetSim.mat'];

% load disease subnetwork similarity for protein pairs interacting with
% protein products of different genes if data file exists
if (load_processed_data == 1) && (exist(diffGeneDisSim_file, 'file') == 2)
    disp('Loading disease subnetwork similarity results for pairs of proteins -');
    disp(['interacting with products of different genes from file ' diffGeneDisSim_file]);
    load(diffGeneDisSim_file)
else
    % calculates disease subnetwork similarity for all protein pairs
    % interacting with protein products of different genes (takes a long time)
    disp('Calculating disease subnetwork similarity for pairs of proteins -');
    disp('interacting with products of different genes');
    netDisProfile = zeros(numGenes,size(geneDisI,2));
    for i = 1:numGenes
        neighbors = find(I(i,:));
        netDisProfile(i,:) = sum(geneDisI(sort([i neighbors]),:),1)>0;
    end
    geneNumDis = sum(netDisProfile,2);
    
    numNeighbors = sum(I,2);
    maxsum = cumsum(1:numGenes-1);
    maxnum = maxsum(numGenes-1);
    diffGene_disSim = zeros(1,maxnum);
    diffGene_pairs = zeros(maxnum,2);
    c = 0;
    for i = 1:numGenes-1
        if ~mod(i,1000)
            disp(['- Calculated disease subnetwork similarity for ' num2str(i) ' genes out of ' num2str(numGenes) ' genes']);
        end
        if geneNumDis(i) > 0
            p1_disProfile = netDisProfile(i,:);
            for j = i+1:numGenes
                if (geneNumDis(j) > 0) && numCommonPartners(i,j) == 0
                    c = c + 1;
                    p2_disProfile = netDisProfile(j,:);
                    diffGene_disSim(c) = 1-sum(p1_disProfile~=p2_disProfile)/sum(p1_disProfile|p2_disProfile);
                    diffGene_pairs(c,:) = [i j];
                end
            end
        end
    end
    diffGene_disSim = diffGene_disSim(1:c);
    diffGene_pairs = diffGene_pairs(1:c,:);
    clear i j c maxsum maxnum
    
    if save_processed_data == 1
        save(diffGeneDisSim_file, 'diffGene_disSim', 'diffGene_pairs');
        disp('Disease subnetwork similarity results for pairs of proteins -');
        disp(['interacting with protein products of different genes were saved in file ' diffGeneDisSim_file]);
    end
end

disp('Mean fraction of disease subnetwoks shared by protein pairs interacting with:');
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
    
    itr = 1000;
    disp('Significance in difference between protein pairs with different isoform interaction profiles');
    disp('and protein pairs interacting with protein products of different genes:');
    bootMeanSim = bootstrapSim(diff_disSim,diffGene_disSim,itr);
    p = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_disSim)-mean(diffGene_disSim)))/size(bootMeanSim,1);
    disp(['bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p)]);
end

% plot results
plotdata = [mean(nodiff_disSim); mean(diff_disSim); mean(diffGene_disSim)];
% standard error of the mean
sderror = [std(nodiff_disSim)/sqrt(length(nodiff_disSim)); std(diff_disSim)/sqrt(length(diff_disSim));...
    std(diffGene_disSim)/sqrt(length(diffGene_disSim))];
figure
hold on
hb = bar(plotdata,0.6,'EdgeColor','none');
set(hb,'facecolor',[1 102/255 102/255]);
ylim([0 0.2]);
set(gca,'XTick',1:3,'XTickLabel',{' Interacting with\newlinethe same subset\newline  of isoforms of\newline the same gene',...
    ' Interacting with\newlinedifferent subsets\newline  of isoforms of\newline the same gene',...
    '   Interacting with\newline  protein products\newline of different genes'},'FontSize',16);
ylabel('Fraction of disease subnetworks\newlineshared by pairs of proteins','FontSize',22);
set(gca,'YGrid','on');
set(gca,'tickDir','out');
box off
pause(0.1);
% plot error bars
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    for j = 1:length(xData)
        line([xData(j) xData(j)], [plotdata(j,ib)-sderror(j,ib) plotdata(j,ib)+sderror(j,ib)],'Color','k','LineWidth',1.5);
    end
end
clear hb ib j xData plotdata sderror
