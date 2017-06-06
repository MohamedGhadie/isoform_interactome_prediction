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

% laod processed data (if exists) from .mat file:
% 1 for yes, 0 to process interactome data from scratch
load_processed_data = 1;

% save interactome processed data to .mat file: 
% 1 for yes, 0 otherwise
save_processed_data = 1;

% Calculate p-values for results: 1 for yes, 0 otherwise
calculate_pvalues = 0;

% processed data directory where interactome processed data will be saved
processed_data_dir = 'interactome_processed/';

%-------------------------------------------------------------------------
% First, create disease subnetwork profiles in the HI-II-14 reference
% interactome
disp('Loading and processing HI-II-14 interactome data');
interactome = 'HI-II-14';
interactomeFile = 'HI-II-14.tsv';
isoformInteractomeFile = [processed_data_dir 'HI-II-14_isoform_interactome'];
spEntrezMapFile = 'HI-II-14_spEntrezMap.tab';
numTimesReported = 1;
removeSelfInteractions = 1;
[rolland_I, rolland_PPIs, rolland_spID, genes] = load_interactome(interactomeFile, [], spEntrezMapFile, numTimesReported, removeSelfInteractions, interactome);
rolland_numGenes = size(rolland_I,1);
rolland_numCommonPartners = double(rolland_I)*double(rolland_I);

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

% Create gene-disease association matrix for HI-II-14 interactome
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

% Create disease subnetwork profiles for HI-II-14 interactome
rolland_netDisProfile = zeros(rolland_numGenes,size(rolland_geneDisI,2));
for i = 1:rolland_numGenes
    neighbors = find(rolland_I(i,:));
    rolland_netDisProfile(i,:) = sum(rolland_geneDisI(sort([i neighbors]),:),1)>0;
end
rolland_geneNumDis = sum(rolland_netDisProfile,2);

%-------------------------------------------------------------------------
% Now process IntAct reference interactome

% select interactome: IntAct
interactome = 'IntAct';

% keep only interactions reported this many times or more
numTimesReported = 2;

% remove self interactions,: 1 for yes, 0 otherwise
removeSelfInteractions = 1;

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

% number of genes in IntAct reference interactome
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
% create gene-disease association matrix for proteins in the IntAct
% reference interactome
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

%-------------------------------------------------------------------------
% calculate disease subnetwork profiles for IntAct reference proteins using
% their interaction profiles in HI-II-14

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

%-------------------------------------------------------------------------
% calculate disease subnetwork similarity for protein pairs interacting
% with different subsets of isoforms of the same gene in the IntAct isoform
% interactome using their interaction profiles in HI-II-14

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
            p1num = sum(p1_disProfile);
            p2num = sum(p2_disProfile);
            if p1num>0 && p2num>0
                diff_disSim = [diff_disSim 1-sum(p1_disProfile~=p2_disProfile)/sum(p1_disProfile|p2_disProfile)];
                diff_numDis = [diff_numDis; p1num p2num];
            end
        end
    end
end

%-------------------------------------------------------------------------
% calculate disease subnetwork similarity for protein pairs interacting
% with the same subset of isoforms of the same gene in the IntAct isoform
% interactome using their interaction profiles in HI-II-14

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
            p1num = sum(p1_disProfile);
            p2num = sum(p2_disProfile);
            if p1num>0 && p2num>0
                nodiff_disSim = [nodiff_disSim 1-sum(p1_disProfile~=p2_disProfile)/sum(p1_disProfile|p2_disProfile)];
                nodiff_numDis = [nodiff_numDis; p1num p2num];
            end
        end
    end
end
clear sp1 sp2 geneNum1 geneNum2 p1num p2num p1_disProfile p2_disProfile

%-------------------------------------------------------------------------
% file directory to load disease subnetwork similarity for protein pairs
% interacting with protein products of different genes
diffGeneDisSim_file = [processed_data_dir 'IntAct_diffGene_disSubNetSim.mat'];

% load disease subnetwork similarity for protein pairs interacting with
% protein products of different genes if data file exists
if (load_processed_data == 1) && (exist(diffGeneDisSim_file, 'file') == 2)
    disp('Loading disease subnetwork similarity results for pairs of proteins -');
    disp(['interacting with products of different genes from file ' diffGeneDisSim_file]);
    load(diffGeneDisSim_file);
else
    % calculates disease subnetwork similarity for all protein pairs
    % interacting with protein products of different genes (takes a long time)
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
        if ~isempty(spID{i}) && (geneNum(i) > 0) && (geneNumDisSubNet(i) > 0)
            for j = i+1:numGenes
                if ~isempty(spID{j}) && (geneNum(j)>0) && (geneNumDisSubNet(j) > 0)
                    if (numCommonPartners(i,j) == 0) && (rolland_numCommonPartners(geneNum(i),geneNum(j)) == 0)
                        c = c + 1;
                        diffGene_disSim(c) = 1-sum(intact_disSubNetProfiles(i,:)~=intact_disSubNetProfiles(j,:))/sum(intact_disSubNetProfiles(i,:)|intact_disSubNetProfiles(j,:));
                        diffGene_pairs(c,:) = [i j];
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
