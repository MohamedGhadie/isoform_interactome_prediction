% This script first predicts the isoform interactome from the reference
% interactome, either HI-II-14 or IntAct interactome. Then it calculates
% tissue specificity of reference proteins which interact, on average, with
% <=50% of the isoforms of their partner and those that interact, on
% average, with >50% of the isoforms of their partner. Tissue specificity
% is calculated using the Tau index.

interactome = 'IntAct';
load_processed_data = 1;
save_processed_data = 1;
processed_data_dir = 'interactome_processed/';
processed_data_exists = 0;
if strcmpi(interactome,'HI-II-14')
    if (load_processed_data == 1) && (exist([processed_data_dir 'HI-II-14_data.mat'], 'file') == 2)
        disp(['Loading interactome processed data from file ' processed_data_dir 'HI-II-14_data.mat']);
        processed_data_exists = 1;
        load([processed_data_dir 'HI-II-14_data.mat']);
    else
        disp('Loading and processing interactome data');
        interactomeFile = 'HI-II-14.xlsx';
        isoformInteractomeFile = [processed_data_dir 'HI-II-14_isoform_interactome.txt'];
        spEntrezMapFile = 'HI-II-14_spEntrezMap.tab';
        numTimesReported = 1;
        removeSelfInteractions = 1;
        [I, PPIs, spID, genes] = load_interactome(interactomeFile, [], spEntrezMapFile, numTimesReported, removeSelfInteractions, interactome);
        outputfile = [processed_data_dir 'HI-II-14_data.mat'];
    end
elseif strcmpi(interactome,'IntAct')
    if (load_processed_data == 1) && (exist([processed_data_dir 'IntAct_data.mat'], 'file') == 2)
        disp(['Loading interactome processed data from file ' processed_data_dir 'IntAct_data.mat']);
        processed_data_exists = 1;
        load([processed_data_dir 'IntAct_data.mat']);
    else
        disp('Loading and processing interactome data');
        interactomeFile = 'intact.txt';
        isoformInteractomeFile = [processed_data_dir 'IntAct_isoform_interactome.txt'];
        spGeneMapFile = 'IntAct_spGeneMap.tab.txt';
        numTimesReported = 2;
        removeSelfInteractions = 1;
        [I, PPIs, spID, genes] = load_interactome(interactomeFile, spGeneMapFile, [], numTimesReported, removeSelfInteractions, interactome);
        outputfile = [processed_data_dir 'IntAct_data.mat'];
    end
else
    disp('Invalid interactome name. Exiting script');
    return
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
        if strcmpi(interactome,'HI-II-14')
            save(outputfile, 'spEntrezMapFile', '-append');
        elseif strcmpi(interactome,'IntAct')
            save(outputfile, 'spGeneMapFile', '-append');
        end
        disp(['Interactome processed data saved to file ' outputfile]);
    end
end

fprintf('\nPredicted isoform interactome statistics:\n');
disp([num2str(length(unique(ref_ref_interactions(:,1:2)))) ' reference proteins and ' num2str(length(unique(alt_ref_interactions(:,1)))) ' alternative isoforms']);
disp([num2str(size(ref_ref_interactions,1)) ' known reference interactions']);
disp([num2str(sum(strcmpi(alt_ref_interactions(:,3),'retained'))) ' predicted retained isoform interactions']);
disp([num2str(sum(strcmpi(alt_ref_interactions(:,3),'lost'))) ' predicted lost isoform interactions']);
disp([num2str(sum(numAltIsoforms(sum(numDDImap)>0)>0)) ' genes have at least two isoforms (including reference)']);
disp([num2str(sum(numLosingIsoforms>0)) ' genes have at least one isoform losing an interaction']);
fprintf('\n');

% Calculate the mean fraction of isoforms of the same gene each protein interacts with
disp('Calculating the mean fraction of isoforms of the same gene each protein interacts with');
numPartners = zeros(numGenes,1);
numDDIpartners = zeros(numGenes,1);
numDDIpartners_2plusIso = zeros(numGenes,1);
meanIsoInterFraction = zeros(numGenes,1);
for i = 1:numGenes
    spid = spID{i};
    % pick genes with at least one alternative isoform
    if ~isempty(spid)
        partners = find(I(i,:));
        partners(partners==i) = [];
        numPartners(i) = length(partners);
        p1domains = find(domPrI(:,i));
        if numPartners(i) > 0
            % pick partners that have a DDI mapping with the canonical of
            % gene i
            for j = numPartners(i):-1:1
                if numDDImap(i,partners(j)) == 0
                    partners(j) = [];
                end
            end
            numDDIpartners(i) = length(partners);
            if numDDIpartners(i) > 0
                for j = numDDIpartners(i):-1:1
                    if numAltIsoforms(partners(j)) == 0
                        partners(j) = [];
                    end
                end
                numDDIpartners_2plusIso(i) = length(partners);
                if numDDIpartners_2plusIso(i) > 0
                    isoInterFraction = zeros(numDDIpartners_2plusIso(i),1);
                    for j = 1:numDDIpartners_2plusIso(i)
                        p2domains = find(domPrI(:,partners(j)));
                        p1domPartners = p2domains(sum(domI(p1domains,p2domains),1)>0);
                        numP1domPartners = length(p1domPartners);
                        if numP1domPartners == 0
                            disp(['Skipping partner ' num2str(j) '. No DDI annotation.']);
                            continue
                        end
                        isoInteractions = 1;
                        isos = altIsoforms{partners(j)};
                        for m = 1:length(isos)
                            if isos(m)+1 > length(isoInterDomains{partners(j)})
                                disp(['Isoform number out of bound for gene ' num2str(partners(j)) '. Skipping']);
                                continue;
                            end
                            isoDom = isoInterDomains{partners(j)}{isos(m)+1};
                            if ~isempty(isoDom)
                                for n = 1:numP1domPartners
                                    if ~isempty(find(isoDom==p1domPartners(n),1))
                                        isoInteractions = isoInteractions + 1;
                                        break
                                    end
                                end
                            end
                        end
                        isoInterFraction(j) = isoInteractions/(numAltIsoforms(partners(j))+1);
                    end
                    meanIsoInterFraction(i) = mean(isoInterFraction);
                end
            end
        end
    end
end
clear i j m n spid partners p1domains p2domains p1domPartners numP1domPartners ...
    isoInteractions isoDom isoInterFraction isos

figure
hist(meanIsoInterFraction(meanIsoInterFraction>0),100);
xlabel('Mean fraction of same-gene alternative isoform partners');
ylabel('Number of proteins in interactome');

% load gene tissue expression data
if exist('E-MTAB-513.tsv.txt', 'file') == 2
    disp('Loading gene tissue-expression data');
    tdfread('E-MTAB-513.tsv.txt');
else
    disp('Human gene tissue expression file (E-MTAB-513.tsv.txt) not found. Exiting script');
    return
end
expr = [adipose adrenal brain breast colon heart kidney leukocyte liver lung lymph_node ovary prostate skeletal_muscle testis thyroid];
expr = log2(expr);
clear adipose adrenal brain breast colon heart kidney leukocyte liver lung lymph_node ovary prostate skeletal_muscle testis thyroid

expr_gene_id = cell(size(Gene_ID,1),1);
expr_gene_names = cell(size(Gene_Name,1),1);
for i = 1:size(Gene_ID,1)
    expr_gene_id{i} = strtrim(Gene_ID(i,:));
    expr_gene_names{i} = strtrim(Gene_Name(i,:));
end
clear i Gene_ID Gene_Name

% identify rows with no expression (-inf) in all tissues
tobeRemoved = [];
for i = 1:size(expr,1)
    if isempty(find(expr(i,:)>-inf,1))
        tobeRemoved = [tobeRemoved i];
    end
end
fprintf('\nThe following rows of the expression matrix have no expression in all tissues:\n');
disp(num2str(tobeRemoved));

% remove rows with no expression in all tissues
expr(tobeRemoved,:) = [];
expr_gene_id(tobeRemoved) = [];
expr_gene_names(tobeRemoved) = [];

% create expression matrix for interacting proteins
disp('Creating tissue expression matrix for all genes in the reference interactome');
numTissues = size(expr,2);
Iexpr = NaN(numGenes,numTissues);
c = [];
for i = 1:numGenes
    if ~isempty(genes{i})
        ind = find(strcmpi(expr_gene_names,genes{i}));
        if length(ind) == 1
            Iexpr(i,:) = expr(ind,:);
        elseif length(ind) > 1
            c = [c i];
        end
    end
end
fprintf('\nThe following genes mapping to more than one expression row were excluded:\n');
disp(num2str(c));
clear i ind tobeRemoved expr_gene_id expr_gene_names c

% calculate tissue specificity using Tau index
tau_ts = NaN(numGenes,1);
for i = 1:numGenes
    notnan = ~isnan(Iexpr(i,:));
    num = sum(notnan);
    if num >= 8
        shiftedExpr = Iexpr(i,notnan) - min(min(Iexpr(i,notnan)),0);
        exprRatio = shiftedExpr/max(shiftedExpr);
        tau_ts(i) = sum(1-exprRatio)/(num-1);
    end
end

figure
hist(tau_ts);
xlabel('Tissue specificity (tau index)');
ylabel('Number of proteins in interactome');

proteomeSize = size(expr,1);
proteome_tau_ts = NaN(proteomeSize,1);
for i = 1:proteomeSize
    notnan = ~isnan(expr(i,:));
    num = sum(notnan);
    if num >= 8
        shiftedExpr = expr(i,notnan) - min(min(expr(i,notnan)),0);
        exprRatio = shiftedExpr/max(shiftedExpr);
        proteome_tau_ts(i) = sum(1-exprRatio)/(num-1);
    end
end

proteome_tau_ts = proteome_tau_ts(~isnan(proteome_tau_ts));

cutoff = 0:0.001:1;
percProteomeTS = zeros(1,length(cutoff));
for i = 1:length(cutoff)
    percProteomeTS(i) = 100*sum(proteome_tau_ts>cutoff(i))/length(proteome_tau_ts);
end
figure
hold on
plot(cutoff,percProteomeTS);
line([0 0.45],[52 52],'LineStyle','--');
line([0.45 0.45],[0 52],'LineStyle','--');
xlabel('tau cutoff','FontSize',16);
ylabel('Percentage of tissue-specific genes in the human genome','FontSize',15);
set(gca,'XGrid','on','YGrid','on');
set(gca,'tickDir','out');
box off

figure
hist(proteome_tau_ts);
xlabel('Tissue specificity (tau index)');
ylabel('Number of proteins in human proteome');

Specific = tau_ts>0.45;
nonSpecific = tau_ts<=0.45;

disp(['Fraction of tissue-specific proteins in the reference interactome: ' num2str(sum(Specific)/(sum(Specific)+sum(nonSpecific)))]);
disp(['Fraction of non-tissue-specific proteins in the reference interactome: ' num2str(sum(nonSpecific)/(sum(Specific)+sum(nonSpecific)))]);

Specific_isoI = Specific & (inIsoInteractome>0);
nonSpecific_isoI = nonSpecific & (inIsoInteractome>0);
disp(['Fraction of tissue-specific proteins in the isoform interactome: ' num2str(sum(Specific_isoI)/(sum(Specific_isoI)+sum(nonSpecific_isoI)))]);
disp(['Fraction of non-tissue-specific proteins in the isoform interactome: ' num2str(sum(nonSpecific_isoI)/(sum(Specific_isoI)+sum(nonSpecific_isoI)))]);

Specific_2plusIso = Specific_isoI & (numDDIpartners_2plusIso>0);
nonSpecific_2plusIso = nonSpecific_isoI & (numDDIpartners_2plusIso>0);

specific_isoFrac = meanIsoInterFraction(Specific_2plusIso);
nonSpecific_isoFrac = meanIsoInterFraction(nonSpecific_2plusIso);

disp(['Mean fraction of isoform partners for tissue-specific proteins: ' num2str(mean(specific_isoFrac))]);
disp(['Mean fraction of isoform partners for non-tissue-specific proteins: ' num2str(mean(nonSpecific_isoFrac))]);

figure
plotdata = [specific_isoFrac' nonSpecific_isoFrac'];
grp = [zeros(1,length(specific_isoFrac)) ones(1,length(nonSpecific_isoFrac))];
boxplot(plotdata,grp);

plotdata = [mean(nonSpecific_isoFrac); mean(specific_isoFrac)];
sderror = [std(nonSpecific_isoFrac)/sqrt(length(nonSpecific_isoFrac)); std(specific_isoFrac)/sqrt(length(specific_isoFrac))];
figure
hold on
hb = bar(plotdata,0.6,'FaceColor','c','LineWidth',1.5,'EdgeColor','k');
pause(0.1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    for j = 1:length(xData)
        line([xData(j) xData(j)], [plotdata(j,ib)-sderror(j,ib) plotdata(j,ib)+sderror(j,ib)],'Color','k','LineWidth',1.5);
    end
end
set(gca,'XTick',1:2,'XTickLabel',{'Non-tissue-specific\newlineproteins','Tissue-specific\newlineproteins'},'FontSize',17);
ylabel('Mean fraction of isoform partners','FontSize',22);
set(gca,'YGrid','on');
set(gca,'tickDir','out');
box off
clear hb ib j xData plotdata sderror

[h, p] = ttest2(specific_isoFrac, nonSpecific_isoFrac, 'vartype', 'unequal');
disp(['Tissue specificity t-test: p = ' num2str(p)]);

%-------------------------------------------------------------------------
% The old way
cutoff = 0.5;
ts1 = tau_ts(meanIsoInterFraction>=cutoff);
ts2 = tau_ts((meanIsoInterFraction>0) & (meanIsoInterFraction<cutoff));
ts1 = ts1(~isnan(ts1));
ts2 = ts2(~isnan(ts2));

disp(['Mean tissue specificity for proteins interactign with >=' num2str(cutoff*100) '% of partner isoforms: ' num2str(mean(ts1))]);
disp(['Mean tissue specificity for proteins interactign with <' num2str(cutoff*100) '% of partner isoforms: ' num2str(mean(ts2))]);

disp(['Mean fraction of tissue-specific proteins among those interacting with >=' num2str(cutoff*100) '% of partner isoforms: ' num2str(sum(ts1>0.45)/length(ts1))]);
disp(['Mean tfraction of tissue-specific proteins among those interacting with <' num2str(cutoff*100) '% of partner isoforms: ' num2str(sum(ts1>0.45)/length(ts1))]);

x = table([sum(ts1>0.45); sum(ts2>0.45)], ...
    [sum(ts1<=0.45); sum(ts2<=0.45)], ...
    'VariableNames',{'Tissue_specific','Non_tissue_specific'}, ...
    'RowNames',{'>50% isoform partners','<=50% isoform partners'})
[h,p,stats] = fishertest(x);
disp(['Tissue specificity Fisher exact test: p = ' num2str(p)]);

% find overall correlation
sel = ~isnan(tau_ts) & (numDDIpartners_2plusIso>0);
[r,~,~,~] = corrcoef([meanIsoInterFraction(sel) tau_ts(sel)],'rows','pairwise')
        
% itr = 100000;
% bootMeanFrac = bootstrapSim(ts2_isoFrac,ts1_isoFrac,itr);
% p = 2*sum((bootMeanFrac(:,1)-bootMeanFrac(:,2))>=abs(mean(ts2_isoFrac)-mean(ts1_isoFrac)))/size(bootMeanFrac,1);
% disp(['Tissue specificity bootstrap test  (' num2str(itr) ' resamplings): p = ' num2str(p)]);
