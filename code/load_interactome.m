function [I, PPIs, spID, genes] = load_interactome(interactomeFile, spGeneMapFile, spEntrezMapFile, numReported, removeSelfPPIs, whichone)

if strcmpi(whichone,'IntAct')
    if ~isempty(dir('interactome_processed/IntAct_human_PPIs.mat'))
        load interactome_processed/IntAct_human_PPIs.mat
    else
        % read IntAct interactions from file line by line since loading the
        % whole file at once may not be possible due to large size
        numLinesGuess = 500000;
        PPIs = cell(numLinesGuess,2);
        fid1 = fopen(interactomeFile);
        fid2 = fopen('interactome_processed/IntAct_human_PPIs.txt','w');
        fprintf(fid2,'Protein_1_UniProt_ID\tProtein_2_UniProt_ID\tProtein_1_IntAct_ID\tProtein_2_IntAct_ID\tInteraction_ID\n');
        headers = strsplit(fgetl(fid1),'\t');
        i = 0;
        c = 0;
        while ~feof(fid1)
            i = i + 1;
            if mod(i,10000) == 0
                disp([num2str(i) ' interactions read from file ' interactomeFile])
            end
            str = strsplit(fgetl(fid1),'\t');
            if strcmpi(str{36},'false')
                if length(str{1}) > 10 && length(str{2}) > 10
                    if strcmpi(str{1}(1:10),'uniprotkb:') && strcmpi(str{2}(1:10),'uniprotkb:')
                        if isempty(strfind(str{1},'-')) && isempty(strfind(str{2},'-'))
                            if ~isempty(strfind(lower(str{10}),'taxid:9606(human)')) && ~isempty(strfind(lower(str{11}),'taxid:9606(human)'))
                                if ~isempty(strfind(lower(str{21}),'(protein)')) && ~isempty(strfind(lower(str{22}),'(protein)'))
                                    c = c + 1;
                                    PPIs{c,1} = str{1}(11:length(str{1}));
                                    PPIs{c,2} = str{2}(11:length(str{2}));
                                    
                                    IDsplit = strtrim(strsplit(str{3},'|'));
                                    ind = find(strncmpi(IDsplit,'intact:EBI-',11),1);
                                    if ~isempty(ind)
                                        PPIs{c,3} = IDsplit{ind}(8:length(IDsplit{ind}));
                                    else
                                        PPIs{c,3} = '-';
                                    end
                                    
                                    IDsplit = strtrim(strsplit(str{4},'|'));
                                    ind = find(strncmpi(IDsplit,'intact:EBI-',11),1);
                                    if ~isempty(ind)
                                        PPIs{c,4} = IDsplit{ind}(8:length(IDsplit{ind}));
                                    else
                                        PPIs{c,4} = '-';
                                    end
                                    
                                    IDsplit = strtrim(strsplit(str{14},'|'));
                                    ind = find(strncmpi(IDsplit,'intact:EBI-',11),1);
                                    if ~isempty(ind)
                                        PPIs{c,5} = IDsplit{ind}(8:length(IDsplit{ind}));
                                    else
                                        PPIs{c,5} = '-';
                                    end
                                    fprintf(fid2,[PPIs{c,1} '\t' PPIs{c,2} '\t' PPIs{c,3} '\t' PPIs{c,4} '\t' PPIs{c,5} '\n']);
                                end
                            end
                        end
                    end
                end
            end
        end
        fclose(fid1);
        fclose(fid2);
        PPIs = PPIs(1:c,:);
        save('interactome_processed/IntAct_human_PPIs.mat', 'PPIs');
    end
    spID1 = [PPIs(:,1); PPIs(:,2)];
    intactID1 = [PPIs(:,3); PPIs(:,4)];
    [spID, ind, ~] = unique(spID1);
    intactID = intactID1(ind);
    numGenes = length(spID);
    
    % load mapping table from Swiss-Prot IDs to gene names
    s = tdfread(spGeneMapFile);
    spGeneMap = cell(size(s.From,1),2);
    for i = 1:size(s.From,1)
        spGeneMap{i,1} = strtrim(s.From(i,:));
        spGeneMap{i,2} = strtrim(s.To(i,:));
    end
    
    % compile list of interacting protein gene names
    genes = cell(numGenes,1);
    for i = 1:length(spID)
        ind = find(strcmpi(spGeneMap(:,1),spID{i}),1);
        if ~isempty(ind)
            genes{i} = spGeneMap{ind,2};
        else
            genes{i} = '';
        end
    end
    
    disp('Creating protein-protein interaction matrix (This will take a while..)');
    I = uint8(zeros(numGenes,numGenes));
    for i = 1:size(PPIs,1)
        ind1 = strcmpi(spID,PPIs{i,1});
        ind2 = strcmpi(spID,PPIs{i,2});
        I(ind1,ind2) = I(ind1,ind2) + 1;
        I(ind2,ind1) = I(ind1,ind2);
    end
    
    if removeSelfPPIs
        disp('Removing self-interactions');
        I = I.*uint8(~diag(ones(1,numGenes)));
    end
    disp(['Removing interactions reported less than ' num2str(numReported) ' times']);
    I = I>=numReported;    % keep only interactions reported numReported times or more
    
    reviewedHumanProteins = {};
    fid = fopen('uniprot_reviewed_human_proteome.list');
    tline = fgetl(fid);
    while ischar(tline)
        reviewedHumanProteins = [reviewedHumanProteins; strtrim(tline)];
        tline = fgetl(fid);
    end
    fclose(fid);
    reviewed = ismember(spID,reviewedHumanProteins);
    
    % update list of sp Ids and gene names for selected interactions
    inI = (sum(I,2)>0) & reviewed;
    I = I(inI,inI);
    spID = spID(inI);
    genes = genes(inI);
    
    % update list of PPIs to keep
    disp('Updating list of PPIs');
    PPIs = cell(sum(sum(triu(I))),4);
    c = 0;
    for i = 1:size(I,1)
        for j = i:size(I,2)
            if I(i,j)>0
                c = c + 1;
                PPIs(c,:) = {spID{i}, spID{j}, intactID{i}, intactID{j}};
            end
        end
    end
elseif strcmpi(whichone,'HI-II-14')
    HI_II_14 = tdfread(interactomeFile);
    PPIs = [strtrim(mat2cell(HI_II_14.Symbol_A, ones(1,size(HI_II_14.Symbol_A,1)), size(HI_II_14.Symbol_A,2))), ...
            strtrim(mat2cell(HI_II_14.Symbol_B, ones(1,size(HI_II_14.Symbol_B,1)), size(HI_II_14.Symbol_B,2)))];
    PPIentrez = [HI_II_14.Entrez_Gene_IDA HI_II_14.Entrez_Gene_IDB];
    
    % compile the list of genes and their entrez IDs from the interactome
    initgenes = [PPIs(:,1); PPIs(:,2)];
    genes = unique(initgenes);
    numGenes = length(genes);
    initentrez = [PPIentrez(:,1); PPIentrez(:,2)];
    entrezID = zeros(numGenes,1);
    for i = 1:numGenes
        entrezID(i) = initentrez(find(strcmpi(genes{i},initgenes),1));
    end
    
    % create mapping table from Swiss-Prot to Entrez IDs
    s = tdfread(spEntrezMapFile);
    sp2entrez = {};
    for i = 1:size(s.From,1)
        sp2entrez{i,1} = strtrim(s.From(i,:));
        sp2entrez{i,2} = num2str(s.To(i));
    end
    % Create list of Swiss-Prot IDs for all genes in the interactome
    disp('Finding Swiss-Prot IDs for all genes. Genes with multiple Swiss-Prot IDs will not be annotated');
    matches = zeros(numGenes,1);
    spID = cell(numGenes,1);
    for i = 1:numGenes
        ind = find(strcmpi(sp2entrez(:,2),num2str(entrezID(i))));
        if ~isempty(ind)
            matches(i) = length(ind);
            if length(ind) > 1
                spID{i} = '';
            else
                spID{i} = sp2entrez{ind,1};
            end
        else
            spID{i} = '';
        end
    end
    disp([num2str(sum(matches>0)) ' out of ' num2str(i)  ' genes with swiss-prot IDs']);
    disp([num2str(sum(matches==1)) ' out of ' num2str(i)  ' genes with one swiss-prot ID']);
    
    disp('Creating protein-protein interaction matrix (This will take a while..)');
    I = uint8(zeros(numGenes,numGenes));
    for i = 1:size(PPIentrez,1)
        ind1 = (entrezID==PPIentrez(i,1));
        ind2 = (entrezID==PPIentrez(i,2));
        I(ind1,ind2) = I(ind1,ind2) + 1;
        I(ind2,ind1) = 1;
    end
    
    if removeSelfPPIs
        disp('Removing self-interactions');
        I = I.*uint8(~diag(ones(1,numGenes)));
    end
    
    disp(['Removing interactions reported less than ' num2str(numReported) ' times']);
    I = I>=numReported;    % keep only interactions reported numReported times or more
    
    % update list of sp Ids and gene names for selected interactions
    inI = (sum(I,2)>0) & (cellfun(@length,spID)>0);
    I = I(inI,inI);
    spID = spID(inI);
    genes = genes(inI);
    entrezID = entrezID(inI);
    
    % update list of PPIs to keep
    disp('Updating list of PPIs');
    PPIs = cell(sum(sum(triu(I))),2);
    c = 0;
    for i = 1:size(I,1)
        for j = i:size(I,2)
            if I(i,j)>0
                c = c + 1;
                PPIs(c,:) = {spID{i},spID{j}};
            end
        end
    end
end
