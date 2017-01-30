function [I, spID, genes] = load_interactome(interactomeFile, spGeneMapFile, spEntrezMapFile, numReported, removeSelfPPIs, whichone)

if strcmpi(whichone,'IntAct')
    if ~isempty(dir('IntAct_human_PPIs.mat'))
        load IntAct_human_PPIs.mat
    else
        % read IntAct interactions from file line by line since loading the
        % whole file at once may not be possible due to large size
        numLinesGuess = 500000;
        PPIs = cell(numLinesGuess,2);
        fid = fopen(interactomeFile);
        headers = strsplit(fgetl(fid),'\t');
        i = 0;
        c = 0;
        while ~feof(fid)
            i = i + 1;
            if mod(i,10000) == 0
                disp([num2str(i) ' interactions read from file ' interactomeFile])
            end
            str = strsplit(fgetl(fid),'\t');
            if strcmpi(str{36},'false')
                if length(str{1}) > 10 && length(str{2}) > 10
                    if strcmpi(str{1}(1:10),'uniprotkb:') && strcmpi(str{2}(1:10),'uniprotkb:')
                        if isempty(strfind(str{1},'-')) && isempty(strfind(str{2},'-'))
                            if ~isempty(strfind(lower(str{10}),'taxid:9606(human)')) && ~isempty(strfind(lower(str{11}),'taxid:9606(human)'))
                                if ~isempty(strfind(lower(str{21}),'(protein)')) && ~isempty(strfind(lower(str{22}),'(protein)'))
                                    c = c + 1;
                                    PPIs{c,1} = str{1}(11:length(str{1}));
                                    PPIs{c,2} = str{2}(11:length(str{2}));
                                end
                            end
                        end
                    end
                end
            end
        end
        fclose(fid);
        PPIs = PPIs(1:c,:);
        save IntAct_human_PPIs.mat PPIs
    end
    
    spID = unique([PPIs(:,1); PPIs(:,2)]);
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
    
    disp('Creating protein-protein interaction matrix');
    I = uint8(zeros(numGenes,numGenes));
    for i = 1:size(PPIs,1)
        ind1 = find(strcmpi(spID,PPIs{i,1}),1);
        ind2 = find(strcmpi(spID,PPIs{i,2}),1);
        if ~isempty(ind1) && ~isempty(ind2)
            I(ind1,ind2) = I(ind1,ind2) + 1;
            I(ind2,ind1) = I(ind1,ind2);
        end
    end
    if removeSelfPPIs
        I = I.*uint8(~diag(ones(1,numGenes)));
    end
    I = I>=numReported;    % keep only interactions reported numReported times or more
elseif strcmpi(whichone,'HI-II-14')
    [~,HI_II_14] = xlsread(interactomeFile);
    PPIs = HI_II_14(:,[2 4]);
    PPIentrez = HI_II_14(:,[1 3]);
    
    % compile the list of genes and their entrez IDs from the interactome
    initgenes = [PPIs(:,1); PPIs(:,2)];
    genes = unique(initgenes);
    numGenes = length(genes);
    initentrez = str2double([PPIentrez(:,1); PPIentrez(:,2)]);
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
    disp('Finding Swiss-Prot IDs for all genes');
    matches = zeros(numGenes,1);
    spID = cell(numGenes,1);
    for i = 1:numGenes
        ind = find(strcmpi(sp2entrez(:,2),num2str(entrezID(i))));
        if ~isempty(ind)
            matches(i) = length(ind);
            if length(ind) > 1
                disp(['Multiple SwissProt IDs found for gene ' num2str(i) '. Gene will not be annotated.']);
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
    
    disp('Creating protein-protein interaction matrix');
    I = zeros(numGenes,numGenes);
    for i = 1:size(PPIs,1)
        id1 = str2double(PPIentrez{i,1});
        id2 = str2double(PPIentrez{i,2});
        ind1 = find(entrezID==id1,1);
        ind2 = find(entrezID==id2,1);
        I(ind1,ind2) = I(ind1,ind2) + 1;
        I(ind2,ind1) = 1;
    end
    if removeSelfPPIs == 1
        I = I.*~diag(ones(1,numGenes));
    end
    I = I>=numReported;
end

% update list of sp Ids and gene names for selected interactions 
inI = sum(I)>0;
I(~inI,:) = [];
I(:,~inI) = [];
spID = spID(inI);
genes = genes(inI);
