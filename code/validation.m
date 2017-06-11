% This script validates the performance of the domain-based isoform
% interaction prediction method on the experimental dataset of Yang et al.
% (2016). Interactions between orfeome proteins and newly cloned reference
% proteins are used as the reference interactome and domain-domain
% interactions are mapped onto these reference interactions. Interactions
% for the newly cloned alternative isoforms are then predicted from the
% DDI-annotated reference interactions using two methods:

% - First method uses reference interactions having full DDI annotations only. 
% A reference interaction has a full DDI annotation if the newly cloned 
% protein contains an interacting domain of a DDI and its orfeome interaction 
% partner contains the other interacting domain of the same DDI. An alternative 
% isoform of the newly cloned reference protein is then predicted to lose 
% the interaction if it loses all DDI annotations with the orfeome partner, 
% otherwise the interaction is retained.

% - Second method uses reference interactions that have either a full DDI 
% annotation or a partial DDI annotation. A reference interaction has a 
% partial DDI annotation if the newly cloned protein contains an interacting 
% domain of a DDI whereas its orfeome interaction partner does not contain 
% the other interacting domain of the DDI. An alternative isoform of the 
% newly cloned reference protein is then predicted to lose the interaction 
% if it loses all of its interacting domains, otherwise the interaction is retained.

% load orfeome sequences
orfData = fastaread('horf81_cloneInfo20120427.fa.txt');

% create list of ORF IDs
orfID = cell(size(orfData,1),1);
for i = 1:size(orfData,1)
    str = strsplit(orfData(i).Header,'|');
    if strcmpi(str{2}(1:8),'ORF_ID: ') && strcmpi(str{5}(1:9),'GENE_ID: ')
        orfID{i} = ['ORF' str{2}(9:length(str{2})) '_GENE' str{5}(10:length(str{5}))];
    else
        disp(['ORF id not readbale for sequence # ' num2str(i)]);
    end
end

%-------------------------------------------------------------------------
% load domain mappings for newly cloned reference isoforms

disp('Loading domain mappings for newly cloned reference isoforms');
numLines =  2154;   % e-value = 1e-3
domCanMap = cell(numLines,2);
domCanPos = zeros(numLines,2);

fid = fopen('hmmscan_ref_output_e10-3.domtab');
tline = fgetl(fid);
c = 0;
i = 0;
while ischar(tline)
    i = i + 1;
    alignment = strsplit(tline);
    dot = find(alignment{2}=='.',1);
    if strcmp(alignment{2}(1:2),'PF') && ~isempty(dot) ...
            && sum(isstrprop(alignment{20},'digit'))==length(alignment{20}) ...
            && sum(isstrprop(alignment{21},'digit'))==length(alignment{21})
        c = c + 1;
        domCanMap{c,1} = alignment{2}(1:dot-1);
        domCanMap{c,2} = alignment{4};
        domCanPos(c,1) = str2double(alignment{20});
        domCanPos(c,2) = str2double(alignment{21});
    else
        disp(['Skipping line ' num2str(i) '. Unexpected line format.']);
        disp(alignment);
    end
    tline = fgetl(fid);
end
fclose(fid);
clear fid tline alignment dot numLines
domCanMap = domCanMap(1:c,:);
domCanPos = domCanPos(1:c,:);
spPFmap = fliplr(domCanMap);
domPos = domCanPos;
clear domCanMap domCanPos

%-------------------------------------------------------------------------
% load domain mappings for newly cloned alternative isoforms

disp('Loading domain mappings for newly cloned alternative isoforms');
numLines = 2560;    % e-value = 1e-3
domIsoMap = cell(numLines,2);
domIsoPos = zeros(numLines,2);

fid = fopen('hmmscan_alt_output_e10-3.domtab');
tline = fgetl(fid);
c = 0;
i = 0;
while ischar(tline)
    i = i + 1;
    alignment = strsplit(tline);
    dot = find(alignment{2}=='.',1);
    if strcmp(alignment{2}(1:2),'PF') && ~isempty(dot) ...
            && sum(isstrprop(alignment{20},'digit'))==length(alignment{20}) ...
            && sum(isstrprop(alignment{21},'digit'))==length(alignment{21})
        c = c + 1;
        domIsoMap{c,1} = alignment{2}(1:dot-1);
        domIsoMap{c,2} = alignment{4};
        domIsoPos(c,1) = str2double(alignment{20});
        domIsoPos(c,2) = str2double(alignment{21});
    else
        disp(['Skipping line ' num2str(i) '. Unexpected line format.']);
        disp(alignment);
    end
    tline = fgetl(fid);
end
fclose(fid);
clear fid tline alignment dot numLines
domIsoMap = domIsoMap(1:c,:);
domIsoPos = domIsoPos(1:c,:);

%-------------------------------------------------------------------------
% load domain mappings for orfeome proteins

disp('Loading domain mappings for orfeome proteins');
numLines = 51380;   % e-value = 1e-3
orfDomMap = cell(numLines,2);
orfDomPos = zeros(numLines,2);

% load domain-protein mapping for ORF sequences
fid = fopen('hmmscan_orf_output_e10-3.domtab');
tline = fgetl(fid);
c = 0;
i = 0;
while ischar(tline)
    i = i + 1;
    alignment = strsplit(tline);
    dot = find(alignment{2}=='.',1);
    if strcmp(alignment{2}(1:2),'PF') && ~isempty(dot) ...
            && sum(isstrprop(alignment{20},'digit'))==length(alignment{20}) ...
            && sum(isstrprop(alignment{21},'digit'))==length(alignment{21})
        c = c + 1;
        orfDomMap{c,1} = alignment{4};
        orfDomMap{c,2} = alignment{2}(1:dot-1);
        orfDomPos(c,1) = str2double(alignment{20});
        orfDomPos(c,2) = str2double(alignment{21});
    else
        disp(['Skipping line ' num2str(i) '. Unexpected line format.']);
        disp(alignment);
    end
    tline = fgetl(fid);
end
fclose(fid);
clear fid tline alignment dot numLines
orfDomMap = orfDomMap(1:c,:);
orfDomPos = orfDomPos(1:c,:);

%-------------------------------------------------------------------------
% load positive and negative PPIs from experimental dataset

disp('Loading experimental dataset');
tdfread('ts2b.tsv');
Isoform_ID = strtrim(mat2cell(Isoform_ID,ones(1,size(Isoform_ID,1)),size(Isoform_ID,2)));
AD_clone_ID = strtrim(mat2cell(AD_clone_ID,ones(1,size(AD_clone_ID,1)),size(AD_clone_ID,2)));
AD_GeneID = strtrim(mat2cell(AD_GeneID,ones(1,size(AD_GeneID,1)),size(AD_GeneID,2)));
Interaction_Found = strtrim(mat2cell(Interaction_Found,ones(1,size(Interaction_Found,1)),size(Interaction_Found,2)));
PPI = cell(length(AD_clone_ID),1);
for i = 1:length(AD_clone_ID)
    orfind = strfind(AD_clone_ID{i},'_ORF');
    PPI{i} = [AD_clone_ID{i}(orfind+1:length(AD_clone_ID{i})) '_GENE' AD_GeneID{i}];
end
PPI = [Isoform_ID PPI Interaction_Found];
PPI(strcmpi(PPI(:,3),'N/A'),:) = [];
PPI(strcmpi(PPI(:,1),PPI(:,2)),:) = [];
clear Entrez_Gene_ID Gene_Symbol Isoform_ID Category AD_clone_ID AD_GeneID AD_symbol Interaction_Found

%-------------------------------------------------------------------------
% create the reference interactome by selecting positive interactions that
% involve reference isoforms.

PisoI = {};
for i = size(PPI,1):-1:1
    len1 = length(PPI{i,1});
    if ~strcmpi(PPI{i,1}(len1-1:len1),'_1')
        PisoI = [{PPI{i,1},PPI{i,2},PPI{i,3}}; PisoI];
        PPI(i,:) = [];
    end
end

%-------------------------------------------------------------------------
% load domain-domain interactions

domineFile = 'domine_interactions.xlsx';
did3File = '3did_flat.txt';
if (exist(domineFile, 'file') == 2) && (exist(did3File, 'file') == 2)
    disp('Loading domain-domain interactions');
    DDIs = load_3did_and_domine_DDIs(did3File, domineFile);
else
    disp('Domain-domain interaction file not found. Exiting script');
    return
end
numDDI = size(DDIs,1);
intDomains = unique([DDIs(:,1); DDIs(:,2)]);

%-------------------------------------------------------------------------
% annotate reference interactome with full domain-domain interactions

disp('Annotating reference interactome with full domain-domain interactions');
PPI_interDomains = cell(size(PPI,1),1);
numFullDDImap = NaN(size(PPI,1),1);
for i = 1:size(PPI,1)
    if strcmpi(PPI{i,3},'Positive')
        ind1 = find(strcmpi(spPFmap(:,1),PPI{i,1}));
        dash = strfind(PPI{i,2},'_GENE');
        orfpart = PPI{i,2}(1:dash);
        ind2 = find(strncmpi(orfDomMap(:,1),orfpart,length(orfpart)));
        if ~isempty(ind1) && ~isempty(ind2)
            domList1 = unique(spPFmap(ind1,2));
            domList2 = unique(orfDomMap(ind2,2));
            ddi1 = ismember(DDIs(:,1),domList1) & ismember(DDIs(:,2),domList2);
            PPI_interDomains{i} = [PPI_interDomains{i}; DDIs(ddi1,:)];
            ddi2 = ismember(DDIs(:,1),domList2) & ismember(DDIs(:,2),domList1);
            ddi2(strcmpi(DDIs(:,1),DDIs(:,2))) = 0;
            PPI_interDomains{i} = [PPI_interDomains{i}; fliplr(DDIs(ddi2,:))];
            numFullDDImap(i) = size(PPI_interDomains{i},1);
        end
    end
end

disp([num2str(sum(numFullDDImap>0)) ' reference interactions with full DDI mapping']);

%-------------------------------------------------------------------------
% predict isoform loss of interaction using interactions with full DDI
% annotations only

disp('Predicting isoform loss of interaction using interactions with full DDI annotations only');
c = 0;
fullDDImap_prediction = (-1)*ones(size(PisoI,1),1);
for i = 1:size(PPI,1)
    if strcmpi(PPI{i,3},'positive');
        if ~isempty(PPI_interDomains{i})
            c = c + 1;
            isoind = find(strncmpi(PisoI(:,1),PPI{i,1},length(PPI{i,1})-1) & strcmpi(PisoI(:,2),PPI{i,2}));
            if ~isempty(isoind)
                domList1 = PPI_interDomains{i}(:,1);
                for j = 1:length(isoind)
                    isoDom = unique(domIsoMap(strcmpi(domIsoMap(:,2),PisoI{isoind(j),1}),1));
                    if isempty(find(ismember(domList1,isoDom),1))
                        fullDDImap_prediction(isoind(j)) = 0;
                    else
                        fullDDImap_prediction(isoind(j)) = 1;
                    end
                end
            end
        end
    end
end
disp(['Predicted ' num2str(sum(fullDDImap_prediction>-1)) ' isoform interactions from ' ...
    num2str(c) ' reference interactions with full DDI mapping']);

%-------------------------------------------------------------------------
% predict an isoform loss of interaction using interactions with full or
% partial DDI annotations

disp('Predicting isoform loss of interaction using interactions with full or partial DDI annotations');
c = 0;
partialOrFullDDImap_prediction = (-1)*ones(size(PisoI,1),1);
for i = 1:size(PPI,1)
    if strcmpi(PPI{i,3},'positive');
        ind1 = find(strcmpi(spPFmap(:,1),PPI{i,1}));
        domList1 = unique(spPFmap(ind1,2));
        p1IntDom = intDomains(ismember(intDomains,domList1));
        if ~isempty(p1IntDom)
            c = c + 1;
            isoind = find(strncmpi(PisoI(:,1),PPI{i,1},length(PPI{i,1})-1) & strcmpi(PisoI(:,2),PPI{i,2}));
            if ~isempty(isoind)
                for j = 1:length(isoind)
                    isoDom = unique(domIsoMap(strcmpi(domIsoMap(:,2),PisoI{isoind(j),1}),1));
                    if isempty(find((ismember(p1IntDom,isoDom)),1))
                        partialOrFullDDImap_prediction(isoind(j)) = 0;
                    else
                        partialOrFullDDImap_prediction(isoind(j)) = 1;
                    end
                end
            end
        end
    end
end
disp(['Predicted ' num2str(sum(partialOrFullDDImap_prediction>-1)) ' isoform interactions from ' ...
    num2str(c) ' reference interactions with full or partial DDI mapping']);

%-------------------------------------------------------------------------
% create contingency table of predictions from reference interactions with
% full DDI mappings
TP = sum(fullDDImap_prediction==1 & strcmpi(PisoI(:,3),'positive'));
FP = sum(fullDDImap_prediction==1 & strcmpi(PisoI(:,3),'negative'));
TN = sum(fullDDImap_prediction==0 & strcmpi(PisoI(:,3),'negative'));
FN = sum(fullDDImap_prediction==0 & strcmpi(PisoI(:,3),'positive'));

disp('Prediction results using reference interactions with full DDI mappings:');
x = table([TP; FP],[FN; TN], ...
    'VariableNames',{'pred_pos','pred_neg'}, ...
    'RowNames',{'actual_pos','actual_neg'})

[h,p,~] = fishertest(x);
disp(['Fisher''s exact test p-value = ' num2str(p)]);

%-------------------------------------------------------------------------
% create contingency table of predictions from reference interactions with
% full or partial DDI mappings
TP = sum(partialOrFullDDImap_prediction==1 & strcmpi(PisoI(:,3),'positive'));
FP = sum(partialOrFullDDImap_prediction==1 & strcmpi(PisoI(:,3),'negative'));
TN = sum(partialOrFullDDImap_prediction==0 & strcmpi(PisoI(:,3),'negative'));
FN = sum(partialOrFullDDImap_prediction==0 & strcmpi(PisoI(:,3),'positive'));

disp('Prediction results using reference interactions with full or partial DDI mappings:');
x = table([TP; FP],[FN; TN], ...
    'VariableNames',{'pred_pos','pred_neg'}, ...
    'RowNames',{'actual_pos','actual_neg'})

[h,p,~] = fishertest(x);
disp(['Fisher''s exact test p-value = ' num2str(p)]);
