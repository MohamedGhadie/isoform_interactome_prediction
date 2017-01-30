function [domIsoMap, domIsoPos] = load_hmmscan_alt(filename)

numLines = 200000;
domIsoMap = cell(numLines,2);
domIsoPos = zeros(numLines,2);

fid = fopen(filename);
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
        spSplit = strsplit(alignment{4},'|');
        domIsoMap{c,2} = spSplit{2};
        domIsoPos(c,1) = str2double(alignment{20});
        domIsoPos(c,2) = str2double(alignment{21});
    else
        disp(['Skipping line ' num2str(i) '. Unexpected line format.']);
    end
    tline = fgetl(fid);
end
fclose(fid);

domIsoMap = domIsoMap(1:c,:);
domIsoPos = domIsoPos(1:c,:);
