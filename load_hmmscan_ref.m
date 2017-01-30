function [spPFmap, domPos] = load_hmmscan_ref(filename)

numLines = 200000;
spPFmap = cell(numLines,2);
domPos = zeros(numLines,2);

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
        spPFmap{c,1} = alignment{4};
        spPFmap{c,2} = alignment{2}(1:dot-1);
        domPos(c,1) = str2double(alignment{20});
        domPos(c,2) = str2double(alignment{21});
    else
        disp(['Skipping line ' num2str(i) '. Unexpected line format.']);
    end
    tline = fgetl(fid);
end
fclose(fid);

spPFmap = spPFmap(1:c,:);
domPos = domPos(1:c,:);
