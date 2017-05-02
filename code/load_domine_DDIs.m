function DDIs = load_domine_DDIs(filename)

[num,domineDDIs] = xlsread(filename);
numDDI = 0;
DDIs = {};
for i = 1:size(num,1)
    if num(i,1) == 1 || num(i,2) == 1
        numDDI = numDDI + 1;
        DDIs{numDDI,1} = domineDDIs{i,1};
        DDIs{numDDI,2} = domineDDIs{i,2};
    end
end