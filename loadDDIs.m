function dominePPI = loadDDIs(filename)

[num,dominePPIinit] = xlsread(filename);
numDDI = 0;
dominePPI = {};
for i = 1:size(num,1)
    if num(i,1) == 1 || num(i,2) == 1
        numDDI = numDDI + 1;
        dominePPI{numDDI,1} = dominePPIinit{i,1};
        dominePPI{numDDI,2} = dominePPIinit{i,2};
    end
end