load GlobalNames.mat

for i = 2:length(allvarnames) 
    if(~(strcmp(allvarnames{i},'Tperp1')||strcmp(allvarnames{i},'Tpar')||strcmp(allvarnames{i},'Tperp2')))
        fluxpcolorhalfgif(['Zoomed ' allfilenames{i} ' 3'], [allvarnames{i} ''''], 1:192,allzeroincs{i})
        [allvarnames{i} ' done.']
    end
end