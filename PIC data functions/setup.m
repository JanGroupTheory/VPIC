global folderpath
global allfilenames
global allvarnames
global allzeroincs

folderpath = 'D:\IslandCoalescenceTimeSlices\';

timeslices = 1:4:473;

scalfile = {'Psi','Ne','Ppar','Pparp1','Pparp2','Pp1p2','Pperp1','Pperp2','Tpar','Tperp1','Tperp2'};
scalvar = {'Psi','ne','Ppar','Pparp1','Pparp2','Pp1p2','Pperp1','Pperp2','Tpar','Tperp1','Tperp2'};
vecfile = {'B','J','E','Ue'};
vecvar = {'b','j','e','ue'};
symtensfile = {'P'};
symtensvar = {'P'};
coordnames = {'x','y','z'};
ndim = length(coordnames);

pcolorconts = [-90:10:0 7:7:49 53:4:65];

%Create temperature data from pressure data

maketempdata(timeslices,folderpath);

%Create magnitude data for all vectors

makemagdata(timeslices,folderpath,vecfile,vecvar,coordnames);

%Create variable name arrays

scalfilenames = scalfile;
scalvarnames = scalvar;
scalzero = repmat('',1,length(scalvarnames));

vecfilenames = cell(length(vecfile),ndim+1);
vecvarnames = cell(length(vecvar),ndim+1);
%magexpress = cell(length(vecvar));
%magtitles = cell(length(vecvar));
veczero = repmat([repmat({'zero'},1,ndim) {''}],length(vecfile),1);

symtensfilenames = cell(length(symtensfile),ndim*(ndim+1)/2);
symtensvarnames = cell(length(symtensvar),ndim*(ndim+1)/2);
symtenszero = repmat('',length(symtensvar),ndim*(ndim+1)/2);

for j = 1:length(vecfile)
    for i = 1:ndim
        vecfilenames{j,i} = [vecfile{j} coordnames{i}];
        vecvarnames{j,i} = [vecvar{j} coordnames{i}];
    end
    
    vecfilenames{j,ndim+1} = ['mag' vecfile{j}];
    vecvarnames{j,ndim+1} = ['mag' vecvar{j} ];
end


for j = 1:length(symtensfile)
    q = 1;
    for i = 1:ndim
        for k = i:ndim
            symtensfilenames(j,q) = {[symtensfile{j} coordnames{i} coordnames{k}]};
            symtensvarnames(j,q) = {[symtensvar{j} coordnames{i} coordnames{k}]};
            q = q + 1;
        end
    end
end

allfilenames = cat(2,scalfilenames,reshape(vecfilenames,1,[]),reshape(symtensfilenames,1,[]));
allvarnames = cat(2,scalvarnames,reshape(vecvarnames,1,[]),reshape(symtensvarnames,1,[]));
allzeroincs = cat(2,scalzero,reshape(veczero,1,[]),reshape(symtenszero,1,[]));

save('GlobalNames.mat','folderpath','allfilenames','allvarnames','allzeroincs')

%Create precalculated flux contours

if(~exist([folderpath 'pcolorCmat.mat'],'file'))
    makefluxCmat(folderpath, 'pcolorCmat.mat',timeslices,pcolorconts)
end
if(~exist([folderpath 'ContAvgCmat.mat'],'file'))
    load([folderpath 'Psi_' timeslices(1) '.mat']);
    avgconts = findcontspacing(Psi);
    makefluxCmat(folderpath, 'FluxAvgCmat.mat',timeslices,avgconts)
end
if(~exist([folderpath 'FluxAvgCmat.mat'],'file'))
    makefluxCmat(folderpath, 'pcolorCmat.mat',timeslices)
end

%Create maximum and minimum arrays

if(~(exist([folderpath 'sliceminmax.mat'],'file')&& exist([folderpath 'fluxavgminmax.mat'],'file')))

    nslice = length(timeslices);
    nfiles = length(allfilenames);

    avgvalscell = cell(1,nfiles);    
    maxvalcell = cell(1,nfiles);
    minvalcell = cell(1,nfiles);
    maxavgcell = cell(1,nfiles);
    minavgcell = cell(1,nfiles);

    maxvals = struct;
    minvals = struct;
    maxavg = struct;
    minavg = struct;

    load FluxAvgCmat.mat

    parfor j = 1:nfiles

        avgvals = struct;
        avgmax = zeros(1,nslice);
        avgmin = zeros(1,nslice);
        filmax = zeros(1,nslice);
        filmin = zeros(1,nslice);
        
        q = 0 ;

        tic

        for is = timeslices

            q = q + 1;

            slicestr = ['_' num2str(is) '.mat'];

            filpath = [folderpath allfilenames{j} slicestr];

            Cmat = PsiCmat.(['Cmat' num2str(is)]);

            if(exist(filpath,'file') ~= 0)
                fil = load(filpath);
                filvals = fil.(allvarnames{j});
                avg = fluxsurfavg(filvals',Cmat);
                avgvals.(['slice_' num2str(is)]) = avg;
                filmax(q) = maxall(filvals);
                filmin(q) = minall(filvals);
                avgmax(q) = maxall(avg);
                avgmin(q) = minall(avg);
            else
                ['Could not find ' allvarnames{j} slicestr]
                continue
            end


            if(mod(q,20) == 19)
               [num2str(q) ': ' allvarnames{j} ' ' num2str(round(100*q/nslice)) '% done.']
               toc
            end
        end

        avgvalscell{j} = avgvals;
        maxavgcell{j} = avgmax;
        minavgcell{j} = avgmin;
        maxvalcell{j} = filmax;
        minvalcell{j} = filmin;
    end

    for j = 1:nfiles
        maxvals.(allvarnames{j}) = maxvalcell{j};
        minvals.(allvarnames{j}) = minvalcell{j};
        maxavg.(allvarnames{j}) = maxavgcell{j};
        minavg.(allvarnames{j}) = minavgcell{j};
        eval([allvarnames{j} '_avg = avgvalscell{j};']);
        save([folderpath allvarnames{j} '_avgs.mat'], [allvarnames{j} '_avg']);
    end

    save([folderpath 'fluxavgminmax.mat'],'maxavg','minavg')
    save([folderpath 'sliceminmax.mat'], 'maxvals', 'minvals')

end
