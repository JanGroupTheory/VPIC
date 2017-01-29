function [] = makefluxCmat(folderpath, filename, timeslices, varargin)
    
    PsiCmat = struct;
    CmatHold = cell(size(timeslices));
        
    
    parfor i = 1:length(timeslices)
        
%         if(mod(slicenum, 10) == 9)
%             slicenum + 1
%             toc
%         end
        
        slicestr = ['_' num2str(timeslices(i)) '.mat'];

        temp = [folderpath 'Psi' slicestr];
        if(exist(temp,'file') ~= 0)
            ps = load(temp);
            Psi = ps.Psi;
            %['Loaded ' char(strcat(sname,'.'))]
        else
            ['Could not find Psi' slicestr '.'] 
            continue
        end
        
        if(~isempty(varargin))
            contvals = varargin{1};
        else
            contvals = findcontspacing(Psi');
        end
        
        [CmatHold{i}, ~] = contour(Psi',contvals);
        
        ['Slice ' num2str(timeslices(i)) ' done.']

    end
    
    for i = 1:length(timeslices)
        PsiCmat.(['Cmat' num2str(slicenum(i))]) = CmatHold{i};
    end
    
    save([folderpath filename],'PsiCmat')
        
        