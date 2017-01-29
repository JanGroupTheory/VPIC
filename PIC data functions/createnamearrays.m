    scalfile = {'Psi','Ne','Ppar','Pparp1','Pparp2','Pp1p2','Pperp1','Pperp2','Tpar','Tperp1','Tperp2'};
    scalvar = {'Psi','ne','Ppar','Pparp1','Pparp2','Pp1p2','Pperp1','Pperp2','Tpar','Tperp1','Tperp2'};
    vecfile = {'B','J','E','Ue'};
    vecvar = {'b','j','e','ue'};
    symtensfile = {'P'};
    symtensvar = {'P'};
    coordnames = {'x','y','z'};
    ndim = length(coordnames);
    
    scalfilenames = scalfile;
    scalvarnames = scalvar;
    
    vecfilenames = cell(length(vecfile),ndim);
    vecvarnames = cell(length(vecvar),ndim);
    magexpress = cell(length(vecvar));
    magtitles = cell(length(vecvar));
    veczero = repmat({'zero'},length(vecvar),ndim);
    
    symtensfilenames = cell(length(symtensfile),ndim*(ndim+1)/2);
    symtensvarnames = cell(length(symtensvar),ndim*(ndim+1)/2);
    
    for j = 1:length(vecfile)
        temp = 'sqrt(0';
        for i = 1:ndim
            vecfilenames(j,i) = {[vecfile{j} coordnames{i}]};
            vecvarnames(j,i) = {[vecvar{j} coordnames{i}]};
            temp = [temp '+' vecvar{j} coordnames{i} '.^2'];
        end
        magexpress(j) = {[temp ')']};
        magtitles(j) = {['mag ' vecfile{j}]};
    end
    
    q = 1;
    
    for j = 1:length(symtensfile)
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