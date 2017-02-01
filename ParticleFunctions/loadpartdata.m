function loadpartdata(directory,tslice)

cd(directory);
loopvar = true;

while loopvar
    dirstruct = dir('data*');
    szvar = size(dirstruct);
    if szvar(1) ~= 0
        loopvar = false;
        datadir = dirstruct.name;
        datadir = [pwd '\' datadir];
    end
    cd ..;
end
cd(directory);
fID = fopen([directory '\info']);
masterstring = fscanf(fID,'%s');
stringcell = strsplit(masterstring,'=');
lookfors = [1;1;1;1;1;1];
for i = 1:(length(stringcell)-1)
    if lookfors(1)
    IIs = strfind(stringcell{i},'wpe/wce');
        if ~isempty(IIs);
            temp = regexp(stringcell{i+1},...
                '[abcdfghijklmnopqrstuvwxyzQWERTYUIOPASDFGHJKLZXCVBNM]', 'split');
            wpewce = str2double(temp{1});
            lookfors(1) = 0;
        end
    end
    if lookfors(2)
        IIs = strfind(stringcell{i},'mi/me');
        if ~isempty(IIs);
            temp = regexp(stringcell{i+1},...
                '[abcdfghijklmnopqrstuvwxyzQWERTYUIOPASDFGHJKLZXCVBNM]', 'split');
            mime = str2double(temp{1});
            lookfors(2) = 0;
        end
    end
    if lookfors(3)
        IIs = strfind(stringcell{i},'Lx/de');
        if ~isempty(IIs);
            temp = regexp(stringcell{i+1},...
                '[abcdfghijklmnopqrstuvwxyzQWERTYUIOPASDFGHJKLZXCVBNM]', 'split');
            Lx = str2double(temp{1});
            lookfors(3) = 0;
        end
    end
    if lookfors(4)
        IIs = strfind(stringcell{i},'Lz/de');
        if ~isempty(IIs);
            temp = regexp(stringcell{i+1},...
                '[abcdfghijklmnopqrstuvwxyzQWERTYUIOPASDFGHJKLZXCVBNM]', 'split');
            Lz = str2double(temp{1});
            lookfors(4) = 0;
        end
    end
    if lookfors(5)
        IIs = strfind(stringcell{i},'nx');
        if ~isempty(IIs);
            temp = regexp(stringcell{i+1},...
                '[abcdfghijklmnopqrstuvwxyzQWERTYUIOPASDFGHJKLZXCVBNM]', 'split');
            nx = str2double(temp{1});
            lookfors(5) = 0;
        end
    end
    if lookfors(6)
        IIs = strfind(stringcell{i},'nz');
        if ~isempty(IIs);
            temp = regexp(stringcell{i+1},...
                '[abcdfghijklmnopqrstuvwxyzQWERTYUIOPASDFGHJKLZXCVBNM]', 'split');
            nz = str2double(temp{1});
            lookfors(6) = 0;
        end   
    end
end
fclose(fID);
% temp = strsplit(stringcell{5},'m');
% wpewce = str2double(temp{1});
% temp = strsplit(stringcell{6},'t');
% mime = str2double(temp{1});
% temp = strsplit(stringcell{10},'L');
% Lx = str2double(temp{1});
% temp = strsplit(stringcell{12},'L');
% Lz = str2double(temp{1});
% temp = strsplit(stringcell{16},'n');
% nx = str2double(temp{1});
% temp = strsplit(stringcell{18},'d');
% nz = str2double(temp{1});

IIx = round(nx/6):round(nx/3); IIz = round(nz/6):round(nz/3);

J0 = 1/wpewce/sqrt(mime);
% tslice =239987; %25th output
% tslice = 25;

% qs = {'ne','uex','uey','uez','uiy','uix','uiz',...
%     'Pexx','Pexy','Pexz','Peyy','Peyz','Pezz','bx','by','bz',...
%     'Pixx','Pixy','Pixz','Piyy','Piyz','Pizz','emix','ex','ey','ez'    };
qs = {'ne','bx','by','bz','ex','ey','ez'};


ne = zeros(nx/2,nz/2);

for ii = 1:length(qs)
    temp = loadslice(qs{ii},tslice,nx,nz,datadir);
    com = [qs{ii} ' = temp'';'];
    eval(com)
end


xv = linspace(0,Lx,nx/2);
zv = linspace(0,Lz,nz/2)-Lz/2;

%make Ay and find X point
Ay = zeros(size(bx));
Ay0 = cumsum(bx(1,:));
Ay = cumsum(bz,1) - repmat(Ay0,[size(bx,1) 1]);
Ay = Ay*(xv(2)-xv(1));
Psi =Ay;

    bxr=bx(IIx,IIz);
    bzr=bz(IIx,IIz);
    
 bzr=sm(bzr,[2,2]);
  bxr=sm(bxr,[2,2]);

Cx=contourc(bxr',[0 0]);
Cz=contourc(bzr',[0 0]);

bzc = interp2(bzr,Cx(2,:),Cx(1,:));

[~,II] = min(abs(bzc));

Xpoint = round([Cx(1,II), Cx(2,II)]+[IIx(1) IIz(1)]-1);

ix = Xpoint(1);
iz = Xpoint(2);

PsiX=Psi(ix,iz);

%

 B = sqrt(bx.^2+by.^2+bz.^2);
% Pepar = (bx.^2.*Pexx + by.^2.*Peyy + bz.^2.*Pezz...
%     +2*bx.*by.*Pexy + 2*bx.*bz.*Pexz + 2*by.*bz.*Peyz)./B.^2;
% Peperp = 0.5*(Pexx+Peyy+Pezz-Pepar);
% 
% Pipar = (bx.^2.*Pixx + by.^2.*Piyy + bz.^2.*Pizz...
%     +2*bx.*by.*Pixy + 2*bx.*bz.*Pixz + 2*by.*bz.*Piyz)./B.^2;
% Piperp = 0.5*(Pixx+Piyy+Pizz-Pipar);

varlist = {'nx';'nz';'xv';'zv';'tslice';'bx';'by';'bz';'mime';'Lx';'Lz'...
    ;'J0';'wpewce'};

save('partsetup.mat',varlist{:});

