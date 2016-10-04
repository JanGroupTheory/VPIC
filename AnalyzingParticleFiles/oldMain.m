% readPartcvaparvparFINALUPDATE
%
% Reading in the particle data
%
% Requires:
% ZSortCreate, ZSortPtsDetFn, ZSortPtsSort, findXpts, findZpts,
% loadVarsloaddir, and readFunctionPumpX2
%
% Optional: varycolor (for making nice plots)


tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializing the filepaths
folderName = '17';
loaddir = ['/ext/VPICRuns/AlfvenWavesInBox/T70ObliqPropAmpUp/'];
savedir = ['/ext/VPICRuns/AlfvenWavesInBox/T70ObliqPropAmpUp/particlesDist'];


% Variables determined from the simulation (can be found in the info file)
[fileName, nlx, nlz, Lx, Lz, mime, Npart, Vlim, Elim, maxZAbs, ...
    minZAbs, Xbloc, Zbloc, maxXAbs, minXAbs, delT] = loadVarsloaddir( loaddir, folderName );

% Computed quantities
xv=(1:nlx)*Lx/nlx-Lx/2;                                                     % x-coordinate when using readFunctionPump
zv=(1:nlz)*Lz/nlz-Lz/2;                                                     % z-coordinate when using readFunctionPump
dx = Lx/nlx;
dz = Lz/nlz;

% Parameters that are used to delineate particle regions
% When modifying these, you must also go into the functions mentioned in
% the first comment and change the size of the struct - as well as the
% sorting portion of the main code
Ndiv = 5;                                                                   % Number of divisions on the z direction
Ndivk = 32;                                                                 % Number of divisions in the x direction

% Internal variables used in this program
nFiles = size(fileName, 2);
Nx=nlz;
Nvpar=201;
Nvperp=Nvpar/2 - 0.5;
NEkin = Nvpar+Nvperp;
Xlim1 = -480;
Xlim2 = 0;

Zlim1 = -Lz/2;
Zlim2 = Lz/2;

clear Nbound1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading in the particle data and storing it in .mat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for i=66:80
    i
    % Going to the correct folder
    eval (['cd ' loaddir]);                                                 % go to the right directory
    fileinp=['electrons-' fileName{i} '.bin'];
    
    % Making the necessary directories
    if exist(savedir, 'dir') == 0
        eval(['mkdir ' savedir ]);
        eval(['mkdir ' savedir '/Run' folderName 'T' fileName{i}]);
    elseif exist( [savedir '/Run' folderName 'T' fileName{i}], 'dir') == 0
        eval(['mkdir ' savedir '/Run' folderName 'T' fileName{i}]);
    end
    
    
    % Loading the particle data from the electrons file
    f = fopen(fileinp);
    dat=zeros(0,[5,3.4e5]);
    Nskip=2e3;                                                              % Take every 40,000th particle for size reasons
    
    X=[];
    Z=[];
    
    for k=1:round(Npart/Nskip)
        fseek(f,(Nskip*k)*5*4,'bof');                                       % There are 5 elements in every particle's data - x, z, vx, vy, vz
        dum=fread(f,5,'single');
        X=[X dum(1)];                                                       % Read in the x coordinate
        Z=[Z dum(2)];                                                       % Read in the z coordinate
    end

    % Finding out the region boundaries of the particles
    IIX=round(X/dx+0.5);                                                    % calculate x-index
    diffII=diff(IIX);                                                       % calculating the difference between each x-index value for all the chosen particles
    diffIISort = sort(diffII);
    
    nZLevels = ((Zlim2 - Zlim1)/Zbloc/dz);
    nXLevels = round(abs(maxXAbs - minXAbs)/Xbloc/dx);
    Nbound0 = zeros(1, (nXLevels-1)*nZLevels);
    zBnds = ones(1, nZLevels + 1);
    lenX = length(X);
    Nbound1 = zeros(1, length(Nbound0) + length(zBnds)-1);
    
    k1 = 1;
    while (k1<(nZLevels))
        tempdiffIIVar = find(diffII == diffIISort(k1));
        
        if (numel(tempdiffIIVar) == 1)
            zBnds(k1+1) = tempdiffIIVar;
            k1 = k1 + 1;
        elseif (numel(tempdiffIIVar) == 2)
            zBnds(k1+1) = tempdiffIIVar(1);
            zBnds(k1+2) = tempdiffIIVar(2);
            k1 = k1+2;
        elseif (numel(tempdiffIIVar) == 3 && k1 == 1)
            zBnds(k1+1) = tempdiffIIVar(1);
            zBnds(k1+2) = tempdiffIIVar(2);
            zBnds(k1+3) = tempdiffIIVar(3);
            k1 = k1+3; 
        else
            disp('Something is wrong with the bound finding part -line 116')
            pause()
        end
        
    end
    zBnds(end) = length(diffII);
    zBnds = sort(zBnds);
    
    for k1 = 1:nZLevels
        Xp = X(zBnds(k1):(zBnds(k1+1)-1));
        for k=1:nXLevels - 1
            
            temp0 = (find((Xp - (Xbloc*dx*k + minXAbs)) >0));
            temp0Sort = sort(temp0);
            if temp0Sort(1) == 1
                temp1 = temp0Sort(2);
            else
                temp1 = temp0Sort(1);
            end
            if (isempty(temp1) ~= 1)
                Nbound0((k1-1)*(nXLevels-1) + k) = temp1+zBnds(k1)-1;
            else
                disp([num2str(k) ' -  first set - cannot find this point'])
            end
        end
    end
    
    Nbound1 = [Nbound0 zBnds(2:(end))];
    Nbound1 = sort(Nbound1);
    
    %     tempLengthNb1 = length(Nbound1); %comment this out when add in halfPt vect
    halfPt = min(find(abs(diffII) == max(abs(diffII))));
    %
    %     Xp = X((halfPt+1):end);
    %
    %     for k=1:round(abs(maxXAbs - minXAbs)/Xbloc/dx) - 1
    %         temp1 = min(find((Xp - (Xbloc*dx*k + minXAbs)) >0));
    %         if (isempty(temp1) ~= 1)
    %             Nbound1(k+tempLengthNb1) = temp1 + (halfPt+1);
    %         else
    %             disp([num2str(k) ' -  second set - cannot find this point'])
    %         end
    %     end
    
    %Getting rid of divisions that are too close
    %     diffNb1 = diff(Nbound1);
    %     tooClosePts = find(abs(diffNb1) <= 10);
    %     if(isempty(tooClosePts) ~= 1)
    %
    %         Nbound1New = [];x1 = -3:.2:3; x2 = -3:.2:3;
    % [X1,X2] = meshgrid(x1,x2);
    % F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    % F = reshape(F,length(x2),length(x1));
    % surf(x1,x2,F);
    % caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
    % axis([-3 3 -3 3 0 .4])
    % xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
    %         oldpt1 = 1;
    %         for k=1:length(tooClosePts)
    %             pt1 = tooClosePts(k);
    %             pt2 = pt1 + 1;
    %             newpt = (Nbound1(pt1) + Nbound1(pt2))/2;
    %             Nbound1(pt1) = newpt;
    %             Nbound1New = [ Nbound1New Nbound1(oldpt1:pt1)];
    %             oldpt1 = pt2 + 1;
    %         end
    %
    %         Nbound1 = Nbound1New;
    %     end
    
    fclose(f);
    
    %Finding the particle numbers at the boundaries
    f = fopen(fileinp);
    kMax = length(Nbound1);
    Nbound1b = [];
    Nnum = kMax - 1;
    NBdk = max(find(min(abs(Nbound1 - halfPt)) == abs(Nbound1 - halfPt))-1);  %CHANGE THIS USING HALFPT VECT
    for k=[1:(NBdk-1) (NBdk+1):Nnum]
        fseek(f,Nskip*(Nbound1(k)-2)*5*4,'bof');%fseek(f,Nskip*(Nbound1(k)-2)*5*4,'bof')
        dum=fread(f,4*Nskip*6,'single');
        if k <= (NBdk-1)
            IIX=find(dum(1:5:end)>((k-kMax)*(((maxXAbs - minXAbs)/kMax))));
        else
            IIX=find(dum(1:5:end)>(((k-NBdk)-kMax)*((maxXAbs - minXAbs)/kMax)));
        end
        Nbound1b(k)=(Nbound1(k)-2)*Nskip+IIX(1);
    end
    if NBdk > 2
        Nbound1b(NBdk) = Nbound1b(NBdk -1) + round(0.95*(Nbound1b(NBdk -1)-Nbound1b(NBdk -2)));
    elseif NBdk == 1
        Nbound1b(NBdk) = Nbound1b(NBdk +1) - round((Nbound1b(NBdk +2) - Nbound1b(NBdk +1)));
    elseif NBdk == 0 
        
    else
        Nbound1b(NBdk) = round(0.5*(Nbound1b(NBdk -1) + Nbound1b(NBdk +1)));
    end
    Nbound1b(Nnum + 1) = Nbound1b(Nnum) + round(0.95*(Nbound1b(Nnum)-Nbound1b(Nnum -1)));
    
    fclose(f);
    
    
    %Initializing the variables
    Nbound1c=[0 Nbound1b];                                                  % Adding 0 to the
    f = fopen(fileinp);
    
    %Reading in the magnetic fields for the particular time slices
    
    tslice = str2num(fileName{i})/delT + 1;
    bx=readFunctionBlueSlice7(tslice, 'Bx', [loaddir 'results/tslice/']);
    by=readFunctionBlueSlice7(tslice, 'By', [loaddir 'results/tslice/']);
    bz=readFunctionBlueSlice7(tslice, 'Bz', [loaddir 'results/tslice/']);
    
    %Creating the structs to sort the particles
    ZSortPts=ZSortCreate( );
    ZSortPts2=ZSortCreate( );
    ZSortPts3=ZSortCreate( );
    
    kcount = 1;                                                             % Keeps track of the k num for increased box size
    
    %Optional Plotting Variables for N89
    NumberOfPlots=2;
    ColorSet=varycolor(NumberOfPlots);
    N89 = zeros(2, round(3*length(Nbound1c)/8)*Ndivk);                      % total number of particles saved - diagnostic tool
    
    kstart = 1;%288;%round(length(Nbound1c)/8);%349;%round(5*length(Nbound1c)/32);
    kend = round(length(Nbound1c)); %864;%round(3*length(Nbound1c)/8);%474;%round(7*length(Nbound1c)/32);
    
    %Choose kstart%3 == 1 so that the sorting works properly
    if mod(kstart, 3) == 2
        kstart = kstart - 1;
    elseif mod(kstart,3) == 0;
        kstart = kstart - 2;
    end
    
    %     %Adding in the midpoint in the z direction - assuming only two regions
    %     NpoiVec = Nbound1c(2:end) - Nbound1c(1:(end-1));
    %     midPt = find(NpoiVec == max(NpoiVec)) - 1;
    %     Nbound1cTemp = zeros(1, length(Nbound1c) + 1);
    %     Nbound1cTemp(1:midPt) = Nbound1c(1:midPt);
    %     Nbound1cTemp(midPt + 1) = (Nbound1c(midPt) + Nbound1c(midPt + 1))/2;
    %     Nbound1cTemp((midPt + 2):length(Nbound1cTemp)) = Nbound1c((midPt + 1):length(Nbound1c));
    %     Nbound1c = Nbound1cTemp;
    
    %Loading, Sorting and Saving the particles
    for k=kstart:kend
        
        if (k ~= kend)
            %Timing
            %if (k == kstart || k == kstart+1 || mod(k, 10) == 0 || k == kend)
                k
            %end
            
            %For each general boundary point determine how many points are in that region
            Npoi=Nbound1c(k+1)-Nbound1c(k);                                     % Number of points in a divison?
            fseek(f,(Nbound1c(k)-1)*5*4,'bof');                                 % find the (Nbound1c(k)-1)*5*4th point from the beginning of the file
            
            %Initializing the variables for the dist stuff
            cntmat=repmat(0,[Nx,Nvpar,Nvperp]);                                 % Create a blank 3 dim mat for the countmap
            cntVmat=repmat(0,[Nx,Nvperp]);                                      % Create a blank 3 dim mat for the countmap
            cntEkinmat=repmat(0,[Nx,Nvpar,Nvperp]);                             % Create a blank 3 dim mat for the countmap
            
            %Creating the struct that will contain the z-sorted points and
            %x-sorted points
            
            
            data1=fread(f,Npoi*5,'single');
            %Reshape the data to the proper size
            data=reshape(data1,[5,length(data1)/5]);
            
            %         %Add in the points on the equivalent half of the data
            %         fseek(f,(Nbound1c(k + midPt + 1)-1)*5*4,'bof');
            %         Npoi=Nbound1c(k + midPt  + 2)-Nbound1c(k + midPt +1);
            %         data2=fread(f,Npoi*5,'single');
            %         data2 = reshape(data2,[5,length(data2)/5]);
            %         data = [data data2];
            
            XZPoints = data(1:2,:);
            ktemp = mod(k,3);
            if ktemp == 0
                ktemp = 3;
            end
            minXTemp = minXAbs + (ktemp-1)*Xbloc*dx;
            maxXTemp = minXTemp + Xbloc*dx;
            diffXTemp = abs(maxXTemp - minXTemp)/Ndivk;
            
            %Spatially sorting the zpoints for a given xpoint range
            kztemp = ((k-ktemp)/3) + 1;
            minZTemp = minZAbs + (kztemp-1)*Zbloc*dz;
            maxZTemp = minZTemp + Zbloc*dz;
            diffZTemp = (maxZTemp - minZTemp)/Ndiv;
            
            %Preallocating some arrays
            [ Xpts0, Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32, ...
                Xpts33 ] = findXpts( data, minXTemp, diffXTemp );
            
            [ Zpts0, Zpts1, Zpts2, Zpts3, Zpts4] = findZpts( data, minZTemp, diffZTemp );
            
            %Sorting the particles - including the particles that belong in
            %adjacent boxes
            if(mod(k,3) == 1)
                ZSortPts = ZSortPtsSort( ZSortPts, data, Zpts0, Zpts1, Zpts2,...
                    Zpts3, Zpts4,...
                    Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                    Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                    Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                    Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                
                ktempm1=mod(k-1,3);
                ktempp1=mod(k+1,3);
                if ktempm1 == 0
                    ktempm1 = 3;
                end
                if ktempp1 == 0
                    ktempp1 = 3;
                end
                kztempm1 = (((k-1)-ktempm1)/3)-1;
                kztempp1 = (((k+1)-ktempm1)/3)-1;
                
                if (isempty(Xpts0) ~= 1 && (kztempm1 == kztemp))
                    minXTemp = minXAbs + (ktemp-2)*Xbloc*dx;
                    maxXTemp = minXTemp + Xbloc*dx*(ktemp-1);
                    
                    [ Xpts0, Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16,...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32, ...
                        Xpts33 ] = findXpts( data, minXTemp, diffXTemp );
                    
                    ZSortPts3 = ZSortPtsSort( ZSortPts3, data, Zpts0, Zpts1, Zpts2,...
                        Zpts3, Zpts4, ...
                        Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16,...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                    
                    if(isempty(Xpts0) ~= 1)
                        disp('Xpts0 is still not empty and has this many elements:')
                        size(Xpts0)
                    end
                    
                end
                
                
                if (isempty(Xpts33) ~= 1 && (kztempp1 == kztemp))
                    minXTemp = minXAbs + (ktemp)*Xbloc*dx;
                    maxXTemp = minXTemp + Xbloc*dx*(ktemp+1);
                    
                    [ Xpts0, Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32, ...
                        Xpts33 ] = findXpts( data, minXTemp, diffXTemp );
                    
                    ZSortPts2 = ZSortPtsSort( ZSortPts2, data, Zpts0, Zpts1, Zpts2,...
                        Zpts3, Zpts4, ...
                        Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                    
                    if(isempty(Xpts33) ~= 1)
                        disp('Xpts33 is still not empty and has this many elements:')
                        size(Xpts33)
                    end
                end
                
                
            elseif (mod(k,3) == 2)
                ZSortPts2 = ZSortPtsSort( ZSortPts2, data, Zpts0, Zpts1, Zpts2,...
                    Zpts3, Zpts4, ...
                    Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                    Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                    Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                    Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                
                ktempm1=mod(k-1,3);
                ktempp1=mod(k+1,3);
                if ktempm1 == 0
                    ktempm1 = 3;
                end
                if ktempp1 == 0
                    ktempp1 = 3;
                end
                kztempm1 = (((k-1)-ktempm1)/3)-1;
                kztempp1 = (((k+1)-ktempm1)/3)-1;
                
                
                if (isempty(Xpts0) ~= 1 && (kztempm1 == kztemp))
                    minXTemp = minXAbs + (ktemp-2)*Xbloc*dx;
                    maxXTemp = minXTemp + Xbloc*dx*(ktemp-1);
                    
                    [ Xpts0, Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32, ...
                        Xpts33 ] = findXpts( data, minXTemp, diffXTemp );
                    
                    ZSortPts = ZSortPtsSort( ZSortPts, data, Zpts0, Zpts1, Zpts2,...
                        Zpts3, Zpts4, ...
                        Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                    
                    if(isempty(Xpts0) ~= 1)
                        disp('Xpts0 is still not empty and has this many elements:')
                        size(Xpts0)
                    end
                    
                end
                
                
                if (isempty(Xpts33) ~= 1 && (kztempp1 == kztemp))
                    minXTemp = minXAbs + (ktemp)*Xbloc*dx;
                    maxXTemp = minXTemp + Xbloc*dx*(ktemp+1);
                    
                    [ Xpts0, Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32, ...
                        Xpts33 ] = findXpts( data, minXTemp, diffXTemp );
                    
                    ZSortPts3 = ZSortPtsSort( ZSortPts3, data, Zpts0, Zpts1, Zpts2,...
                        Zpts3, Zpts4, ...
                        Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                    
                    if(isempty(Xpts33) ~= 1)
                        disp('Xpts33 is still not empty and has this many elements:')
                        size(Xpts33)
                    end
                end
                
            else
                ZSortPts3 = ZSortPtsSort( ZSortPts3, data, Zpts0, Zpts1, Zpts2,...
                    Zpts3, Zpts4, ...
                    Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                    Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                    Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                    Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                
                ktempm1=mod(k-1,3);
                ktempp1=mod(k+1,3);
                if ktempm1 == 0
                    ktempm1 = 3;
                end
                if ktempp1 == 0
                    ktempp1 = 3;
                end
                kztempm1 = (((k-1)-ktempm1)/3)-1;
                kztempp1 = (((k+1)-ktempm1)/3)-1;
                
                if (isempty(Xpts0) ~= 1 && (kztempm1 == kztemp))
                    minXTemp = minXAbs + (ktemp-2)*Xbloc*dx;
                    maxXTemp = minXTemp + Xbloc*dx*(ktemp-1);
                    
                    [ Xpts0, Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32, ...
                        Xpts33 ] = findXpts( data, minXTemp, diffXTemp );
                    
                    ZSortPts2 = ZSortPtsSort( ZSortPts2, data, Zpts0, Zpts1, Zpts2,...
                        Zpts3, Zpts4, ...
                        Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                    
                    if(isempty(Xpts0) ~= 1)
                        disp('Xpts0 is still not empty and has this many elements:')
                        size(Xpts0)
                    end
                    
                end
                
                
                if (isempty(Xpts33) ~= 1 && (kztempp1 == kztemp))
                    minXTemp = minXAbs + (ktemp)*Xbloc*dx;
                    maxXTemp = minXTemp + Xbloc*dx*(ktemp+1);
                    
                    [ Xpts0, Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32, ...
                        Xpts33 ] = findXpts( data, minXTemp, diffXTemp );
                    
                    ZSortPts = ZSortPtsSort( ZSortPts, data, Zpts0, Zpts1, Zpts2,...
                        Zpts3, Zpts4, ...
                        Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
                        Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, ...
                        Xpts17, Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24,...
                        Xpts25, Xpts26, Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32);
                    
                    if(isempty(Xpts33) ~= 1)
                        disp('Xpts33 is still not empty and has this many elements:')
                        size(Xpts33)
                    end
                end
                
            end
        end
        
        % Saving the particles and making them into a countmap/distribution
        if (k ~= kstart)
            for n1 = 1:Ndivk
                for l = 0:(Ndiv-1)
                    
                    % Loading in the correct set of spatially distributed points
                    data = [ ];
                    if (mod(k,3) == 1 && k ~= kstart)
                        data = ZSortPtsDetFn(ZSortPts3, n1,l, data);
                    elseif (mod(k,3) == 2)
                        data = ZSortPtsDetFn(ZSortPts, n1,l, data);
                    else
                        data = ZSortPtsDetFn(ZSortPts2, n1,l, data);
                    end
                    
                    kcount = (k-1)*Ndivk + n1;
                    if l == (Ndiv - 2)
                        N89(1, kcount) = size(data, 2);
                    elseif l==(Ndiv - 1)
                        N89(2, kcount) = size(data, 2);
                    end
                    
                    cntVmat=repmat(0,[Nx,Nvpar,Nvperp]);                            % Create a blank 3 dim mat for the countmap
                    cntEkinmat=repmat(0,[Nx,Nvpar,Nvperp]);                         % Create a blank 3 dim mat for the countmap
                    cntAngmat=repmat(0,[Nx,Nvpar,Nvperp]);
                    
                    nlxp = nlx/2;
                    nlzp = nlz/2;
                    
                    xvp=linspace(-Lx/2 , Lx/2, size(bx,1));%(1:nlxp)*Lx/nlxp-Lx/2;                                                     % x-coordinate when using readFunctionPump
                    zvp=linspace(-Lz/2 , Lz/2, size(bx,2));%(1:nlzp)*Lz/nlzp-Lz/2;                                                     % z-coordinate when using readFunctionPump
                    
                    
                    %Finding the magnetic fields at regular intervals
                    bxcut=interp2(xvp,zvp,bx',data(1,:),data(2,:));
                    bycut=interp2(xvp,zvp,by',data(1,:),data(2,:));
                    bzcut=interp2(xvp,zvp,bz',data(1,:),data(2,:));
                    
                    %Normalizing the values
                    Btot=sqrt(bxcut.^2+bycut.^2+bzcut.^2);
                    bxcut=bxcut./Btot;
                    bycut=bycut./Btot;
                    bzcut=bzcut./Btot;
                    
                    %Finding the v_parallel and v_perp along B for each particle
                    vpar=data(3,:).*bxcut+data(4,:).*bycut+data(5,:).*bzcut;
                    vperp=sqrt(data(3,:).^2+data(4,:).^2+data(5,:).^2-vpar.^2);
                    
                    %Data now has vpar and vperp
                    data=[data; vpar; vperp];
                    
                    %Starting to do the dist stuff
                    IImat=data([2,6,7],:);
                    Nz = nlz;
                    IImat(1,:)=round((abs(IImat(1,:))-Zlim1)/(Zlim2-Zlim1)*Nz +0.5);     % index for z coordinate
                    IImat(2,:)=round((IImat(2,:)+Vlim)/(2*Vlim)*Nvpar +0.5);        % index for vpar coordinate
                    IImat(3,:)=round((IImat(3,:))/Vlim*Nvperp +0.5);                % index for vperp coordinate
                    
                    %Energy dist
                    Ekin = [];
                    Vtot=sqrt(data([6],:).^2+data([7],:).^2);
                    Ekin(1,:) = sqrt(1 - IImat(2,:).^2) - 1;
                    Ekin(2,:) = sqrt(1 - IImat(3,:).^2) - 1;
                    IIVtot=round((Vtot)/Vlim*Nvperp +0.5);                          % index for vperp coordinate
                    IIEkin = round(abs((Ekin)/Elim*NEkin) +0.5);
                    
                    %Ksi and V dist
                    AngMat = [];
                    AngMat(1,:) = (data(6,:)./Vtot);
                    AngMat(2,:) = Vtot;
                    IIAngMat = zeros(size(AngMat,1), size(AngMat,2));
                    IIAngMat(1,:) = round((AngMat(1,:) + 1)/2*Nvpar + 0.5);
                    IIAngMat(2,:) = round(AngMat(2,:)/(Vlim)*Nvperp + 0.5);
                    
                    %If the vpar and vperp are the right values, increment the correct box by 1 in the countmap
                    for m=1:size(data,2)
                        
%                         if IImat(1,m)>0 & IImat(1,m)< (Nz+1) & ...
%                                 IImat(2,m)>0 & IImat(2,m)< (Nvpar+1) & ...
%                                 IImat(3,m)>0 & IImat(3,m)< (Nvperp+1)
%                             
%                             cntVmat(IImat(1,m),IImat(2,m),IImat(3,m)) = ...
%                                 cntVmat(IImat(1,m),IImat(2,m),IImat(3,m))+1;
%                         end
%                         
%                         if IImat(1,m)>0 & IImat(1,m)< (Nz+1) & ...
%                                 IIEkin(1,m)>0 & IIEkin(1,m)< (NEkin+1) & ...
%                                 IIEkin(2,m)>0 & IIEkin(2,m)< (NEkin+1)
%                             cntEkinmat(IImat(1,m),IIEkin(2,m), IIEkin(3,m)) = ...
%                                 cntEkinmat(IImat(1,m),IIEkin(2,m), IIEkin(3,m))+1;
%                         end
                        
                        if IImat(1,m)>0 & IImat(1,m)< (Nz+1) & ...
                                IIAngMat(1,m)>0 & IIAngMat(1,m)< (Nvpar+1) & ...
                                IIAngMat(2,m)>0 & IIAngMat(2,m)< (Nvperp+1)
                            cntAngmat(IImat(1,m),IIAngMat(1,m), IIAngMat(2,m)) = ...
                                cntAngmat(IImat(1,m),IIAngMat(1,m), IIAngMat(2,m))+1;
                        end
                        
                    end
                    
%                     if (k == 8)
%                         disp('Paused here so could look at distribution')
%                         pause()
%                     end
                    
                    %Finding the Jacobian
                    vpar=linspace(-Vlim,Vlim,Nvpar);                                % vpar at center of intevals RENAME
                    vperp=linspace(0,Vlim,Nvperp+1);                                % vperp for boundaries of intervals RENAME
                    JackV=vperp(2:end).^2-vperp(1:(end-1)).^2;
                    
                    epar = linspace(0, Elim, Nvpar);
                    eperp = linspace(0, Elim, Nvperp+1);
                    JackEkin = eperp(2:end) - eperp(1:(end-1));
                    
                    vVec = linspace(0, Vlim, Nvperp+1); 
                    JackAng =vVec(2:end).^3-vVec(1:(end-1)).^3;
                    
                    %Finding the distribution from the countmap by dividing by the Jacobian
                    distV = cntVmat;
                    distEkin = cntEkinmat;
                    
                    distAng = cntAngmat;
                    
                    for m=1:length(JackV)
                        distV(:,:,m)=distV(:,:,m)/JackV(m);
                        distEkin(:,:,m)=distEkin(:,:,m)/JackEkin(m);
                        distAng(:,:, m) = distAng(:,:,m)/JackAng(m); 
                    end
                    
                    
                    %Saving the distribution file (see extra code for
                    %saving Edist file - omitted for time
%                     com=['save ' savedir '/Run' folderName 'T' fileName{i} '/distV_k'...
%                         num2str(k) 'n1' num2str(n1) 'l' num2str(l) 'v1' ' distV cntVmat'];
%                     eval(com);
%                     
                    %                     com=['save ' savedir '/Run' folderName 'T' fileName{i} '/distEkin_k'...
                    %                         num2str(k) 'n1' num2str(n1) 'l' num2str(l) 'v1' ' distEkin cntEkinmat'];
                    %                     eval(com);
                    
                    com=['save ' savedir '/Run' folderName 'T' fileName{i} '/distAng_k'...
                        num2str(k) 'n1' num2str(n1) 'l' num2str(l) 'v1' ' distAng cntAngmat'];
                    eval(com);                    
                    
                end
            end
            
            % Emptying the ZSortPts struct after its data has been saved
            if (mod(k,3) == 1 && k ~= kstart)
                ZSortPts3=ZSortCreate( );
            elseif (mod(k,3) == 2)
                ZSortPts=ZSortCreate( );
            else
                ZSortPts2=ZSortCreate( );
            end
            
        end
        
        
    end
    
    fclose(f);
    
end
toc
%%%%%%%%%%%%%%
% Extra Code %
%%%%%%%%%%%%%%

% loaddir = ['/ext/VPICRuns/Pumping2DSinglePump/T2Nppc200'];
% savedir = ['/ext/VPICRuns/ChgGeo/Pumping2DSinglePump/T2Nppc200/particlesDist'];
% loaddir = ['/ext/VPICRuns/ChgGeo/T12B305V2/'];
% savedir = ['/ext/VPICRuns/ChgGeo/T12B305V2/particlesDist'];
% loaddir = ['/ext/VPICRuns/ChgGeo/T11B305V05T' folderName '/'];
% savedir = ['/ext/VPICRuns/ChgGeo/T11B305V05T' folderName '/particlesDist'];
% loaddir = ['/ext/VPICRuns/Pumping2DSinglePump/Test' folderName '/particles'];
% savedir = ['/ext/VPICRuns/Pumping2DSinglePump/Test' folderName '/particlesDist'];

%kcount = (k-1)*Ndivk + n1;

%                 com=['save ' savedir '/Run' folderName 'T' fileName{i} '/distEkin_k'...
%                     num2str(kcount) 'l' num2str(l) 'v1'  ' distEkin cntEkinmat'];
%                 eval(com);

%         figure(1), clf
%         plot(N89(1,:), 'Color',ColorSet(1,:))
%         hold on;
%         plot(N89(2,:), 'Color',ColorSet(2,:))
%         hold on;
%         plot(N89(3,:), 'Color',ColorSet(3,:))
%         hold on;
%         plot(N89(4,:), 'Color',ColorSet(4,:))
%         hold on;
%         plot(N89(5,:), 'Color',ColorSet(5,:))
%         hold on;
%         xlim([1 k*Ndivk])
%         legend('5', '6', '7', '8', '9');
%
% for k = 1:length(Nbound1)
% plot([1 1]*Nbound1(k), [-720, 720], 'g')
% hold on;
% end
% for k = 1:length(Nbound1)
% plot([1 Nbound1(end)], [(-720 + Xbloc*dx*k) (-720 + Xbloc*dx*k)], 'r')
% hold on;
% end