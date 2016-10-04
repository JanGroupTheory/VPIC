function [fileName, nlx, nlz, Lx, Lz, mime, Npart, Vlim, Elim, maxZAbs, ...
    minZAbs, Xbloc, Zbloc, maxXAbs, minXAbs, delT] = loadVarsloaddir( loaddir, folderName )
%loadVarsloaddir Summary of this function goes here
%   By Emily Lichko

if strcmp(loaddir,['/ext/VPICRuns/Pumping2DSinglePump/Test' folderName '/particles'])
    nlx = 12096;
    nlz = 800;
    Lx = 1200;
    Lz = 80;
    mime = 100;
    Npart = 242660000;
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 40.0;
    minZAbs = 30.0;
    
%     maxZAbs1 = 10.0;
%     minZAbs2 = 0.0;
    
    maxXAbs = 600.0;
    minXAbs = -600.0;
    
    nlx = nlx/2;
    
    %X size of read in block
    Xbloc = 8;
    
elseif strcmp(loaddir,['/ext/VPICRuns/Pumping2DSinglePump/T2Nppc200'])
    fileName = {'0' '45050' '90100'};
    nlx = 12096;
    nlz = 800;
    Lx = 1200;
    Lz = 80;
    mime = 100;
    Npart = 242660000;
    
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 40.0;
    minZAbs = 30.0;
    
%     maxZAbs1 = 10.0;
%     minZAbs2 = 0.0;
    
    maxXAbs = 600.0;
    minXAbs = -600.0;
    
    nlx = nlx/2;
    
    %X size of read in block
    Xbloc = 8;
    
elseif strcmp(loaddir,['/ext/VPICRuns/ChgGeo/T11B305V05T' folderName '/'])
    fileName = {'0'};
    nlx = 18432;
    nlz = 576;
    Lx = 1440;
    Lz = 60;
    mime = 10000;
    Npart = 637730000;
    delT = 253;
    
    Vlim=0.5;%1;%1.4;                                                       % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 0.1340;                                                          % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 21.0;
    maxXAbs = 720.0;
    minXAbs = -720.0;
    
    %X size of read in block
    Xbloc = 16;
    
elseif strcmp(loaddir,['/ext/VPICRuns/ChgGeo/T12B305V2/'])
    fileName = {'0' '8000' '16000'};
    nlx = 4608;
    nlz = 576;
    Lx = 1440;
    Lz = 60;
    mime = 10000;
    Npart = 130395000;%1592525000;
    delT = 160;
    
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 21.0;
    maxXAbs = 720.0;
    minXAbs = -720.0;
    
    %X size of read in block
    Xbloc = 16; %[8 16 16]; %for 16000 - 8 for 0 and
    
elseif strcmp(loaddir,['/ext/VPICRuns/ChgGeo/T12B305V1/'])
    fileName = {'0' '3660' '7320' '10980' '14640'};
    nlx = 9216;
    nlz = 576;
    Lx = 1440;
    Lz = 60;
    mime = 10000;
    Npart = 321148000;
    delT = 183;
    
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 21.0;
    maxXAbs = 720.0;
    minXAbs = -720.0;
    
    %X size of read in block
    Xbloc = 16;%[8 8 8 8];
    
elseif strcmp(loaddir,['/ext/VPICRuns/ChgGeo/T12B305V05/'])
    fileName = {'0' '5060' '10120' '15180' '20240'};
    nlx = 18432;
    nlz = 576;
    Lx = 1440;
    Lz = 60;
    mime = 10000;
    Npart = 640428000;%642786000;%650486000;%665180000;
    delT = 253;
    
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 21.0;
    maxXAbs = 720.0;
    minXAbs = -720.0;
    
    %X size of read in block
    Xbloc = 16;%[8 8 8 8];
    
elseif strcmp(loaddir,['/ext/VPICRuns/ChgGeo/T13B305V05RealMassRatio/'])
    fileName = {'0' '5060' '10120' '15180' '20240'};
    nlx = 18432;
    nlz = 576;
    Lx = 1440;
    Lz = 60;
    mime = 1836;
    Npart = 640428000;%647292000;%703478000;
    delT = 253;
    
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 21.0;
    maxXAbs = 720.0;
    minXAbs = -720.0;
    
    %X size of read in block
    Xbloc = 16;%[8 8 8 8];
    
elseif strcmp(loaddir,['/ext/VPICRuns/ChgGeo/T15B305V05ChgPumpRegion/'])
    fileName = {'0' '5060' '10120' '15180' '20240'};
    nlx = 18432;
    nlz = 576;
    Lx = 1440;
    Lz = 60;
    mime = 1836;
    Npart = 640428000;%647292000;%703478000;
    delT = 253;
    
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 21.0;
    maxXAbs = 720.0;
    minXAbs = -720.0;
    
    %X size of read in block
    Xbloc = 16;%[8 8 8 8];
    
elseif strcmp(loaddir,['/ext/VPICRuns/ChgGeo/T16B305V05Tsteps25/'])
    fileName = {'0' '500' '1000' '1500' '2000' '2500' '3000' '3500' '4000'...
        '4500' '5000' '5500' '6000' '6500' '7000' '7500' '8000' '8500' ...
        '9000' '9500' '10000' '10500' '11000' '11500' '12000' '12500' ...
        '13000' '13500' '14000' '14500' '15000' '15500' '16000' '16500' ...
        '17000' '17500' '18000' '18500' '19000' '19500' '20000' '20500'...
        '21000' '21500' '22000' '22500'};
    nlx = 18432;
    nlz = 576;
    Lx = 1440;
    Lz = 60;
    mime = 10000;
    Npart = 47416000;%48274000;%48460000;%50192000;
    delT = 25;
    
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 28.0;
    maxXAbs = 0.0;
    minXAbs = -480.0;
    
    %X size of read in block
    Xbloc = 16;%[8 8 8 8];
elseif strcmp(loaddir,['/ext/VPICRuns/ChgGeo/T17B305V05Tsteps25RMR/'])
    fileName = {'0' '500' '1000' '1500' '2000' '2500' '3000' '3500' '4000'...
        '4500' '5000' '5500' '6000' '6500' '7000' '7500' '8000' '8500' ...
        '9000' '9500' '10000' '10500' '11000' '11500' '12000' '12500' ...
        '13000' '13500' '14000' '14500' '15000' '15500' '16000' '16500' ...
        '17000' '17500' '18000' '18500' '19000' '19500' '20000' '20500'...
        '21000' '21500' '22000' '22500'};
    nlx = 18432;
    nlz = 576;
    Lx = 1440;
    Lz = 60;
    mime = 1836;
    Npart = 940428000;%647292000;%703478000;
    delT = 25;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 28.0;
    maxXAbs = 0.0;
    minXAbs = -480.0;
    
    %X size of read in block
    Xbloc = 16;%[8 8 8 8];
    
elseif strcmp(loaddir,['/ext/VPICRuns/Pumping2DSinglePump/T5Nppc50'])
    fileName = {'0' '4500' '9000' '13500' '18000' '22500' '27000' '31500'};
    nlx = 12096;
    nlz = 800;
    Lx = 1200;
    Lz = 80;
    mime = 10000;
    Npart = 967680000;
    delT = 225;
    
    Vlim=1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 30.0;
    minZAbs = 21.0;
    maxXAbs = 720.0;
    minXAbs = -720.0;
    
    %X size of read in block
    Xbloc = 16;
    
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T59703StdStandCollOff/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960' '95040' '96120'...
        '97200' '98280' '99360' '100440' '101520' '102600' '103680'  '104760'...
        '105840' '106920' '108000' '109080' '110160' '111240' '112320' '113400'...
        '114480' '115560' '116640' '117720' '118800' '119880'};
    nlx = 2688;
    nlz = 536;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart = 81044000;%81054000;%81102000;%83296000;%83428000;%474160000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 112;%[8 8 8 8];
    Zbloc = 134;
    
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T60703ObliqStandAng45/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960'};
%     
%         fileName = {'0' '3240' '5400' '8640' '10800' '12960' '15120' '17280' '19440' ...
%         '22680' '24840' '28080' '30240' '32400' '35640' '37800' '39960' '41040' ...
%         '43200' '46440' '49680' '52920' '55080' '57240' '59400' '61560' '63720' ...
%         '65880' '69120' '71280' '74520' '76680' '78840' '79920' '82080' '85320' ...
%         '88560' '91800' '93960'};
    nlx = 2688;
    nlz = 536;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart = 81046000;%81118000;%81188000;%81876000;%82704000;%83048000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 112;%[8 8 8 8];
    Zbloc = 134;
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T61StdStandVthe2/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960' '95040' '96120'...
        '97200' '98280' '99360' '100440' '101520' '102600' '103680'  '104760'...
        '105840' '106920' '108000' '109080' '110160' '111240' '112320' '113400'...
        '114480' '115560' '116640' '117720' '118800' '119880' '120960' '122040'...
        '123120' '124200' '125280' '126360' '127440' '128520' '129600' '130680'};
    
    nlx = 1344;
    nlz = 268;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart = 20296000;%20318000;%20360000;%20452000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 56;%112;%[8 8 8 8];
    Zbloc = 67;%134;
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T62ObliqStandAng45Vthe2/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960' '95040' '96120'...
        '97200' '98280' '99360' '100440' '101520' '102600' '103680'  '104760'...
        '105840' '106920' '108000' '109080' '110160' '111240' '112320' '113400'...
        '114480' '115560' '116640' '117720' '118800' '119880' '120960' '122040'...
        '123120' '124200' '125280' '126360' '127440' '128520' '129600' '130680'};
    
    nlx = 1344;
    nlz = 268;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart = 20276000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 56;%112;%[8 8 8 8];
    Zbloc = 67;%134;
    
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T63StdProp/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960' '95040' '96120'...
        '97200' '98280' '99360' '100440' '101520' '102600' '103680'  '104760'...
        '105840' '106920' '108000' '109080' '110160' '111240' '112320' '113400'...
        '114480' '115560'};
    nlx = 2688;
    nlz = 536;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart =  81062000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 112;%[8 8 8 8];
    Zbloc = 134;
    
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T65StdPropAmpUp/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960' '95040' '96120'...
        '97200' '98280' '99360' '100440' '101520' '102600' '103680'  '104760'...
        '105840' '106920' '108000' '109080' '110160' '111240' '112320' '113400'...
        '114480' '115560'};
    nlx = 2688;
    nlz = 536;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart =  81110000;%81138000;%81062000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 112;%[8 8 8 8];
    Zbloc = 134;
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T66ObliqPropAmpUp/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960' '95040' '96120'...
        '97200' '98280' '99360' '100440' '101520' '102600' '103680'  '104760'...
        '105840' '106920' '108000' '109080' '110160' '111240' '112320' '113400'...
        '114480' '115560'};
    nlx = 2688;
    nlz = 536;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart =  81062000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 112;%[8 8 8 8];
    Zbloc = 134;
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T69StdPropAmpUp/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960' '95040' '96120'...
        '97200' '98280' '99360' '100440' '101520' '102600' '103680'  '104760'...
        '105840' '106920' '108000' '109080' '110160' '111240' '112320' '113400'...
        '114480' '115560'};
    nlx = 2688;
    nlz = 536;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart =  81076000;%83062000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 112;%[8 8 8 8];
    Zbloc = 134;
elseif strcmp(loaddir,['/ext/VPICRuns/AlfvenWavesInBox/T70ObliqPropAmpUp/'])
    fileName = {'0' '1080' '2160' '3240' '4320' '5400' '6480' '7560' '8640'...
        '9720' '10800' '11880' '12960' '14040' '15120' '16200' '17280' '18360' ...
        '19440' '20520' '21600' '22680' '23760' '24840' '25920' '27000' '28080'...
        '29160' '30240' '31320' '32400' '33480' '34560' '35640' '36720' '37800'...
        '38880' '39960' '41040' '42120' '43200' '44280' '45360' '46440' '47520'...
        '48600' '49680' '50760' '51840' '52920' '54000' '55080' '56160' '57240'...
        '58320' '59400' '60480' '61560' '62640' '63720' '64800' '65880' '66960'...
        '68040' '69120' '70200' '71280' '72360' '73440' '74520' '75600' '76680'...
        '77760' '78840' '79920' '81000' '82080' '83160' '84240' '85320' '86400'...
        '87480' '88560' '89640' '90720' '91800' '92880' '93960' '95040' '96120'...
        '97200' '98280' '99360' '100440' '101520' '102600' '103680'  '104760'...
        '105840' '106920' '108000' '109080' '110160' '111240' '112320' '113400'...
        '114480' '115560'};
    nlx = 2688;
    nlz = 536;
    Lx = 320;
    Lz = 320;
    mime = 25;
    Npart =  81078000;
    delT = 54;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 160.0;
    minZAbs = -160.0;
    maxXAbs = 160.0;
    minXAbs = 120.0;
    
    %X size of read in block
    Xbloc = 112;%[8 8 8 8];
    Zbloc = 134;
    
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1D/1Dnu0/'])
    fileName = {'181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart =  51230000;
    delT = 453;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
    
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1D/1Dnu30/'])
    fileName = {'181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart =  506280000;
    delT = 453;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1D/1Dnu10/'])
    fileName = {'181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart =  50766000;
    delT = 453;
    
    Vlim = 1.4;                                                             % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1D/1Dnu3/'])
    fileName = {'181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart =  51916000;%5062800000;
    delT = 453;
    
    Vlim = 1.4;                                                             % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
    
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1DCollChg/1Dnu30Big/'])
    fileName = {'0' '90600' '181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart = 36628000;
    delT = 453;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1DCollChg/1Dnu10Big/'])
    fileName = {'0' '90600' '181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart = 36628000;
    delT = 453;
    
    Vlim = 1.4;                                                             % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1DCollChg/1Dnu3Big/'])
    fileName = {'0' '90600' '181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart =  519160000;%5062800000;
    delT = 453;
    
    Vlim = 1.4;                                                             % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
    
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1DCollChgAmpDecr/Big/1Dnu30/'])
    fileName = {'0' '90600' '181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart = 366280000;
    delT = 453;
    
    Vlim = 1.4;                                                               % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1DCollChgAmpDecr/Big/1Dnu10/'])
    fileName = {'0' '90600' '181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart =  5076600000;
    delT = 453;
    
    Vlim = 1.4;                                                             % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1DCollChgAmpDecr/Big/1Dnu3/'])
    fileName = {'0' '90600' '181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart = 366280000;%519160000;%5062800000;
    delT = 453;
    
    Vlim = 1.4;                                                             % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
elseif strcmp(loaddir,['/ext/VPICRuns/MagneticPumping/1DCollChgAmpDecr/Big/1Dnu0/'])
    fileName = {'0' '90600' '181200' '271800' '362400' '453000' '543600' '634200' '724800' ...
        '815400' '906000' '996600' '1087200' '1177800' '1268400' '1359000' ...
        '1449600' '1540200' '1630800' '1721400' '1812000'};
    nlx = 80;
    nlz = 1632;
    Lx = 4;
    Lz = 80;
    mime = 100;
    Npart = 366280000;%519160000;%5062800000;
    delT = 453;
    
    Vlim = 1.4;                                                             % (relativistic gamma)*v/c (units of all Bill's sim) - max velocity
    Elim = 1.0;                                                             % sqrt(1 - Vlim^2) - 1;
    
    % Input Variables
    maxZAbs = 10.0;
    minZAbs = -10.0;
    maxXAbs = 2.0;
    minXAbs = -2.0;
    
    %X size of read in block
    Xbloc = 10;%[8 8 8 8];
    Zbloc = 34;
    
end

end

