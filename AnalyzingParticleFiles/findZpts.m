function [ Zpts0, Zpts1, Zpts2, Zpts3, Zpts4, Zpts5, Zpts6, Zpts7, Zpts8, ...
    Zpts9] = findZpts( data, minZAbs, diffZAbs )
%findZpts Summary of this function goes here
%   Detailed explanation goes here

        Zpts0 = [];
        Zpts1 = [];
        Zpts2 = [];
        Zpts3 = [];
        Zpts4 = [];
%         Zpts5 = [];
%         Zpts6 = [];
%         Zpts7 = [];
%         Zpts8 = [];
%         Zpts9 = [];
        [foo, Zpts0] = find((data(2,:)) <= (minZAbs + diffZAbs));
        [foo, Zpts1] = find((data(2,:)) <= (minZAbs + 2*diffZAbs));
        [foo, Zpts2] = find((data(2,:)) <= (minZAbs + 3*diffZAbs));
        [foo, Zpts3] = find((data(2,:)) <= (minZAbs + 4*diffZAbs));
        [foo, Zpts4] = find((data(2,:)) <= (minZAbs + 5*diffZAbs));
%         [foo, Zpts5] = find(abs(data(2,:)) <= (minZAbs + 6*diffZAbs));
%         [foo, Zpts6] = find(abs(data(2,:)) <= (minZAbs + 7*diffZAbs));
%         [foo, Zpts7] = find(abs(data(2,:)) <= (minZAbs + 8*diffZAbs));
%         [foo, Zpts8] = find(abs(data(2,:)) <= (minZAbs + 9*diffZAbs));
%         [foo, Zpts9] = find(abs(data(2,:)) <= (minZAbs + 10*diffZAbs));
%         Zpts9 = setdiff(Zpts9, Zpts8);
%         Zpts8 = setdiff(Zpts8, Zpts7);
%         Zpts7 = setdiff(Zpts7, Zpts6);
%         Zpts6 = setdiff(Zpts6, Zpts5);
%         Zpts5 = setdiff(Zpts5, Zpts4);
        Zpts4 = setdiff(Zpts4, Zpts3);
        Zpts3 = setdiff(Zpts3, Zpts2);
        Zpts2 = setdiff(Zpts2, Zpts1);
        Zpts1 = setdiff(Zpts1, Zpts0);


end

