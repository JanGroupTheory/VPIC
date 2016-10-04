function [ Xpts0, Xpts1, Xpts2, Xpts3, Xpts4, Xpts5, Xpts6, Xpts7, Xpts8, ...
    Xpts9, Xpts10, Xpts11, Xpts12, Xpts13, Xpts14, Xpts15, Xpts16, Xpts17, ...
    Xpts18, Xpts19, Xpts20, Xpts21, Xpts22, Xpts23, Xpts24, Xpts25, Xpts26, ...
    Xpts27, Xpts28, Xpts29, Xpts30, Xpts31, Xpts32, Xpts33]...
    = findXpts( data, minXTemp, diffXTemp )
%findXpts Summary of this function goes here
%   Detailed explanation goes here

Xpts0 = [];
Xpts1 = [];
Xpts2 = [];
Xpts3 = [];
Xpts4 = [];
Xpts5 = [];
Xpts6 = [];
Xpts7 = [];
Xpts8 = [];
Xpts9 = [];
Xpts10 = [];
Xpts11 = [];
Xpts12 = [];
Xpts13 = [];
Xpts14 = [];
Xpts15 = [];
Xpts16 = [];
Xpts17 = [];
Xpts18 = [];
Xpts19 = [];
Xpts20 = [];
Xpts21 = [];
Xpts22 = [];
Xpts23 = [];
Xpts24 = [];
Xpts25 = [];
Xpts26 = [];
Xpts27 = [];
Xpts28 = [];
Xpts29 = [];
Xpts30 = [];
Xpts31 = [];
Xpts32 = [];
Xpts33 = [];
[foo, Xpts0] = find(data(1,:) < minXTemp);
[foo, Xpts1] = find(data(1,:) <= minXTemp + diffXTemp);
[foo, Xpts2] = find(data(1,:) <= minXTemp + 2*diffXTemp);
[foo, Xpts3] = find(data(1,:) <= minXTemp + 3*diffXTemp);
[foo, Xpts4] = find(data(1,:) <= minXTemp + 4*diffXTemp);
[foo, Xpts5] = find(data(1,:) <= minXTemp + 5*diffXTemp);
[foo, Xpts6] = find(data(1,:) <= minXTemp + 6*diffXTemp);
[foo, Xpts7] = find(data(1,:) <= minXTemp + 7*diffXTemp);
[foo, Xpts8] = find(data(1,:) <= minXTemp + 8*diffXTemp);
[foo, Xpts9] = find(data(1,:) <= minXTemp + 9*diffXTemp);
[foo, Xpts10] = find(data(1,:) <= minXTemp + 10*diffXTemp);
[foo, Xpts11] = find(data(1,:) <= minXTemp + 11*diffXTemp);
[foo, Xpts12] = find(data(1,:) <= minXTemp + 12*diffXTemp);
[foo, Xpts13] = find(data(1,:) <= minXTemp + 13*diffXTemp);
[foo, Xpts14] = find(data(1,:) <= minXTemp + 14*diffXTemp);
[foo, Xpts15] = find(data(1,:) <= minXTemp + 15*diffXTemp);
[foo, Xpts16] = find(data(1,:) <= minXTemp + 16*diffXTemp);
[foo, Xpts17] = find(data(1,:) <= minXTemp + 17*diffXTemp);
[foo, Xpts18] = find(data(1,:) <= minXTemp + 18*diffXTemp);
[foo, Xpts19] = find(data(1,:) <= minXTemp + 19*diffXTemp);
[foo, Xpts20] = find(data(1,:) <= minXTemp + 20*diffXTemp);
[foo, Xpts21] = find(data(1,:) <= minXTemp + 21*diffXTemp);
[foo, Xpts22] = find(data(1,:) <= minXTemp + 22*diffXTemp);
[foo, Xpts23] = find(data(1,:) <= minXTemp + 23*diffXTemp);
[foo, Xpts24] = find(data(1,:) <= minXTemp + 24*diffXTemp);
[foo, Xpts25] = find(data(1,:) <= minXTemp + 25*diffXTemp);
[foo, Xpts26] = find(data(1,:) <= minXTemp + 26*diffXTemp);
[foo, Xpts27] = find(data(1,:) <= minXTemp + 27*diffXTemp);
[foo, Xpts28] = find(data(1,:) <= minXTemp + 28*diffXTemp);
[foo, Xpts29] = find(data(1,:) <= minXTemp + 29*diffXTemp);
[foo, Xpts30] = find(data(1,:) <= minXTemp + 30*diffXTemp);
[foo, Xpts31] = find(data(1,:) <= minXTemp + 31*diffXTemp);
[foo, Xpts32] = find(data(1,:) <= minXTemp + 32*diffXTemp);
[foo, Xpts33] = find(data(1,:) > minXTemp + 32*diffXTemp);
Xpts33 = setdiff(Xpts33, Xpts32);
Xpts32 = setdiff(Xpts32, Xpts31);
Xpts31 = setdiff(Xpts31, Xpts30);
Xpts30 = setdiff(Xpts30, Xpts29);
Xpts29 = setdiff(Xpts29, Xpts28);
Xpts28 = setdiff(Xpts28, Xpts27);
Xpts27 = setdiff(Xpts27, Xpts26);
Xpts26 = setdiff(Xpts26, Xpts25);
Xpts25 = setdiff(Xpts25, Xpts24);
Xpts24 = setdiff(Xpts24, Xpts23);
Xpts23 = setdiff(Xpts23, Xpts22);
Xpts22 = setdiff(Xpts22, Xpts21);
Xpts21 = setdiff(Xpts21, Xpts20);
Xpts20 = setdiff(Xpts20, Xpts19);
Xpts19 = setdiff(Xpts19, Xpts18);
Xpts18 = setdiff(Xpts18, Xpts17);
Xpts17 = setdiff(Xpts17, Xpts16);
Xpts16 = setdiff(Xpts16, Xpts15);
Xpts15 = setdiff(Xpts15, Xpts14);
Xpts14 = setdiff(Xpts14, Xpts13);
Xpts13 = setdiff(Xpts13, Xpts12);
Xpts12 = setdiff(Xpts12, Xpts11);
Xpts11 = setdiff(Xpts11, Xpts10);
Xpts10 = setdiff(Xpts10, Xpts9);
Xpts9 = setdiff(Xpts9, Xpts8);
Xpts8 = setdiff(Xpts8, Xpts7);
Xpts7 = setdiff(Xpts7, Xpts6);
Xpts6 = setdiff(Xpts6, Xpts5);
Xpts5 = setdiff(Xpts5, Xpts4);
Xpts4 = setdiff(Xpts4, Xpts3);
Xpts3 = setdiff(Xpts3, Xpts2);
Xpts2 = setdiff(Xpts2, Xpts1);
Xpts1 = setdiff(Xpts1, Xpts0);

end

