function [ data ] = ZSortPtsDetFn(ZSortPts, n1 ,l, data )
%ZSortPtsDetFn Determines which subset of ZSortPtsDetFn to use
%   Used to make the code look cleaner
%   Takes m [1, 32]
%   And   l [0, 11]
%   By Emily Lichko

if n1 == 1
    if l == 0
        data = ZSortPts.l00k01;
    elseif l == 1
        data = ZSortPts.l01k01;
    elseif l == 2
        data = ZSortPts.l02k01;
    elseif l == 3
        data = ZSortPts.l03k01;
    elseif l == 4
        data = ZSortPts.l04k01;
    elseif l == 5
        data = ZSortPts.l05k01;
    elseif l == 6
        data = ZSortPts.l06k01;
    elseif l == 7
        data = ZSortPts.l07k01;
    elseif l == 8
        data = ZSortPts.l08k01;
    elseif l == 9
        data = ZSortPts.l09k01;
    elseif l == 10
        data = ZSortPts.l10k01;
    elseif l == 11
        data = ZSortPts.l11k01;
    end
elseif n1 == 2
    if l == 0
        data = ZSortPts.l00k02;
    elseif l == 1
        data = ZSortPts.l01k02;
    elseif l == 2
        data = ZSortPts.l02k02;
    elseif l == 3
        data = ZSortPts.l03k02;
    elseif l == 4
        data = ZSortPts.l04k02;
    elseif l == 5
        data = ZSortPts.l05k02;
    elseif l == 6
        data = ZSortPts.l06k02;
    elseif l == 7
        data = ZSortPts.l07k02;
    elseif l == 8
        data = ZSortPts.l08k02;
    elseif l == 9
        data = ZSortPts.l09k02;
    elseif l == 10
        data = ZSortPts.l10k02;
    elseif l == 11
        data = ZSortPts.l11k02;
    end
elseif n1 == 3
    if l == 0
        data = ZSortPts.l00k03;
    elseif l == 1
        data = ZSortPts.l01k03;
    elseif l == 2
        data = ZSortPts.l02k03;
    elseif l == 3
        data = ZSortPts.l03k03;
    elseif l == 4
        data = ZSortPts.l04k03;
    elseif l == 5
        data = ZSortPts.l05k03;
    elseif l == 6
        data = ZSortPts.l06k03;
    elseif l == 7
        data = ZSortPts.l07k03;
    elseif l == 8
        data = ZSortPts.l08k03;
    elseif l == 9
        data = ZSortPts.l09k03;
    elseif l == 10
        data = ZSortPts.l10k03;
    elseif l == 11
        data = ZSortPts.l11k03;
    end
elseif n1 == 4
    if l == 0
        data = ZSortPts.l00k04;
    elseif l == 1
        data = ZSortPts.l01k04;
    elseif l == 2
        data = ZSortPts.l02k04;
    elseif l == 3
        data = ZSortPts.l03k04;
    elseif l == 4
        data = ZSortPts.l04k04;
    elseif l == 5
        data = ZSortPts.l05k04;
    elseif l == 6
        data = ZSortPts.l06k04;
    elseif l == 7
        data = ZSortPts.l07k04;
    elseif l == 8
        data = ZSortPts.l08k04;
    elseif l == 9
        data = ZSortPts.l09k04;
    elseif l == 10
        data = ZSortPts.l10k04;
    elseif l == 11
        data = ZSortPts.l11k04;
    end
elseif n1 == 5
    if l == 0
        data = ZSortPts.l00k05;
    elseif l == 1
        data = ZSortPts.l01k05;
    elseif l == 2
        data = ZSortPts.l02k05;
    elseif l == 3
        data = ZSortPts.l03k05;
    elseif l == 4
        data = ZSortPts.l04k05;
    elseif l == 5
        data = ZSortPts.l05k05;
    elseif l == 6
        data = ZSortPts.l06k05;
    elseif l == 7
        data = ZSortPts.l07k05;
    elseif l == 8
        data = ZSortPts.l08k05;
    elseif l == 9
        data = ZSortPts.l09k05;
    elseif l == 10
        data = ZSortPts.l10k05;
    elseif l == 11
        data = ZSortPts.l11k05;
    end
elseif n1 == 6
    if l == 0
        data = ZSortPts.l00k06;
    elseif l == 1
        data = ZSortPts.l01k06;
    elseif l == 2
        data = ZSortPts.l02k06;
    elseif l == 3
        data = ZSortPts.l03k06;
    elseif l == 4
        data = ZSortPts.l04k06;
    elseif l == 5
        data = ZSortPts.l05k06;
    elseif l == 6
        data = ZSortPts.l06k06;
    elseif l == 7
        data = ZSortPts.l07k06;
    elseif l == 8
        data = ZSortPts.l08k06;
    elseif l == 9
        data = ZSortPts.l09k06;
    elseif l == 10
        data = ZSortPts.l10k06;
    elseif l == 11
        data = ZSortPts.l11k06;
    end
elseif n1 == 7
    if l == 0
        data = ZSortPts.l00k07;
    elseif l == 1
        data = ZSortPts.l01k07;
    elseif l == 2
        data = ZSortPts.l02k07;
    elseif l == 3
        data = ZSortPts.l03k07;
    elseif l == 4
        data = ZSortPts.l04k07;
    elseif l == 5
        data = ZSortPts.l05k07;
    elseif l == 6
        data = ZSortPts.l06k07;
    elseif l == 7
        data = ZSortPts.l07k07;
    elseif l == 8
        data = ZSortPts.l08k07;
    elseif l == 9
        data = ZSortPts.l09k07;
    elseif l == 10
        data = ZSortPts.l10k07;
    elseif l == 11
        data = ZSortPts.l11k07;
    end
elseif n1 == 8
    if l == 0
        data = ZSortPts.l00k08;
    elseif l == 1
        data = ZSortPts.l01k08;
    elseif l == 2
        data = ZSortPts.l02k08;
    elseif l == 3
        data = ZSortPts.l03k08;
    elseif l == 4
        data = ZSortPts.l04k08;
    elseif l == 5
        data = ZSortPts.l05k08;
    elseif l == 6
        data = ZSortPts.l06k08;
    elseif l == 7
        data = ZSortPts.l07k08;
    elseif l == 8
        data = ZSortPts.l08k08;
    elseif l == 9
        data = ZSortPts.l09k08;
    elseif l == 10
        data = ZSortPts.l10k08;
    elseif l == 11
        data = ZSortPts.l11k08;
    end
elseif n1 == 9
    if l == 0
        data = ZSortPts.l00k09;
    elseif l == 1
        data = ZSortPts.l01k09;
    elseif l == 2
        data = ZSortPts.l02k09;
    elseif l == 3
        data = ZSortPts.l03k09;
    elseif l == 4
        data = ZSortPts.l04k09;
    elseif l == 5
        data = ZSortPts.l05k09;
    elseif l == 6
        data = ZSortPts.l06k09;
    elseif l == 7
        data = ZSortPts.l07k09;
    elseif l == 8
        data = ZSortPts.l08k09;
    elseif l == 9
        data = ZSortPts.l09k09;
    elseif l == 10
        data = ZSortPts.l10k09;
    elseif l == 11
        data = ZSortPts.l11k09;
    end
elseif n1 == 10
    if l == 0
        data = ZSortPts.l00k10;
    elseif l == 1
        data = ZSortPts.l01k10;
    elseif l == 2
        data = ZSortPts.l02k10;
    elseif l == 3
        data = ZSortPts.l03k10;
    elseif l == 4
        data = ZSortPts.l04k10;
    elseif l == 5
        data = ZSortPts.l05k10;
    elseif l == 6
        data = ZSortPts.l06k10;
    elseif l == 7
        data = ZSortPts.l07k10;
    elseif l == 8
        data = ZSortPts.l08k10;
    elseif l == 9
        data = ZSortPts.l09k10;
    elseif l == 10
        data = ZSortPts.l10k10;
    elseif l == 11
        data = ZSortPts.l11k10;
    end
elseif n1 == 11
    if l == 0
        data = ZSortPts.l00k11;
    elseif l == 1
        data = ZSortPts.l01k11;
    elseif l == 2
        data = ZSortPts.l02k11;
    elseif l == 3
        data = ZSortPts.l03k11;
    elseif l == 4
        data = ZSortPts.l04k11;
    elseif l == 5
        data = ZSortPts.l05k11;
    elseif l == 6
        data = ZSortPts.l06k11;
    elseif l == 7
        data = ZSortPts.l07k11;
    elseif l == 8
        data = ZSortPts.l08k11;
    elseif l == 9
        data = ZSortPts.l09k11;
    elseif l == 10
        data = ZSortPts.l10k11;
    elseif l == 11
        data = ZSortPts.l11k11;
    end
elseif n1 == 12
    if l == 0
        data = ZSortPts.l00k12;
    elseif l == 1
        data = ZSortPts.l01k12;
    elseif l == 2
        data = ZSortPts.l02k12;
    elseif l == 3
        data = ZSortPts.l03k12;
    elseif l == 4
        data = ZSortPts.l04k12;
    elseif l == 5
        data = ZSortPts.l05k12;
    elseif l == 6
        data = ZSortPts.l06k12;
    elseif l == 7
        data = ZSortPts.l07k12;
    elseif l == 8
        data = ZSortPts.l08k12;
    elseif l == 9
        data = ZSortPts.l09k12;
    elseif l == 10
        data = ZSortPts.l10k12;
    elseif l == 11
        data = ZSortPts.l11k12;
    end
elseif n1 == 13
    if l == 0
        data = ZSortPts.l00k13;
    elseif l == 1
        data = ZSortPts.l01k13;
    elseif l == 2
        data = ZSortPts.l02k13;
    elseif l == 3
        data = ZSortPts.l03k13;
    elseif l == 4
        data = ZSortPts.l04k13;
    elseif l == 5
        data = ZSortPts.l05k13;
    elseif l == 6
        data = ZSortPts.l06k13;
    elseif l == 7
        data = ZSortPts.l07k13;
    elseif l == 8
        data = ZSortPts.l08k13;
    elseif l == 9
        data = ZSortPts.l09k13;
    elseif l == 10
        data = ZSortPts.l10k13;
    elseif l == 11
        data = ZSortPts.l11k13;
    end
elseif n1 == 14
    if l == 0
        data = ZSortPts.l00k14;
    elseif l == 1
        data = ZSortPts.l01k14;
    elseif l == 2
        data = ZSortPts.l02k14;
    elseif l == 3
        data = ZSortPts.l03k14;
    elseif l == 4
        data = ZSortPts.l04k14;
    elseif l == 5
        data = ZSortPts.l05k14;
    elseif l == 6
        data = ZSortPts.l06k14;
    elseif l == 7
        data = ZSortPts.l07k14;
    elseif l == 8
        data = ZSortPts.l08k14;
    elseif l == 9
        data = ZSortPts.l09k14;
    elseif l == 10
        data = ZSortPts.l10k14;
    elseif l == 11
        data = ZSortPts.l11k14;
    end
elseif n1 == 15
    if l == 0
        data = ZSortPts.l00k15;
    elseif l == 1
        data = ZSortPts.l01k15;
    elseif l == 2
        data = ZSortPts.l02k15;
    elseif l == 3
        data = ZSortPts.l03k15;
    elseif l == 4
        data = ZSortPts.l04k15;
    elseif l == 5
        data = ZSortPts.l05k15;
    elseif l == 6
        data = ZSortPts.l06k15;
    elseif l == 7
        data = ZSortPts.l07k15;
    elseif l == 8
        data = ZSortPts.l08k15;
    elseif l == 9
        data = ZSortPts.l09k15;
    elseif l == 10
        data = ZSortPts.l10k15;
    elseif l == 11
        data = ZSortPts.l11k15;
    end
elseif n1 == 16
    if l == 0
        data = ZSortPts.l00k16;
    elseif l == 1
        data = ZSortPts.l01k16;
    elseif l == 2
        data = ZSortPts.l02k16;
    elseif l == 3
        data = ZSortPts.l03k16;
    elseif l == 4
        data = ZSortPts.l04k16;
    elseif l == 5
        data = ZSortPts.l05k16;
    elseif l == 6
        data = ZSortPts.l06k16;
    elseif l == 7
        data = ZSortPts.l07k16;
    elseif l == 8
        data = ZSortPts.l08k16;
    elseif l == 9
        data = ZSortPts.l09k16;
    elseif l == 10
        data = ZSortPts.l10k16;
    elseif l == 11
        data = ZSortPts.l11k16;
    end
elseif n1 == 17
    if l == 0
        data = ZSortPts.l00k17;
    elseif l == 1
        data = ZSortPts.l01k17;
    elseif l == 2
        data = ZSortPts.l02k17;
    elseif l == 3
        data = ZSortPts.l03k17;
    elseif l == 4
        data = ZSortPts.l04k17;
    elseif l == 5
        data = ZSortPts.l05k17;
    elseif l == 6
        data = ZSortPts.l06k17;
    elseif l == 7
        data = ZSortPts.l07k17;
    elseif l == 8
        data = ZSortPts.l08k17;
    elseif l == 9
        data = ZSortPts.l09k17;
            elseif l == 10
        data = ZSortPts.l10k17;
    elseif l == 11
        data = ZSortPts.l11k17;
    end
elseif n1 == 18
    if l == 0
        data = ZSortPts.l00k18;
    elseif l == 1
        data = ZSortPts.l01k18;
    elseif l == 2
        data = ZSortPts.l02k18;
    elseif l == 3
        data = ZSortPts.l03k18;
    elseif l == 4
        data = ZSortPts.l04k18;
    elseif l == 5
        data = ZSortPts.l05k18;
    elseif l == 6
        data = ZSortPts.l06k18;
    elseif l == 7
        data = ZSortPts.l07k18;
    elseif l == 8
        data = ZSortPts.l08k18;
    elseif l == 9
        data = ZSortPts.l09k18;
            elseif l == 10
        data = ZSortPts.l10k18;
    elseif l == 11
        data = ZSortPts.l11k18;
    end
elseif n1 == 19
    if l == 0
        data = ZSortPts.l00k19;
    elseif l == 1
        data = ZSortPts.l01k19;
    elseif l == 2
        data = ZSortPts.l02k19;
    elseif l == 3
        data = ZSortPts.l03k19;
    elseif l == 4
        data = ZSortPts.l04k19;
    elseif l == 5
        data = ZSortPts.l05k19;
    elseif l == 6
        data = ZSortPts.l06k19;
    elseif l == 7
        data = ZSortPts.l07k19;
    elseif l == 8
        data = ZSortPts.l08k19;
    elseif l == 9
        data = ZSortPts.l09k19;
            elseif l == 10
        data = ZSortPts.l10k19;
    elseif l == 11
        data = ZSortPts.l11k19;
    end
elseif n1 == 20
    if l == 0
        data = ZSortPts.l00k20;
    elseif l == 1
        data = ZSortPts.l01k20;
    elseif l == 2
        data = ZSortPts.l02k20;
    elseif l == 3
        data = ZSortPts.l03k20;
    elseif l == 4
        data = ZSortPts.l04k20;
    elseif l == 5
        data = ZSortPts.l05k20;
    elseif l == 6
        data = ZSortPts.l06k20;
    elseif l == 7
        data = ZSortPts.l07k20;
    elseif l == 8
        data = ZSortPts.l08k20;
    elseif l == 9
        data = ZSortPts.l09k20;
    elseif l == 10
        data = ZSortPts.l10k20;
    elseif l == 11
        data = ZSortPts.l11k20;
    end
elseif n1 == 21
    if l == 0
        data = ZSortPts.l00k21;
    elseif l == 1
        data = ZSortPts.l01k21;
    elseif l == 2
        data = ZSortPts.l02k21;
    elseif l == 3
        data = ZSortPts.l03k21;
    elseif l == 4
        data = ZSortPts.l04k21;
    elseif l == 5
        data = ZSortPts.l05k21;
    elseif l == 6
        data = ZSortPts.l06k21;
    elseif l == 7
        data = ZSortPts.l07k21;
    elseif l == 8
        data = ZSortPts.l08k21;
    elseif l == 9
        data = ZSortPts.l09k21;
    elseif l == 10
        data = ZSortPts.l10k21;
    elseif l == 11
        data = ZSortPts.l11k21;
    end
elseif n1 == 22
    if l == 0
        data = ZSortPts.l00k22;
    elseif l == 1
        data = ZSortPts.l01k22;
    elseif l == 2
        data = ZSortPts.l02k22;
    elseif l == 3
        data = ZSortPts.l03k22;
    elseif l == 4
        data = ZSortPts.l04k22;
    elseif l == 5
        data = ZSortPts.l05k22;
    elseif l == 6
        data = ZSortPts.l06k22;
    elseif l == 7
        data = ZSortPts.l07k22;
    elseif l == 8
        data = ZSortPts.l08k22;
    elseif l == 9
        data = ZSortPts.l09k22;
    elseif l == 10
        data = ZSortPts.l10k22;
    elseif l == 11
        data = ZSortPts.l11k22;
    end
elseif n1 == 23
    if l == 0
        data = ZSortPts.l00k23;
    elseif l == 1
        data = ZSortPts.l01k23;
    elseif l == 2
        data = ZSortPts.l02k23;
    elseif l == 3
        data = ZSortPts.l03k23;
    elseif l == 4
        data = ZSortPts.l04k23;
    elseif l == 5
        data = ZSortPts.l05k23;
    elseif l == 6
        data = ZSortPts.l06k23;
    elseif l == 7
        data = ZSortPts.l07k23;
    elseif l == 8
        data = ZSortPts.l08k23;
    elseif l == 9
        data = ZSortPts.l09k23;
    elseif l == 10
        data = ZSortPts.l10k23;
    elseif l == 11
        data = ZSortPts.l11k23;
    end
elseif n1 == 24
    if l == 0
        data = ZSortPts.l00k24;
    elseif l == 1
        data = ZSortPts.l01k24;
    elseif l == 2
        data = ZSortPts.l02k24;
    elseif l == 3
        data = ZSortPts.l03k24;
    elseif l == 4
        data = ZSortPts.l04k24;
    elseif l == 5
        data = ZSortPts.l05k24;
    elseif l == 6
        data = ZSortPts.l06k24;
    elseif l == 7
        data = ZSortPts.l07k24;
    elseif l == 8
        data = ZSortPts.l08k24;
    elseif l == 9
        data = ZSortPts.l09k24;
    elseif l == 10
        data = ZSortPts.l10k24;
    elseif l == 11
        data = ZSortPts.l11k24;
    end
elseif n1 == 25
    if l == 0
        data = ZSortPts.l00k25;
    elseif l == 1
        data = ZSortPts.l01k25;
    elseif l == 2
        data = ZSortPts.l02k25;
    elseif l == 3
        data = ZSortPts.l03k25;
    elseif l == 4
        data = ZSortPts.l04k25;
    elseif l == 5
        data = ZSortPts.l05k25;
    elseif l == 6
        data = ZSortPts.l06k25;
    elseif l == 7
        data = ZSortPts.l07k25;
    elseif l == 8
        data = ZSortPts.l08k25;
    elseif l == 9
        data = ZSortPts.l09k25;
    elseif l == 10
        data = ZSortPts.l10k25;
    elseif l == 11
        data = ZSortPts.l11k25;
    end
elseif n1 == 26
    if l == 0
        data = ZSortPts.l00k26;
    elseif l == 1
        data = ZSortPts.l01k26;
    elseif l == 2
        data = ZSortPts.l02k26;
    elseif l == 3
        data = ZSortPts.l03k26;
    elseif l == 4
        data = ZSortPts.l04k26;
    elseif l == 5
        data = ZSortPts.l05k26;
    elseif l == 6
        data = ZSortPts.l06k26;
    elseif l == 7
        data = ZSortPts.l07k26;
    elseif l == 8
        data = ZSortPts.l08k26;
    elseif l == 9
        data = ZSortPts.l09k26;
    elseif l == 10
        data = ZSortPts.l10k26;
    elseif l == 11
        data = ZSortPts.l11k26;
    end
elseif n1 == 27
    if l == 0
        data = ZSortPts.l00k27;
    elseif l == 1
        data = ZSortPts.l01k27;
    elseif l == 2
        data = ZSortPts.l02k27;
    elseif l == 3
        data = ZSortPts.l03k27;
    elseif l == 4
        data = ZSortPts.l04k27;
    elseif l == 5
        data = ZSortPts.l05k27;
    elseif l == 6
        data = ZSortPts.l06k27;
    elseif l == 7
        data = ZSortPts.l07k27;
    elseif l == 8
        data = ZSortPts.l08k27;
    elseif l == 9
        data = ZSortPts.l09k27;
    elseif l == 10
        data = ZSortPts.l10k27;
    elseif l == 11
        data = ZSortPts.l11k27;
    end
elseif n1 == 28
    if l == 0
        data = ZSortPts.l00k28;
    elseif l == 1
        data = ZSortPts.l01k28;
    elseif l == 2
        data = ZSortPts.l02k28;
    elseif l == 3
        data = ZSortPts.l03k28;
    elseif l == 4
        data = ZSortPts.l04k28;
    elseif l == 5
        data = ZSortPts.l05k28;
    elseif l == 6
        data = ZSortPts.l06k28;
    elseif l == 7
        data = ZSortPts.l07k28;
    elseif l == 8
        data = ZSortPts.l08k28;
    elseif l == 9
        data = ZSortPts.l09k28;
    elseif l == 10
        data = ZSortPts.l10k28;
    elseif l == 11
        data = ZSortPts.l11k28;
    end
elseif n1 == 29
    if l == 0
        data = ZSortPts.l00k29;
    elseif l == 1
        data = ZSortPts.l01k29;
    elseif l == 2
        data = ZSortPts.l02k29;
    elseif l == 3
        data = ZSortPts.l03k29;
    elseif l == 4
        data = ZSortPts.l04k29;
    elseif l == 5
        data = ZSortPts.l05k29;
    elseif l == 6
        data = ZSortPts.l06k29;
    elseif l == 7
        data = ZSortPts.l07k29;
    elseif l == 8
        data = ZSortPts.l08k29;
    elseif l == 9
        data = ZSortPts.l09k29;
    elseif l == 10
        data = ZSortPts.l10k29;
    elseif l == 11
        data = ZSortPts.l11k29;
    end
elseif n1 == 30
    if l == 0
        data = ZSortPts.l00k30;
    elseif l == 1
        data = ZSortPts.l01k30;
    elseif l == 2
        data = ZSortPts.l02k30;
    elseif l == 3
        data = ZSortPts.l03k30;
    elseif l == 4
        data = ZSortPts.l04k30;
    elseif l == 5
        data = ZSortPts.l05k30;
    elseif l == 6
        data = ZSortPts.l06k30;
    elseif l == 7
        data = ZSortPts.l07k30;
    elseif l == 8
        data = ZSortPts.l08k30;
    elseif l == 9
        data = ZSortPts.l09k30;
    elseif l == 10
        data = ZSortPts.l10k30;
    elseif l == 11
        data = ZSortPts.l11k30;
    end
elseif n1 == 31
    if l == 0
        data = ZSortPts.l00k31;
    elseif l == 1
        data = ZSortPts.l01k31;
    elseif l == 2
        data = ZSortPts.l02k31;
    elseif l == 3
        data = ZSortPts.l03k31;
    elseif l == 4
        data = ZSortPts.l04k31;
    elseif l == 5
        data = ZSortPts.l05k31;
    elseif l == 6
        data = ZSortPts.l06k31;
    elseif l == 7
        data = ZSortPts.l07k31;
    elseif l == 8
        data = ZSortPts.l08k31;
    elseif l == 9
        data = ZSortPts.l09k31;
    elseif l == 10
        data = ZSortPts.l10k31;
    elseif l == 11
        data = ZSortPts.l11k31;
    end
elseif n1 == 32
    if l == 0
        data = ZSortPts.l00k32;
    elseif l == 1
        data = ZSortPts.l01k32;
    elseif l == 2
        data = ZSortPts.l02k32;
    elseif l == 3
        data = ZSortPts.l03k32;
    elseif l == 4
        data = ZSortPts.l04k32;
    elseif l == 5
        data = ZSortPts.l05k32;
    elseif l == 6
        data = ZSortPts.l06k32;
    elseif l == 7
        data = ZSortPts.l07k32;
    elseif l == 8
        data = ZSortPts.l08k32;
    elseif l == 9
        data = ZSortPts.l09k32;
    elseif l == 10
        data = ZSortPts.l10k32;
    elseif l == 11
        data = ZSortPts.l11k32;
    end
    
    elseif n1 == 33
    if l == 0
        data = ZSortPts.l00k33;
    elseif l == 1
        data = ZSortPts.l01k33;
    elseif l == 2
        data = ZSortPts.l02k33;
    elseif l == 3
        data = ZSortPts.l03k33;
    elseif l == 4
        data = ZSortPts.l04k33;
    elseif l == 5
        data = ZSortPts.l05k33;
    elseif l == 6
        data = ZSortPts.l06k33;
    elseif l == 7
        data = ZSortPts.l07k33;
    elseif l == 8
        data = ZSortPts.l08k33;
    elseif l == 9
        data = ZSortPts.l09k33;
    elseif l == 10
        data = ZSortPts.l10k33;
    elseif l == 11
        data = ZSortPts.l11k33;
    end
elseif n1 == 34
    if l == 0
        data = ZSortPts.l00k34;
    elseif l == 1
        data = ZSortPts.l01k34;
    elseif l == 2
        data = ZSortPts.l02k34;
    elseif l == 3
        data = ZSortPts.l03k34;
    elseif l == 4
        data = ZSortPts.l04k34;
    elseif l == 5
        data = ZSortPts.l05k34;
    elseif l == 6
        data = ZSortPts.l06k34;
    elseif l == 7
        data = ZSortPts.l07k34;
    elseif l == 8
        data = ZSortPts.l08k34;
    elseif l == 9
        data = ZSortPts.l09k34;
    elseif l == 10
        data = ZSortPts.l10k34;
    elseif l == 11
        data = ZSortPts.l11k34;
    end
elseif n1 == 35
    if l == 0
        data = ZSortPts.l00k35;
    elseif l == 1
        data = ZSortPts.l01k35;
    elseif l == 2
        data = ZSortPts.l02k35;
    elseif l == 3
        data = ZSortPts.l03k35;
    elseif l == 4
        data = ZSortPts.l04k35;
    elseif l == 5
        data = ZSortPts.l05k35;
    elseif l == 6
        data = ZSortPts.l06k35;
    elseif l == 7
        data = ZSortPts.l07k35;
    elseif l == 8
        data = ZSortPts.l08k35;
    elseif l == 9
        data = ZSortPts.l09k35;
    elseif l == 10
        data = ZSortPts.l10k35;
    elseif l == 11
        data = ZSortPts.l11k35;
    end
elseif n1 == 36
    if l == 0
        data = ZSortPts.l00k36;
    elseif l == 1
        data = ZSortPts.l01k36;
    elseif l == 2
        data = ZSortPts.l02k36;
    elseif l == 3
        data = ZSortPts.l03k36;
    elseif l == 4
        data = ZSortPts.l04k36;
    elseif l == 5
        data = ZSortPts.l05k36;
    elseif l == 6
        data = ZSortPts.l06k36;
    elseif l == 7
        data = ZSortPts.l07k36;
    elseif l == 8
        data = ZSortPts.l08k36;
    elseif l == 9
        data = ZSortPts.l09k36;
    elseif l == 10
        data = ZSortPts.l10k36;
    elseif l == 11
        data = ZSortPts.l11k36;
    end
elseif n1 == 37
    if l == 0
        data = ZSortPts.l00k37;
    elseif l == 1
        data = ZSortPts.l01k37;
    elseif l == 2
        data = ZSortPts.l02k37;
    elseif l == 3
        data = ZSortPts.l03k37;
    elseif l == 4
        data = ZSortPts.l04k37;
    elseif l == 5
        data = ZSortPts.l05k37;
    elseif l == 6
        data = ZSortPts.l06k37;
    elseif l == 7
        data = ZSortPts.l07k37;
    elseif l == 8
        data = ZSortPts.l08k37;
    elseif l == 9
        data = ZSortPts.l09k37;
    elseif l == 10
        data = ZSortPts.l10k37;
    elseif l == 11
        data = ZSortPts.l11k37;
    end
elseif n1 == 38
    if l == 0
        data = ZSortPts.l00k38;
    elseif l == 1
        data = ZSortPts.l01k38;
    elseif l == 2
        data = ZSortPts.l02k38;
    elseif l == 3
        data = ZSortPts.l03k38;
    elseif l == 4
        data = ZSortPts.l04k38;
    elseif l == 5
        data = ZSortPts.l05k38;
    elseif l == 6
        data = ZSortPts.l06k38;
    elseif l == 7
        data = ZSortPts.l07k38;
    elseif l == 8
        data = ZSortPts.l08k38;
    elseif l == 9
        data = ZSortPts.l09k38;
    elseif l == 10
        data = ZSortPts.l10k38;
    elseif l == 11
        data = ZSortPts.l11k38;
    end
elseif n1 == 39
    if l == 0
        data = ZSortPts.l00k39;
    elseif l == 1
        data = ZSortPts.l01k39;
    elseif l == 2
        data = ZSortPts.l02k39;
    elseif l == 3
        data = ZSortPts.l03k39;
    elseif l == 4
        data = ZSortPts.l04k39;
    elseif l == 5
        data = ZSortPts.l05k39;
    elseif l == 6
        data = ZSortPts.l06k39;
    elseif l == 7
        data = ZSortPts.l07k39;
    elseif l == 8
        data = ZSortPts.l08k39;
    elseif l == 9
        data = ZSortPts.l09k39;
    elseif l == 10
        data = ZSortPts.l10k39;
    elseif l == 11
        data = ZSortPts.l11k39;
    end
elseif n1 == 40
    if l == 0
        data = ZSortPts.l00k40;
    elseif l == 1
        data = ZSortPts.l01k40;
    elseif l == 2
        data = ZSortPts.l02k40;
    elseif l == 3
        data = ZSortPts.l03k40;
    elseif l == 4
        data = ZSortPts.l04k40;
    elseif l == 5
        data = ZSortPts.l05k40;
    elseif l == 6
        data = ZSortPts.l06k40;
    elseif l == 7
        data = ZSortPts.l07k40;
    elseif l == 8
        data = ZSortPts.l08k40;
    elseif l == 9
        data = ZSortPts.l09k40;
    elseif l == 10
        data = ZSortPts.l10k40;
    elseif l == 11
        data = ZSortPts.l11k40;
    end
elseif n1 == 41
    if l == 0
        data = ZSortPts.l00k41;
    elseif l == 1
        data = ZSortPts.l01k41;
    elseif l == 2
        data = ZSortPts.l02k41;
    elseif l == 3
        data = ZSortPts.l03k41;
    elseif l == 4
        data = ZSortPts.l04k41;
    elseif l == 5
        data = ZSortPts.l05k41;
    elseif l == 6
        data = ZSortPts.l06k41;
    elseif l == 7
        data = ZSortPts.l07k41;
    elseif l == 8
        data = ZSortPts.l08k41;
    elseif l == 9
        data = ZSortPts.l09k41;
    elseif l == 10
        data = ZSortPts.l10k41;
    elseif l == 11
        data = ZSortPts.l11k41;
    end
elseif n1 == 42
    if l == 0
        data = ZSortPts.l00k42;
    elseif l == 1
        data = ZSortPts.l01k42;
    elseif l == 2
        data = ZSortPts.l02k42;
    elseif l == 3
        data = ZSortPts.l03k42;
    elseif l == 4
        data = ZSortPts.l04k42;
    elseif l == 5
        data = ZSortPts.l05k42;
    elseif l == 6
        data = ZSortPts.l06k42;
    elseif l == 7
        data = ZSortPts.l07k42;
    elseif l == 8
        data = ZSortPts.l08k42;
    elseif l == 9
        data = ZSortPts.l09k42;
    elseif l == 10
        data = ZSortPts.l10k42;
    elseif l == 11
        data = ZSortPts.l11k42;
    end
elseif n1 == 43
    if l == 0
        data = ZSortPts.l00k43;
    elseif l == 1
        data = ZSortPts.l01k43;
    elseif l == 2
        data = ZSortPts.l02k43;
    elseif l == 3
        data = ZSortPts.l03k43;
    elseif l == 4
        data = ZSortPts.l04k43;
    elseif l == 5
        data = ZSortPts.l05k43;
    elseif l == 6
        data = ZSortPts.l06k43;
    elseif l == 7
        data = ZSortPts.l07k43;
    elseif l == 8
        data = ZSortPts.l08k43;
    elseif l == 9
        data = ZSortPts.l09k43;
    elseif l == 10
        data = ZSortPts.l10k43;
    elseif l == 11
        data = ZSortPts.l11k43;
    end
elseif n1 == 44
    if l == 0
        data = ZSortPts.l00k44;
    elseif l == 1
        data = ZSortPts.l01k44;
    elseif l == 2
        data = ZSortPts.l02k44;
    elseif l == 3
        data = ZSortPts.l03k44;
    elseif l == 4
        data = ZSortPts.l04k44;
    elseif l == 5
        data = ZSortPts.l05k44;
    elseif l == 6
        data = ZSortPts.l06k44;
    elseif l == 7
        data = ZSortPts.l07k44;
    elseif l == 8
        data = ZSortPts.l08k44;
    elseif l == 9
        data = ZSortPts.l09k44;
    elseif l == 10
        data = ZSortPts.l10k44;
    elseif l == 11
        data = ZSortPts.l11k44;
    end
elseif n1 == 45
    if l == 0
        data = ZSortPts.l00k45;
    elseif l == 1
        data = ZSortPts.l01k45;
    elseif l == 2
        data = ZSortPts.l02k45;
    elseif l == 3
        data = ZSortPts.l03k45;
    elseif l == 4
        data = ZSortPts.l04k45;
    elseif l == 5
        data = ZSortPts.l05k45;
    elseif l == 6
        data = ZSortPts.l06k45;
    elseif l == 7
        data = ZSortPts.l07k45;
    elseif l == 8
        data = ZSortPts.l08k45;
    elseif l == 9
        data = ZSortPts.l09k45;
    elseif l == 10
        data = ZSortPts.l10k45;
    elseif l == 11
        data = ZSortPts.l11k45;
    end
elseif n1 == 46
    if l == 0
        data = ZSortPts.l00k46;
    elseif l == 1
        data = ZSortPts.l01k46;
    elseif l == 2
        data = ZSortPts.l02k46;
    elseif l == 3
        data = ZSortPts.l03k46;
    elseif l == 4
        data = ZSortPts.l04k46;
    elseif l == 5
        data = ZSortPts.l05k46;
    elseif l == 6
        data = ZSortPts.l06k46;
    elseif l == 7
        data = ZSortPts.l07k46;
    elseif l == 8
        data = ZSortPts.l08k46;
    elseif l == 9
        data = ZSortPts.l09k46;
    elseif l == 10
        data = ZSortPts.l10k46;
    elseif l == 11
        data = ZSortPts.l11k46;
    end
elseif n1 == 47
    if l == 0
        data = ZSortPts.l00k47;
    elseif l == 1
        data = ZSortPts.l01k47;
    elseif l == 2
        data = ZSortPts.l02k47;
    elseif l == 3
        data = ZSortPts.l03k47;
    elseif l == 4
        data = ZSortPts.l04k47;
    elseif l == 5
        data = ZSortPts.l05k47;
    elseif l == 6
        data = ZSortPts.l06k47;
    elseif l == 7
        data = ZSortPts.l07k47;
    elseif l == 8
        data = ZSortPts.l08k47;
    elseif l == 9
        data = ZSortPts.l09k47;
    elseif l == 10
        data = ZSortPts.l10k47;
    elseif l == 11
        data = ZSortPts.l11k47;
    end
elseif n1 == 48
    if l == 0
        data = ZSortPts.l00k48;
    elseif l == 1
        data = ZSortPts.l01k48;
    elseif l == 2
        data = ZSortPts.l02k48;
    elseif l == 3
        data = ZSortPts.l03k48;
    elseif l == 4
        data = ZSortPts.l04k48;
    elseif l == 5
        data = ZSortPts.l05k48;
    elseif l == 6
        data = ZSortPts.l06k48;
    elseif l == 7
        data = ZSortPts.l07k48;
    elseif l == 8
        data = ZSortPts.l08k48;
    elseif l == 9
        data = ZSortPts.l09k48;
    elseif l == 10
        data = ZSortPts.l10k48;
    elseif l == 11
        data = ZSortPts.l11k48;
    end
elseif n1 == 49
    if l == 0
        data = ZSortPts.l00k49;
    elseif l == 1
        data = ZSortPts.l01k49;
    elseif l == 2
        data = ZSortPts.l02k49;
    elseif l == 3
        data = ZSortPts.l03k49;
    elseif l == 4
        data = ZSortPts.l04k49;
    elseif l == 5
        data = ZSortPts.l05k49;
    elseif l == 6
        data = ZSortPts.l06k49;
    elseif l == 7
        data = ZSortPts.l07k49;
    elseif l == 8
        data = ZSortPts.l08k49;
    elseif l == 9
        data = ZSortPts.l09k49;
    elseif l == 10
        data = ZSortPts.l10k49;
    elseif l == 11
        data = ZSortPts.l11k49;
    end
elseif n1 == 50
    if l == 0
        data = ZSortPts.l00k50;
    elseif l == 1
        data = ZSortPts.l01k50;
    elseif l == 2
        data = ZSortPts.l02k50;
    elseif l == 3
        data = ZSortPts.l03k50;
    elseif l == 4
        data = ZSortPts.l04k50;
    elseif l == 5
        data = ZSortPts.l05k50;
    elseif l == 6
        data = ZSortPts.l06k50;
    elseif l == 7
        data = ZSortPts.l07k50;
    elseif l == 8
        data = ZSortPts.l08k50;
    elseif l == 9
        data = ZSortPts.l09k50;
    elseif l == 10
        data = ZSortPts.l10k50;
    elseif l == 11
        data = ZSortPts.l11k50;
    end
    
end


end

