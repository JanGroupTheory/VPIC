function [ ZSortPts ] = ZSortCreate( )
%ZSortCreate Creates a struct with Ndivk*Ndivl subsets
%   By Emily Lichko

 ZSortPts=struct('l00k01', [ ], 'l00k02', [ ], 'l00k03', [ ], 'l00k04', [ ],...
        'l00k05', [ ], 'l00k06', [ ], 'l00k07', [ ], 'l00k08', [ ], 'l00k09', [ ], ...
        'l00k10', [ ], 'l00k11', [ ], 'l00k12', [ ], 'l00k13', [ ], 'l00k14', [ ],...
        'l00k15', [ ], 'l00k16', [ ], 'l00k17', [ ], 'l00k18', [ ], 'l00k19', [ ], ...
        'l00k20', [ ], 'l00k21', [ ], 'l00k22', [ ], 'l00k23', [ ], 'l00k24', [ ],...
        'l00k25', [ ], 'l00k26', [ ], 'l00k27', [ ], 'l00k28', [ ], 'l00k29', [ ], ...
        'l00k30', [ ], 'l00k31', [ ], 'l00k32', [ ],...    
        'l01k01', [ ], 'l01k02', [ ], 'l01k03', [ ], 'l01k04', [ ],...
        'l01k05', [ ], 'l01k06', [ ], 'l01k07', [ ], 'l01k08', [ ], 'l01k09', [ ], ...
        'l01k10', [ ], 'l01k11', [ ], 'l01k12', [ ], 'l01k13', [ ], 'l01k14', [ ],...
        'l01k15', [ ], 'l01k16', [ ], 'l01k17', [ ], 'l01k18', [ ], 'l01k19', [ ], ...
        'l01k20', [ ], 'l01k21', [ ], 'l01k22', [ ], 'l01k23', [ ], 'l01k24', [ ],...
        'l01k25', [ ], 'l01k26', [ ], 'l01k27', [ ], 'l01k28', [ ], 'l01k29', [ ], ...
        'l01k30', [ ], 'l01k31', [ ], 'l01k32', [ ],...
        'l02k01', [ ], 'l02k02', [ ], 'l02k03', [ ], 'l02k04', [ ],...
        'l02k05', [ ], 'l02k06', [ ], 'l02k07', [ ], 'l02k08', [ ], 'l02k09', [ ], ...
        'l02k10', [ ], 'l02k11', [ ], 'l02k12', [ ], 'l02k13', [ ], 'l02k14', [ ],...
        'l02k15', [ ], 'l02k16', [ ], 'l02k17', [ ], 'l02k18', [ ], 'l02k19', [ ], ...
        'l02k20', [ ], 'l02k21', [ ], 'l02k22', [ ], 'l02k23', [ ], 'l02k24', [ ],...
        'l02k25', [ ], 'l02k26', [ ], 'l02k27', [ ], 'l02k28', [ ], 'l02k29', [ ], ...
        'l02k30', [ ], 'l02k31', [ ], 'l02k32', [ ],...
        'l03k01', [ ], 'l03k02', [ ], 'l03k03', [ ], 'l03k04', [ ],...
        'l03k05', [ ], 'l03k06', [ ], 'l03k07', [ ], 'l03k08', [ ], 'l03k09', [ ], ...
        'l03k10', [ ], 'l03k11', [ ], 'l03k12', [ ], 'l03k13', [ ], 'l03k14', [ ],...
        'l03k15', [ ], 'l03k16', [ ], 'l03k17', [ ], 'l03k18', [ ], 'l03k19', [ ], ...
        'l03k20', [ ], 'l03k21', [ ], 'l03k22', [ ], 'l03k23', [ ], 'l03k24', [ ],...
        'l03k25', [ ], 'l03k26', [ ], 'l03k27', [ ], 'l03k28', [ ], 'l03k29', [ ], ...
        'l03k30', [ ], 'l03k31', [ ], 'l03k32', [ ],...
        'l04k01', [ ], 'l04k02', [ ], 'l04k03', [ ], 'l04k04', [ ],...
        'l04k05', [ ], 'l04k06', [ ], 'l04k07', [ ], 'l04k08', [ ], 'l04k09', [ ], ...
        'l04k10', [ ], 'l04k11', [ ], 'l04k12', [ ], 'l04k13', [ ], 'l04k14', [ ],...
        'l04k15', [ ], 'l04k16', [ ], 'l04k17', [ ], 'l04k18', [ ], 'l04k19', [ ], ...
        'l04k20', [ ], 'l04k21', [ ], 'l04k22', [ ], 'l04k23', [ ], 'l04k24', [ ],...
        'l04k25', [ ], 'l04k26', [ ], 'l04k27', [ ], 'l04k28', [ ], 'l04k29', [ ], ...
        'l04k30', [ ], 'l04k31', [ ], 'l04k32', [ ]); %,...
    
    
%         'l05k01', [ ], 'l05k02', [ ], 'l05k03', [ ], 'l05k04', [ ],...
%         'l05k05', [ ], 'l05k06', [ ],'l05k07', [ ], 'l05k08', [ ],'l05k09', [ ], ...
%         'l05k10', [ ], 'l05k11', [ ], 'l05k12', [ ],'l05k13', [ ], 'l05k14', [ ],...
%         'l05k15', [ ], 'l05k16', [ ],'l05k17', [ ], 'l05k18', [ ],'l05k19', [ ], ...
%         'l05k20', [ ],...
%         'l06k01', [ ], 'l06k02', [ ], 'l06k03', [ ], 'l06k04', [ ],...
%         'l06k05', [ ], 'l06k06', [ ],'l06k07', [ ], 'l06k08', [ ],'l06k09', [ ], ...
%         'l06k10', [ ], 'l06k11', [ ], 'l06k12', [ ],'l06k13', [ ], 'l06k14', [ ],...
%         'l06k15', [ ], 'l06k16', [ ],'l06k17', [ ], 'l06k18', [ ],'l06k19', [ ], ...
%         'l06k20', [ ],...
%         'l07k01', [ ], 'l07k02', [ ], 'l07k03', [ ], 'l07k04', [ ],...
%         'l07k05', [ ], 'l07k06', [ ],'l07k07', [ ], 'l07k08', [ ],'l07k09', [ ], ...
%         'l07k10', [ ], 'l07k11', [ ], 'l07k12', [ ],'l07k13', [ ], 'l07k14', [ ],...
%         'l07k15', [ ], 'l07k16', [ ],'l07k17', [ ], 'l07k18', [ ],'l07k19', [ ], ...
%         'l07k20', [ ],...
%         'l08k01', [ ], 'l08k02', [ ], 'l08k03', [ ], 'l08k04', [ ],...
%         'l08k05', [ ], 'l08k06', [ ],'l08k07', [ ], 'l08k08', [ ],'l08k09', [ ], ...
%         'l08k10', [ ], 'l08k11', [ ], 'l08k12', [ ],'l08k13', [ ], 'l08k14', [ ],...
%         'l08k15', [ ], 'l08k16', [ ],'l08k17', [ ], 'l08k18', [ ],'l08k19', [ ], ...
%         'l08k20', [ ],...
%         'l09k01', [ ], 'l09k02', [ ], 'l09k03', [ ], 'l09k04', [ ],...
%         'l09k05', [ ], 'l09k06', [ ],'l09k07', [ ], 'l09k08', [ ],'l09k09', [ ], ...
%         'l09k10', [ ], 'l09k11', [ ], 'l09k12', [ ],'l09k13', [ ], 'l09k14', [ ],...
%         'l09k15', [ ], 'l09k16', [ ],'l09k17', [ ], 'l09k18', [ ],'l09k19', [ ], ...
%         'l09k20', [ ]);
end

