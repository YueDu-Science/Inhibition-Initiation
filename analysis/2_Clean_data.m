%%
%%% this code clean up the raw data
%%% e.g., deal with bad trails such as responding way too early even before
%%% stop signal appears

clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));
%DATA = readtable('tr-training.txt');
load Init_Inhb_Raw.mat;  % D
%head(DATA,3)  % check if data look right

RAW = STOPSIG_RAW.EXP;
sub = RAW.id;
sub_name = unique(sub);
STOPSIG_CLEAN = table;

% deal with individual subject's data 
for s = 1:length(sub_name)
        data = [];
        ind_sub = RAW.id == sub_name(s);
        data = RAW(ind_sub == 1,:);

        % response before 300ms considered as immatural response
        %index = data.t_choice < 0.35 & (data.choice == 0 | data.choice == 1);
        index = data.t_choice < 0.35 & (data.choice == 0 | data.choice == 1);
        data(index,:) = [];
        STOPSIG_CLEAN = [STOPSIG_CLEAN; data];
end
STOPSIG_CLEAN.EXP = STOPSIG_CLEAN;
datafname = ['Init_Inhb_Clean.mat'];
save(datafname, 'STOPSIG_CLEAN');
