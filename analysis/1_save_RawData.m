% This code extract experimental data from orignial files (.json & .csv)
% And create a raw data file that includes all data informaiton.
% this code takes a while to complete.

clear; clc; close all;
% the code path
code_path = ('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis');
addpath(code_path);

path = ('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\data'); % the data folder path
folder = dir(path);
cd(path);

STOPSIG = table; % initialize table for data
PRAC = table;    % initialize table for prac trials

for i = 3:size(char(folder.name),1) % first two are not folders
    
        stopsig = [];
        prac = [];
        
        sub_folder = folder(i).name; % specify subject
        
        block_folder = dir(sub_folder); % specifiy task folder within participant
        cd(folder(i).name)
        study_folder = dir(block_folder(end).name);
        
        for j = 3:size(char(study_folder.name),1)
            task_folder = study_folder(j).name;
            cd([block_folder(end).name]);
            if isfolder(task_folder)
               blk_folder = dir(task_folder);
               cd(task_folder); 
               for n = 3:size(char(blk_folder.name),1)
                   files_folder = dir(blk_folder(n).name);
                   cd(blk_folder(n).name);
                   setting = jsondecode(fileread('block_settings.json')); % this file has experimental info for each block
                    % read the data files
                    csv = dir('*.*csv'); % this is the file with experimental data
                   if ~isempty(csv)
                       data = readtable(csv.name);
                       if ismember('color',data.Properties.VariableNames)
                            data = removevars(data,{'color'});
                       end
                       
                       if length(setting.name) >= 7 && strcmp(setting.name(end-7:end),'practice') % extract practice data
                            data.practice(:) = 1;
                            
                            for m = 3:size(char(files_folder.name),1)
                                trace_folder = files_folder(m).name;
                                if isfolder(trace_folder)
                                    force_folder = dir(trace_folder);
                                    cd(trace_folder);
                                    list_file = dir('*trace*'); % find out trials
                                    % this has the continuous force data for each trial; not used for this study
                                    file_num = size(list_file,1);
                                    if file_num ~= 0
                                          for nn = 1:file_num
                                              temp = jsondecode(fileread(list_file(nn).name));
                                              time = temp.forces.time;    %n*1
                                              button = temp.buttons.data; %n*2
                                              if max(temp.forces.data(:,1)) >=  max(temp.forces.data(:,2))
                                                    force = temp.forces.data(:,1); % n*1
                                              else
                                                    force = temp.forces.data(:,2); % n*1
                                              end
                                              
                                              ind_button1 = find(button(:,1) == 1);
                                              ind_button2 = find(button(:,2) == 1);
                                              ind_button = min([ind_button1;ind_button2]);
                                              %mean_force = nanmean(force(1:600));
                                              %std_force = nanstd(force(1:600));
                                              %ind_force = find(force >= mean_force + 5*std_force, 1);
                                               if length(ind_button) == 1
                                                   %tmp = (ind_button - ind_force)/1000;
                                                   button_time = (ind_button - 300)/1000;
                                               else
                                                   %tmp = -99;
                                                   button_time = -99;
                                               end
                                              data.force(nn) = {force'};
                                              data.time(nn) = {time'};
                                              data.button(nn) = button_time;
                                          end
                                    end
                                end
                            end
                            cd ..;
                            prac = [prac; data];
                       else % extract experimental data
                            
                            data.practice(:) = 0;
                            data.t_prep2 = data.t_prep;
                            data.t_prep2(data.choice > -1) = data.t_choice(data.choice > -1) - data.t_max(data.choice > -1) + data.t_prep(data.choice > -1);
                            
                            data.t_prep3 = data.t_prep; % for button click data
                            
                            data.correct_nolate = data.correct_choice;
                            data.t_prep_nolate = data.t_prep;
                            if data.initial == 0 % nogo to go
                                data.correct_nolate(data.t_choice > 0.53) = 0; % timing tolerance was 30ms; 60 ms in total
                                data.t_prep_nolate(data.correct_nolate == 1) ...
                                    = data.t_choice(data.correct_nolate == 1)...
                                    - data.t_max(data.correct_nolate == 1)...
                                    + data.t_prep(data.correct_nolate == 1);
                                
                            elseif data.initial == 1
                                data.t_prep_nolate(data.correct_nolate == 0) ...
                                    = data.t_choice(data.correct_nolate == 0)...
                                    - data.t_max(data.correct_nolate == 0)...
                                    + data.t_prep(data.correct_nolate == 0);
                            end
                            
                            
                            % using different timing tolerance to determine
                            % whether a response trial is correct or not
                            data.correct_loose20 = data.correct_choice;
                            data.t_prep_loose20 = data.t_prep;
                            if data.initial == 0 % nogo to go
                                
                                data.correct_loose20(data.t_choice > 0.55) = 0; % timing tolerance was 30ms; 60 ms in total
                                data.t_prep_loose20(data.correct_loose20 == 1) ...
                                    = data.t_choice(data.correct_loose20 == 1)...
                                    - data.t_max(data.correct_loose20 == 1)...
                                    + data.t_prep(data.correct_loose20 == 1);
                                
                            elseif data.initial == 1
                                data.t_prep_loose20(data.correct_loose20 == 0) ...
                                    = data.t_choice(data.correct_loose20 == 0)...
                                    - data.t_max(data.correct_loose20 == 0)...
                                    + data.t_prep(data.correct_loose20 == 0);
                            end
                            
                            % using different timing tolerance to determine
                            % whether a response trial is correct or not
                            data.correct_loose10 = data.correct_choice;
                            data.t_prep_loose10 = data.t_prep;
                            if data.initial == 0 % nogo to go
                                
                                data.correct_loose10(data.t_choice > 0.54) = 0; % timing tolerance was 30ms; 60 ms in total
                                data.t_prep_loose10(data.correct_loose10 == 1) ...
                                    = data.t_choice(data.correct_loose10 == 1)...
                                    - data.t_max(data.correct_loose10 == 1)...
                                    + data.t_prep(data.correct_loose10 == 1);
                                
                            elseif data.initial == 1
                                data.t_prep_loose10(data.correct_loose10 == 0) ...
                                    = data.t_choice(data.correct_loose10 == 0)...
                                    - data.t_max(data.correct_loose10 == 0)...
                                    + data.t_prep(data.correct_loose10 == 0);
                            end
                            
                            % using different timing tolerance to determine
                            % whether a response trial is correct or not
                            data.correct_loose40 = data.correct_choice;
                            data.t_prep_loose40 = data.t_prep;
                            if data.initial == 0 % nogo to go
                                
                                data.correct_loose40(data.t_choice > 0.57) = 0; % timing tolerance was 30ms; 60 ms in total
                                data.t_prep_loose40(data.correct_loose40 == 1) ...
                                    = data.t_choice(data.correct_loose40 == 1)...
                                    - data.t_max(data.correct_loose40 == 1)...
                                    + data.t_prep(data.correct_loose40 == 1);
                                
                            elseif data.initial == 1
                                data.t_prep_loose40(data.correct_loose40 == 0) ...
                                    = data.t_choice(data.correct_loose40 == 0)...
                                    - data.t_max(data.correct_loose40 == 0)...
                                    + data.t_prep(data.correct_loose40 == 0);
                            end
                            
                            if strcmp(sub_folder,'101')
                                data.correct_correct(:) = data.correct_choice(:);
                                data.correct_choice(:) = 0;
                                data.correct_choice(data.final == 0 & data.choice == -1) = 1;
                                data.correct_choice(data.final == 1 & data.choice > -1) = 1;
                            else
                                data.correct_correct(:) = data.correct_choice(:) .* data.good_timing(:);
                            end
                            
                            for m = 3:size(char(files_folder.name),1)
                                trace_folder = files_folder(m).name;
                                if strcmp(trace_folder, 'traces')
                                    force_folder = dir(trace_folder);
                                    cd(trace_folder);
                                    list_file = dir('*trace*'); % find out trials
                                    file_num = size(list_file,1);
                                    if file_num ~= 0
                                          for nn = 1:file_num
                                              force = [];
                                              time = [];
                                              file_name = ['traces_trial' num2str(nn-1) '.json'];
                                              %temp = jsondecode(fileread(list_file(nn).name));
                                              temp = jsondecode(fileread(file_name));
                                              time = temp.forces.time;    %n*1
                                              button = temp.buttons.data; %n*2
                                              if max(temp.forces.data(:,1)) >=  max(temp.forces.data(:,2))
                                                    force = temp.forces.data(:,1); % n*1
                                              else
                                                    force = temp.forces.data(:,2); % n*1
                                              end
                                              ind_button1 = find(button(:,1) == 1);
                                              ind_button2 = find(button(:,2) == 1);
                                              ind_button = min([ind_button1;ind_button2]);
                                              %mean_force = nanmean(force(1:600));
                                              %std_force = nanstd(force(1:600));
                                              %ind_force = find(force >= mean_force + 5*std_force, 1);
                                               if length(ind_button) == 1
                                                   %tmp = (ind_button - ind_force)/1000;
                                                   button_time = (ind_button - 300)/1000;
                                                   data.t_prep3(nn) = button_time - data.t_max(nn) + data.t_prep(nn);
                                               else
                                                   %tmp = -99;
                                                   button_time = -99;
                                               end
                                              data.force(nn) = {force'};
                                              data.time(nn) = {time'};
                                              data.button(nn) = button_time;
                                          end
                                    end
                                 end
                            end
                            cd ..;
                            stopsig = [stopsig; data];
                       end
                   end
                   cd ..;
               end
               cd ..; 
            end
            cd ..;
        end
        STOPSIG = [STOPSIG; stopsig];
        PRAC = [PRAC; prac];
        cd ..;
end

STOPSIG_RAW.EXP = STOPSIG;
STOPSIG_RAW.PRAC = PRAC;
%%
cd(code_path);       
datafname = ['Init_Inhb_Raw.mat'];
save(datafname, 'STOPSIG_RAW');
%save('gono_data_raw.mat', 'gono_data_raw');