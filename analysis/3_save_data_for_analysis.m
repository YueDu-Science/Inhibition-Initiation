%%% this code restructualize the raw data and make it ready for analysis
%%% e.g., create SAT; model and estimate the speeds etc.
%% WARNING: This code takes a long time (A FEW HOURS) to run because of bootstrapping/simulation processes
%%% the output file "Init_Inhb_Analysis.mat" is already in the github
%%% folder.
%%
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'))
load Init_Inhb_Clean.mat;  % load raw data
%head(DATA,3)  % check if data look right
x_size_gono = 0.05;  % sliding window for speed-accuracy tradeoff visualization
x_size_nogo = 0.05;        

x0 = [0.3, 0.1, 0.001, 0.99];

DATA = STOPSIG_CLEAN.EXP;
sub = DATA.id;
sub_name = unique(sub);

for s = 1:length(sub_name)
%%  get data for different trial types       
    data = [];
    ind_sub = DATA.id == sub_name(s);
    data = DATA(ind_sub == 1,:);
    D.sub_name{s} = sub_name(s); % save subject name

    % find out R-to-NR and NR-to-NR blocks and trials
    % initial & final: the stimulus state at the begining and end of the trial
    % 1: go     0: nogo

    ind_gono = (data.initial == 1)';  % column array
    ind_nogo = (data.initial == 0)';

    gono_all = data(ind_gono == 1,:); % all trials starting from go
    nogo_all = data(ind_nogo == 1,:); % all trials starting from nogo

    nogo = nogo_all(nogo_all.final == 1,:); % switch trials starting from nogo
    gono = gono_all(gono_all.final == 0,:); % switch trials startging from go

    % determine xplot
    xplot = sort(unique(gono.t_prep)); % xplot is used to constructed speed-accuracy tradeoff (for visualization)
    %xplot = 0.05:0.001:0.5 - delay;

    nono = nogo_all(nogo_all.final == 0,:); % non-switch trials starting from nogo
    gogo = gono_all(gono_all.final == 1,:); % non-switch trials starting from go
%% the same fitting proceudure applies multiple times;
%%% each time the t_prep is a bit different.
%%% how accuracy is defined is different.
%%% see paper for more details

%% first
%%% when a response is needed, it is considered being accurate when a
%%% response is generated, no matter when.
%%% e.g., a response can be generated 100 ms later than the target line,
%%% but still considered as a correct response

    %%%    fitting model to estimate speed   
    % NOGO trials
    % t_prep2 is the acutal RT; see 1_save_RawData 
    %(from stimulus onset to the time of response if there is a response)
    % otherwise, t_prep2 = t_prep, the original preset rt.
    y_time = nogo.t_prep2; 
    hit = nogo.correct_choice;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.nogo{s}  = para;
    ycdf_nogo = ycdf;

    % GONO trials
    y_time = gono.t_prep2;
    hit = gono.correct_choice;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.gono{s}  = para;
    ycdf_gono = ycdf;

    %%% sliding window to get SAT (for visualization only)
    [f N] = sliding_window(gono.t_prep2, gono.correct_choice, xplot,x_size_gono);
    p_gono = f; clear f;

    [f N] = sliding_window(nogo.t_prep2, nogo.correct_choice, xplot,x_size_nogo);
    p_nogo = f; clear f;
          
          
%% Second
%%% basically, we measured the force trace for each trial
%%% this can be used to verify our key press data.
%%% resulst from this measure is consistent with key press data
%%% This is generally the same as the First data set above.
          
    % NOGO data
    nogo.button_choice(:) = 0;
    nogo.button_choice(nogo.button ~= -99) = 1;
    y_time = nogo.t_prep3; % t_prep3 = t_prep, the preset rt
    hit = nogo.button_choice;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.nogo_button{s}  = para;
    ycdf_nogo_button = ycdf;

    % GONO data
    gono.button_choice(:) = 1;
    gono.button_choice(gono.button ~= -99) = 0;
    y_time = gono.t_prep3;
    hit = gono.button_choice;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.gono_button{s}  = para;
    ycdf_gono_button = ycdf;

    %%% sliding window to get SAT (for visualization only)
    [f N] = sliding_window(gono.t_prep3, gono.button_choice, xplot,x_size_gono);
    p_gono_button = f; clear f;
    [f N] = sliding_window(nogo.t_prep3, nogo.button_choice, xplot,x_size_nogo);
    p_nogo_button = f; clear f;

%% Third
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% button data for simulation %%%%%%%%%%%%%%

    % NOGO data
    nogo.button_choice(:) = 0;
    nogo.button_choice(nogo.button ~= -99) = 1;
    y_time = nogo.t_prep;
    hit = nogo.button_choice;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.nogo_button_sim{s} = para;
    ycdf_nogo_button_sim = ycdf;

    % input data
    gono.button_choice(:) = 1;
    gono.button_choice(gono.button ~= -99) = 0;
    y_time = gono.t_prep;
    hit = gono.button_choice;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.gono_button_sim{s}  = para;
    ycdf_gono_button_sim = ycdf;

    %%% sliding window to get SAT (for visualization only)
    [f N] = sliding_window(gono.t_prep, gono.button_choice, xplot,x_size_gono);
    p_gono_button_sim = f; clear f;
        
    [f N] = sliding_window(nogo.t_prep, nogo.button_choice, xplot,x_size_nogo);
    p_nogo_button_sim = f; clear f;
    
%% Fourth
%%% The difference between Fourth and First is the criterion to determine
%%% whether a response is correct or not
%%% Here, a response is made and only considered as correct if it is not
%%% later than the target response time; that is, the stimulus circle still
%%% overlaps with the target line when a response is made.

    % NOGO data
    y_time = nogo.t_prep_nolate; % similar to t_prep2, the difference is included in manuscript
    hit = nogo.correct_nolate;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.nogo_nolate{s}  = para;
    ycdf_nogo_nolate = ycdf;

    % GONO data
    y_time = gono.t_prep2;
    hit = gono.correct_choice;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.gono_nolate{s}  = para;
    ycdf_gono_nolate = ycdf;
    
    %%% sliding window to get SAT (for visualization only)
    [f N] = sliding_window(gono.t_prep2, gono.correct_nolate, xplot,x_size_gono);
    p_gono_nolate = f; clear f;
    [f N] = sliding_window(nogo.t_prep_nolate, nogo.correct_nolate, xplot,x_size_nogo);
    p_nogo_nolate = f; clear f;
      

%% Fifth
%%% same as no late data
%%% but use different timing tolerance to define whether a response is
%%% correct or not
%%% The original is 30 ms; here is 40 ms

    % NOGO
    y_time = nogo.t_prep_loose10;
    hit = nogo.correct_loose10;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.nogo_loose10{s}  = para;
    ycdf_nogo_loose10 = ycdf;

    % GONO
    y_time = gono.t_prep_loose10;
    hit = gono.correct_loose10;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.gono_loose10{s}  = para;
    ycdf_gono_loose10 = ycdf;

    [f N] = sliding_window(gono.t_prep_loose10, gono.correct_loose10, xplot,x_size_gono);
    p_gono_loose10 = f; clear f;
    [f N] = sliding_window(nogo.t_prep_loose10, nogo.correct_loose10, xplot,x_size_nogo);
    p_nogo_loose10 = f; clear f;

%% Sixth
%%% same as no late data
%%% but use different timing tolerance to define whether a response is
%%% correct or not
%%% The original is 30 ms; here is 50 ms

    % NOGO
    y_time = nogo.t_prep_loose20;
    hit = nogo.correct_loose20;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.nogo_loose20{s}  = para;
    ycdf_nogo_loose20 = ycdf;

    % GONO
    y_time = gono.t_prep_loose20;
    hit = gono.correct_loose20;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.gono_loose20{s}  = para;
    ycdf_gono_loose20 = ycdf;

    [f N] = sliding_window(gono.t_prep_loose20, gono.correct_loose20, xplot,x_size_gono);
    p_gono_loose20 = f; clear f;
    [f N] = sliding_window(nogo.t_prep_loose20, nogo.correct_loose20, xplot,x_size_nogo);
    p_nogo_loose20 = f; clear f;    
    
%% Seventh
%%% same as no late data
%%% but use different timing tolerance to define whether a response is
%%% correct or not
%%% The original is 30 ms; here is 70 ms

    % NOGO
    y_time = nogo.t_prep_loose40;
    hit = nogo.correct_loose40;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.nogo_loose40{s}  = para;
    ycdf_nogo_loose40 = ycdf;

    % GONO
    y_time = gono.t_prep_loose40;
    hit = gono.correct_loose40;
    [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
    model.gono_loose40{s}  = para;
    ycdf_gono_loose40 = ycdf;

    [f N] = sliding_window(gono.t_prep_loose40, gono.correct_loose40, xplot,x_size_gono);
    p_gono_loose40 = f; clear f;
    [f N] = sliding_window(nogo.t_prep_loose40, nogo.correct_loose40, xplot,x_size_nogo);
    p_nogo_loose40 = f; clear f;    
%% save all data
    % delay causes some DT to negative
    % not useful here; because all t_prep were preset.
    ind = find(xplot > 0);
    xplot = xplot(ind);

    D.p_nogo{s} = p_nogo(ind);
    D.p_gono{s} = p_gono(ind);
    D.p_nogo_button{s} = p_nogo_button(ind);
    D.p_gono_button{s} = p_gono_button(ind);
    D.p_nogo_button_sim{s} = p_nogo_button_sim(ind);
    D.p_gono_button_sim{s} = p_gono_button_sim(ind);
    D.p_nogo_nolate{s} = p_nogo_nolate(ind);
    D.p_gono_nolate{s} = p_gono_nolate(ind);
    D.p_nogo_loose20{s} = p_nogo_loose20(ind);
    D.p_gono_loose20{s} = p_gono_loose20(ind);
    D.p_nogo_loose10{s} = p_nogo_loose10(ind);
    D.p_gono_loose10{s} = p_gono_loose10(ind);
    D.p_nogo_loose40{s} = p_nogo_loose40(ind);
    D.p_gono_loose40{s} = p_gono_loose40(ind);

    D.nono{s} = nono;
    D.nogo{s} = nogo;
    D.gono{s} = gono;
    D.gogo{s} = gogo;

    D.gono_all{s} = gono_all;
    D.nogo_all{s} = nogo_all;

    D.ycdf_gono{s} = ycdf_gono(ind);
    D.ycdf_nogo{s} = ycdf_nogo(ind);
    D.ycdf_gono_button{s} = ycdf_gono_button(ind);
    D.ycdf_nogo_button{s} = ycdf_nogo_button(ind);
    D.ycdf_gono_button_sim{s} = ycdf_gono_button_sim(ind);
    D.ycdf_nogo_button_sim{s} = ycdf_nogo_button_sim(ind);
    D.ycdf_gono_nolate{s} = ycdf_gono_nolate(ind);
    D.ycdf_nogo_nolate{s} = ycdf_nogo_nolate(ind);
    D.ycdf_gono_loose20{s} = ycdf_gono_loose20(ind);
    D.ycdf_nogo_loose20{s} = ycdf_nogo_loose20(ind);
    D.ycdf_gono_loose10{s} = ycdf_gono_loose10(ind);
    D.ycdf_nogo_loose10{s} = ycdf_nogo_loose10(ind);
    D.ycdf_gono_loose40{s} = ycdf_gono_loose40(ind);
    D.ycdf_nogo_loose40{s} = ycdf_nogo_loose40(ind);
end
D.xplot= xplot;
D.x_size_nogo = x_size_nogo;
D.x_size_gono = x_size_gono;
D.model = model;

clearvars -except DATA D sub_name;

%% So far, D data file is generally for sanity check because how to define accuracy and which RT to use may affect our results
%%% In fact, none of these different things qualitatively change our
%%% results

%% Repeat above procedure to create new SAT
%%% this time, the RT for nogo trials were adjusted based on the potential
%%% response time if a response that is withheld was produced.
%%% This is the major data we used for our paper. It provides more previce
%%% estimation of the speeds.
%%% But again, our conclusion does not rely on this manipulation.
rn = rng(01092022);  % this is just for replicate our figures
                     % the random seed does not change our results.
rn = RandStream('mlfg6331_64'); 
iter = 1000;

x_size = 0.05;

for sub = 1:length(sub_name)
%%  get data for different trial types       
    data = [];
    ind_sub = DATA.id == sub_name(sub);
    data = DATA(ind_sub == 1,:);
    D.sub_name{sub} = sub_name(sub); % save subject name

    % find out R-to-NR and NR-to-NR blocks and trials
    % initial & final: the stimulus state at the begining and end of the trial
    % 1: go     0: nogo

    ind_gono = (data.initial == 1)';  % column array
    ind_nogo = (data.initial == 0)';

    gono_all = data(ind_gono == 1,:); % all trials starting from go
    nogo_all = data(ind_nogo == 1,:); % all trials starting from nogo

    nogo = nogo_all(nogo_all.final == 1,:); % switch trials starting from nogo
    gono = gono_all(gono_all.final == 0,:); % switch trials startging from go

    % determine xplot
    xplot = sort(unique(gono.t_prep)); % xplot is used to constructed speed-accuracy tradeoff (for visualization)
    %xplot = 0.05:0.001:0.5 - delay;

    nono = nogo_all(nogo_all.final == 0,:); % non-switch trials starting from nogo
    gogo = gono_all(gono_all.final == 1,:); % non-switch trials starting from go

%% first
%%% 
    % extract trials with a response 
    % to get the respone time (the time a respone was made) distriution
    gogo_resp = gogo(gogo.correct_choice == 1,:); % block started from 2 to 7 

    P_gono = [];
    P_nogo = [];
    PARA_gono = [];
    PARA_nogo = [];

    count = 0;
    
    while count < iter % just add mean to data
        % GONO
        % could create a matrix T row(gogo_resp.t_choice) x col(iter)
        % then datasample(rn, T, ,length(gono.t_prep_nolate),'Replace',false)
        % this will be much faster but computing speed is not important
        % here
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        gono.sim_t_prep = gono.t_prep_nolate;
        ind = find(gono.correct_nolate == 1);
        gono.sim_t_prep(ind) = gono.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(gono.sim_t_prep, gono.correct_nolate,xplot,x_size);
        p_gono = f; clear f;
        P_gono = [P_gono p_gono];
        
        y_time = gono.sim_t_prep;
        hit = gono.correct_nolate;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_gono = [PARA_gono para'];
        
        % NOGO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        nogo.sim_t_prep = nogo.t_prep_nolate;
        ind = find(nogo.correct_nolate == 0);
        nogo.sim_t_prep(ind) = nogo.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(nogo.sim_t_prep, nogo.correct_nolate,xplot,x_size);
        p_nogo = f; clear f;
        P_nogo = [P_nogo p_nogo];
        
        y_time = nogo.sim_t_prep;
        hit = nogo.correct_nolate;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_nogo = [PARA_nogo para'];

        count = count + 1;
    end
    Sim.p_gono = nansum(P_gono,2)/iter;
    Sim.p_nogo = nansum(P_nogo,2)/iter;

    model_gono  = nanmean(PARA_gono,2);
    para = nanmean(PARA_gono,2);
    ycdf_gono = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));  

    model_nogo  = nanmean(PARA_nogo,2);
    para = nanmean(PARA_nogo,2);  
    ycdf_nogo = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));
         
    ycdf = para(4)*normcdf(0.001:0.001:0.5,para(1),para(2))+para(3)*(1-normcdf(0.001:0.001:0.5,para(1),para(2)));
    % get t_min which is the minimum time participants started to make
    % choice instead of guessing
    % we set it as the time point where the SAT exceed the chance level by 5%
    t_min(sub) = find(ycdf >= para(3) + 0.05, 1);    
    
    t_max(sub) = find(ycdf >= para(4) - 0.05, 1);
   
    A.p_nogo_nolate{sub} = Sim.p_nogo;
    A.p_gono_nolate{sub} = Sim.p_gono;

    A.ycdf_gono_nolate{sub} = ycdf_gono;
    A.ycdf_nogo_nolate{sub} = ycdf_nogo;

    A.model_nogo_nolate{sub}= model_nogo;
    A.model_gono_nolate{sub}= model_gono;

%% Second
%%% Response accuracy tolerance loose 10
    P_gono = [];
    P_nogo = [];
    PARA_gono = [];
    PARA_nogo = [];

    count = 0;
    
    while count < iter % just add mean to data
        % GONO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        gono.sim_t_prep = gono.t_prep_loose10;
        ind = find(gono.correct_loose10 == 1);
        gono.sim_t_prep(ind) = gono.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(gono.sim_t_prep, gono.correct_loose10,xplot,x_size);
        p_gono = f; clear f;
        P_gono = [P_gono p_gono];
        
        y_time = gono.sim_t_prep;
        hit = gono.correct_loose10;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_gono = [PARA_gono para'];
        
        % NOGO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        nogo.sim_t_prep = nogo.t_prep_loose10;
        ind = find(nogo.correct_loose10 == 0);
        nogo.sim_t_prep(ind) = nogo.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(nogo.sim_t_prep, nogo.correct_loose10,xplot,x_size);
        p_nogo = f; clear f;
        P_nogo = [P_nogo p_nogo];
        
        y_time = nogo.sim_t_prep;
        hit = nogo.correct_loose10;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_nogo = [PARA_nogo para'];
        
        count = count + 1;
    end
    Sim.p_gono = nansum(P_gono,2)/iter;
    Sim.p_nogo = nansum(P_nogo,2)/iter;

    model_gono  = nanmean(PARA_gono,2);
    para = nanmean(PARA_gono,2);
    ycdf_gono = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));  

    model_nogo  = nanmean(PARA_nogo,2);
    para = nanmean(PARA_nogo,2);  
    ycdf_nogo = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));

    A.p_nogo_loose10{sub} = Sim.p_nogo;
    A.p_gono_loose10{sub} = Sim.p_gono;

    A.ycdf_gono_loose10{sub} = ycdf_gono;
    A.ycdf_nogo_loose10{sub} = ycdf_nogo;

    A.model_nogo_loose10{sub}= model_nogo;
    A.model_gono_loose10{sub}= model_gono;
    
%% Third
%%% Response accuracy tolerance loose 20
    P_gono = [];
    P_nogo = [];
    PARA_gono = [];
    PARA_nogo = [];

    count = 0;
    
    while count < iter % just add mean to data
        % GONO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        gono.sim_t_prep = gono.t_prep_loose20;
        ind = find(gono.correct_loose20 == 1);
        gono.sim_t_prep(ind) = gono.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(gono.sim_t_prep, gono.correct_loose20,xplot,x_size);
        p_gono = f; clear f;
        P_gono = [P_gono p_gono];
        
        y_time = gono.sim_t_prep;
        hit = gono.correct_loose20;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_gono = [PARA_gono para'];
        
        % NOGO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        nogo.sim_t_prep = nogo.t_prep_loose20;
        ind = find(nogo.correct_loose20 == 0);
        nogo.sim_t_prep(ind) = nogo.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(nogo.sim_t_prep, nogo.correct_loose20,xplot,x_size);
        p_nogo = f; clear f;
        P_nogo = [P_nogo p_nogo];
        
        y_time = nogo.sim_t_prep;
        hit = nogo.correct_loose20;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_nogo = [PARA_nogo para'];
        
        count = count + 1;
    end
    Sim.p_gono = nansum(P_gono,2)/iter;
    Sim.p_nogo = nansum(P_nogo,2)/iter;

    model_gono  = nanmean(PARA_gono,2);
    para = nanmean(PARA_gono,2);
    ycdf_gono = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));  

    model_nogo  = nanmean(PARA_nogo,2);
    para = nanmean(PARA_nogo,2);  
    ycdf_nogo = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));

    A.p_nogo_loose20{sub} = Sim.p_nogo;
    A.p_gono_loose20{sub} = Sim.p_gono;

    A.ycdf_gono_loose20{sub} = ycdf_gono;
    A.ycdf_nogo_loose20{sub} = ycdf_nogo;

    A.model_nogo_loose20{sub}= model_nogo;
    A.model_gono_loose20{sub}= model_gono;
    
%% Fourth
%%% Response accuracy tolerance loose 40
    P_gono = [];
    P_nogo = [];
    PARA_gono = [];
    PARA_nogo = [];

    count = 0;
    
    while count < iter % just add mean to data
        % GONO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        gono.sim_t_prep = gono.t_prep_loose40;
        ind = find(gono.correct_loose40 == 1);
        gono.sim_t_prep(ind) = gono.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(gono.sim_t_prep, gono.correct_loose40,xplot,x_size);
        p_gono = f; clear f;
        P_gono = [P_gono p_gono];
        
        y_time = gono.sim_t_prep;
        hit = gono.correct_loose40;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_gono = [PARA_gono para'];
        
        % NOGO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        nogo.sim_t_prep = nogo.t_prep_loose40;
        ind = find(nogo.correct_loose40 == 0);
        nogo.sim_t_prep(ind) = nogo.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(nogo.sim_t_prep, nogo.correct_loose40,xplot,x_size);
        p_nogo = f; clear f;
        P_nogo = [P_nogo p_nogo];
        
        y_time = nogo.sim_t_prep;
        hit = nogo.correct_loose40;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_nogo = [PARA_nogo para'];
        
        count = count + 1;
    end
    Sim.p_gono = nansum(P_gono,2)/iter;
    Sim.p_nogo = nansum(P_nogo,2)/iter;

    model_gono = nanmean(PARA_gono,2);
    para = nanmean(PARA_gono,2);
    ycdf_gono = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));  

    model_nogo = nanmean(PARA_nogo,2);
    para = nanmean(PARA_nogo,2);  
    ycdf_nogo = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));

    A.p_nogo_loose40{sub} = Sim.p_nogo;
    A.p_gono_loose40{sub} = Sim.p_gono;

    A.ycdf_gono_loose40{sub} = ycdf_gono;
    A.ycdf_nogo_loose40{sub} = ycdf_nogo;

    A.model_nogo_loose40{sub}= model_nogo;
    A.model_gono_loose40{sub}= model_gono;
    
%% Fifth
%%% Response accuracy tolerance: wide tolerance, as long as responding, it
%%% is correct; use force trace data
    P_gono = [];
    P_nogo = [];
    PARA_gono = [];
    PARA_nogo = [];

    count = 0;
    
    while count < iter % just add mean to data
        % GONO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        
        gono.button_choice(:) = 1;
        gono.button_choice(gono.button ~= -99) = 0;
        gono.sim_t_prep = gono.t_prep3;
        ind = find(gono.button_choice == 1);
        gono.sim_t_prep(ind) = gono.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(gono.sim_t_prep, gono.button_choice,xplot,x_size);
        p_gono = f; clear f;
        P_gono = [P_gono p_gono];
        
        y_time = gono.sim_t_prep;
        hit = gono.button_choice;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_gono = [PARA_gono para'];
        
        % NOGO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        nogo.button_choice(:) = 0;
        nogo.button_choice(nogo.button ~= -99) = 1;
        nogo.sim_t_prep = nogo.t_prep3;
        ind = find(nogo.button_choice == 0);
        nogo.sim_t_prep(ind) = nogo.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(nogo.sim_t_prep, nogo.button_choice,xplot,x_size);
        p_nogo = f; clear f;
        P_nogo = [P_nogo p_nogo];
        
        y_time = nogo.sim_t_prep;
        hit = nogo.button_choice;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_nogo = [PARA_nogo para'];
        
        count = count + 1;
    end
    Sim.p_gono = nansum(P_gono,2)/iter;
    Sim.p_nogo = nansum(P_nogo,2)/iter;

    model_gono  = nanmean(PARA_gono,2);
    para = nanmean(PARA_gono,2);
    ycdf_gono = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));  

    model_nogo  = nanmean(PARA_nogo,2);
    para = nanmean(PARA_nogo,2);  
    ycdf_nogo = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));

    A.p_nogo_button{sub} = Sim.p_nogo;
    A.p_gono_button{sub} = Sim.p_gono;

    A.ycdf_gono_button{sub} = ycdf_gono;
    A.ycdf_nogo_button{sub} = ycdf_nogo;

    A.model_nogo_button{sub}= model_nogo;
    A.model_gono_button{sub}= model_gono;    
    
%% Sixth
%%% Response accuracy tolerance: wide tolerance, as long as responding, it
%%% is correct; use key data
    P_gono = [];
    P_nogo = [];
    PARA_gono = [];
    PARA_nogo = [];

    count = 0;
    
    while count < iter % just add mean to data
        % GONO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        gono.sim_t_prep = gono.t_prep2;
        ind = find(gono.correct_choice == 1);
        gono.sim_t_prep(ind) = gono.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(gono.sim_t_prep, gono.correct_choice,xplot,x_size);
        p_gono = f; clear f;
        P_gono = [P_gono p_gono];
        
        y_time = gono.sim_t_prep;
        hit = gono.correct_choice;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_gono = [PARA_gono para'];
        
        % NOGO
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        nogo.sim_t_prep = nogo.t_prep2;
        ind = find(nogo.correct_choice == 0);
        nogo.sim_t_prep(ind) = nogo.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(nogo.sim_t_prep, nogo.correct_choice,xplot,x_size);
        p_nogo = f; clear f;
        P_nogo = [P_nogo p_nogo];
        
        y_time = nogo.sim_t_prep;
        hit = nogo.correct_choice;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_nogo = [PARA_nogo para'];
        
        count = count + 1;
    end
    Sim.p_gono = nansum(P_gono,2)/iter;
    Sim.p_nogo = nansum(P_nogo,2)/iter;

    model_gono = nanmean(PARA_gono,2);
    para = nanmean(PARA_gono,2);
    ycdf_gono = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));  

    model_nogo = nanmean(PARA_nogo,2);
    para = nanmean(PARA_nogo,2);  
    ycdf_nogo = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));

    A.p_nogo_key{sub} = Sim.p_nogo
    A.p_gono_key{sub} = Sim.p_gono;

    A.ycdf_gono_key{sub} = ycdf_gono;
    A.ycdf_nogo_key{sub} = ycdf_nogo;

    A.model_nogo_key{sub}= model_nogo;
    A.model_gono_key{sub}= model_gono;    
end

%% Repeat above procedure to create new SAT for each 3 blocks
%%% to examine whether the expectation of switch trials change differently
%%% between conditions
for sub = 1:length(sub_name)
% get data for different blocks
%% first
%%% first three blocks
    data = [];
    ind_sub = DATA.id == sub_name(sub);
    ind_blk = (DATA.block_key >= 2 & DATA.block_key <= 4) | (DATA.block_key >= 8 & DATA.block_key <= 10); % block from 2 to 7 and 8 to 13
    data = DATA(ind_sub == 1 & ind_blk == 1,:);
    D.sub_name{sub} = sub_name(sub); % save subject name

    % find out R-to-NR and NR-to-NR blocks and trials
    % initial & final: the stimulus state at the begining and end of the trial
    % 1: go     0: nogo

    ind_gono = (data.initial == 1)';  % column array
    ind_nogo = (data.initial == 0)';

    gono_all = data(ind_gono == 1,:); % all trials starting from go
    nogo_all = data(ind_nogo == 1,:); % all trials starting from nogo

    nogo = nogo_all(nogo_all.final == 1,:); % switch trials starting from nogo
    gono = gono_all(gono_all.final == 0,:); % switch trials startging from go

    % determine xplot
    xplot = sort(unique(gono.t_prep)); % xplot is used to constructed speed-accuracy tradeoff (for visualization)
    %xplot = 0.05:0.001:0.5 - delay;

    nono = nogo_all(nogo_all.final == 0,:); % non-switch trials starting from nogo
    gogo = gono_all(gono_all.final == 1,:); % non-switch trials starting from go


    % extract trials with a response 
    % to get the respone time (the time a respone was made) distriution
    gogo_resp = gogo(gogo.correct_choice == 1,:); % block started from 2 to 7 

    P_gono = [];
    P_nogo = [];
    PARA_gono = [];
    PARA_nogo = [];

    count = 0;
    
    while count < iter % just add mean to data
        % GONO
        % could create a matrix T row(gogo_resp.t_choice) x col(iter)
        % then datasample(rn, T, ,length(gono.t_prep_nolate),'Replace',false)
        % this will be much faster but computing speed is not important
        % here
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        gono.sim_t_prep = gono.t_prep_nolate;
        ind = find(gono.correct_nolate == 1);
        gono.sim_t_prep(ind) = gono.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(gono.sim_t_prep, gono.correct_nolate,xplot,x_size);
        p_gono = f; clear f;
        P_gono = [P_gono p_gono];
        
        y_time = gono.sim_t_prep;
        hit = gono.correct_nolate;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_gono = [PARA_gono para'];
        
        % NOGO
        tmp = datasample(rn,gogo_resp.t_choice,length(nogo.t_prep_nolate),'Replace',false) - 0.5;
        nogo.sim_t_prep = nogo.t_prep_nolate;
        ind = find(nogo.correct_nolate == 0);
        nogo.sim_t_prep(ind) = nogo.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(nogo.sim_t_prep, nogo.correct_nolate,xplot,x_size);
        p_nogo = f; clear f;
        P_nogo = [P_nogo p_nogo];
        
        y_time = nogo.sim_t_prep;
        hit = nogo.correct_nolate;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_nogo = [PARA_nogo para'];

        count = count + 1;
    end
    Sim.p_gono = nansum(P_gono,2)/iter;
    Sim.p_nogo = nansum(P_nogo,2)/iter;

    model_gono  = nanmean(PARA_gono,2);
    para = nanmean(PARA_gono,2);
    ycdf_gono = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));  

    model_nogo  = nanmean(PARA_nogo,2);
    para = nanmean(PARA_nogo,2);  
    ycdf_nogo = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));
         
    A.p_nogo_nolate_early{sub} = Sim.p_nogo;
    A.p_gono_nolate_early{sub} = Sim.p_gono;

    A.ycdf_gono_nolate_early{sub} = ycdf_gono;
    A.ycdf_nogo_nolate_early{sub} = ycdf_nogo;

    A.model_nogo_nolate_early{sub}= model_nogo;
    A.model_gono_nolate_early{sub}= model_gono;

%% second
%%% last three blocks
    data = [];
    ind_sub = DATA.id == sub_name(sub);
    ind_blk = (DATA.block_key >= 2 & DATA.block_key <= 4) | (DATA.block_key >= 8 & DATA.block_key <= 10); % block from 2 to 7 and 8 to 13
    data = DATA(ind_sub == 1 & ind_blk == 0,:);
    D.sub_name{sub} = sub_name(sub); % save subject name

    % find out R-to-NR and NR-to-NR blocks and trials
    % initial & final: the stimulus state at the begining and end of the trial
    % 1: go     0: nogo

    ind_gono = (data.initial == 1)';  % column array
    ind_nogo = (data.initial == 0)';

    gono_all = data(ind_gono == 1,:); % all trials starting from go
    nogo_all = data(ind_nogo == 1,:); % all trials starting from nogo

    nogo = nogo_all(nogo_all.final == 1,:); % switch trials starting from nogo
    gono = gono_all(gono_all.final == 0,:); % switch trials startging from go

    % determine xplot
    xplot = sort(unique(gono.t_prep)); % xplot is used to constructed speed-accuracy tradeoff (for visualization)
    %xplot = 0.05:0.001:0.5 - delay;

    nono = nogo_all(nogo_all.final == 0,:); % non-switch trials starting from nogo
    gogo = gono_all(gono_all.final == 1,:); % non-switch trials starting from go


    % extract trials with a response 
    % to get the respone time (the time a respone was made) distriution
    gogo_resp = gogo(gogo.correct_choice == 1,:); % block started from 2 to 7 

    P_gono = [];
    P_nogo = [];
    PARA_gono = [];
    PARA_nogo = [];

    count = 0;
    
    while count < iter % just add mean to data
        % GONO
        % could create a matrix T row(gogo_resp.t_choice) x col(iter)
        % then datasample(rn, T, ,length(gono.t_prep_nolate),'Replace',false)
        % this will be much faster but computing speed is not important
        % here
        tmp = datasample(rn,gogo_resp.t_choice,length(gono.t_prep_nolate),'Replace',false) - 0.5;
        
        gono.sim_t_prep = gono.t_prep_nolate;
        ind = find(gono.correct_nolate == 1);
        gono.sim_t_prep(ind) = gono.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(gono.sim_t_prep, gono.correct_nolate,xplot,x_size);
        p_gono = f; clear f;
        P_gono = [P_gono p_gono];
        
        y_time = gono.sim_t_prep;
        hit = gono.correct_nolate;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_gono = [PARA_gono para'];
        
        % NOGO
        tmp = datasample(rn,gogo_resp.t_choice,length(nogo.t_prep_nolate),'Replace',false) - 0.5;
        nogo.sim_t_prep = nogo.t_prep_nolate;
        ind = find(nogo.correct_nolate == 0);
        nogo.sim_t_prep(ind) = nogo.sim_t_prep(ind) + tmp(ind);
        
        [f N] = sliding_window(nogo.sim_t_prep, nogo.correct_nolate,xplot,x_size);
        p_nogo = f; clear f;
        P_nogo = [P_nogo p_nogo];
        
        y_time = nogo.sim_t_prep;
        hit = nogo.correct_nolate;
        [para, ycdf, fval] = fit_SAT(y_time, hit, x0, xplot);
        PARA_nogo = [PARA_nogo para'];

        count = count + 1;
    end
    Sim.p_gono = nansum(P_gono,2)/iter;
    Sim.p_nogo = nansum(P_nogo,2)/iter;

    model_gono  = nanmean(PARA_gono,2);
    para = nanmean(PARA_gono,2);
    ycdf_gono = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));  

    model_nogo  = nanmean(PARA_nogo,2);
    para = nanmean(PARA_nogo,2);  
    ycdf_nogo = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));
         
    A.p_nogo_nolate_late{sub} = Sim.p_nogo;
    A.p_gono_nolate_late{sub} = Sim.p_gono;

    A.ycdf_gono_nolate_late{sub} = ycdf_gono;
    A.ycdf_nogo_nolate_late{sub} = ycdf_nogo;

    A.model_nogo_nolate_late{sub}= model_nogo;
    A.model_gono_nolate_late{sub}= model_gono;
end
%% save all data
    % delay causes some DT to negative
    % not useful here; because all t_prep were preset.
ind = find(xplot > 0);
xplot = xplot(ind);
A.xplot = xplot;

A.t_min = t_min;
A.t_max = t_max;

Analysis.D = D;
Analysis.A = A;

%% save data file
datafname = ['Init_Inhb_Analysis.mat'];
save(datafname, 'Analysis');
