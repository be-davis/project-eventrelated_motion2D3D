function [sub_dirs] = event_related_pphys(run_exp)
    close all;
    if nargin < 1
        run_exp = 3;
    else
    end
    addpath(genpath('/Volumes/Denali_DATA1/Brandon/code/git/mrC'));
    exclude_subs = {'empty'};
    if run_exp == 1
        top_path = '/Volumes/Denali_DATA1/Brandon/eventrelated_motion2D3D/pphys_v1';
        exp_time = 2; % two second experiment
    elseif run_exp == 2
        top_path = '/Volumes/Denali_DATA1/Brandon/eventrelated_motion2D3D/pphys_v2';
        exp_time = 4; % four second experiment
    elseif run_exp == 3
        top_path = '/Volumes/Denali_DATA1/Brandon/eventrelated_motion2D3D/EEG_exp1';
        exp_time = 4; % four second experiment
        exclude_subs = {'20180713_nl-0014','20180719_nl-0045','20180720_nl-0037','20180720_nl-0037_DELETE' };
    else
        msg = sprintf('\n unknown experiment: %d',run_exp);
        error(msg);
    end
    lower_cutoff = exp_time/2+.2; % lowest allowed RT, in seconds
    %% PREPARE ANALYSIS
    cond_names = {'right','left','near','far'};
    if run_exp == 3
        sub_dirs = subfolders([top_path,'/2018*'],1);% list of subject directories
    else
        sub_dirs = subfolders([top_path,'/nl-*'],1);
    end
    % get rid of excluded subjects
    exclude_subs = cellfun(@(x) fullfile(top_path,x), exclude_subs,'uni',false);
    sub_dirs = sub_dirs(~ismember(sub_dirs,exclude_subs));
    
    num_subs = length(sub_dirs);
    %% RUN ANALYSIS
    total_trials = [];
    for s = 1:length(sub_dirs) % loop over subject directories
        % check if subject is in exclude_subs
        % organize behavioral data for pphys and EEG_exp1 experiments
        if run_exp == 3
            cur_dir = subfolders(sub_dirs{s},1);
            cur_dir = cur_dir{2};
            cur_dir = [cur_dir, '/Exp_MATL_HCN_128_Avg'];
            cur_files = subfiles([cur_dir, '/RT*'],1);
            axx_files = subfiles([cur_dir, '/Axx*trials.mat'],1);
        else
            cur_dir = [sub_dirs{s},'/Exp_MATL'];
            cur_files = subfiles([cur_dir,'/RT*'],1);
        end
        
        cur_rt = []; cur_resp = []; cur_cond = [];
        % loop over files/blocks, and concanetate trial data
        for f = 1:length(cur_files)
            cur_data = load(cur_files{f});
            if ~isempty(cur_data.TimeLine)
                cur_rt = cat(1,cur_rt,cur_data.TimeLine.respTimeSec);
                cur_resp = cat(1,cur_resp,cat(2,{cur_data.TimeLine.respString})');
                cur_cond = cat(1,cur_cond, cur_data.TimeLine.cndNmb);
            else
            end
        end
        cur_rt = (cur_rt-exp_time/2) * 1000;
        % converting string response values to numbers
        [~,cur_resp] = ismember(cur_resp,{'Ra','La','Da','Ua'});
        % count missed responses
        proportion_misses(s) = size(find(cur_resp == 0),1) / size(cur_resp,1) * 100;
        % gather behavioral data into cell variable, on cell per subject
        beh_data{s}(:,1) = cur_cond;
        beh_data{s}(:,2) = cur_rt;
        beh_data{s}(:,3) = cur_resp;
        if run_exp == 3
            % load in eeg
            for c = 1:length(axx_files)
                cur_axx = load(axx_files{c});
                eeg_raw{s}(:,:,cur_cond==c) = cur_axx.Wave;
            end
        else
        end
    end
    beh_data = cellfun(@(x) beh_preproc(x,lower_cutoff), beh_data,'uni',false);
    [rt_mean,p_correct,trial_rts,conf_mat] = cellfun(@(x) beh_average(x), beh_data,'uni',false);
    all_trials = cell(1,4); for s = 1:4; for c=1:4 all_trials{c} = cat(1,all_trials{c},trial_rts{s}{c}); end; end;    
    
    % PLOT EEG
    plot(mean(eeg_raw{1}(:,75,beh_data{s}(:,1)==1 & beh_data{s}(:,4)==1 & beh_data{s}(:,5)==1),3))
    hold on
    plot(mean(eeg_raw{1}(:,75,beh_data{s}(:,1)==2 & beh_data{s}(:,4)==1 & beh_data{s}(:,5)==1),3))
    plot(mean(eeg_raw{1}(:,75,beh_data{s}(:,1)==3 & beh_data{s}(:,4)==1 & beh_data{s}(:,5)==1),3))
    plot(mean(eeg_raw{1}(:,75,beh_data{s}(:,1)==4 & beh_data{s}(:,4)==1 & beh_data{s}(:,5)==1),3))
    
    %% MAKE FIGURES 
    % displaying histogram
    subplot(5,1,1);
    c1 = histogram(all_trials{1},20, 'FaceColor', 'y');
    hold on
    subplot(5,1,2)
    c2 = histogram(all_trials{2}, 20, 'FaceColor', 'm');
    hold on
    subplot(5,1,3)
    c3 = histogram (all_trials{3}, 20, 'FaceColor', 'c');
    hold on
    subplot(5,1,4)
    c4 = histogram (all_trials{4}, 20, 'FaceColor', 'r');
    hold on
end

function beh_data = beh_preproc(beh_data,lower_cutoff)
    % get correct responses
    beh_data(:,4) = beh_data(:,1) == beh_data(:,3);
    % make list of responses to keep
    beh_data(:,5) = beh_data{s}(:,3) ~= 0;
    % indicate early responses
    beh_data(:,5) = beh_data{s}(:,5) & (beh_data{s}(:,2) >=lower_cutoff);
end

function [rt_mean,percent_correct,correct_rts,conf_mat] = beh_average(beh_data)
    % get rid of missing trials
    beh_data = beh_data(beh_data(:,5)==1,:);
    % get list of conditions
    cond_list = unique(beh_data(:,1));
    % grab trials where the first column == cond_list(c)
    for c = 1:length(cond_list)
        correct_rts{c} = beh_data( beh_data(:,1)==cond_list(c) & beh_data(:,4)==1 ,2);
        % compute average RT for correct trials
        rt_mean(c) = mean( correct_rts{c} );
        % compute percent correct trials
        percent_correct(c) = sum(beh_data(beh_data(:,1)==cond_list(c),4))./length(beh_data(beh_data(:,1)==cond_list(c),4))*100;
        % compute confusion
        conf_mat(c,:) = hist(beh_data(beh_data(:,1)==cond_list(c),3),cond_list)./length(beh_data(beh_data(:,1)==cond_list(c))); 
    end
end
function [] = stim_averaging(beh_data,eeg_data)
    % STIMULUS LOCKED AVERAGING
    % RETURN TIME X ELECTRODE MATRIX X CONDITION
    mean_eeg(c) = mean(all_eeg(c),3); %average over trials of each condition
    % AVERAGE OVER TRIALS
    mean_all_eeg = mean(all_eeg, 3); %all trials
    % IGNORE MISSING AND INCORRECT TRIALS
    eeg_raw{s} = eeg_raw{s}(:,:,(beh_data{s}(:,4)==1 & beh_data{s}(:,5)==1));
    for t = 1:360 %ignore for incorrect
        if beh_data{s}(t,4)==0
            all_eeg(:,:,t)=[];
        else
            all_eeg(:,:,t) = all_eeg(:,:,t);
        end
    end
   for u = 1:size(all_eeg,3) %ignore for missed
        if beh_data{s}(t,5)==0
            all_eeg(:,:,t)=[];
        else
            all_eeg(:,:,t) = all_eeg(:,:,t);
        end
   end   
end
%function [] resp_averaging(beh_data,eeg_data)
    % RESPONSE-LOCKED AVERAGING
    
    % RETURN TIME X ELECTRODE MATRIX X CONDITION
    
    % AVERAGE OVER TRIALS 
    
    % IGNORE MISSING AND INCORRECT TRIALS
%end
    
    




