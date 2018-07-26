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
        exp_dur = 2000; % two second experiment
    elseif run_exp == 2
        top_path = '/Volumes/Denali_DATA1/Brandon/eventrelated_motion2D3D/pphys_v2';
        exp_dur = 4000; % four second experiment
    elseif run_exp == 3
        top_path = '/Volumes/Denali_DATA1/Brandon/eventrelated_motion2D3D/EEG_exp1';
        exp_dur = 4000; % four second experiment
        exclude_subs = {'20180713_nl-0014','20180720_nl-0037_DELETE'};
    else
        msg = sprintf('\n unknown experiment: %d',run_exp);
        error(msg);
    end
    %% USER VARIABLES
    lower_cutoff = 200;  % lowest allowed RT, in milliseconds
    upper_cutoff = 1800; % highest allowed RT, in milliseconds
    time_res = 420;
    comp_channel = 75;
    pre_time = 1000; post_time = 200; % pre and post response time to include in response-locking
    rca_time_roi = (1+2*time_res):(2.5*time_res); % temporal ROI for RCA
    run_analysis = false;

    %% RUN ANALYSIS
    % prepare variables
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
    if ~run_analysis
        load(sprintf('%s/analyzed_data.mat',top_path),'rca_data','beh_data','egi_data');
    else
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
            cur_rt = cur_rt*1000 - exp_dur/2; % convert to ms, and subtract half of exp_dur
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
            msg = sprintf('\n finished subject %s',sub_dirs{s});
            disp(msg);
        end
        % preprocess behavioral data
        beh_data = cellfun(@(x) beh_preproc(x,lower_cutoff,upper_cutoff), beh_data,'uni',false);
        
        % eeg data
        [egi_data.stim,stim_trials] = cellfun(@(x,y) stim_averaging(x,y), beh_data,eeg_raw,'uni',false);
        [egi_data.resp,resp_trials] = cellfun(@(x,y) resp_averaging(x,y,time_res,exp_dur,pre_time,post_time), beh_data,eeg_raw,'uni',false);
        
        %% DO RCA
        for c = 1: length(cond_names)
            new_stim(:,c) = arrayfun(@(x,y) stim_trials{x}{c}, 1:length(stim_trials),'uni',false);
            new_resp(:,c) = arrayfun(@(x,y) resp_trials{x}{c}, 1:length(resp_trials),'uni',false);
        end
        stim_trials = new_stim'; clear new_stim;
        resp_trials = new_resp'; clear new_resp;

        rca_trials = cellfun(@(x) x(rca_time_roi,:,:),stim_trials,'uni',false);

        [rca_data.Data,rca_data.W,rca_data.A,rca_data.RXX,rca_data.RYY,rca_data.RXY]=rcaRun(rca_trials,11,7);
        rca_data.stim = rcaProject(egi_data.stim_trials,rca_data.W);
        rca_data.resp = rcaProject(egi_data.resp_trials,rca_data.W);
        rca_data = rmfield(rca_data,'Data');
        
        % average cross trials
        % cross-trials
        rca_temp_stim = cellfun(@(x) nanmean(x,3), rca_data.stim,'uni',false);
        rca_temp_resp = cellfun(@(x) nanmean(x,3), rca_data.resp,'uni',false);
        rca_data = rmfield(rca_data,'stim'); rca_data = rmfield(rca_data,'resp');
        egi_temp_stim = egi_data.stim; 
        egi_temp_resp = egi_data.resp;
        egi_data = rmfield(egi_data,'stim'); egi_data = rmfield(egi_data,'resp');
        for c = 1: length(cond_names)
            rca_data.stim(:,:,:,c) = cat(3,rca_temp_stim{c,:});
            rca_data.resp(:,:,:,c) = cat(3,rca_temp_resp{c,:});
            for s = 1:length(egi_temp_stim)
                egi_data.stim(:,:,s,c) = egi_temp_stim{s}(:,:,c);
                egi_data.resp(:,:,s,c) = egi_temp_resp{s}(:,:,c);
            end
        end
        clear *temp*;
        % and then save
        save(sprintf('%s/analyzed_data.mat',top_path),'rca_data','beh_data','egi_data','-v7.3');
    end
    %% AVERAGING
    % average behavioral data
    [rt_mean,p_correct,trial_rts,conf_mat] = cellfun(@(x) beh_average(x), beh_data,'uni',false);
    % compute trials distributions
    all_rts = cell(1,4); for s = 1:4; for c=1:4 all_rts{c} = cat(1,all_rts{c},trial_rts{s}{c}); end; end;

    stim_mean = cat(2,rca_data.stim,egi_data.stim(:,comp_channel,:,:));    
    stim_mean(:,:,:,length(cond_names)+1) = stim_mean(:,:,:,3) - mean(stim_mean(:,:,:,1:2),4);
    stim_mean(:,:,:,length(cond_names)+2) = stim_mean(:,:,:,4) - mean(stim_mean(:,:,:,1:2),4);
    resp_mean = cat(2,rca_data.resp,egi_data.resp(:,comp_channel,:,:));    
    resp_mean(:,:,:,length(cond_names)+1) = resp_mean(:,:,:,3) - mean(resp_mean(:,:,:,1:2),4);
    resp_mean(:,:,:,length(cond_names)+2) = resp_mean(:,:,:,4) - mean(resp_mean(:,:,:,1:2),4);

    % cross-subjects
    plot_stim_mean = squeeze(mean(stim_mean,3));
    plot_stim_err = squeeze(std(stim_mean,0,3)./sqrt(num_subs));
    plot_resp_mean = squeeze(mean(resp_mean,3));
    plot_resp_err = squeeze(std(resp_mean,0,3)./sqrt(num_subs));
    
    %% PLOT BEHAVIOR
        
    %% MAKE FIGURES 
    % displaying histogram
    subplot(5,1,1);
    c1 = histogram(all_trials{1},20, 'FaceColor', 'y');
    xlim([0,2000]);
    hold on
    subplot(5,1,2)
    c2 = histogram(all_trials{2}, 20, 'FaceColor', 'm');
    xlim([0,2000]);
    hold on
    subplot(5,1,3)
    c3 = histogram (all_trials{3}, 20, 'FaceColor', 'c');
    xlim([0,2000]);
    hold on
    subplot(5,1,4)
    c4 = histogram (all_trials{4}, 20, 'FaceColor', 'r');
    xlim([0,2000]);
    hold on
    
    
    %% PLOT EEG
    plot_diff = false; % plot difference waveforms (true/false)
    close all;
    f_size = 12;
    l_width = 2;
    gcaOpts = {'tickdir','out','ticklength',[0.0200,0.0200],'box','off','fontsize',f_size,'fontname','Helvetica','linewidth',l_width,'clipping','on'};
    plot_comps = [1,2,3,4,size(plot_stim_mean,2)];
    for e = 1:length(plot_comps)
        for q = 1:2
            if e < length(plot_comps)
                % plot topography
                egi_h(e) = subplot(length(plot_comps),5,(3)+5*(e-1)); 
                hold on
                mrC.plotOnEgi(rca.A(:,e));
                hold off
            else
            end
            for c = 1:(length(cond_names)+2)
                if ~plot_diff && c > length(cond_names)
                    continue;
                else
                end
                if q == 1
                    title_str = 'stimulus-locked';
                    eeg_h(e) = subplot(length(plot_comps),5,(1:2)+5*(e-1));
                    y_vals = plot_stim_mean(:,plot_comps(e),c);
                    err_vals = plot_stim_err(:,plot_comps(e),c);
                    x_min = 2000;
                    x_max = 3500;
                    x_vals = 1000/time_res:1000/time_res:exp_dur;
                else
                    title_str = 'response-locked';
                    eeg_h(e) = subplot(length(plot_comps),5,(4:5)+5*(e-1));
                    y_vals = plot_resp_mean(:,plot_comps(e),c);
                    err_vals = plot_resp_err(:,plot_comps(e),c);
                    x_min = -pre_time;
                    x_max = post_time;
                    x_vals = -pre_time:1000/time_res:post_time;
                    plot(zeros(2,1),[x_min,x_max],'-k','linewidth',l_width)
                end
                eeg_pos(:,e) = get(eeg_h(e),'position');
                hold on;
                p_h(c) = plot(x_vals,y_vals,'linewidth',l_width);
                if c > length(cond_names)
                    fill([x_vals;flipud(x_vals)],[y_vals-err_vals;flipud(y_vals+err_vals)],get(p_h(c),'color'),'linestyle','none','facealpha',.25);
                else
                end
            end
            arrayfun(@(x) uistack(x), p_h,'uni',false);
            if e == 1
                title(title_str,'fontsize',f_size,'fontname','Helvetica')
            elseif e == length(plot_comps) && q == 1
                if plot_diff
                    l_h = legend(p_h,[cond_names,'near-diff','far-diff'],'location','northeast');
                else
                    l_h = legend(p_h,cond_names,'location','northeast');
                end
                legend boxoff;
                l_pos = get(l_h,'position');
                l_pos(1) = l_pos(1) + l_pos(3) * 1.5;
                l_pos(2) = l_pos(2) + l_pos(4) * 0.2;
                set(l_h,'position',l_pos);
                xlabel('time (ms)','fontsize',f_size,'fontname','Helvetica')
                ylabel('amplitude (\muV)','fontsize',f_size,'fontname','Helvetica')
            else
            end
            xlim([x_min;x_max]);
            switch e
                case 1
                    if q == 1
                        y_max = 20; y_min = -15; y_unit = 5;
                    else
                        y_max = 10; y_min = -10; y_unit = 5;
                    end
                case 2
                    if q == 1
                        y_max = 10; y_min = -10; y_unit = 5;
                    else
                        y_max = 10; y_min = -10; y_unit = 5;
                    end
                otherwise
                    if q == 1
                        y_max = 4; y_min = -6; y_unit = 2;
                    else
                        y_max = 4; y_min = -4; y_unit = 2;
                    end
            end
            
            ylim([y_min,y_max]);
            set(gca,gcaOpts{:},'ytick',(y_min:y_unit:y_max));
            hold off;
        end
    end
    set(gcf,'units','centimeters');
    fig_pos = get(gcf,'position');
    fig_pos(3) = 20;
    fig_pos(4) = 30;
    % adjust eeg pos
    for e = 1:length(egi_h)
        egi_pos = get(egi_h(e),'position');
        egi_pos(3:4) = egi_pos(3:4)*1.5;
        egi_pos(1) = egi_pos(1)-egi_pos(1)*.1;
        egi_pos(2) = eeg_pos(2,e)-eeg_pos(2,1)*.05;
        set(egi_h(e),'position',egi_pos);
    end
    set(gcf,'position',fig_pos);
    export_fig(sprintf('%s/eeg_data.pdf',top_path),'-pdf',gcf);
end

function beh_data = beh_preproc(beh_data,lower_cutoff,upper_cutoff)
    % get correct responses
    beh_data(:,4) = beh_data(:,1) == beh_data(:,3);
    % make list of responses to keep
    beh_data(:,5) = beh_data(:,3) ~= 0;
    % indicate early responses
    beh_data(:,5) = beh_data(:,5) & (beh_data(:,2) >=lower_cutoff);
    beh_data(:,5) = beh_data(:,5) & (beh_data(:,2) <=upper_cutoff);
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
function [eeg_mean,eeg_trials] = stim_averaging(beh_data,eeg_raw)
    % STIMULUS LOCKED AVERAGING
    % RETURN TIME X ELECTRODE MATRIX X CONDITION
    % get list of conditions
    cond_list = unique(beh_data(:,1));
    % grab trials where the first column == cond_list(c)
    for c = 1:length(cond_list)
        % get rid of missing and incorrect trials
        eeg_trials{c} = eeg_raw(:,:, beh_data(:,1)==c & beh_data(:,4)==1 & beh_data(:,5)==1 );
        eeg_mean(:,:,c) = squeeze(mean(eeg_trials{c},3));
    end
end
function [resp_mean,resp_trials] = resp_averaging(beh_data,eeg_raw,time_res,exp_dur,pre_time,post_time)
    % RESPONSE-LOCKED AVERAGING
    % RETURN TIME X ELECTRODE MATRIX X CONDITION
    % get list of conditions
    cond_list = unique(beh_data(:,1));
    % grab trials where the first column == cond_list(c)
    for c = 1:length(cond_list)
        % get rid of missing and incorrect trials
        trial_idx = beh_data(:,1)==c & beh_data(:,4)==1 & beh_data(:,5)==1;
        get_trials = eeg_raw(:,:, trial_idx );
        get_timing = beh_data(trial_idx,2) + exp_dur/2;
        timing_idx = round(get_timing/1000*time_res);
        pre_samples = time_res*(pre_time/1000);
        post_samples = time_res*(post_time/1000);
        resp_temp = arrayfun(@(x) ...
                get_trials(timing_idx(x)-pre_samples:timing_idx(x)+post_samples,:,x),1:length(timing_idx),'uni',false);
        resp_trials{c} = cat(3,resp_temp{:});
        resp_mean(:,:,c) = squeeze(mean(resp_trials{c},3));
    end
end
    
    




