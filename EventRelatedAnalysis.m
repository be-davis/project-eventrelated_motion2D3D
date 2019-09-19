function [sub_dirs] = EventRelatedAnalysis(run_exp,run_analysis,comp_channel,time_res)
    close all;
    if nargin < 1 || isempty(run_exp)
        run_exp = 3;
    else
    end
    if nargin < 2 || isempty(run_analysis)
        run_analysis = false;
    else
    end
    if nargin < 3 || isempty(comp_channel)
        comp_channel = 75;
    else
    end
    if nargin < 4 || isempty(time_res)
        time_res = 420;
    else
    end
        
    exclude_subs = {'empty'};
    if run_exp == 1
        top_path = '/Volumes/Denali_DATA1/Brandon/eventrelated_motion2D3D/pphys_v1';
        fig_path = '/Users/kohler/Dropbox/WRITING/Articles/2019_KohlerMotion2D3D_eventrelated/figures/pphys_v1';
        exp_dur = 2000; % two second experiment
    elseif run_exp == 2
        top_path = '/Volumes/Denali_DATA1/Brandon/eventrelated_motion2D3D/pphys_v2';
        fig_path = '/Users/kohler/Dropbox/WRITING/Articles/2019_KohlerMotion2D3D_eventrelated/figures/pphys_v2';
        exp_dur = 4000; % four second experiment
    elseif run_exp == 3
        top_path = '/Volumes/Denali_DATA1/Brandon/eventrelated_motion2D3D/EEG_exp1';
        fig_path = '/Users/kohler/Dropbox/WRITING/Articles/2019_KohlerMotion2D3D_eventrelated/figures/EEG';
        exp_dur = 4000; % four second experiment
        exclude_subs = {'20180713_nl-0014','20180720_nl-0037_DELETE', '20180803_nl-1668','20180815_nl-1681','20180815_nl-1676'};
    else
        msg = sprintf('\n unknown experiment: %d',run_exp);
        error(msg);
    end
    %% USER VARIABLES
    lower_cutoff = 200;  % lowest allowed RT, in milliseconds
    upper_cutoff = 1800; % highest allowed RT, in milliseconds
    pre_time = 1800; post_time = 200; % pre and post response time to include in response-locking
    rca_time_roi = (1+2*time_res):(2.5*time_res); % temporal ROI for RCA
    num_bins = 3; split_by_cond = true; resp_rca = false;

    %% RUN ANALYSIS
    % prepare variables
    cond_names = {'right','left','towards','away'};
    
    if ~run_analysis
        fprintf('\n ... LOADING DATA ... \n');
        load(sprintf('%s/analyzed_data.mat',fig_path));
    else
        if run_exp == 3
            sub_dirs = subfolders([top_path,'/2018*'],1);% list of subject directories
        else
            sub_dirs = subfolders([top_path,'/nl-*'],1);
        end
        sub_dirs = sub_dirs(cellfun(@(x) isempty(strfind(x,'_noblink')), sub_dirs));
        % get rid of excluded subjects
        exclude_subs = cellfun(@(x) fullfile(top_path,x), exclude_subs,'uni',false);
        sub_dirs = sub_dirs(~ismember(sub_dirs,exclude_subs)); 
        sub_dirs = cellfun(@(x) [x,'_noblink'], sub_dirs,'uni',false);
        
        sub_idx = 0;
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
            rt_correction = 120;
            cur_rt = cur_rt*1000 - exp_dur/2 - rt_correction; % convert to ms, and subtract half of exp_dur
            % converting string response values to numbers
            [~,cur_resp] = ismember(cur_resp,{'Ra','La','Da','Ua'});
            % count missed responses
            proportion_misses(s) = size(find(cur_resp == 0),1) / size(cur_resp,1) * 100;
            % gather behavioral data into cell variable, one cell per subject
            beh_temp(:,1) = cur_cond;
            beh_temp(:,2) = cur_rt;
            beh_temp(:,3) = cur_resp;
            if run_exp == 3
                % load in eeg
                for c = 1:length(axx_files)
                    cur_axx = load(axx_files{c});
                    eeg_raw(:,:,cur_cond==c) = cur_axx.Wave;
                    clear cur_axx;
                end
            else
            end
            % preprocess behavioral data
            [beh_temp, proportion_rejected(s,:),proportion_incorrect(s,:) ] = beh_preproc(beh_temp,lower_cutoff,upper_cutoff, false);
            % how many trials are missing
            proportion_missing(s) = sum(squeeze(sum(all(eeg_raw == 0,1),2)) == size(eeg_raw,2))./size(eeg_raw,3);
            
            skip_sub = false;
            if proportion_missing(s) > .2
                fprintf('\n skipping subject %s, too many missing EEG trials \n',sub_dirs{s});
                skip_sub = true;
            else
            end
            if any(proportion_rejected(s,:) > 1/4, 2)
                fprintf('\n skipping subject %s, too many rejected trials \n',sub_dirs{s});
                skip_sub = true;
            else
            end
            if any(mean(proportion_incorrect(s,[1,2]),2) > .4, 2) || ...
               any(mean(proportion_incorrect(s,[3,4]),2) > .4, 2)     
                fprintf('\n skipping subject %s, too many incorrect trials \n',sub_dirs{s});
                skip_sub = true;
            else
            end
            if ~skip_sub
                % no more than 20 % of data missing from any one condition
                sub_idx = sub_idx + 1;
                beh_data{sub_idx} = beh_temp;
                
                % create shuffled data
                beh_data_shuf{sub_idx} = beh_preproc(beh_temp,lower_cutoff,upper_cutoff, true);
                
                % split on response time, add as additional column
                beh_data{sub_idx} = beh_split(beh_data{sub_idx}, num_bins, split_by_cond);
                beh_data_shuf{sub_idx} = beh_split(beh_data_shuf{sub_idx}, num_bins, split_by_cond);
                         
                % get rid of zeros in raw eeg
                missing_trials(sub_idx) = sum(squeeze(sum(all(eeg_raw == 0,1),2)) == size(eeg_raw,2));
                eeg_raw(repmat(all(eeg_raw == 0,1),size(eeg_raw,1),1,1)) = NaN;
                egi_data{sub_idx} = eeg_raw;
                
                fprintf('\n finished subject %s \n',sub_dirs{s});
            else  
            end
            clear eeg_raw; 
        end
        
        % make rca trials
        delete(gcp('nocreate'))
        for q = 1:2
            for z = 1:length(egi_data)
                cur_data = egi_data{z};
                if q == 1
                    rca_trials(:,z) = stim_averaging(beh_data{z}, cur_data(rca_time_roi,:,:));
                else
                    % create response-locked dataset for rca
                    rca_trials(:,z) =  resp_averaging(beh_data{z}, egi_data{z}, time_res, exp_dur, pre_time, post_time);
                end
            end
            % run RCA
            [rca_struct.Data,rca_struct.W,rca_struct.A,Rxx,Ryy,Rxy,dGen]=rcaRun(rca_trials,11,7);
            covData.Rxx = Rxx;
            covData.Ryy = Ryy;
            covData.Rxy = Rxy;
            covData.sortedGeneralizedEigenValues = dGen;
            rca_struct.covData = covData;
            clear rca_trials;
            rca_struct = rmfield(rca_struct,'Data');
            rca_raw = rcaProject(egi_data, rca_struct.W);
            comp_channel = 75;
            rca_raw = arrayfun(@(x) cat(2,rca_raw{x},egi_data{x}(:,comp_channel,:)), 1:length(rca_raw), 'uni', false);
            if q == 1
                stim_rca_raw = rca_raw;
                stim_rca_struct = rca_struct;
            else
                resp_rca_raw = rca_raw;
                resp_rca_struct = rca_struct;
            end
        end

        % and then save
        fprintf('\n ... SAVING DATA ... \n');
        save(sprintf('%s/analyzed_data.mat',fig_path),'stim_rca_*','resp_rca_*','beh_data*','-v7.3');
    end
    
     % GET RCA EXPLAINED
    [ rca_rel_expl(:,1), rca_var_expl(:,1) ,pca_var_expl(:,1) ] = rcaExplained(stim_rca_struct,2);
    [ rca_rel_expl(:,2), rca_var_expl(:,2) ,pca_var_expl(:,2) ] = rcaExplained(resp_rca_struct,2);
    
    if resp_rca
        rca_raw = resp_rca_raw;
        rca_struct = resp_rca_struct;
        rca_str = '_resp_rca';
    else
        rca_raw = stim_rca_raw;
        rca_struct = stim_rca_struct;
        rca_str = '';
    end
    clear stim_rca_*; clear resp_rca_*;
    
    % create various data sets
    num_subs = size(beh_data,2);
    rca_sorted = cell(size(rca_raw,1),2);
    for s = 1:length(rca_raw)
        [rca_struct.intact.stim.trials(:,s), rca_struct.intact.stim.ave(:,:,:,s)] = ...
            stim_averaging(beh_data{s}, rca_raw{s});
        [rca_struct.intact.resp.trials(:,s), rca_struct.intact.resp.ave(:,:,:,s)] = ...
            resp_averaging(beh_data{s}, rca_raw{s}, time_res, exp_dur, pre_time, post_time);
        [rca_struct.intact.stim_sort.trials(:,s), rca_struct.intact.stim_sort.ave(:,:,:,s)] = ...
            stim_averaging(beh_data{s}, rca_raw{s}, 6);
        [rca_struct.intact.resp_sort.trials(:,s), rca_struct.intact.resp_sort.ave(:,:,:,s)] = ...
            resp_averaging(beh_data{s}, rca_raw{s}, time_res, exp_dur, pre_time, post_time, 6);
        [rca_struct.shuffle.stim.trials(:,s), rca_struct.shuffle.stim.ave(:,:,:,s)] = ... 
            stim_averaging(beh_data_shuf{s}, rca_raw{s});
        [rca_struct.shuffle.resp.trials(:,s), rca_struct.shuffle.resp.ave(:,:,:,s)] = ...
            resp_averaging(beh_data_shuf{s}, rca_raw{s}, time_res, exp_dur, pre_time, post_time);
        [rca_struct.shuffle.stim_sort.trials(:,s), rca_struct.shuffle.stim_sort.ave(:,:,:,s)] = ...
            stim_averaging(beh_data_shuf{s}, rca_raw{s}, 6);
        [rca_struct.shuffle.resp_sort.trials(:,s), rca_struct.shuffle.resp_sort.ave(:,:,:,s)] = ... 
            resp_averaging(beh_data_shuf{s}, rca_raw{s}, time_res, exp_dur, pre_time, post_time, 6);
        [~, sort_idx] = sort(beh_data{s}(:,2), 'descend');
        temp_beh = beh_data{s}(sort_idx,:);
        temp_rca = rca_raw{s}(:,:,sort_idx);
        keep_trials = ~squeeze(sum(all(isnan(temp_rca),1),2)==size(temp_rca,2));
        for q = 1:2
            if q == 1
                trial_idx = any(temp_beh(:,1)==[1,2],2) & temp_beh(:,4)==1 & temp_beh(:,5)==1 & keep_trials==1;
            else
                trial_idx = any(temp_beh(:,1)==[3,4],2) & temp_beh(:,4)==1 & temp_beh(:,5)==1 & keep_trials==1;
            end
            beh_sorted{s,q} = temp_beh(trial_idx,2);
            rca_sorted_sl{s,q} = temp_rca(:,:,trial_idx);
            % do response locking
            pre_samples = time_res*(pre_time/1000);
            post_samples = time_res*(post_time/1000);
            get_timing = temp_beh(trial_idx,2) + exp_dur/2;
            timing_idx = round(get_timing/1000*time_res);
            for t = 1:length(timing_idx)
                rca_sorted_rl{s,q}(:,:,t) = rca_sorted_sl{s,q}(timing_idx(t)-pre_samples:timing_idx(t)+post_samples,:,t);
            end
        end
        clear temp_*
    end
    
    %% PLOT PARAMS
    close all;
    f_size = 12;
    l_width = 2;
    gcaOpts = {'tickdir','out','ticklength',[0.0300,0.0300],'box','off','fontsize',f_size,'fontname','Helvetica','linewidth',l_width,'clipping','off'};
    text_params = {'fontweight','normal','fontname','Helvetica','fontsize', f_size};
    cBrewer = load('colorBrewer_new.mat');
    cond_colors = [cBrewer.rgb20(3,:); cBrewer.rgb20(7,:); cBrewer.rgb20(5,:); cBrewer.rgb20(17,:)];
    im_average = false;
    plot_corr = false;
% %% PLOT IMAGES
%     % average behavioral data
%     
%     for i = 1:2
%         for rc = [1:3,8]
%             figure;
%             for q = 1:2
%                 if im_average
%                     im_reshaped = zeros(t_dim,100);
%                     for s = 1:length(rca_sorted)
%                         im_reshaped = im_reshaped + imresize(squeeze(rca_sorted{s,q}(:,1,:)),[t_dim,100]);
%                     end
%                     im_reshaped = im_reshaped'./length(rca_sorted);
%                 else
%                     if i == 2
%                         all_image = cat(3,rca_sorted_rl{:,q});
%                         im_time_roi = 1:size(all_image,1);
%                         x_vals = -pre_time:1000/time_res:post_time;
%                         x_unit = 200;
%                     else
%                         all_image = cat(3,rca_sorted_sl{:,q});
%                         im_time_roi = (1+2*time_res):(4*time_res); % temporal ROI for RCA
%                         x_vals = (1000/time_res:1000/time_res:exp_dur)-2000;
%                         x_unit = 500;
%                     end
%                     x_vals = x_vals(im_time_roi);
%                     t_dim = size(all_image,1);
%                     not_nan = all(squeeze(~isnan(all_image(:,rc,:))),1);
%                     all_image = squeeze(all_image(:,rc, not_nan));
%                     all_rt = cat(1,beh_sorted{:,q});
%                     all_rt = all_rt(not_nan);
%                     [all_sort, all_idx] = sort(all_rt, 'descend');
%                     all_image = all_image(:,all_idx);
%                     im_reshaped = imresize(all_image,[t_dim,100])';
%                     rt_reshaped = imresize(all_sort, [100,1]);
%                     [rt_corr, corr_p] = corr(rt_reshaped, abs(im_reshaped));
%                 end
% 
%                 subplot(2,1,q);
%                 im_h = imagesc(im_reshaped(:,im_time_roi),[-15 10]);
%                 hold on
% 
%                 set(im_h, 'XData', x_vals);
%                 x_lims = [floor(min(x_vals)/100)*100, ceil(max(x_vals)/100)*100];
%                 if i == 2
%                     [ax_h, rt_h, corr_h] = plotyy(gca, zeros(1,length(rt_reshaped)), 1:length(rt_reshaped), x_vals, rt_corr(im_time_roi));
%                 else
%                     [ax_h, rt_h, corr_h] = plotyy(gca, rt_reshaped, 1:length(rt_reshaped), x_vals, rt_corr(im_time_roi));
%                 end
%                 set(rt_h, 'linewidth',l_width, 'color', [1,1,1]);
%                 set(corr_h, 'linewidth',l_width, 'color', cond_colors(2,:)) 
%                 set(ax_h(1), 'xlim', [x_lims(1), x_lims(2)], 'xtick', x_lims(1):x_unit:x_lims(2), gcaOpts{:}, 'ylim', [0,100], 'ytick', 0:50:100, 'ycolor', 'k')
%                 set(ax_h(2), 'xlim', [x_lims(1), x_lims(2)], 'xtick', x_lims(1):x_unit:x_lims(2), gcaOpts{:}, 'ylim', [-.5,.5], 'ytick', -1:.5:1, 'ycolor', 'k')
%                 %plot([min(x_vals), max(x_vals)], ones(2,1) , 'k-', 'linewidth',l_width)
%                 ylabel(ax_h(1), 'trial bins (fastest \rightarrow slowest)', text_params{:});
%                 label_pos = max(get(ax_h(1), 'xlim'))+diff(get(ax_h(1), 'xlim'))*0.05;
%                 ylabel(ax_h(2), 'correlation (r)', text_params{:}, 'position', [label_pos,0,0]);
%                 if plot_corr
%                     leg_h = legend([rt_h,corr_h], {'reaction times','RT \times amplitude'});
%                     corr_str = '_corr';
%                 else
%                     set(corr_h,'visible','off')
%                     set(ax_h(2),'visible','off')
%                     leg_h = legend(rt_h, {'reaction times'});
%                     corr_str = '';
%                 end
%                 set(leg_h, 'TextColor', [1,1,1], 'box', 'off', text_params{:});
%                 if q == 1
%                     title('2D conditions', text_params{:}, 'color', cond_colors(1,:))
%                     im_pos(:,q) = get(gca, 'position');
%                     im_pos(2,q) = im_pos(2,q)+im_pos(4,q)*0.1;
%                     im_pos(3,q) = im_pos(3,q)*.9;
%                     set(gca, 'position', im_pos(:,q), 'xticklabels','');
%                 else
%                     title('3D conditions', text_params{:}, 'color', cond_colors(3,:))
%                     c_bar = colorbar(text_params{:});
%                     c_bar.Label.String = 'amplitude (\muV)';
%                     im_pos(:,q) = get(gca, 'position');
%                     im_pos([1,3,4],q) = im_pos([1,3,4],1);
%                     im_pos(2,q) = im_pos(2,1) - im_pos(4,1)*1.2;
%                     set(gca, 'position', im_pos(:,q));
%                     c_pos = get(c_bar, 'position');
%                     if plot_corr
%                         c_pos(1) = im_pos(1,q)+im_pos(3,q)*1.1;
%                     else
%                         c_pos(1) = im_pos(1,q)+im_pos(3,q)*1.02;
%                     end
%                     c_pos(2) = c_pos(2) + c_pos(4)/2; 
%                     c_pos(3) = 0.02; c_pos(4) = im_pos(4,q);
%                     set(c_bar, 'position', c_pos);
%                     xlabel('time (ms)', text_params{:});
%                 end
%                 if i == 2
%                     resp_subplot(q,rc) = gca;
%                 else
%                     stim_subplot(q,rc) = gca;
%                 end
%             end
%             set(gcf,'units','centimeters');
%             fig_pos = get(gcf,'position');
%             fig_pos(3) = 25;
%             fig_pos(4) = 15;
%             set(gcf,'position',fig_pos);
%             if im_average 
%                 %set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[fig_pos(3), fig_pos(4)])
%                 %print(gcf, '-dpdf', '-bestfit', sprintf('%s/motion2D3D_image_ave.pdf',fig_path), '-r0');
%                 if i == 2
%                     export_fig(sprintf('%s/motion2D3D_respimage_ave_rc%d_%d-%d%s.tif',fig_path, x_lims(1), x_lims(2), corr_str),'-tif','-painters','-r300', '-transparent',gcf);
%                 else
%                     export_fig(sprintf('%s/motion2D3D_stimimage_ave_rc%d_%d-%d%s.tif',fig_path, x_lims(1), x_lims(2), corr_str),'-tif','-painters','-r300', '-transparent',gcf);
%                 end
%             else
%                 %set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[fig_pos(3), fig_pos(4)])
%                 %print(gcf, '-dpdf', '-bestfit', sprintf('%s/motion2D3D_image_all.pdf',fig_path), '-r0');
%                 if i == 2
%                     export_fig(sprintf('%s/motion2D3D_respimage_all_rc%d_%d-%d%s.tif',fig_path, rc, x_lims(1), x_lims(2), corr_str),'-tif','-painters','-r300', '-transparent',gcf);
%                 else
%                     export_fig(sprintf('%s/motion2D3D_stimimage_all_rc%d_%d-%d%s.tif',fig_path, rc, x_lims(1), x_lims(2), corr_str),'-tif','-painters','-r300', '-transparent',gcf);
%                 end
%             end
%         end
%     end
    
    %% AVERAGING
    [rt_mean,p_correct,trial_rts,conf_mat] = cellfun(@(x) beh_average(x), beh_data,'uni',false);
    % compute trials distributions
    for c=1:4
        temp_rts = cellfun(@(x) x{c}, trial_rts,'uni',false);
        all_rts{c} = cat(1,temp_rts{:});
    end
    all_rts{5} = cat(1,all_rts{1:2}); all_rts{6} = cat(1,all_rts{3:4});
    
    data_types = {'stim','resp','stim_sort','resp_sort'};
    shuffle_types = {'intact','shuffle'};
    for z = 1:length(shuffle_types)
        for d = 1:length(data_types)
            %cur_data = cat(2,rca_struct.(data_types{d}) ,egi_data.(data_types{d})(:,comp_channel,:,:));
            cur_data = rca_struct.(shuffle_types{z}).(data_types{d}).ave;
            if isempty(strfind(data_types{d},'sort'))
                % average 2D and 3D conditions
                cur_data(:,:,length(cond_names)+1,:) = nanmean(cur_data(:,:,1:2,:),3);
                cur_data(:,:,length(cond_names)+2,:) = nanmean(cur_data(:,:,3:4,:),3);
                % subtract 2D conditions from each of the 3D conditions
                cur_data(:,:,length(cond_names)+3,:) = cur_data(:,:,3,:) - nanmean(cur_data(:,:,1:2,:),3);
                cur_data(:,:,length(cond_names)+4,:) = cur_data(:,:,4,:) - nanmean(cur_data(:,:,1:2,:),3);
            else
            end
            if z == 1
                % save into output struct
                intact_data(:,d).subs = cur_data;
                intact_data(:,d).ave = squeeze(nanmean(cur_data,4));
                intact_data(:,d).err = squeeze(nanstd(cur_data,0,4)./sqrt(num_subs));
            else
                % save into output struct
                shuffle_data(:,d).subs = cur_data;
                shuffle_data(:,d).ave = squeeze(nanmean(cur_data,4));
                shuffle_data(:,d).err = squeeze(nanstd(cur_data,0,4)./sqrt(num_subs));
            end
        end
    end
    
    %% do comparison
%     shuf_ready(:,:,1) = shuffle_data(:,2).subs(:,1,5,:);
%     shuf_ready(:,:,2) = shuffle_data(:,2).subs(:,1,6,:);
%     real_ready(:,:,1) = intact_data(:,2).subs(:,1,5,:);
%     real_ready(:,:,2) = intact_data(:,2).subs(:,1,6,:);
%     
%     diff_ready = real_ready-shuf_ready;
%     [sig_uncorr(:,1), p_uncorr(:,1), sig_corr(:,1)] = ttest_permute(diff_ready(:,:,1) , 10000);
%     [sig_uncorr(:,2), p_uncorr(:,2), sig_corr(:,2)] = ttest_permute(diff_ready(:,:,2) , 10000);
%     
%     % plot
%     x_min = -pre_time;
%     x_max = post_time;
%     x_vals = -pre_time:1000/time_res:post_time;
%     plot(zeros(2,1),[-3,3],'-k','linewidth',2)
%     hold on
%     plot([x_min,x_max], zeros(2,1),'-k','linewidth',2)
%     plot(x_vals,squeeze(mean(real_ready(:,:,1),2)),'r');
%     plot(x_vals,squeeze(mean(shuf_ready(:,:,1),2)),'r--');
%     plot(x_vals(sig_corr(:,1)==1),squeeze(mean(real_ready(sig_corr(:,1)==1,:,1),2)),'r','linewidth',3);
%     plot(x_vals,squeeze(mean(real_ready(:,:,2),2)),'b');
%     plot(x_vals,squeeze(mean(shuf_ready(:,:,2),2)),'b--');
%     plot(x_vals(sig_corr(:,2)==1),squeeze(mean(real_ready(sig_corr(:,2)==1,:,2),2)),'b','linewidth',3);
    
    %% PLOT BEHAVIOR
    % plot histogram
    figure;
    bin_size = 100;
    hist_ymax = 300;
    hist_unit = 50;
    bin_edges = lower_cutoff:bin_size:upper_cutoff;
    hold on
    for c = 1:length(cond_names)
        subplot(2,4,c+4);        
        hist_y = histc(all_rts{c},bin_edges);
        hist_y = hist_y(1:end-1);
        hist_x = bin_edges+bin_size/2;
        hist_x = hist_x(1:end-1);
        bar(hist_x,hist_y, 'facecolor', cond_colors((c>2)*2+1,:), 'BarWidth', .8, 'linewidth',l_width, 'edgecolor','none', 'facealpha', 1, 'linestyle','none')
        xlim([200,1800]);
        ylim([0,hist_ymax]);
        set(gca,gcaOpts{:},'ytick',0:hist_unit:hist_ymax,'xtick',0:500:2000,'xticklabels',0:500:2000)
        text(250, hist_ymax*.9, cond_names{c}, text_params{:})
        if c == 1
            xlabel('time (ms)', text_params{:});
            ylabel('# responses \times RT', text_params{:});
        else
        end
    end
    subplot(2,4,3:4);
    hold on
    mean_rt_beh = mean(cat(3,rt_mean{:}),3);
    SEM_rt_beh = std(cat(3,rt_mean{:}),[],3)./sqrt(num_subs);
    mean_corr_beh = mean(cat(3,p_correct{:}),3);
    std_corr_beh = std(cat(3,p_correct{:}),[],3);
    for b = 1:4
        bar(b, mean_rt_beh(b), 'facecolor', cond_colors((b>2)*2+1,:),'edgecolor','none','linestyle','none');
        text(b,200,[sprintf('%.01f', mean_corr_beh(b)),'\pm', sprintf('%.01f%%', std_corr_beh(b))], ... 
            text_params{:}, 'color','w','horizontalalignment','center', 'fontsize',8);
        text(b,100,'correct', text_params{:}, 'color','w','horizontalalignment','center');
    end
    errorb(1:4, mean_rt_beh, SEM_rt_beh, 'color','k','linewidth',2);
    x_min = .25; x_max = 4.75;
    xlim([x_min,x_max])
    plot([x_min, x_max],zeros(2,1),'k','linewidth', l_width);
    ylabel('reaction times', text_params{:});
    set(gca,gcaOpts{:},'ytick',0:200:1000,'xticklabel',cond_names, 'xTick', 1:numel(cond_names)) 
    hold off
    set(gcf,'units','centimeters');
    hist_pos = get(gcf,'position');
    hist_pos(3) = 30;
    hist_pos(4) = 15;
    set(gcf,'position',hist_pos);
    export_fig(sprintf('%s/motion2D3D_beh.pdf',fig_path),'-png','-r1000','-transparent',gcf);
     
%     %plot pdfs: Normalized 
%     figure,
%     for c = 1:length(cond_names)
%         subplot(1,4,c);hold on
%         [f,xi] = ksdensity(all_rts{c});
%         plot(xi,f,'color',cond_colors(c,:),'linewidth',2)
%         xlim([0 2000]);
%         ylim([0,.003]);
%         fill([xi flip(xi)],[f zeros(size(f))],cond_colors(c,:),'LineStyle','none');
%         alpha(0.25)
%         set(gca,gcaOpts{:},'xtick',0:500:2000)
%     end
%     set(gcf,'units','centimeters');
%     hist_pos = get(gcf,'position');
%     hist_pos(3) = 40;
%     hist_pos(4) = 10;
%     set(gcf,'position',hist_pos);
%     
%     
%     % plot average RT and percent correct
%     % bar plots! 
%     mean_rt_bar = mean(cat(3,rt_mean{:}),3);
%     SEM_rt_bar = std(cat(3,rt_mean{:}),[],3)./sqrt(num_subs);% Standard error of mean
%     rt_mean_temp = cat(1,rt_mean{:});
%     rt_mean2D3D = [mean(rt_mean_temp(:,1:2),2) mean(rt_mean_temp(:,3:4),2)];
%     %ttests for significance of rt data
%     [~,p_rt_comp,~,stats_rt_comp] = ttest(rt_mean2D3D(:,1),rt_mean2D3D(:,2));
% 
%     %ttests for significance of percent correct data
%     mean_corr_bar = mean(cat(3,p_correct{:}),3);
%     SEM_corr_bar = std(cat(3,p_correct{:}),[],3)./sqrt(num_subs);
%     corr_mean_temp = cat(1,p_correct{:});
%     corr_mean_2D3D = [mean(corr_mean_temp(:,1:2),2) mean(corr_mean_temp(:,3:4),2)];
%     [~,p_corr_comp,~,stats_corr_comp] = ttest(corr_mean_2D3D(:,1),corr_mean_2D3D(:,2));
%  
%     figure;
%     %bar graph of average reaction time
%     subplot(2,1,1);
%     hold on
%     for i = 1:length(cond_names)
%         h = bar(i, mean_rt_bar(i), 'facecolor', cond_colors((c>2)*2+1,:));
%     end
%     set(gca, gcaOpts{:},'xticklabel',cond_names, 'XTick', 1:numel(cond_names));
%     ylabel('Reaction time (ms)')
%    
%     errorb(1:4, mean_rt_bar, SEM_rt_bar,'.','color','k','linewidth',2);
    
%     d = 50;
%     line([1.5 3.5], [1200 1200],'color','k');
%     text(2.5, 1220, '*','fontsize',25)
%     

%     
%     %bar graph of percent correct
%     subplot(2,1,2);
%     hold on;
%     for i = 1:length(cond_names)
%         b = bar(i, mean_corr_bar(i), 'facecolor', cond_colors(i,:), 'basevalue', 0);
%     end
%     errorbar(1:4, mean_corr_bar, SEM_corr_bar,'.','color','k','linewidth',2);
%     set(gca, gcaOpts{:},'xticklabel', cond_names, 'XTick', 1:numel(cond_names));
%     
%     ylabel('Percent correct')
    
   
%     % confusion matrices, averaged over subjects
%     figure;
%     mean_conf_mat = mean(cat(3, conf_mat{:}),3);
%     imagesc(mean_conf_mat)
%     colormap(flipud(gray));
%     colormap;
%     textStrings = num2str(mean_conf_mat(:), '%0.2f');% Create strings from the matrix values
%     textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
%     [x, y] = meshgrid(1:4);  % Create x and y coordinates for the strings
%     hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center', 'FontSize', 20);
%     midValue = mean(get(gca, 'CLim'));  % Get the 3Ddle value of the color range
%     textColors = repmat(mean_conf_mat(:) > midValue, 1, 3);  % Choose white or black for the
%                                                %   text color of the strings so
%                                                %   they can be easily seen over
%                                                %   the background color
%     set(hStrings, {'Color'}, num2cell(textColors, 2));
%     xlabel('response');
%     ylabel('condition');
%     
%     set(gca,gcaOpts{:}, 'XTick', [1:4], 'YTick', [1:4], 'YTickLabel', {'right', 'left', 'near', 'far'}, 'XTickLabel', {'right', 'left', 'down', 'up'});
    %print('/Users/babylab/Desktop/','-r300','-dpng') % prints a high
    %quality image of the figure
    %% PLOT SORTED BY CONDITION
    close all;
    % significance colors
    p_colormap = jmaColors('pval');
    p_colormap(end,:) = [1 1 1];
    plot_type = 'ave'; % 'diff','conds'
    rcs_to_plot = [1,2];
    for q = 1:2
        figure;
        set(gcf,'units','centimeters');
        fig_pos = get(gcf,'position');
        fig_pos(3) = 20;
        fig_pos(4) = 20;
        set(gcf,'position',fig_pos);
        for e = 1:length(rcs_to_plot)
            if e < size(rca_struct.A,2)
                % plot topography
                egi_h(e) = subplot(6,length(rcs_to_plot),e); 
                mrC.plotOnEgi(rca_struct.A(:,e));
            else
            end
             switch e
                case 1
                    if q == 1
                        y_max = 10; y_min = -15; y_unit = 5;
                    else
                        y_max = 10; y_min = -5; y_unit = 5;
                    end
                case 2
                    if q == 1
                        y_max = 10; y_min = -15; y_unit = 5;
                    else
                        y_max = 10; y_min = -5; y_unit = 5;
                    end
                otherwise
                    if q == 1
                        y_max = 4; y_min = -6; y_unit = 2;
                    else
                        y_max = 4; y_min = -4; y_unit = 2;
                    end
            end
            plot_hist = true;
            switch plot_type
                case 'ave'
                    conds_to_plot = [5,6];
                    leg_names = {'2D movement', '3D movement'};
                    cur_colors = cond_colors([1,3],:);
                case 'diff'
                    conds_to_plot = [7,8]; 
                    leg_names = {'near-diff','far-diff'};
                    plot_hist = false;
                case 'all'
                    conds_to_plot = [1:4,7:8];
                    leg_names = [cond_names,'near-diff','far-diff'];
                otherwise
                    conds_to_plot = 1:4; hist_to_plot = 1:4;
                    leg_names = cond_names;
                    cur_colors = cond_colors;
            end
            clear p_h;   
            if q == 1
                title_str = 'stimulus-locked';
                eeg_h(e) = subplot(6,length(rcs_to_plot),[9,11]+(e-1));
                x_min = 0;
                x_max = 1800;
                time_roi = (1+2*time_res):(3.8*time_res); 
                x_vals = (1000/time_res:1000/time_res:exp_dur)-2000;
                x_vals = x_vals(time_roi);
            else
                title_str = 'response-locked';
                eeg_h(e) = subplot(6,length(rcs_to_plot),[9,11]+(e-1));
                x_min = -pre_time;
                x_max = post_time;
                x_vals = -pre_time:1000/time_res:post_time;
                time_roi = 1:length(x_vals);
            end
            for c = 1:(length(conds_to_plot))
                y_vals = intact_data(:,q).ave(time_roi,rcs_to_plot(e), conds_to_plot(c));
                err_vals = intact_data(:,q).err(time_roi,rcs_to_plot(e), conds_to_plot(c));
                eeg_pos(:,e) = get(eeg_h(e),'position');
                if plot_hist && (conds_to_plot(c) <= size(all_rts,2)) && (q == 1) && (e == 1)
                    hist_y = histc(all_rts{conds_to_plot(c)},bin_edges);
                    hist_y = hist_y(1:end-1);
                    hist_x = bin_edges+bin_size/2;
                    hist_x = hist_x(1:end-1);
                    if c == 1
                        [ax_h,p_h(c),h_h(c)] = plotyy(x_vals,y_vals,hist_x,hist_y,'plot','bar');
                        set(p_h(c), 'linewidth', l_width,'color',cur_colors(c,:));
                        hist_params = {'BarWidth', .8, 'linewidth',l_width, 'edgecolor','w', 'facealpha', .25, 'linestyle','none'};
                        set(h_h(c), hist_params{:}, 'FaceColor',cur_colors(c,:));
                        hold(ax_h(1),'on');
                        hold(ax_h(2),'on');
                    else
                        p_h(c) = plot(ax_h(1), x_vals,y_vals,'linewidth',l_width,'color',cur_colors(c,:));
                        h_h(c) = bar(ax_h(2), hist_x, hist_y, hist_params{:}, 'FaceColor',cur_colors(c,:));
                    end
                    fill(ax_h(1),[x_vals,fliplr(x_vals)],[(y_vals-err_vals)',fliplr((y_vals+err_vals)')],cur_colors(c,:),'linestyle','none','facealpha',.25);
                    arrayfun(@(x) uistack(x,'top'), p_h(c),'uni',false);
                    arrayfun(@(x) uistack(x,'bottom'), h_h(c),'uni',false);
                else
                    ax_h = gca;
                    hold(ax_h(1),'on');
                    p_h(c) = plot(ax_h(1), x_vals,y_vals,'linewidth',l_width,'color',cur_colors(c,:));
                    fill(ax_h(1), [x_vals,fliplr(x_vals)],[(y_vals-err_vals)',fliplr((y_vals+err_vals)')],cur_colors(c,:),'linestyle','none','facealpha',.25);
                    arrayfun(@(x) uistack(x,'top'), p_h(c),'uni',false);
                end
                if c == length(conds_to_plot)   
                    temp_h(1) = plot(ones(2,1)*x_min,[y_min,y_max],'-k','linewidth',l_width);
                    
                    if q == 2
                        temp_h(2) = plot(zeros(2,1),[y_min,y_max],'-','linewidth',l_width, 'color', cond_colors(2,:));
                    else
                        if e == 1
                            temp_h(2) = plot(ones(2,1)*x_max,[y_min,y_max],'-k','linewidth',l_width);
                        else
                        end
                    end
                    arrayfun(@(x) uistack(x,'top'), temp_h,'uni',false);
                else
                end
                clear temp_h;
                % do shuffle comparison
                diff_ready = zeros(length(time_roi),15);
                if q == 2
                    diff_ready = squeeze(intact_data(:,q).subs(time_roi,rcs_to_plot(e), conds_to_plot(c),:) - shuffle_data(:,q).subs(time_roi,rcs_to_plot(e), conds_to_plot(c),:));
                else
                    if length(conds_to_plot) == 2 && c == 2
                        diff_ready = squeeze(intact_data(:,q).subs(time_roi,rcs_to_plot(e), conds_to_plot(1),:) - intact_data(:,q).subs(time_roi,rcs_to_plot(e), conds_to_plot(2),:));
                    else
                    end
                end
                if ~all(diff_ready(:) == 0)
                    [sig_uncorr, p_uncorr, t_uncorr, sig_corr] = ttest_permute_sstats(diff_ready, 10000, 'mass');
                
                    if any(sig_uncorr)
                        % plot corrected t-values
                        if q == 1
                            sig_pos(1) = y_max;
                            sig_pos(2) = y_max - (y_max-y_min) *0.06;
                            sig_pos = sig_pos + (y_max-y_min) *0.08;
                        else
                            sig_pos(1) = y_max - (y_max-y_min) *0.03 * (c-1)*2;
                            sig_pos(2) = y_max - (y_max-y_min) *0.03 * (c)*2;
                            sig_pos = sig_pos + (y_max-y_min) *0.16;
                        end
                        
                        regionIdx = bwlabel(sig_corr);
                        for m=1:max(regionIdx)
                            tmp = regionprops(regionIdx == m,'centroid');
                            idx = round(tmp.Centroid(2));
                            hTxt = text(x_vals(idx),sig_pos(1),'*','fontsize',18,'fontname','Helvetica','horizontalalignment','center','verticalalignment','cap');
                        end

                        % plot uncorrected t-values
                        cur_p = repmat( p_uncorr', 20,1 );
                        h_img = image([min(x_vals),max(x_vals)],[sig_pos(1),sig_pos(2)], cur_p, 'CDataMapping', 'scaled','Parent',gca);
                        colormap( gca, p_colormap );   
                        c_mapmax = .05+2*.05/(size(p_colormap,1));
                        set( gca, 'CLim', [ 0 c_mapmax ] ); % set range for color scale
                        set(gca, gcaOpts{:});
                        uistack(h_img,'bottom')
                    else
                    end
                else
                end
            end

            if e == 1
                text(x_max+(x_max-x_min)*0.15, y_max+(y_max-y_min)*1.75, title_str, text_params{:}, 'fontsize', 18, 'horizontalalignment','center');
                ylabel('amplitude (\muV)', text_params{:})
            else
            end
            if q == 1
                if e == 1
                    l_h = legend(p_h,leg_names,'location','northeast');
                else
                    l_h = legend(p_h,leg_names,'location','northeast');
                end
                legend boxoff;
                l_pos = get(l_h,'position');
                l_pos(1) = l_pos(1) + l_pos(3) * 0.03;
                %l_pos(2) = l_pos(2) + l_pos(4) * 0.4;
                set(l_h,'position',l_pos, text_params{:});
                xlabel('time (ms)',text_params{:})
            else
                xlabel('time before button-press (ms)',text_params{:})
            end
            set(ax_h(1),gcaOpts{:},'ytick',(y_min:y_unit:y_max),'ylim',[y_min,y_max],'xlim', [x_min;x_max], 'ycolor','k','clipping','off');
            set(ax_h(1), 'color', 'none')
            ax_pos = get(ax_h(1), 'position');
            if numel(ax_h) > 1
                set(ax_h(2),gcaOpts{:},'ytick',0:(hist_unit*2):hist_ymax,'xlim', [x_min;x_max],'ylim',[0,hist_ymax*4],'clipping','off')
                uistack(ax_h(2),'bottom')
                set(ax_h(2), 'color', 'none')
            else
            end
            hold off;
            for t = 1:2
                if q == 2
                    all_image = cat(3,rca_sorted_rl{:,t});
                    im_time_roi = 1:size(all_image,1);
                    x_vals = -pre_time:1000/time_res:post_time;
                    x_unit = 200;
                else
                    all_image = cat(3,rca_sorted_sl{:,t});
                    im_time_roi = (1+2*time_res):(3.8*time_res); % temporal ROI for RCA
                    x_vals = (1000/time_res:1000/time_res:exp_dur)-2000;
                    x_unit = 500;
                end
                x_vals = x_vals(im_time_roi);
                t_dim = size(all_image,1);
                not_nan = all(squeeze(~isnan(all_image(:,e,:))),1);
                all_image = squeeze(all_image(:,e, not_nan));
                all_rt = cat(1,beh_sorted{:,t});
                all_rt = all_rt(not_nan);
                [all_sort, all_idx] = sort(all_rt, 'descend');
                all_image = all_image(:,all_idx);
                im_reshaped = imresize(all_image,[t_dim,100])';
                rt_reshaped = imresize(all_sort, [100,1]);
                [rt_corr, corr_p] = corr(rt_reshaped, abs(im_reshaped));

                subplot(6,length(rcs_to_plot),5+(e-1)+(t-1)*2);
                im_h = imagesc(im_reshaped(:,im_time_roi),[y_min y_max]);
                hold on
                set(im_h, 'XData', x_vals);
                x_lims = [floor(min(x_vals)/100)*100, ceil(max(x_vals)/100)*100];
                if q == 2
                    [ax_h, rt_h, corr_h] = plotyy(gca, zeros(1,length(rt_reshaped)), 1:length(rt_reshaped), x_vals, rt_corr(im_time_roi));
                else
                    [ax_h, rt_h, corr_h] = plotyy(gca, rt_reshaped, 1:length(rt_reshaped), x_vals, rt_corr(im_time_roi));
                end
                set(rt_h, 'linewidth',l_width, 'color', cond_colors(2,:));
                set(corr_h, 'linewidth',l_width, 'color', cond_colors(2,:)) 
                if q == 1
                    x_ticks = x_lims(1):x_unit:x_lims(2);
                else
                    x_ticks = [-1500:500:0];
                end
                set(ax_h(1), 'xlim', [x_lims(1), x_lims(2)], 'xtick', x_ticks, gcaOpts{:}, 'ylim', [0,100], 'ytick', 0:50:100, 'ycolor', 'k')
                set(ax_h(2), 'xlim', [x_lims(1), x_lims(2)], 'xtick', x_ticks, gcaOpts{:}, 'ylim', [-.5,.5], 'ytick', -1:.5:1, 'ycolor', 'k')

                label_pos = max(get(ax_h(1), 'xlim'))+diff(get(ax_h(1), 'xlim'))*0.05;
                ylabel(ax_h(2), 'correlation (r)', text_params{:}, 'position', [label_pos,0,0]);
                if plot_corr
                    leg_h = legend([rt_h,corr_h], {'RT','RT \times amplitude'}, 'location', 'southeast');
                else
                    set(corr_h,'visible','off')
                    set(ax_h(2),'visible','off')
                    leg_h = legend(rt_h, {'RT'}, 'location', 'southeast');
                end
                set(leg_h, 'TextColor', [0,0,0], 'box', 'off', text_params{:});
                if q == 2
                    set(leg_h, 'location', 'southwest')
                else
                end
                title(leg_names{t}, text_params{:}, 'color', cond_colors(t+(t-1),:))
                if t == 1
                    im_pos(:,t) = get(gca, 'position');
                    im_pos(4,t) = im_pos(4,t)*1.2;
                    im_pos(2,t) = im_pos(2,t)+im_pos(4,t)*0.4;
                    im_pos(3,t) = im_pos(3,t);
                    set(gca, 'position', im_pos(:,t), 'xticklabels','');
                else
                    if e == 1
                        text(x_min-(x_max-x_min)*.25, y_min, ...
                            {'trial bins'; '(fastest \rightarrow slowest)'}, text_params{:}, 'rotation', 90, 'horizontalalignment','center');
                    else
                    end
                    
                    im_pos(:,t) = get(gca, 'position');
                    im_pos([1,3,4],t) = im_pos([1,3,4],1);
                    im_pos(2,t) = im_pos(2,1) - im_pos(4,1)*1.3;
                    set(gca, 'position', im_pos(:,t), 'xticklabels','');
                    if e == 2
                        c_bar = colorbar(text_params{:});
                        c_bar.Label.String = 'amplitude (\muV)';
                        c_pos = get(c_bar, 'position');
                        if plot_corr
                            c_pos(1) = im_pos(1,t)+im_pos(3,t)*1.1;
                        else
                            c_pos(1) = im_pos(1,t)+im_pos(3,t)*1.02;
                        end
                        c_pos(2) = im_pos(2,t); 
                        c_pos(3) = 0.01; c_pos(4) = abs(diff(im_pos(2,:)))+im_pos(4,1);
                        set(c_bar, 'position', c_pos, text_params{:});
                    else
                    end
                end
            end
            hor_pos(e,:) = im_pos([1,3],t);
        end
        for e = 1:length(egi_h)
            egi_pos = get(egi_h(e),'position');
            egi_pos(1) = egi_pos(1)-egi_pos(3)*2/4;
            egi_pos(2) = egi_pos(2)-egi_pos(4)*.75;
            egi_pos(3:4) = egi_pos(3:4)*2;
            set(egi_h(e),'position',egi_pos);
        end
        set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[fig_pos(3), fig_pos(4)])
        if q == 2
            print(gcf, '-dpdf', '-fillpage', sprintf('%s/motion2D3D_erp_condsort_resp%s.pdf',fig_path, rca_str), '-r0');
            export_fig(gcf, '-png', sprintf('%s/motion2D3D_erp_condsort_resp%s.png',fig_path, rca_str), '-r1000', '-transparent');
        else
            print(gcf, '-dpdf', '-fillpage', sprintf('%s/motion2D3D_erp_condsort_stim%s.pdf',fig_path, rca_str), '-r0');
            export_fig(gcf, '-png', sprintf('%s/motion2D3D_erp_condsort_stim%s.png',fig_path, rca_str), '-r1000', '-transparent');
        end
    end
    
    
    %% PLOT SORTED BY RT
    close all;
    rcs_to_plot = [1,2];
    for q = 1:2
        figure;
        set(gcf,'units','centimeters');
        fig_pos = get(gcf,'position');
        fig_pos(3) = 20;
        fig_pos(4) = 20;
        set(gcf,'position',fig_pos);
        for e = 1:length(rcs_to_plot)
            switch e
            case 1
                if q == 1
                    y_max = 10; y_min = -15; y_unit = 5;
                else
                    y_max = 10; y_min = -5; y_unit = 5;
                end
            case 2
                if q == 1
                    y_max = 10; y_min = -15; y_unit = 5;
                else
                    y_max = 10; y_min = -5; y_unit = 5;
                end
            otherwise
                if q == 1
                    y_max = 4; y_min = -6; y_unit = 2;
                else
                    y_max = 4; y_min = -4; y_unit = 2;
                end
            end
            if e < size(rca_struct.A,2)
                % plot topography
                egi_h(e) = subplot(6,length(rcs_to_plot),e); 
                hold on
                mrC.plotOnEgi(rca_struct.A(:,e));
                hold off
            else
            end
            for t = 1:2
                if t == 1
                    conds_to_plot = [1,3];
                    leg_names = {'fast-2D', 'slow-2D'};
                    cur_colors = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
                else
                    conds_to_plot = [4,6];
                    leg_names = {'fast-3D','slow-3D'};
                    cur_colors = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
                end
                for c = 1:(length(conds_to_plot))
                    if q == 1
                        title_str = 'stimulus-locked';
                        eeg_h(e) = subplot(6,length(rcs_to_plot),[5,7]+(e-1)+(t-1)*4);
                        x_min = 0;
                        x_max = 1800;
                        time_roi = (1+2*time_res):(3.8*time_res); 
                        x_vals = (1000/time_res:1000/time_res:exp_dur)-2000;
                        x_vals = x_vals(time_roi);
                    else
                        title_str = 'response-locked';
                        eeg_h(e) = subplot(6,length(rcs_to_plot),[5,7]+(e-1)+(t-1)*4);
                        x_min = -pre_time;
                        x_max = post_time;
                        x_vals = -pre_time:1000/time_res:post_time;
                        time_roi = 1:length(x_vals);
                    end
                    y_vals = intact_data(:,q+2).ave(time_roi, rcs_to_plot(e), conds_to_plot(c));
                    err_vals = intact_data(:,q+2).err(time_roi, rcs_to_plot(e), conds_to_plot(c));
                    hold on;
                    p_h(c) = plot(x_vals,y_vals,'linewidth',l_width,'color',cur_colors(c,:));
                    fill([x_vals,fliplr(x_vals)],[(y_vals-err_vals)',fliplr((y_vals+err_vals)')],cur_colors(c,:),'linestyle','none','facealpha',.25);
                    
                    clear temp_h
                    if c == length(conds_to_plot)   
                        temp_h(1) = plot(ones(2,1)*x_min,[y_min,y_max],'-k','linewidth',l_width);
                        if q == 2
                            temp_h(2) = plot(zeros(2,1),[y_min,y_max],'-','linewidth',l_width, 'color',cond_colors(2,:));
                        else
                        end
                        arrayfun(@(x) uistack(x,'top'), temp_h,'uni',false);
                    else
                    end

                    % do shuffle comparison
                    diff_ready = zeros(length(time_roi),15);
                    if q == 2
                        diff_ready = squeeze(intact_data(:,q+2).subs(time_roi,rcs_to_plot(e), conds_to_plot(c),:) - shuffle_data(:,q+2).subs(time_roi,rcs_to_plot(e), conds_to_plot(c),:));
                    else
                        if length(conds_to_plot) == 2 && c == 2
                            diff_ready = squeeze(intact_data(:,q+2).subs(time_roi,rcs_to_plot(e), conds_to_plot(1),:) - intact_data(:,q+2).subs(time_roi,rcs_to_plot(e), conds_to_plot(2),:));
                        else
                        end
                    end
                    if ~all(diff_ready(:) == 0)
                        [sig_uncorr, p_uncorr, t_uncorr, sig_corr] = ttest_permute_sstats(diff_ready, 10000, 'mass');

                        if any(sig_uncorr)
                            % plot corrected t-values
                            if q == 1
                                sig_pos(1) = y_max;
                                sig_pos(2) = y_max - (y_max-y_min) *0.06;
                                sig_pos = sig_pos + (y_max-y_min) *0.08;
                            else
                                sig_pos(1) = y_max - (y_max-y_min) *0.03 * (c-1)*2;
                                sig_pos(2) = y_max - (y_max-y_min) *0.03 * (c)*2;
                                sig_pos = sig_pos + (y_max-y_min) *0.16;
                            end
                            regionIdx = bwlabel(sig_corr);
                            for m=1:max(regionIdx)
                                tmp = regionprops(regionIdx == m,'centroid');
                                idx = round(tmp.Centroid(2));
                                hTxt = text(x_vals(idx),sig_pos(1),'*','fontsize',18,'fontname','Helvetica','horizontalalignment','center','verticalalignment','cap');
                            end

                            % plot uncorrected t-values
                            cur_p = repmat( p_uncorr', 20,1 );
                            h_img = image([min(x_vals),max(x_vals)],[sig_pos(1),sig_pos(2)], cur_p, 'CDataMapping', 'scaled','Parent',gca);
                            colormap( gca, p_colormap );   
                            c_mapmax = .05+2*.05/(size(p_colormap,1));
                            set( gca, 'CLim', [ 0 c_mapmax ] ); % set range for color scale
                            set(gca, gcaOpts{:});
                            uistack(h_img,'bottom')
                        else
                        end
                    else
                    end
                end
                arrayfun(@(x) uistack(x), p_h,'uni',false);
                if q == 1
                    if t == 2
                        xlabel('time (ms)',text_params{:});
                    else
                    end
                    leg_loc = 'northeast';
                else
                    if t == 2
                        xlabel('time before button-press (ms)', text_params{:})
                    else
                    end
                    leg_loc = 'northwest';
                end
                if e == 1
                    l_h = legend(p_h,leg_names,'location',leg_loc);
                    legend boxoff;
                    l_pos = get(l_h,'position');
                    l_pos(1) = l_pos(1) + l_pos(3) * 0.2;
                    l_pos(2) = l_pos(2) + l_pos(4) * 0.2;
                    %set(l_h,'position',l_pos);
                    ylabel('amplitude (\muV)', text_params{:})
                else
                end
                xlim([x_min;x_max]);
                ylim([y_min,y_max]);
                eeg_pos(:,e) = get(eeg_h(e),'position');
                if t == 1               
                    set(gca, gcaOpts{:}, 'xticklabels',[], 'ytick',(y_min:y_unit:y_max), 'clipping', 'off', 'color', 'none');
                else
                    eeg_pos(2,e) = eeg_pos(2,e)+eeg_pos(4,e)*.1;
                    set(gca,gcaOpts{:},'ytick',(y_min:y_unit:y_max), 'clipping', 'off', 'color', 'none');
                end
                
                eeg_pos([1,3],e) = hor_pos(e,:);
                eeg_pos(4,e) = im_pos(4,t)*1.5;
                set(eeg_h(e),'position', eeg_pos(:,e));
                hold off;
            end
        end
        % adjust eeg pos
%         for e = 1:length(egi_h)
%             egi_pos = get(egi_h(e),'position');
%             egi_pos(3:4) = egi_pos(3:4)*1.5;
%             egi_pos(1) = egi_pos(1)-egi_pos(1)*.1;
%             egi_pos(2) = eeg_pos(2,e)-eeg_pos(2,1)*.05;
%             set(egi_h(e),'position',egi_pos);
%         end
        set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[fig_pos(3), fig_pos(4)])
        if q == 2
            %print(gcf, '-dpdf', '-fillpage', sprintf('%s/motion2D3D_erp_rtsort_resp_%s.pdf',fig_path, rca_str), '-r0');
            export_fig(gcf, '-png', sprintf('%s/motion2D3D_erp_rtsort_resp%s.png',fig_path, rca_str), '-r1000', '-transparent');
        else
            %print(gcf, '-dpdf', '-fillpage', sprintf('%s/motion2D3D_erp_rtsort_stim_%s.pdf',fig_path, rca_str), '-r0');
            export_fig(gcf, '-png', sprintf('%s/motion2D3D_erp_rtsort_stim%s.png',fig_path, rca_str), '-r1000', '-transparent');
        end
    end
end

function [beh_out, prop_rej, prop_inc] = beh_preproc(beh_in,lower_cutoff,upper_cutoff, rand_split)
    % get correct responses
    beh_in(:,4) = beh_in(:,1) == beh_in(:,3);
    % make list of responses to keep
    beh_in(:,5) = beh_in(:,3) ~= 0;
    % indicate early responses
    beh_in(:,5) = beh_in(:,5) & (beh_in(:,2) >=lower_cutoff);
    beh_in(:,5) = beh_in(:,5) & (beh_in(:,2) <=upper_cutoff);
    cond_list = unique(beh_in(:,1));
    
    beh_out = beh_in;
    for c = 1:length(cond_list)
        cond_idx = beh_in(:,1) == cond_list(c);
        prop_rej(c) = sum(beh_in(cond_idx,5)==0,1)./sum(cond_idx);
        prop_inc(c) = sum(beh_in(cond_idx,4)==0,1)./sum(cond_idx);
        % randomize rts within condition
        if rand_split
            rand_idx = all([cond_idx,beh_in(:,5)],2);
            rand_rts = beh_in(rand_idx, 2);
            rand_rts = rand_rts(randperm(length(rand_rts)));
            beh_out(rand_idx, 2) = rand_rts;
            if ~all(ismember(beh_out(cond_idx,2),beh_in(cond_idx,2)))
                print('hello');
            else
            end
        else
        end
        
    end
    

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
function [eeg_trials , eeg_mean] = stim_averaging(beh_data,eeg_raw, cond_idx)
    % STIMULUS LOCKED AVERAGING
    % RETURN TIME X ELECTRODE MATRIX X CONDITION
    % get list of conditions
    if nargin < 3
        % index of condition column
        cond_idx = 1;
    else
    end
    cond_list = unique(beh_data(:,cond_idx));
    cond_list = cond_list(cond_list > 0);
    % grab trials where the first column == cond_list(c)
    for c = 1:length(cond_list)
        % get rid of missing and incorrect trials
        eeg_trials{c} = eeg_raw(:,:, beh_data(:,cond_idx)==cond_list(c) & beh_data(:,4)==1 & beh_data(:,5)==1 );
        eeg_mean(:,:,c) = squeeze(nanmean(eeg_trials{c},3));
    end
end

function [resp_trials, resp_mean] = resp_averaging(beh_data,eeg_raw,time_res,exp_dur,pre_time,post_time,cond_idx)
    % RESPONSE-LOCKED AVERAGING
    % RETURN TIME X ELECTRODE MATRIX X CONDITION
    % get list of conditions
    if nargin < 7
        % index of condition column
        cond_idx = 1;
    else
    end
    cond_list = unique(beh_data(:,cond_idx));
    cond_list = cond_list(cond_list > 0);
    % pre and post samples
    pre_samples = time_res*(pre_time/1000);
    post_samples = time_res*(post_time/1000);
    % grab trials where the first column == cond_list(c)
    for c = 1:length(cond_list)
        % get rid of missing and incorrect trials
        trial_idx = beh_data(:,cond_idx)==cond_list(c) & beh_data(:,4)==1 & beh_data(:,5)==1;
        get_trials = eeg_raw(:,:, trial_idx );
        get_timing = beh_data(trial_idx,2) + exp_dur/2;
        timing_idx = round(get_timing/1000*time_res);
        resp_temp = arrayfun(@(x) ...
                get_trials(timing_idx(x)-pre_samples:timing_idx(x)+post_samples,:,x),1:length(timing_idx),'uni',false);
        resp_trials{c} = cat(3,resp_temp{:});
        resp_mean(:,:,c) = squeeze(nanmean(resp_trials{c},3));
    end
end

function beh_sort = beh_split(beh_data, num_bins, by_cond)
    % do trial-wise variability
    if by_cond
        cond_sets = [1,2; 3,4];
    else
        cond_sets = [1,2,3,4];
    end
    % make new version of beh input
    beh_sort = beh_data;
    beh_sort(:,end+1) = 0;
    for c = 1:size(cond_sets,1)
        trial_idx = any(beh_data(:,1) == cond_sets(c,:),2) & beh_data(:,4)==1 & beh_data(:,5)==1;
        num_trials = sum(trial_idx);
        [y,sort_idx] = sort(beh_data(trial_idx,2)); % sorting results
        num_vals = floor(num_trials/num_bins);
        def_bins = num_bins*num_vals;
        extra_trials = num_trials-def_bins;
        extra_bins = randperm(num_bins,extra_trials);
        bin_idx = reshape(repmat(1:num_bins,num_vals,1),def_bins,1); % getting ordered bin numbers
        bin_idx = [bin_idx; extra_bins']; % accounting for results above calculated bin number
        bin_idx = sort(bin_idx) + (c-1) * num_bins; % make sure bin_idx is ordered, and multiply by c
        bin_idx(sort_idx) = bin_idx; % now sort bin_idx according to trial order
        beh_sort(trial_idx,end) = bin_idx; % assigning bin numbers to RTs
    end
end
    
function plot_rca_behavior(rca_data,stim_mean,all_rts,Conds,cond_names2,cond_colors,Smoothparam,stat_type,time_res,exp_dur,top_path)
    dostats = true;
    plot_stim_mean = squeeze(mean(stim_mean,4));
    l_width = 2;
    f_size = 12;

    diff2D3D = squeeze(mean(stim_mean(:,:,1:2,:),3)-mean(stim_mean(:,:,3:4,:),3));%2D vs. 3D difference

    figure;
    plot_diff = false; % plot difference waveforms (true/false)
    rcs_to_plot = [1,2];
    for e = 1:length(rcs_to_plot)
        if e <= length(rcs_to_plot)
            % plot topography
            egi_h(e) = subplot(length(rcs_to_plot)+1,3,(3)+3*(e-1)); 
            hold on
            mrC.plotOnEgi(rca_data.A(:,e));
            hold off
        else
        end
        for c = 1:(length(cond_names2))
            title_str = 'stimulus-locked';
            eeg_h(e) = subplot(length(rcs_to_plot)+1,3,(1:2)+3*(e-1));
            % ----------smoothing?-----------------------------------
            %Smoothparam = 2;
            FP = ones(1,Smoothparam)/Smoothparam;
            fDelay = (length(FP)-1)/2;
            %--------------------------------------------------------
            y_vals = filter(FP,1,mean(plot_stim_mean(:,rcs_to_plot(e),Conds{c}),3));
            %y_vals = mean(plot_stim_mean(:,rcs_to_plot(e),Conds{c}),3);
            x_min = 2000;
            x_max = 4000;
            x_vals = 1000/time_res:1000/time_res:exp_dur;

            eeg_pos(:,e) = get(eeg_h(e),'position');
            hold on;
            p_h(c) = plot(x_vals-fDelay,y_vals,'linewidth',l_width,'color',mean(cond_colors(Conds{c},:),1));
            
        end
        
        arrayfun(@(x) uistack(x), p_h,'uni',false);
        if e == 1
            title(title_str,'fontsize',f_size,'fontname','Helvetica')
        elseif e == length(rcs_to_plot)
            ylabel('amplitude (\muV)','fontsize',f_size,'fontname','Helvetica')
        end
        xlim([x_min;x_max]);
        set(gca,'xtick',2000:500:4000,'xticklabel',0:500:2000)
        switch e
            case 1
                    y_max = 20; y_min = -15; y_unit = 5;
            case 2
                    y_max = 10; y_min = -10; y_unit = 5;
            otherwise
                    y_max = 4; y_min = -6; y_unit = 2;   
        end
        %---------------------Stats on RCs-------------------------
        if dostats ,
            %p_condcomp(:,e) = arrayfun(@(x) signrank(squeeze(mean(stim_mean(x,e,1:2,:))),squeeze(mean(stim_mean(x,e,3:4,:)))),1:size(stim_mean)); % uncorrected for multiple comparision        

            [realT,p_condcomp(:,e),corrT,critVal,clustDistrib]= ttest_permute_sstats(squeeze(diff2D3D(:,e,:)));
            sig_ind1 = [find(diff(p_condcomp(:,e)<0.05)==1) find(diff(p_condcomp(:,e)<0.05)==-1)];
            sig_ind2 = [find(diff(corrT)==1) find(diff(corrT)==-1)];
            arrayfun(@(x) fill([x_vals(sig_ind1(x,:)) flip(x_vals(sig_ind1(x,:)))],[ones(1,2)*(y_min+(y_max-y_min)*.08) ones(1,2)*(y_min)],'y','LineStyle','none'),1:size(sig_ind1,1));
            alpha(.45);
            arrayfun(@(x) fill([x_vals(sig_ind2(x,:)) flip(x_vals(sig_ind2(x,:)))],[ones(1,2)*(y_min+(y_max-y_min)*.08) ones(1,2)*(y_min)],'r','LineStyle','none'),1:size(sig_ind2,1));
            
        end
        
        %-----------------------------------------------------------    
        ylim([y_min,y_max]);
        set(gca,'ytick',(y_min:y_unit:y_max));
        hold off;
    end
    
    
    % plot histograms
    eeg_h(e) = subplot(length(rcs_to_plot)+1,3,(1:2)+3*(length(rcs_to_plot)));
    
    for c = 1:length(cond_names2)
        hold on
        all_rts_conds = all_rts(Conds{c});
        [f,xi] = ksdensity(cat(1,all_rts_conds{:}),1:1:2000); 
        if strcmp(stat_type,'cdf'), f = cumsum(f); end% if CDF
        p_h(c) = plot(xi,f,'color',mean(cond_colors(Conds{c},:),1),'linewidth',2);
        xlim([0 2000]);
        if strcmp(stat_type,'cdf'),
            ylim([0,1]);
        else
            ylim([0,.0025]);
        end
        %fill([xi flip(xi)],[f zeros(size(f))],mean(cond_colors(Conds{c},:),1),'LineStyle','none');
        alpha(0.25)
        set(gca,'xtick',0:500:2000)
    end
    
     
    l_h = legend(p_h,cond_names2,'location','northeast');
    legend boxoff;
    l_pos = get(l_h,'position');
    l_pos(1) = l_pos(1) + l_pos(3) * 1.5;
    l_pos(2) = l_pos(2) + l_pos(4) * 0.2;
    set(l_h,'position',l_pos);
    xlabel('time (ms)','fontsize',f_size,'fontname','Helvetica')
    ylabel(stat_type,'fontsize',f_size,'fontname','Helvetica')
     
    
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
    %export_fig(sprintf('%s/eeg_data.pdf',top_path),'-pdf',gcf);
end    




