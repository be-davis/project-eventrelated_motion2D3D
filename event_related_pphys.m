function event_related_pphys(run_exp)
    close all;
    if nargin < 1
        run_exp = 2;
    else
    end
    addpath(genpath('/Volumes/Denali_4D2/Brandon/code/git/mrC'));
    if run_exp == 1
        top_path = '/Volumes/Denali_4D2/Brandon/eventrelated_motion2D3D/pphys_v1';
        exp_time = 2; % two second experiment
    elseif run_exp == 2
        top_path = '/Volumes/Denali_4D2/Brandon/eventrelated_motion2D3D/pphys_v2';
        exp_time = 4; % four second experiment
    else
        msg = sprintf('\n unknown experiment: %d',run_exp);
        error(msg);
    end
    lower_cutoff = exp_time/2+.2; % lowest allowed RT, in seconds
    %% PREPARE ANALYSIS
    cond_names = {'right','left','near','far'};
    sub_dirs = subfolders([top_path,'/nl-*'],1); % list of subject directories
    num_subs = length(sub_dirs);
    %% RUN ANALYSIS
    for s = 1:length(sub_dirs) % loop over subject directories
        cur_dir = [sub_dirs{s},'/Exp_MATL'];
        cur_files = subfiles([cur_dir,'/RT*'],1);
        cur_rt = [];
        cur_resp = [];
        cur_cond = [];
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
        % converting string response values to numbers
        [~,cur_resp] = ismember(cur_resp,{'Ra','La','Da','Ua'});
        % count missed responses
        proportion_misses(s) = size(find(cur_resp == 0),1) / size(cur_resp,1) * 100;

        % get rid of missed responses
        cur_rt = cur_rt(cur_resp ~= 0);
        cur_cond = cur_cond(cur_resp ~= 0);
        cur_resp = cur_resp(cur_resp ~= 0);

        % get rid of early responses
       if any(cur_rt < lower_cutoff)
           msg = sprintf('subject %s has %d responses below cutoff',sub_dirs{s},sum(cur_rt<lower_cutoff));
           warning(msg);
            cur_rt = cur_rt(cur_rt >= lower_cutoff);
            cur_cond = cur_cond(cur_rt >= lower_cutoff);
            cur_resp = cur_resp(cur_rt >= lower_cutoff);
            cur_rt = cur_rt(cur_rt >= lower_cutoff);
        else
       end

        cur_correct = cur_resp == cur_cond;

        %mean rt, percent correct, and confusion matrix for correct trials
        for c = 1:4
            mean_rt(s,c) = (mean(cur_rt(cur_cond == c & cur_correct==true)-exp_time/2)*1000);
            percent_correct(s,c) = length(cur_rt(cur_cond == c & cur_correct==1))./length(cur_rt(cur_cond == c))*100;
            conf_mat(c,:,s) = arrayfun(@(x) numel(find(cur_resp(cur_cond==c)==x))./numel(cur_resp(cur_cond==c)),1:4);
        end

    end
    %% MAKE FIGURES 
    disp(conf_mat); % Prints the confusion matrix for each subject in the command window
    mean_conf_mat = mean(conf_mat, 3); % Averages responses across all subjects
    %Displaying heat map based on mean_conf_mat

    imagesc(mean_conf_mat);
    fig = gcf;
    ax = fig.CurrentAxes;
    ax.XLabel.String = 'response';
    ax.YLabel.String = 'condition';
    set(ax, 'YTick', [1 2 3 4], 'XTick', [1 2 3 4], 'YTickLabel', {'right', 'left', 'near', 'far'}, 'XTickLabel', {'right', 'left', 'down', 'up'});
    colormap('summer')
    colorbar;
    %Setting window for displaying percent correct and reaction time
    fig_h = 8; % 8 inches
    fig_w = 5; % 4 inches
    % define a cell variable, containing axis options
    fSize = 12;
    lWidth = 2;
    axOpts = {'tickdir','out','ticklength',[0.0100,0.0100],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth,'xlim',[.5,4.5]};

    figure;
    subplot(2,1,1);
    rt_grandmean = (mean(mean_rt,1));
    rt_granderr = std(mean_rt,0,1)./sqrt(num_subs);
    %Making bars different colors for reaction time bar graph
    hold on
    bar_colors = {'y','m','c','r'};
    for i = 1:length(rt_grandmean)
      if i == 1
        h=bar(i, rt_grandmean(i),'facecolor',bar_colors{i});
        elseif i == 2
            h=bar(i, rt_grandmean(i),'facecolor',bar_colors{i});
        elseif i == 3
            h=bar(i, rt_grandmean(i),'facecolor',bar_colors{i});
        elseif i == 4
            h=bar(i, rt_grandmean(i),'facecolor',bar_colors{i});
        else
        end 
    end
    %Creating x and y axis labels    
    hold off        
    set(gca,axOpts{:},'XTickLabel',cond_names, 'XTick',1:numel(cond_names));
    ylabel('reaction time (ms)');
    hold on;
    errorb(1:4,rt_grandmean,rt_granderr);
    hold off;
    subplot(2,1,2);
    correct_grandmean = (mean(percent_correct,1));
    correct_granderr = std(percent_correct,0,1)./sqrt(num_subs);
    %Making bars different colors for percent correct bar graph
    hold on
    for i = 1:length(correct_grandmean)
        if i == 1
        h=bar(i, correct_grandmean(i),'facecolor',bar_colors{i});
        elseif i == 2
            h=bar(i, correct_grandmean(i),'facecolor',bar_colors{i});
        elseif i == 3
            h=bar(i, correct_grandmean(i),'facecolor',bar_colors{i});
        elseif i == 4
            h=bar(i, correct_grandmean(i),'facecolor',bar_colors{i});
        else
        end 
    end
    hold off
    set(gca, axOpts{:},'XTickLabel',cond_names, 'XTick',1:numel(cond_names));
    ylabel('percent correct');
    hold on;
    errorb(1:4,correct_grandmean,correct_granderr);
    set(gcf,'units','inches');
    fig_pos = get(gcf,'position');
    fig_pos(3) = fig_w;
    fig_pos(4) = fig_h;
    set(gcf,'position',fig_pos)
end

