results_mat_dir = '/Users/maullz/Desktop/Confidence_Sets_Manuscript/Figures/3D_Figures';

results_mat_file = fullfile(results_mat_dir,'3DResults.mat');

load(results_mat_file)

nominal_vec = ["nom_80_results","nom_90_results","nom_95_results"];
signal_vec  = ["met_1_sig_1_std_1_results","met_1_sig_1_std_2_results", ...
               "met_2_sig_1_std_1_results","met_2_sig_1_std_2_results", ...
               "met_3_sig_1_std_1_results","met_3_sig_1_std_2_results", ...
               "met_1_sig_2_std_1_results","met_1_sig_2_std_2_results", ...
               "met_2_sig_2_std_1_results","met_2_sig_2_std_2_results", ...
               "met_3_sig_2_std_1_results","met_3_sig_2_std_2_results", ...
               "met_1_sig_3_std_1_results","met_1_sig_3_std_2_results", ...
               "met_2_sig_3_std_1_results","met_2_sig_3_std_2_results", ...
               "met_3_sig_3_std_1_results","met_3_sig_3_std_2_results", ...
               "met_1_sig_4_results","met_2_sig_4_results","met_3_sig_4_results" ...
               ];
color_vec   = 'rbrbrbrbrbrbrbrbrbrrr';

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 1:2
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 1 Small Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
        lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_1_sig_1_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 3:4
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 2 Small Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_2_sig_1_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 5:6
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 3 Small Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_3_sig_1_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 7:8
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 1 Large Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_1_sig_2_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 9:10
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 2 Large Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_2_sig_2_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 11:12
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 3 Large Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_3_sig_2_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 13:14
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.0 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 1 Multi. Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_1_sig_3_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 14:15
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.0 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 2 Multi. Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_2_sig_3_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 17:18
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.0 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 3 Multi. Sphere (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_3_sig_3_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 19
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 1 Biobank (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_1_sig_4_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 20
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 2 Biobank (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_2_sig_4_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 21
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(3,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.7 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Met. 3 Biobank (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'met_3_sig_4_coverage_results.pdf'))