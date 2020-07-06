results_mat_dir = '/Users/maullz/Desktop/Confidence_Sets_Manuscript/Figures/3D_Figures';

results_mat_file = fullfile(results_mat_dir,'3DResults.mat');

load(results_mat_file)

nominal_vec      = ["nom_80_results","nom_90_results","nom_95_results"];
signal_vec       = ["1","2","3","4"];
signal_title_vec = ["Small Sphere", "Large Sphere","Multi. Sphere","UK Biobank"];
std_vec          = ["1","2"];
std_title_vec    = ["Homo. Variance", "Hetero. Variance"];
color_vec        = 'rbm';
color_title_vec  = ["Algorithm 1", "Algorithm 2", "Algorithm 3"];

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

for i = 1:length(signal_vec)
    signal = signal_vec(i);
    if i ~= length(signal_vec)
        for j = 1:length(std_vec)
            std = std_vec(j);
            y_lim   = 0.85; 
            figure, clf
            for k = 1:length(nominal_vec)
                subplot(2,3,k);
                hold on
                percent = nominal_vec(k);
                for l = 1:length(color_vec)
                    results        = strcat("met_",num2str(l),"_sig_",num2str(i),"_std_",num2str(j),"_results");
                    [est_bdry_cov] = results_params.(percent).(results)(2,:);
                    [subs]         = results_params.(percent).(results)(3,:);
                    plot(subs,est_bdry_cov,[color_vec(l) '-' 'x'],'linewidth', 1.5);
                    y_lim = min(y_lim,min(est_bdry_cov));
                end
                % plot the nominal lvl
                nominal_prct_level = results_params.(percent).nominal_prct_level;
                plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
                plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(percent).std_error nominal_prct_level-1.96*results_params.(percent).std_error], 'k--')
                plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(percent).std_error nominal_prct_level+1.96*results_params.(percent).std_error], 'k--')
                % set the range for the axis
                xlim([subs(1)-5, subs(end)+5])
                ylim([y_lim - 0.05 1]);
                % specify tiks
                xticks(subs)
                % put label onto the axis
                xlabel('Sample Size [N]');
                ylabel('Emp. Covering Rate');

                if j == 1 
                    titlename = {sprintf(signal_title_vec(i) + ' (%d%% Nom.)', results_params.(percent).nominal_level), std_title_vec(j)};
                else
                    titlename =  std_title_vec(j); 
                end
                title(titlename);

                set(gca, 'fontsize', 18);
                axis square;
                hold off

                if [k,j] == [length(nominal_vec),length(std_vec)] 
                    % create legend
                    lgd = legend('Algorithm 1', ...
                                 'Algorithm 2', ...
                                 'Algorithm 3', ...
                                 'Nominal Coverage Level', ...
                                 '\pm 1.96 \times Simulation Std. Error');
                    lgd.FontSize = 18;
                end   

            end

            if j == length(std_vec)
                lgd_plot = subplot(2,3,5);
                axis square;
                pos_lgd  = get(lgd_plot,'position');
                lgd.FontWeight = 'bold';
                set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
                set(lgd, 'interpreter', 'tex');
                axis(lgd_plot,'off');
            end

            set(gcf,'position', [-21,120,1195,682]);
            fh = gcf;
            set(fh,'color','w');
            export_fig(fh,fullfile(results_mat_dir,strcat('sig_',num2str(i),'_std_',num2str(j),'_results.pdf')))
        end
    else
        y_lim   = 0.85; 
        figure, clf
            for k = 1:length(nominal_vec)
                subplot(2,3,k);
                hold on
                percent = nominal_vec(k);
                for l = 1:length(color_vec)
                    results        = strcat("met_",num2str(l),"_sig_",num2str(i),"_results");
                    [est_bdry_cov] = results_params.(percent).(results)(2,:);
                    [subs]         = results_params.(percent).(results)(3,:);
                    plot(subs,est_bdry_cov,[color_vec(l) '-' 'x'],'linewidth', 1.5);
                    y_lim = min(y_lim,min(est_bdry_cov));
                end
                % plot the nominal lvl
                nominal_prct_level = results_params.(percent).nominal_prct_level;
                plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
                plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(percent).std_error nominal_prct_level-1.96*results_params.(percent).std_error], 'k--')
                plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(percent).std_error nominal_prct_level+1.96*results_params.(percent).std_error], 'k--')
                % set the range for the axis
                xlim([subs(1)-5, subs(end)+5])
                ylim([y_lim - 0.05 1]);
                % specify tiks
                xticks(subs)
                % put label onto the axis
                xlabel('Sample Size [N]');
                ylabel('Emp. Covering Rate');
                
                titlename = sprintf(signal_title_vec(i) + ' (%d%% Nom.)', results_params.(percent).nominal_level);
                title(titlename);

                set(gca, 'fontsize', 18);
                axis square;
                hold off

                if k == length(nominal_vec)
                    % create legend
                    lgd = legend('Algorithm 1', ...
                                 'Algorithm 2', ...
                                 'Algorithm 3', ...
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
            export_fig(fh,fullfile(results_mat_dir,strcat('sig_',num2str(i),'_results.pdf')))
        
    end
        
end