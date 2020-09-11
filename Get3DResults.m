Basedir = '/storage/maullz/Confidence_Sets_Manuscript/';

outdir  = fullfile(Basedir,'Figures','3D_Figures');

if ~isdir(outdir)
    mkdir(outdir)
end

sim_results_mat = {fullfile(Basedir,'Sim_13_results','30_subjects','Sim_13_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_13_results','60_subjects','Sim_13_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_13_results','120_subjects','Sim_13_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_13_results','240_subjects','Sim_13_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_13_results','480_subjects','Sim_13_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_16_results','30_subjects','Sim_16_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_16_results','60_subjects','Sim_16_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_16_results','120_subjects','Sim_16_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_16_results','240_subjects','Sim_16_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_16_results','480_subjects','Sim_16_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_14_results','30_subjects','Sim_14_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_14_results','60_subjects','Sim_14_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_14_results','120_subjects','Sim_14_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_14_results','240_subjects','Sim_14_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_14_results','480_subjects','Sim_14_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_17_results','30_subjects','Sim_17_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_17_results','60_subjects','Sim_17_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_17_results','120_subjects','Sim_17_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_17_results','240_subjects','Sim_17_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_17_results','480_subjects','Sim_17_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_15_results','30_subjects','Sim_15_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_15_results','60_subjects','Sim_15_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_15_results','120_subjects','Sim_15_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_15_results','240_subjects','Sim_15_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_15_results','480_subjects','Sim_15_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_18_results','30_subjects','Sim_18_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_18_results','60_subjects','Sim_18_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_18_results','120_subjects','Sim_18_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_18_results','240_subjects','Sim_18_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_18_results','480_subjects','Sim_18_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_19_results','30_subjects','Sim_19_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_19_results','60_subjects','Sim_19_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_19_results','120_subjects','Sim_19_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_19_results','240_subjects','Sim_19_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_19_results','480_subjects','Sim_19_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_22_results','30_subjects','Sim_22_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_22_results','60_subjects','Sim_22_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_22_results','120_subjects','Sim_22_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_22_results','240_subjects','Sim_22_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_22_results','480_subjects','Sim_22_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_20_results','30_subjects','Sim_20_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_20_results','60_subjects','Sim_20_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_20_results','120_subjects','Sim_20_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_20_results','240_subjects','Sim_20_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_20_results','480_subjects','Sim_20_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_23_results','30_subjects','Sim_23_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_23_results','60_subjects','Sim_23_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_23_results','120_subjects','Sim_23_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_23_results','240_subjects','Sim_23_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_23_results','480_subjects','Sim_23_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_21_results','30_subjects','Sim_21_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_21_results','60_subjects','Sim_21_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_21_results','120_subjects','Sim_21_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_21_results','240_subjects','Sim_21_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_21_results','480_subjects','Sim_21_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_24_results','30_subjects','Sim_24_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_24_results','60_subjects','Sim_24_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_24_results','120_subjects','Sim_24_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_24_results','240_subjects','Sim_24_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_24_results','480_subjects','Sim_24_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_25_results','30_subjects','Sim_25_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_25_results','60_subjects','Sim_25_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_25_results','120_subjects','Sim_25_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_25_results','240_subjects','Sim_25_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_25_results','480_subjects','Sim_25_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_28_results','30_subjects','Sim_28_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_28_results','60_subjects','Sim_28_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_28_results','120_subjects','Sim_28_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_28_results','240_subjects','Sim_28_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_28_results','480_subjects','Sim_28_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_26_results','30_subjects','Sim_26_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_26_results','60_subjects','Sim_26_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_26_results','120_subjects','Sim_26_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_26_results','240_subjects','Sim_26_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_26_results','480_subjects','Sim_26_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_29_results','30_subjects','Sim_29_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_29_results','60_subjects','Sim_29_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_29_results','120_subjects','Sim_29_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_29_results','240_subjects','Sim_29_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_29_results','480_subjects','Sim_29_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_27_results','30_subjects','Sim_27_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_27_results','60_subjects','Sim_27_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_27_results','120_subjects','Sim_27_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_27_results','240_subjects','Sim_27_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_27_results','480_subjects','Sim_27_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_30_results','30_subjects','Sim_30_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_30_results','60_subjects','Sim_30_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_30_results','120_subjects','Sim_30_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_30_results','240_subjects','Sim_30_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_30_results','480_subjects','Sim_30_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_31_results','30_subjects','Sim_31_30_subjects.mat'), ...
		       fullfile(Basedir,'Sim_31_results','60_subjects','Sim_31_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_31_results','120_subjects','Sim_31_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_31_results','240_subjects','Sim_31_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_31_results','480_subjects','Sim_31_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_32_results','30_subjects','Sim_32_30_subjects.mat'), ...
		       fullfile(Basedir,'Sim_32_results','60_subjects','Sim_32_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_32_results','120_subjects','Sim_32_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_32_results','240_subjects','Sim_32_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_32_results','480_subjects','Sim_32_480_subjects.mat'); ...
                   fullfile(Basedir,'Sim_33_results','30_subjects','Sim_33_30_subjects.mat'), ...
		       fullfile(Basedir,'Sim_33_results','60_subjects','Sim_33_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_33_results','120_subjects','Sim_33_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_33_results','240_subjects','Sim_33_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_33_results','480_subjects','Sim_33_480_subjects.mat'), ...
                   fullfile(Basedir,'Sim_34_results','30_subjects','Sim_34_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_34_results','60_subjects','Sim_34_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_34_results','120_subjects','Sim_34_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_34_results','240_subjects','Sim_34_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_34_results','480_subjects','Sim_34_480_subjects.mat'), ...
                   fullfile(Basedir,'Sim_35_results','30_subjects','Sim_35_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_35_results','60_subjects','Sim_35_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_35_results','120_subjects','Sim_35_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_35_results','240_subjects','Sim_35_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_35_results','480_subjects','Sim_35_480_subjects.mat'), ...
                   fullfile(Basedir,'Sim_36_results','30_subjects','Sim_36_30_subjects.mat'), ...
                   fullfile(Basedir,'Sim_36_results','60_subjects','Sim_36_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_36_results','120_subjects','Sim_36_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_36_results','240_subjects','Sim_36_240_subjects.mat'), ...
                   fullfile(Basedir,'Sim_36_results','480_subjects','Sim_36_480_subjects.mat')};
             
sim_results_mat_dim = size(sim_results_mat);

nominal_levels      = [80, 90, 95];
nominal_prct_levels = nominal_levels/100;
nominal_levels_dim  = size(nominal_levels);

result = load(sim_results_mat{1});
nRlz          = result.nRlz;

std_error_vector = zeros(nominal_levels_dim);
sim_covering_levels_and_subs = zeros([3 5 24 3]);

for i = 1:nominal_levels_dim(2)
    i
    std_error = sqrt((nominal_prct_levels(i)*(1 - nominal_prct_levels(i)))/nRlz);
    std_error_vector(i) = std_error; 
 
    for j = 1:sim_results_mat_dim(1)
        for k = 1:sim_results_mat_dim(2)
            temp = load(char(sim_results_mat(j,k)));
            raw_result                          = sprintf('percentage_success_vector_raw_%d_alternate', nominal_levels(i));
            observed_result                     = sprintf('percentage_success_vector_observed_%d_alternate', nominal_levels(i));
            nSubj                               = sprintf('nSubj');
            sim_covering_levels_and_subs(1,k,j,i) = temp.(raw_result);
            sim_covering_levels_and_subs(2,k,j,i) = temp.(observed_result);
            sim_covering_levels_and_subs(3,k,j,i) = temp.(nSubj);
        end
    end
end

results_params = struct('nom_80_results', struct('nominal_level', nominal_levels(1), ...
                                                 'nominal_prct_level', nominal_prct_levels(1), ...
                                                 'std_error', std_error_vector(1), ...
                                                 'met_1_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,1,1),  ...
                                                 'met_1_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,2,1),  ...
                                                 'met_2_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,3,1),  ...
                                                 'met_2_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,4,1),  ...
                                                 'met_3_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,5,1),  ...
                                                 'met_3_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,6,1),  ...
                                                 'met_1_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,7,1),  ...
                                                 'met_1_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,8,1),  ...
                                                 'met_2_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,9,1),  ...
                                                 'met_2_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,10,1),  ...
                                                 'met_3_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,11,1),  ...
                                                 'met_3_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,12,1),  ...
                                                 'met_1_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,13,1),  ...
                                                 'met_1_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,14,1),  ...
                                                 'met_2_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,15,1),  ...
                                                 'met_2_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,16,1),  ...
                                                 'met_3_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,17,1),  ...
                                                 'met_3_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,18,1),  ...
            						 'met_1_sig_4_results', sim_covering_levels_and_subs(:,:,19,1),  ...
            						 'met_2_sig_4_results', sim_covering_levels_and_subs(:,:,20,1),  ...
            						 'met_3_sig_4_results', sim_covering_levels_and_subs(:,:,21,1),  ...
                                                 'met_1_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,22,1),  ...
                                                 'met_2_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,23,1),  ...
                                                 'met_3_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,24,1)), ...
                        'nom_90_results', struct('nominal_level', nominal_levels(2), ...
                                                 'nominal_prct_level', nominal_prct_levels(2), ...
                                                 'std_error', std_error_vector(2), ...
                                                 'met_1_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,1,2),  ...
                                                 'met_1_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,2,2),  ...
                                                 'met_2_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,3,2),  ...
                                                 'met_2_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,4,2),  ...
                                                 'met_3_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,5,2),  ...
                                                 'met_3_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,6,2),  ...
                                                 'met_1_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,7,2),  ...
                                                 'met_1_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,8,2),  ...
                                                 'met_2_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,9,2),  ...
                                                 'met_2_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,10,2),  ...
                                                 'met_3_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,11,2),  ...
                                                 'met_3_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,12,2),  ...
                                                 'met_1_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,13,2),  ...
                                                 'met_1_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,14,2),  ...
                                                 'met_2_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,15,2),  ...
                                                 'met_2_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,16,2),  ...
                                                 'met_3_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,17,2),  ...
                                                 'met_3_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,18,2),  ...
            						 'met_1_sig_4_results', sim_covering_levels_and_subs(:,:,19,2),  ...
            						 'met_2_sig_4_results', sim_covering_levels_and_subs(:,:,20,2),  ...
            						 'met_3_sig_4_results', sim_covering_levels_and_subs(:,:,21,2),  ...
                                                 'met_1_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,22,2),  ...
                                                 'met_2_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,23,2),  ...
                                                 'met_3_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,24,2)),  ...
                        'nom_95_results', struct('nominal_level', nominal_levels(3), ...
                                                 'nominal_prct_level', nominal_prct_levels(3), ...
                                                 'std_error', std_error_vector(3), ...
                                                 'met_1_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,1,3),  ...
                                                 'met_1_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,2,3),  ...
                                                 'met_2_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,3,3),  ...
                                                 'met_2_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,4,3),  ...
                                                 'met_3_sig_1_std_1_results', sim_covering_levels_and_subs(:,:,5,3),  ...
                                                 'met_3_sig_1_std_2_results', sim_covering_levels_and_subs(:,:,6,3),  ...
                                                 'met_1_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,7,3),  ...
                                                 'met_1_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,8,3),  ...
                                                 'met_2_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,9,3),  ...
                                                 'met_2_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,10,3),  ...
                                                 'met_3_sig_2_std_1_results', sim_covering_levels_and_subs(:,:,11,3),  ...
                                                 'met_3_sig_2_std_2_results', sim_covering_levels_and_subs(:,:,12,3),  ... 
                                                 'met_1_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,13,3),  ...
                                                 'met_1_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,14,3),  ...
                                                 'met_2_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,15,3),  ...
                                                 'met_2_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,16,3),  ...
                                                 'met_3_sig_3_std_1_results', sim_covering_levels_and_subs(:,:,17,3),  ...
                                                 'met_3_sig_3_std_2_results', sim_covering_levels_and_subs(:,:,18,3),  ... 
            						 'met_1_sig_4_results', sim_covering_levels_and_subs(:,:,19,3),  ...
            						 'met_2_sig_4_results', sim_covering_levels_and_subs(:,:,20,3),  ...
            						 'met_3_sig_4_results', sim_covering_levels_and_subs(:,:,21,3),  ...
                                                 'met_1_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,22,3),  ...
                                                 'met_2_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,23,3),  ...
                                                 'met_3_sig_4_thr_05_results', sim_covering_levels_and_subs(:,:,24,3))  ...
                        );
                    
filename = fullfile(outdir, '3DResults.mat');
save(filename,'results_params');
