function Algorithm_3(nSubj, SvNm, nRlz, thr)

%------------Starting Up initialization
if (nargin<1)
  nSubj  = 60;  % Number of subjects
end
if (nargin<2)
  SvNm  = 'Normsim';  % Save name
end
if (nargin<3)
  nRlz = 5000;
end  
if (nargin<4)
  thr = 0.5;
end 
if exist([SvNm '.mat'], 'file')
  error('Will not overwrite sim result')
end

% Define parameters 
ground_truth_dir = '/gpfs2/well/nichols/users/bas627/Confidence_Sets_Manuscript/Biobank_simulation/ground_truth_files';
Biobank_dir = '/well/nichols/projects/UKB/IMAGING/ContourInf/MNI';
MNI_mask_file = fullfile(ground_truth_dir, 'MNI152_T1_2mm_brain_mask.nii');
MNI_mask_file = spm_vol(MNI_mask_file);
MNI_mask_file = spm_read_vols(MNI_mask_file);
tau     = 1/sqrt(nSubj);
nBoot   = 5000;
dim     = [91 109 91];

% Variables for the transformation
a       = sqrt((nSubj-1)/(nSubj -3 ));
b       = sqrt(nSubj)*sqrt((8*nSubj^2 - 17*nSubj + 11)/((5-4*nSubj)^2*(nSubj-3)));
alpha   = b^-1;
beta    = b/a;

transformed_thr = alpha*asinh((1 - (3/(4*nSubj - 5)))^(-1)*thr*beta);
second_deriv_term = (1/(2*nSubj))*b^2*((1 - (3/(4*nSubj - 5)))^(-1)*thr)*(( ((nSubj-1)/(nSubj-3)) + nSubj*(thr^2)*( (8*nSubj^2 - 17*nSubj + 11)/(16*(nSubj-3)*((nSubj - 2)^2)) ))^(-1/2)); 

observed_data  = zeros([prod(dim) nSubj]);

% This stores the vector SupG for each run
% This vector stores the result for each realisation on whether AC^+ < AC < AC^ for each level of smoothing (1 if true, 0 if false) 
subset_success_vector_raw_80           = zeros(nRlz, 1); 
subset_success_vector_raw_90           = zeros(nRlz, 1);
subset_success_vector_raw_95           = zeros(nRlz, 1);
subset_success_vector_observed_80           = zeros(nRlz, 1); 
subset_success_vector_observed_90           = zeros(nRlz, 1);
subset_success_vector_observed_95           = zeros(nRlz, 1);
subset_success_vector_raw_80_alternate = zeros(nRlz, 1); 
subset_success_vector_raw_90_alternate = zeros(nRlz, 1);
subset_success_vector_raw_95_alternate = zeros(nRlz, 1);
subset_success_vector_observed_80_alternate = zeros(nRlz, 1); 
subset_success_vector_observed_90_alternate = zeros(nRlz, 1);
subset_success_vector_observed_95_alternate = zeros(nRlz, 1);

%- This vector stores the threshold value 'c' for each run
threshold_raw_80_store                  = zeros(nRlz, 1);
threshold_raw_90_store                  = zeros(nRlz, 1);
threshold_raw_95_store                  = zeros(nRlz, 1);

threshold_observed_80_store                  = zeros(nRlz, 1);
threshold_observed_90_store                  = zeros(nRlz, 1);
threshold_observed_95_store                  = zeros(nRlz, 1);

%- This vector stores the percentage volumes A^+_c/A_c, A^_c/A_c, A^-_c/A_c
lower_contour_raw_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_raw_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_raw_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_95_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_95_volume_prct_store                     = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store                   = zeros(nBoot, nRlz);
supG_observed_store              = zeros(nBoot, nRlz);

supG_raw                         = zeros(nBoot,1);
supG_observed                    = zeros(nBoot,1);

cohen_d_unmasked = spm_vol(fullfile(ground_truth_dir, 'Biobank_4000_cohens_d.nii'));
cohen_d_unmasked = spm_read_vols(cohen_d_unmasked);

% Removing subjects that were used to get ground truth
ground_truth_subjects = load(fullfile(ground_truth_dir,'4000_random_sample.mat'));
% All the biobank cope_files
cope_files = cellstr(spm_select('FPList', Biobank_dir, ['.*\_cope5_MNI.nii.gz']));
% All the biobank mask_files
mask_files = cellstr(spm_select('FPList', Biobank_dir, ['.*\_mask_MNI.nii.gz']));
% We have a total of 8945 Biobank subject-level cope files...
total_subjects = 1:8945;
% but importantly, we remove the files we used to create the ground truth!
total_subjects = setdiff(total_subjects, ground_truth_subjects.shuffle_ids);

for t=1:nRlz
    fprintf('.');
    observed_mean = zeros(prod(dim),1);
    observed_std  = zeros(prod(dim),1);
    intersection_mask = zeros(prod(dim),1);
    % Create a random subset of nSubj subjects from the Biobank data
    subset_of_subjects = total_subjects(randperm(size(total_subjects,2), nSubj));
      for i=1:nSubj
        % Load in Biobank subject-level copes and mask
        if tID <= floor(4945/(nSubj*nRlz))
            subject_cope = cope_files{total_subjects(i + (t-1)*nSubj + (tID - 1)*nRlz*nSubj)};
            subject_mask = mask_files{total_subjects(i + (t-1)*nSubj)};
        elseif tID == floor(4945/(nsubj*nRlz) + 1)
            k = floor((4945 - tID*nRlz*nSubj)/nSubj);
            if t <= k
              subject_cope = cope_files{total_subjects(i + (t-1)*nSubj + (tID - 1)*nRlz*nSubj)};
              subject_mask = mask_files{total_subjects(i + (t-1)*nSubj)};
            else
              subject_cope = cope_files{subset_of_subjects(i)};
              subject_mask = mask_files{subset_of_subjects(i)};
            end 
        else
            subject_cope = cope_files{subset_of_subjects(i)};
            subject_mask = mask_files{subset_of_subjects(i)};
        end
        
        % We have to copy and unzip tImgs because octave cant deal with .gz
        [~, subject_cope_name, subject_cope_ext] = fileparts(subject_cope);
        % Copying it to the groun_truth_dir where I can unzip and delete it
        copyfile(subject_cope, pwd);
        gunzip(fullfile(pwd, [subject_cope_name subject_cope_ext]));
        
        
        tImgs = spm_vol(fullfile(pwd, subject_cope_name));
        tImgs = spm_read_vols(tImgs);
        tImgs = reshape(tImgs, [prod(dim), 1]);
        
        observed_data(:,i) = tImgs; 
        observed_mean = observed_mean + tImgs;
        observed_std  = observed_std + tImgs.^2;
        
        % We have to copy and unzip tImgs because octave cant deal with .gz
        [~, subject_mask_name, subject_mask_ext] = fileparts(subject_mask);
        % Copying it to the groun_truth_dir where I can unzip and delete it
        copyfile(subject_mask, pwd);
        gunzip(fullfile(pwd, [subject_mask_name subject_mask_ext]));
        
        subject_mask = spm_vol(fullfile(pwd, subject_mask_name));
        subject_mask = spm_read_vols(subject_mask);
        subject_mask = reshape(subject_mask, [prod(dim), 1]);
        
        intersection_mask = intersection_mask + subject_mask; 
        
        % Now we're done with the cope and mask files, we delete them
        delete(fullfile(pwd, subject_cope_name));
        delete(fullfile(pwd, subject_mask_name));
        
      end %========== Loop i (subjects)
      
      observed_mean = observed_mean/nSubj;

      observed_std = sqrt(observed_std/nSubj - observed_mean.^2);
      
      observed_cohen_d = observed_mean./observed_std;

      % We only consider the observed data in the intersection mask
      intersection_mask = intersection_mask > nSubj - 1;
      observed_mean = observed_mean.*intersection_mask;
      observed_std = observed_std.*intersection_mask;
      observed_cohen_d = observed_cohen_d.*intersection_mask;

      % ... and we only want to consider the ground truth in the mask too
      cohen_d = cohen_d_unmasked.*reshape(intersection_mask, dim);
      AC                      = cohen_d >= thr;
      middle_contour          = AC;
      middle_contour_volume   = sum(middle_contour(:));
      cohen_d_boundary_edges   = getBdryparams(cohen_d, thr);
      n_cohen_d_boundary_edges = size(getBdryvalues(cohen_d, cohen_d_boundary_edges),1);

      transformed_observed_cohen_d = alpha*asinh(observed_cohen_d*beta);
            
      observed_boundary_edges   = getBdryparams(reshape(observed_cohen_d, dim), (1 - (3/(4*nSubj - 5)))^(-1)*thr);
      n_observed_boundary_edges = size(getBdryvalues(cohen_d, observed_boundary_edges),1);
      
      resid_boundary_values = zeros([n_cohen_d_boundary_edges nSubj]);
      observed_resid_boundary_values = zeros([n_observed_boundary_edges nSubj]);
      
      for i=1:nSubj
          % We mask each subjects data with the intersection mask
          observed_data(:,i) = observed_data(:,i).*intersection_mask;    
          cohen_resid = create_resid(observed_data(:,i), observed_mean, observed_std, 2); 
          subject_resid_field                 = reshape(cohen_resid, [dim 1]);
          resid_boundary_values(:,i)          = getBdryvalues(subject_resid_field, cohen_d_boundary_edges);
          observed_resid_boundary_values(:,i) = getBdryvalues(subject_resid_field, observed_boundary_edges);
      end
      
      % Standardizing the residuals on the boundary after the interpolation
      observed_boundary_std = sqrt(1 + thr^2*(beta^2));
      resid_boundary_values = (alpha*beta./observed_boundary_std)*resid_boundary_values;
      observed_resid_boundary_values = (alpha*beta./observed_boundary_std)*observed_resid_boundary_values;
      
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          % Applying the bootstrap using Rademacher variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;

          % True boundary
          boundary_bootstrap                = resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
          boundary_resid_field              = sum(boundary_bootstrap, 2)/sqrt(nSubj);
          % Re-standardizing by bootstrap standard deviation
          boot_std                          = std(boundary_bootstrap, 0, 2);
          boundary_resid_field              = boundary_resid_field./boot_std;
          supG_raw(k)                       = max(abs(boundary_resid_field));
          
          % Estimated boundary
          observed_boundary_bootstrap       = observed_resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
          observed_boundary_resid_field     = sum(observed_boundary_bootstrap, 2)/sqrt(nSubj);
          % Re-standardizing by bootstrap standard deviation
          observed_boot_std                 = std(observed_boundary_bootstrap, 0, 2);
          observed_boundary_resid_field     = observed_boundary_resid_field./observed_boot_std;
          supG_observed(k)                  = max(abs(observed_boundary_resid_field));
      end 
    
    transformed_observed_cohen_d     = reshape(transformed_observed_cohen_d, dim);
        
    supGa_raw_80                     = prctile(supG_raw, 80);
    supGa_raw_90                     = prctile(supG_raw, 90);
    supGa_raw_95                     = prctile(supG_raw, 95);
    
    % Gaussian random variable results for the true and estimated boundary
    % True boundary
    lower_contour_raw_80             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_raw_80*tau;
    upper_contour_raw_80             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_raw_80*tau;
    lower_contour_raw_80_volume_prct = sum(lower_contour_raw_80(:))/middle_contour_volume;
    upper_contour_raw_80_volume_prct = sum(upper_contour_raw_80(:))/middle_contour_volume;
    mid_on_upper_raw_80              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80              = middle_contour.*lower_contour_raw_80;
    upper_subset_mid_raw_80          = upper_contour_raw_80 - mid_on_upper_raw_80;
    mid_subset_lower_raw_80          = middle_contour - lower_on_mid_raw_80;
    
    lower_contour_raw_90             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_raw_90*tau;
    upper_contour_raw_90             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_raw_90*tau;
    lower_contour_raw_90_volume_prct = sum(lower_contour_raw_90(:))/middle_contour_volume;
    upper_contour_raw_90_volume_prct = sum(upper_contour_raw_90(:))/middle_contour_volume;
    mid_on_upper_raw_90              = upper_contour_raw_90.*middle_contour;
    lower_on_mid_raw_90              = middle_contour.*lower_contour_raw_90;
    upper_subset_mid_raw_90          = upper_contour_raw_90 - mid_on_upper_raw_90;
    mid_subset_lower_raw_90          = middle_contour - lower_on_mid_raw_90;    
    
    lower_contour_raw_95             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_raw_95*tau;
    upper_contour_raw_95             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_raw_95*tau;
    lower_contour_raw_95_volume_prct = sum(lower_contour_raw_95(:))/middle_contour_volume;
    upper_contour_raw_95_volume_prct = sum(upper_contour_raw_95(:))/middle_contour_volume;
    mid_on_upper_raw_95              = upper_contour_raw_95.*middle_contour;
    lower_on_mid_raw_95              = middle_contour.*lower_contour_raw_95;
    upper_subset_mid_raw_95          = upper_contour_raw_95 - mid_on_upper_raw_95;
    mid_subset_lower_raw_95          = middle_contour - lower_on_mid_raw_95;
    
    % Observed boundary
    supGa_observed_80                     = prctile(supG_observed, 80);
    supGa_observed_90                     = prctile(supG_observed, 90);
    supGa_observed_95                     = prctile(supG_observed, 95);
       
    lower_contour_observed_80             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_observed_80*tau;
    upper_contour_observed_80             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_observed_80*tau;
    lower_contour_observed_80_volume_prct = sum(lower_contour_observed_80(:))/middle_contour_volume;
    upper_contour_observed_80_volume_prct = sum(upper_contour_observed_80(:))/middle_contour_volume;
    mid_on_upper_observed_80              = upper_contour_observed_80.*middle_contour;
    lower_on_mid_observed_80              = middle_contour.*lower_contour_observed_80;
    upper_subset_mid_observed_80          = upper_contour_observed_80 - mid_on_upper_observed_80;
    mid_subset_lower_observed_80          = middle_contour - lower_on_mid_observed_80;
    
    lower_contour_observed_90             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_observed_90*tau;
    upper_contour_observed_90             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_observed_90*tau;
    lower_contour_observed_90_volume_prct = sum(lower_contour_observed_90(:))/middle_contour_volume;
    upper_contour_observed_90_volume_prct = sum(upper_contour_observed_90(:))/middle_contour_volume;
    mid_on_upper_observed_90              = upper_contour_observed_90.*middle_contour;
    lower_on_mid_observed_90              = middle_contour.*lower_contour_observed_90;
    upper_subset_mid_observed_90          = upper_contour_observed_90 - mid_on_upper_observed_90;
    mid_subset_lower_observed_90          = middle_contour - lower_on_mid_observed_90;    
    
    lower_contour_observed_95             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_observed_95*tau;
    upper_contour_observed_95             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_observed_95*tau;
    lower_contour_observed_95_volume_prct = sum(lower_contour_observed_95(:))/middle_contour_volume;
    upper_contour_observed_95_volume_prct = sum(upper_contour_observed_95(:))/middle_contour_volume;
    mid_on_upper_observed_95              = upper_contour_observed_95.*middle_contour;
    lower_on_mid_observed_95              = middle_contour.*lower_contour_observed_95;
    upper_subset_mid_observed_95          = upper_contour_observed_95 - mid_on_upper_observed_95;
    mid_subset_lower_observed_95          = middle_contour - lower_on_mid_observed_95;

    %
    % Storing all variables of interest
    %
    % True boundary variables
    supG_raw_store(:,t)                                    = supG_raw;
    threshold_raw_80_store(t)                              = supGa_raw_80;
    lower_contour_raw_80_volume_prct_store(t)              = lower_contour_raw_80_volume_prct;
    upper_contour_raw_80_volume_prct_store(t)              = upper_contour_raw_80_volume_prct;
 
    threshold_raw_90_store(t)                              = supGa_raw_90;
    lower_contour_raw_90_volume_prct_store(t)              = lower_contour_raw_90_volume_prct;
    upper_contour_raw_90_volume_prct_store(t)              = upper_contour_raw_90_volume_prct;

    threshold_raw_95_store(t)                              = supGa_raw_95;
    lower_contour_raw_95_volume_prct_store(t)              = lower_contour_raw_95_volume_prct;
    upper_contour_raw_95_volume_prct_store(t)              = upper_contour_raw_95_volume_prct;
    
    % Observed boundary variables
    supG_observed_store(:,t)                                    = supG_observed;
    threshold_observed_80_store(t)                              = supGa_observed_80;
    lower_contour_observed_80_volume_prct_store(t)              = lower_contour_observed_80_volume_prct;
    upper_contour_observed_80_volume_prct_store(t)              = upper_contour_observed_80_volume_prct;
 
    threshold_observed_90_store(t)                              = supGa_observed_90;
    lower_contour_observed_90_volume_prct_store(t)              = lower_contour_observed_90_volume_prct;
    upper_contour_observed_90_volume_prct_store(t)              = upper_contour_observed_90_volume_prct;

    threshold_observed_95_store(t)                              = supGa_observed_95;
    lower_contour_observed_95_volume_prct_store(t)              = lower_contour_observed_95_volume_prct;
    upper_contour_observed_95_volume_prct_store(t)              = upper_contour_observed_95_volume_prct;
    
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the true boundary
    lower_condition_80 = transformed_thr - second_deriv_term - supGa_raw_80*tau;
    upper_condition_80 = transformed_thr - second_deriv_term + supGa_raw_80*tau;
    lower_condition_90 = transformed_thr - second_deriv_term - supGa_raw_90*tau;
    upper_condition_90 = transformed_thr - second_deriv_term + supGa_raw_90*tau;
    lower_condition_95 = transformed_thr - second_deriv_term - supGa_raw_95*tau;
    upper_condition_95 = transformed_thr - second_deriv_term + supGa_raw_95*tau;
    
    transformed_observed_cohen_d_true_boundary_values = getBdryvalues(transformed_observed_cohen_d, cohen_d_boundary_edges);
    
    lower_condition_80_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_80;
    upper_condition_80_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_80;
    
    lower_condition_90_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_90;
    upper_condition_90_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_90;
    
    lower_condition_95_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_95;
    upper_condition_95_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_95;
                              
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the observed boundary
    lower_condition_80_observed = transformed_thr - second_deriv_term - supGa_observed_80*tau;
    upper_condition_80_observed = transformed_thr - second_deriv_term + supGa_observed_80*tau;
    lower_condition_90_observed = transformed_thr - second_deriv_term - supGa_observed_90*tau;
    upper_condition_90_observed = transformed_thr - second_deriv_term + supGa_observed_90*tau;
    lower_condition_95_observed = transformed_thr - second_deriv_term - supGa_observed_95*tau;
    upper_condition_95_observed = transformed_thr - second_deriv_term + supGa_observed_95*tau;
    
    lower_condition_80_observed_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_80_observed;
    upper_condition_80_observed_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_80_observed;
    
    lower_condition_90_observed_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_90_observed;
    upper_condition_90_observed_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_90_observed;
    
    lower_condition_95_observed_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_95_observed;
    upper_condition_95_observed_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_95_observed;
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
    % binarized sets for residuals on the true boundary in mult. bootstrap
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:))==0
      subset_success_vector_raw_80(t) = 1;
      fprintf('raw nominal 80 success! \n');
    else 
      subset_success_vector_raw_80(t) = 0; 
      fprintf('raw nominal 80 failure! \n');
    end
    
    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:))==0
      subset_success_vector_raw_90(t) = 1;
      fprintf('raw nominal 90 success! \n');
    else 
      subset_success_vector_raw_90(t) = 0; 
      fprintf('raw nominal 90 failure! \n');
    end

    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:))==0
      subset_success_vector_raw_95(t) = 1;
      fprintf('raw nominal 95 success! \n');
    else 
      subset_success_vector_raw_95(t) = 0; 
      fprintf('raw nominal 95 failure! \n');
    end
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
    % binarized sets for residuals on the observed boundary in mult. bootstrap
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:))==0
      subset_success_vector_observed_80(t) = 1;
      fprintf('observed nominal 80 success! \n');
    else 
      subset_success_vector_observed_80(t) = 0; 
      fprintf('observed nominal 80 failure! \n');
    end
    
    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:))==0
      subset_success_vector_observed_90(t) = 1;
      fprintf('observed nominal 90 success! \n');
    else 
      subset_success_vector_observed_90(t) = 0; 
      fprintf('observed nominal 90 failure! \n');
    end

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:))==0
      subset_success_vector_observed_95(t) = 1;
      fprintf('observed nominal 95 success! \n');
    else 
      subset_success_vector_observed_95(t) = 0; 
      fprintf('observed nominal 95 failure! \n');
    end
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the true boundary
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:))+sum(upper_condition_80_success)+sum(lower_condition_80_success)==0
      subset_success_vector_raw_80_alternate(t) = 1;
      fprintf('raw nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_80_alternate(t) = 0; 
      fprintf('raw nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:))+sum(upper_condition_90_success)+sum(lower_condition_90_success)==0
      subset_success_vector_raw_90_alternate(t) = 1; 
      fprintf('raw nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_90_alternate(t) = 0; 
      fprintf('raw nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:))+sum(upper_condition_95_success)+sum(lower_condition_95_success)==0
      subset_success_vector_raw_95_alternate(t) = 1; 
      fprintf('raw nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_95_alternate(t) = 0; 
      fprintf('raw nominal 95 alternate true boundary failure! \n');
    end 
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the observed boundary
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:))+sum(upper_condition_80_observed_success)+sum(lower_condition_80_observed_success)==0
      subset_success_vector_observed_80_alternate(t) = 1;
      fprintf('observed nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_80_alternate(t) = 0; 
      fprintf('observed nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:))+sum(upper_condition_90_observed_success)+sum(lower_condition_90_observed_success)==0
      subset_success_vector_observed_90_alternate(t) = 1; 
      fprintf('observed nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_90_alternate(t) = 0; 
      fprintf('observed nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:))+sum(upper_condition_95_observed_success)+sum(lower_condition_95_observed_success)==0
      subset_success_vector_observed_95_alternate(t) = 1; 
      fprintf('observed nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_95_alternate(t) = 0; 
      fprintf('observed nominal 95 alternate true boundary failure! \n');
    end 
                              
end

percentage_success_vector_raw_80                         = mean(subset_success_vector_raw_80, 1);
percentage_success_vector_raw_90                         = mean(subset_success_vector_raw_90, 1);
percentage_success_vector_raw_95                         = mean(subset_success_vector_raw_95, 1);

percentage_success_vector_observed_80                    = mean(subset_success_vector_observed_80, 1);
percentage_success_vector_observed_90                    = mean(subset_success_vector_observed_90, 1);
percentage_success_vector_observed_95                    = mean(subset_success_vector_observed_95, 1);

percentage_success_vector_raw_80_alternate               = mean(subset_success_vector_raw_80_alternate, 1);
percentage_success_vector_raw_90_alternate               = mean(subset_success_vector_raw_90_alternate, 1);
percentage_success_vector_raw_95_alternate               = mean(subset_success_vector_raw_95_alternate, 1);

percentage_success_vector_observed_80_alternate          = mean(subset_success_vector_observed_80_alternate, 1);
percentage_success_vector_observed_90_alternate          = mean(subset_success_vector_observed_90_alternate, 1);
percentage_success_vector_observed_95_alternate          = mean(subset_success_vector_observed_95_alternate, 1);

eval(['save ' SvNm ' nSubj nRlz rimFWHM thr nBoot '... 
      'threshold_raw_80_store threshold_raw_90_store threshold_raw_95_store threshold_observed_80_store threshold_observed_90_store threshold_observed_95_store '...
      'subset_success_vector_raw_80 subset_success_vector_raw_90 subset_success_vector_raw_95 subset_success_vector_observed_80 subset_success_vector_observed_90 subset_success_vector_observed_95 subset_success_vector_raw_80_alternate subset_success_vector_raw_90_alternate subset_success_vector_raw_95_alternate subset_success_vector_observed_80_alternate subset_success_vector_observed_90_alternate subset_success_vector_observed_95_alternate '...
      'percentage_success_vector_raw_80 percentage_success_vector_raw_90 percentage_success_vector_raw_95 percentage_success_vector_observed_80 percentage_success_vector_observed_90 percentage_success_vector_observed_95 percentage_success_vector_raw_80_alternate percentage_success_vector_raw_90_alternate percentage_success_vector_raw_95_alternate percentage_success_vector_observed_80_alternate percentage_success_vector_observed_90_alternate percentage_success_vector_observed_95_alternate '...
      'supG_raw_store supG_observed_store '...
      'middle_contour_volume '...
      'lower_contour_raw_80_volume_prct_store lower_contour_raw_90_volume_prct_store lower_contour_raw_95_volume_prct_store lower_contour_observed_80_volume_prct_store lower_contour_observed_90_volume_prct_store lower_contour_observed_95_volume_prct_store '...
      'upper_contour_raw_80_volume_prct_store upper_contour_raw_90_volume_prct_store upper_contour_raw_95_volume_prct_store upper_contour_observed_80_volume_prct_store upper_contour_observed_90_volume_prct_store upper_contour_observed_95_volume_prct_store'])
