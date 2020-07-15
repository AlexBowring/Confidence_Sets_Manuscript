% Directory containing the Biobank copes
basedir = '/well/nichols/projects/UKB/IMAGING/ContourInf/MNI';

% Directory to save images to
outdir = '/Users/maullz/Desktop/test';

% Path to MNI mask image;
MNI_mask_file = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz';


% List all copefiles in the directory
%%cope_files = cellstr(spm_select('FPList', basedir, ['.*\_cope5_MNI.nii.gz']));
%%mask_files = cellstr(spm_select('FPList', basedir, ['.*\_mask_MNI.nii.gz']));
cope_files = cellstr(spm_select('FPList', basedir, ['.*\cope.nii.gz']));

% Select 4000 random copes from the total 8945 available
shuffle_ids = randperm(length(cope_files)); 
%%shuffle_ids = randperm(8945, 4000);


dim = [91, 109, 91];
% Empty datamat, datamat_squared and maskmat for loading in the copes
datamat = zeros([prod(dim) length(shuffle_ids)]);
datamat_squared = zeros([prod(dim) length(shuffle_ids)]);
maskmat = zeros([prod(dim) length(shuffle_ids)]);

for i=1:length(shuffle_ids)
    VY = spm_vol(cope_files{shuffle_ids(i)});
    VM = (reshape(spm_read_vols(VY), [prod(dim) 1]) ~= 0);
    %%VM = spm_vol(mask_files{shuffle_ids(i)});
    datamat(:,i) = reshape(spm_read_vols(VY), [prod(dim) 1]);
    datamat_squared(:,i) = reshape(spm_read_vols(VY), [prod(dim) 1]).^2; 
    maskmat(:,i) = VM;
    %%maskmat(:,i) = reshape(spm_read_vols(VM), [prod(dim) 1]);
end

sum_datamat = sum(datamat,2);
sum_datamat = reshape(sum_datamat, dim);

sum_datamat_squared = sum(datamat_squared,2);
sum_datamat_squared = reshape(sum_datamat_squared, dim);

sum_maskmat = sum(maskmat,2);
sum_maskmat = reshape(sum_maskmat, dim);

sum_maskmat_minus_one = sum_maskmat - ones(dim); 

more_than_100 = (sum(maskmat,2) > 99);
more_than_100 = reshape(more_than_100, dim);

sample_mean = sum_datamat./sum_maskmat;
sample_mean(isnan(sample_mean))=0;
sample_mean = sample_mean.*more_than_100;

sample_variance = sum_datamat_squared./sum_maskmat_minus_one - (sum_datamat.^2)./(sum_maskmat.*sum_maskmat_minus_one);
sample_variance(isnan(sample_variance)) = 0;
sample_variance = sample_variance.*more_than_100;
sample_sd = sqrt(sample_variance);

sample_cohens_d = sample_mean./sample_sd;
sample_cohens_d(isnan(sample_cohens_d)) = 0;

sample_cohens_d_threshold_02 = sample_cohens_d > 0.2;
sample_cohens_d_threshold_05 = sample_cohens_d > 0.5;
sample_cohens_d_threshold_08 = sample_cohens_d > 0.8;

% Finally, mask all important images with the MNI152 mask
MNI_mask = spm_vol(MNI_mask_file);
MNI_mask = spm_read_vols(MNI_mask);

sample_mean = sample_mean.*MNI_mask;
sample_sd = sample_sd.*MNI_mask;
sample_cohens_d = sample_cohens_d.*MNI_mask;
more_than_100 = more_than_100.*MNI_mask;

% Save images of interest
% Cohen's d image
cd(outdir);
Vout = VY;
Vout.fname = 'Biobank_4000_cohens_d.nii'; % crucially, change the file name!
Vout.descrip = 'Cohens d image from 4000 Biobank subjects'; % Actually, put something more
                                        % informative here
spm_write_vol(Vout, sample_cohens_d);

% 0.2 Thresholded Cohen's d 
Vout = VY;
Vout.fname = 'Biobank_4000_cohens_d_02.nii'; % crucially, change the file name!
Vout.descrip = 'Cohens d thresholded at c = 0.2'; % Actually, put something more
                                        % informative here
spm_write_vol(Vout, sample_cohens_d_threshold_02);

% 0.5 Thresholded Cohen's d 
Vout = VY;
Vout.fname = 'Biobank_4000_cohens_d_05.nii'; % crucially, change the file name!
Vout.descrip = 'Cohens d thresholded at c = 0.5'; % Actually, put something more
                                        % informative here
spm_write_vol(Vout, sample_cohens_d_threshold_05);

% 0.8 Thresholded Cohen's d 
Vout = VY;
Vout.fname = 'Biobank_4000_cohens_d_08.nii'; % crucially, change the file name!
Vout.descrip = 'Cohens d thresholded at c = 0.8'; % Actually, put something more
                                        % informative here
spm_write_vol(Vout, sample_cohens_d_threshold_08);

% Sample mean image
Vout = VY;
Vout.fname = 'Biobank_4000_sample_mean.nii'; % crucially, change the file name!
Vout.descrip = 'sample mean image from 4000 Biobank subjects'; % Actually, put something more
                                        % informative here
spm_write_vol(Vout, sample_mean);

% Sample standard deviation image
Vout = VY;
Vout.fname = 'Biobank_4000_sample_sd.nii'; % crucially, change the file name!
Vout.descrip = 'sample standard deviation from 4000 Biobank subjects'; % Actually, put something more

spm_write_vol(Vout, sample_sd);

% Saving the random sample we choose so this sample is excluded for testing
save('4000_random_sample.mat', 'shuffle_ids')


