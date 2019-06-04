nSubj=240;
nRlz=30;

cd /storage/maullz/Confidence_Sets_Manuscript/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Confidence_Sets_Manuscript/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_05(nSubj,['Sim_05_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
