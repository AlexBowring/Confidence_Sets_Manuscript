nSubj=60;
nRlz=30;

cd /storage/maullz/Confidence_Sets_Manuscript/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Confidence_Sets_Manuscript/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_07(nSubj,['Sim_07_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
