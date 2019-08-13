nSubj=60;
nRlz=5;

cd /storage/maullz/Confidence_Sets_Manuscript/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Confidence_Sets_Manuscript/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_32(nSubj,['Sim_32_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
