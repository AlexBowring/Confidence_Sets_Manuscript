rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'));
get_ground_truth(tID)