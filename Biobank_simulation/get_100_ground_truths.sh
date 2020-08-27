#!/bin/bash
#$ -N get_100_ground_truths
#$ -q short.qe
#$ -t 1:100
#$ -o /well/nichols/users/bas627/Confidence_Sets_Manuscript/Biobank_simulation/logs
#$ -e /well/nichols/users/bas627/Confidence_Sets_Manuscript/Biobank_simulation/logs
#$ -pe shmem 2

. /etc/profile
. ~/.bashrc

module add Octave/4.4.1-foss-2018b

cd /well/nichols/users/bas627/Confidence_Sets_Manuscript/Biobank_simulation

octave get_100_ground_truths.m
