# Recording date and time of inference
INF=$LSB_JOBNAME
echo "Performing inference $INF"
START=$( date '+%F_%H:%M:%S' )
echo $START

Rscript ./code/analyses.R $1 $2 $3 $4 $5 $6

# Recoding end date and time of inference
echo "Finished inference $INF"
END=$( date '+%F_%H:%M:%S' )
echo $END

T=$(printf '\t')
echo "$INF$T$START$T$END" >> /storage1/fs1/michael.landis/Active/SSBcoding/joblogs/run_log.txt

# OPTIONS:
# NUMBER                    > 001 : 010
# USERANDOM                 > "TRUE", "FALSE"
# ANALYSIS                  > "NONE", "MCCTREE", "SPECMATS", "MCMCGLMM", "PAGEL"
# SSBTYPE                   > "NONE", "SSB", "FSSB", "MSSB"
# HYPOTHESIS (for MCMCglmm) > "NONE", "MI", "MGS", "SDT", "RDM", "AM", "IUCN"
# VARIABLE (for Pagel)      > "None", "NAC", "SPI", "TSP", "SDT", "RDM"