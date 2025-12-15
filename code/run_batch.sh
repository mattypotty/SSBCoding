RUN_LIST=$(seq -w 001 001)

for i in ${RUN_LIST[@]}
do
	source run_job.sh ${i} "FALSE" "NONE" "NONE" "NONE" "NONE"
done

# OPTIONS:
# NUMBER                    > 001 : 010
# USERANDOM                 > "TRUE", "FALSE"
# ANALYSIS                  > "NONE", "MCCTREE", "SPECMATS", "MCMCGLMM", "PAGEL"
# SSBTYPE                   > "NONE", "SSB", "FSSB", "MSSB"
# HYPOTHESIS (for MCMCglmm) > "NONE", "MI", "MGS", "SDT", "RDM", "AM", "IUCN"
# VARIABLE (for Pagel)      > "None", "NAC", "SPI", "TSP", "SDT", "RDM"