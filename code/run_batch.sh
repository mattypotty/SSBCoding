NUMBER_LIST=$(seq -w 001 001)
SSBTYPE_LIST=("SSB" "FSSB" "MSSB")
HYPOTHESIS_LIST=("MI" "MGS" "SDT" "SOC" "RDM" "AM" "IUCN")
VARIABLE_LIST=("NAC" "SPI" "TSP" "SDT" "RDM")

USERANDOM="FALSE"

for NUMBER in ${NUMBER_LIST[@]}
do
	for SSBTYPE in ${SSBTYPE_LIST[@]}
	do
		for HYPOTHESIS in ${HYPOTHESIS_LIST[@]}
		do
			source run_job.sh ${NUMBER} ${USERANDOM} "MCMCGLMM" ${SSBTYPE} ${HYPOTHESIS} "NONE"
		done
		for VARIABLE in ${VARIABLE_LIST[@]}
		do
			source run_job.sh ${NUMBER} ${USERANDOM} "PAGEL" ${SSBTYPE} "NONE" ${VARIABLE}
		done
	done
done

# OPTIONS:
# NUMBER                    > 001 : 010
# USERANDOM                 > "TRUE", "FALSE"
# ANALYSIS                  > "NONE", "MCCTREE", "SPECMATS", "MCMCGLMM", "PAGEL"
# SSBTYPE                   > "NONE", "SSB", "FSSB", "MSSB"
# HYPOTHESIS (for MCMCglmm) > "NONE", "MI", "MGS", "SDT", "SOC", "RDM", "AM", "IUCN"
# VARIABLE (for Pagel)      > "None", "NAC", "SPI", "TSP", "SDT", "RDM"
