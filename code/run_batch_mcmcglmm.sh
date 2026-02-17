NUMBER_LIST=$(seq -w 001 100)
NUMBER_LIST_2=$(seq -w 001 010)
SSBTYPE_LIST=("SSB" "FSSB" "MSSB")
HYPOTHESIS_LIST=("MI" "MGS" "SDT" "SOC" "RDM" "AM" "IUCN")

USERANDOM="FALSE"

for NUMBER in ${NUMBER_LIST_2[@]}
do
	for SSBTYPE in ${SSBTYPE_LIST[@]}
	do
		for HYPOTHESIS in ${HYPOTHESIS_LIST[@]}
		do
			source run_job.sh ${NUMBER} ${USERANDOM} "MCMCGLMM" ${SSBTYPE} ${HYPOTHESIS} "NONE"
		done
	done
done

USERANDOM="TRUE"

for NUMBER in ${NUMBER_LIST[@]}
do
	for SSBTYPE in ${SSBTYPE_LIST[@]}
	do
		for HYPOTHESIS in ${HYPOTHESIS_LIST[@]}
		do
			source run_job.sh ${NUMBER} ${USERANDOM} "MCMCGLMM" ${SSBTYPE} ${HYPOTHESIS} "NONE"
		done
	done
done

# OPTIONS:
# NUMBER                    > 001 : 100
# USERANDOM                 > "TRUE", "FALSE"
# ANALYSIS                  > "NONE", "MCCTREE", "SPECMATS", "MCMCGLMM", "PAGEL"
# SSBTYPE                   > "NONE", "SSB", "FSSB", "MSSB"
# HYPOTHESIS (for MCMCglmm) > "NONE", "MI", "MGS", "SDT", "SOC", "RDM", "AM", "IUCN", "NAC"
# NEW HYPOTHESES: "NAC", "SOC", "HAR"
