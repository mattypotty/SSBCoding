NUMBER_LIST=$(seq -w 001 100)
SSBTYPE_LIST=("SSB" "FSSB" "MSSB")
VARIABLE_LIST=("NAC" "SPI" "TSP" "SDT" "RDM")

for SSBTYPE in ${SSBTYPE_LIST[@]}
do
	for VARIABLE in ${VARIABLE_LIST[@]}
	do
		source run_job.sh 001 "FALSE" "PAGEL" ${SSBTYPE} "NONE" ${VARIABLE}
	done
done

for NUMBER in ${NUMBER_LIST[@]}
do
	for SSBTYPE in ${SSBTYPE_LIST[@]}
	do
		for VARIABLE in ${VARIABLE_LIST[@]}
		do
			source run_job.sh ${NUMBER} "TRUE" "PAGEL" ${SSBTYPE} "NONE" ${VARIABLE}
		done
	done
done

# OPTIONS:
# NUMBER                    > 001 : 100
# USERANDOM                 > "TRUE", "FALSE"
# ANALYSIS                  > "NONE", "MCCTREE", "SPECMATS", "MCMCGLMM", "PAGEL"
# SSBTYPE                   > "NONE", "SSB", "FSSB", "MSSB"
# VARIABLE (for Pagel)      > "None", "NAC", "SPI", "TSP", "SDT", "RDM", "HAR"
# NEW VARIABLES: "HAR"
