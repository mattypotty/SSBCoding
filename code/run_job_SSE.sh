# Set up volumes
LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active/SSBCoding:/storage1/fs1/michael.landis/Active/SSBCoding"
JOBDIR="/storage1/fs1/michael.landis/Active/SSBCoding/joblogs"

DATASET=$1

# Submit job
bsub -G compute-michael.landis \
-g /k.swiston/SSE \
-cwd /storage1/fs1/michael.landis/Active/SSBCoding/ \
-o $JOBDIR/$DATASET.stdout.txt \
-J $DATASET \
-q general \
-n 8 -M 4GB -R "rusage [mem=4GB] span[hosts=1]" \
-a 'docker(sswiston/phylo_docker:slim_amd64)' /bin/bash /storage1/fs1/michael.landis/Active/SSBCoding/code/run_SSE.sh $DATASET

# OPTIONS:
# DATASET > Afrotheria Artiodactyla_Aquatic Artiodactyla_Terrestrial Basal_Euarchontoglires Carnivora_Caniformia Carnivora_Feliformia Chiroptera_Pteropodiformes Chiroptera_Vespertilioniformes Eulipotyphla Lagomorpha Marsupalia Panperissodactyla Primates Rodentia_Hystricomorpha Rodentia_Sciuromorpha Rodentia_Supramyomorpha_Myomorphi Rodentia_Supramyomorpha_non-Myomorphi Xenarthra
