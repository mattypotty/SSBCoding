# Recording date and time of inference
INF=$LSB_JOBNAME
echo "Performing inference $INF"
START=$( date '+%F_%H:%M:%S' )
echo $START

echo "dataset = \"$1\"; source(\"/storage1/fs1/michael.landis/Active/SSBCoding/code/SSE.Rev\")" | rb

# Recoding end date and time of inference
echo "Finished inference $INF"
END=$( date '+%F_%H:%M:%S' )
echo $END

T=$(printf '\t')
echo "$INF$T$START$T$END" >> /storage1/fs1/michael.landis/Active/SSBCoding/joblogs/run_log.txt

# OPTIONS:
# DATASET > Afrotheria Artiodactyla_Aquatic Artiodactyla_Terrestrial Basal_Euarchontoglires Carnivora_Caniformia Carnivora_Feliformia Chiroptera_Pteropodiformes Chiroptera_Vespertilioniformes Eulipotyphla Lagomorpha Marsupalia Panperissodactyla Primates Rodentia_Hystricomorpha Rodentia_Sciuromorpha Rodentia_Supramyomorpha_Myomorphi Rodentia_Supramyomorpha_non-Myomorphi Xenarthra
