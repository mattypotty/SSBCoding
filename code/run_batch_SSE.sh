#DATASET_LIST=("Afrotheria" "Artiodactyla_Aquatic" "Artiodactyla_Terrestrial" "Basal_Euarchontoglires" "Carnivora_Caniformia" "Carnivora_Feliformia" "Chiroptera_Pteropodiformes" "Chiroptera_Vespertilioniformes" "Eulipotyphla" "Lagomorpha" "Marsupalia" "Panperissodactyla" "Primates" "Rodentia_Hystricomorpha" "Rodentia_Sciuromorpha" "Rodentia_Supramyomorpha_Myomorphi" "Rodentia_Supramyomorpha_non-Myomorphi" "Xenarthra")
DATASET_LIST=("Primates")

for DATASET in ${DATASET_LIST[@]}
do
source run_job_SSE.sh ${DATASET}
done

# OPTIONS:
# DATASET > Afrotheria Artiodactyla_Aquatic Artiodactyla_Terrestrial Basal_Euarchontoglires Carnivora_Caniformia Carnivora_Feliformia Chiroptera_Pteropodiformes Chiroptera_Vespertilioniformes Eulipotyphla Lagomorpha Marsupalia Panperissodactyla Primates Rodentia_Hystricomorpha Rodentia_Sciuromorpha Rodentia_Supramyomorpha_Myomorphi Rodentia_Supramyomorpha_non-Myomorphi Xenarthra
