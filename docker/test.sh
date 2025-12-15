# R
echo ""; echo "Checking for R..."; echo ""
R --version | head -n 1

# R PACKAGES
echo ""; echo "Testing for R packages..."; echo ""
TR="TRUE"
for package in "this.path" "ape" "phangorn" "castor" "caper" "phytools" "diversitree" "MCMCglmm" "coda" "MCMCvis" "mcmcr" "styler"
do
    output=`echo $(R -e "\"${package}\" %in% rownames(installed.packages())") | tail -c 7 | head -c 4`
    if [ "$output" = "$TR" ]; then echo "$package installed"; else echo "$package not found"; fi
done
echo ""; echo "Full List:"; echo ""
echo $(R -e "cat(rownames(installed.packages()))") | cut -d ">" -f2- | cut -c 38- | rev | cut -c 4- | rev