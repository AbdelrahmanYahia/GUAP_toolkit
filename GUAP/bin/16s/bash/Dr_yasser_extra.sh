#export the both absoulte and relative abundance 
mkdir -p exports
mkdir -p exports/abundance
###############################
echo -e "${YEL}###########      STARTED    ###########${NC}"

absolutes (){
    i=$1
    taxaa=$2
    echo -e "${YEL}##### running absoulutes on $i $taxaa ######${NC}"

    #absolute
    qiime taxa collapse  \
    --i-table table.qza \
    --i-taxonomy classify/taxonomy.qza  \
    --p-level $i  \
    --o-collapsed-table exports/abundance/level-$i-$taxaa-table.qza 
    error_cheker $?
    qiime tools export \
    --input-path exports/abundance/level-$i-$taxaa-table.qza \
    --output-path taxonomy
    error_cheker $?
    biom convert -i taxonomy/feature-table.biom -o taxonomy/abs-level-$i-$taxaa-table.tsv \
    --to-tsv
    #relative
    qiime feature-table relative-frequency \
    --i-table exports/abundance/level-$i-$taxaa-table.qza \
    --o-relative-frequency-table exports/abundance/rel-level-$i-$taxaa-table.qza
    error_cheker $?
    qiime tools export \
    --input-path exports/abundance/rel-level-$i-$taxaa-table.qza \
    --output-path taxonomy
    error_cheker $?
    biom convert -i taxonomy/feature-table.biom -o taxonomy/rel-level-$i-$taxaa-table.tsv \
    --to-tsv
    error_cheker $?
}

for i in {1..7} 
do 
    if [ $i == "1" ]; then taxaa="Kingdom"; fi
    if [ $i == "2" ]; then taxaa="Phylum"; fi
    if [ $i == "3" ]; then taxaa="Class"; fi
    if [ $i == "4" ]; then taxaa="Order"; fi
    if [ $i == "5" ]; then taxaa="Family"; fi
    if [ $i == "6" ]; then taxaa="Genus"; fi
    if [ $i == "7" ]; then taxaa="Species"; fi

    absolutes $i $taxaa &
    

done
wait


mkdir -p composition_tables
mkdir -p Differential_abundance
#Ancom Differential expression
for i in {1..7} 
do 
    if [ $i == "1" ]; then taxaa="Kingdom"; fi
    if [ $i == "2" ]; then taxaa="Phylum"; fi
    if [ $i == "3" ]; then taxaa="Class"; fi
    if [ $i == "4" ]; then taxaa="Order"; fi
    if [ $i == "5" ]; then taxaa="Family"; fi
    if [ $i == "6" ]; then taxaa="Genus"; fi
    if [ $i == "7" ]; then taxaa="Species"; fi
    echo -e "${YEL}##### Generating table $i $taxaa ######${NC}"

    qiime composition add-pseudocount  \
    --i-table table.qza   \
    --o-composition-table composition_tables/table_$taxa.qza 
    
    
    echo -e "${YEL}##### running RDA on $i $taxaa ######${NC}"

    for comparison in condition p_name 
    do
        mkdir -p Differential_abundance/$comparison
        qiime composition ancom  \
        --i-table composition_tables/table_$taxa.qza  \
        --m-metadata-file ../sample-metadata.tsv  \
        --m-metadata-column $comparison \
        --o-visualization Differential_abundance/$comparison/level-$i-$taxaa-ancom-$comparison.qzv &
    done
done
wait
echo -e "${YEL}##### DONE ######${NC}"
