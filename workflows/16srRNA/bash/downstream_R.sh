export_ps_figs(){
	echo -e "${YEL}########    Phyloseq Rscript running   ########${NC}"
	mkdir -p "${OUTPUT}/Phyloseq_Figures"
	Rscript ${GUAP_DIR}/bin/16s/scripts/Export_phyloseq_figs.R \
		-i "${OUTPUT}/ps.rds" \
		-n ${name} \
		-c ${condition_name} \
		-w "${OUTPUT}/Phyloseq_Figures"
}
