
# script to generate txt file containing all donor ID's for TOB
# used to run vireo against all donor genotypes for TOB samples 
ALL_IDS_OUT="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/all_drews_ids_tob.txt";
touch $ALL_IDS_OUT;
for file in /share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/samples_in_pools_drews_names/*.tsv;
do 
    cat "$file" >> $ALL_IDS_OUT;
done; 

sort $ALL_IDS_OUT | uniq > $ALL_IDS_OUT.temp && mv $ALL_IDS_OUT.temp $ALL_IDS_OUT  # remove duplicate donor ids 