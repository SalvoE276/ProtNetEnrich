#!/bin/bash

working_dir=tests/ATGHub

queryfile=$working_dir/query_proteins.txt
raw_interaction_data=$working_dir/ppi_data.json

echo Starting PPI-Network Analyzer pipeline...
# Step 1. Retrieval of Protein-to-protein interaction data from STRING-DB.org
echo Colleting data from STRING-DB.org...
python3 src/ppi_data_retrieval.py $queryfile $working_dir -il 15
# Step 2. Build PPI network using Networkx
echo Building network...
python3 src/network_generator.py $working_dir/ppi_data.json $working_dir $queryfile
# Step 3. PPI network Analysis
echo Starting network analysis...
python3 src/network_analysis.py $working_dir/raw_network.gml $working_dir -stds 2 --all_metrics --remove_dead_ends
# Step 4. Enrichment analysis
echo Starting hubs enrichment analysis...
python3 src/enrichment_analysis.py $working_dir/hubs.cand.json $working_dir
# Step 5. GSEA of candidate hubs
echo Starting hubs GSEA...
python3 src/GSEA.py $working_dir/hubs.cand.json $working_dir

echo "*** Process completed ***"