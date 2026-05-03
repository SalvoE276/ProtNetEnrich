#!/bin/bash

queryfile="tests/ATGHub_annotated_only_gene_symbol.txt"

echo Starting PPI-Network Analyzer pipeline...
# Step 1. Retrieval of Protein-to-protein interaction data from STRING-DB.org
echo Colleting data from STRING-DB.org...
python3 src/ppi_data_retrieval.py $queryfile tests/ppi_data.json -il 15
# Step 2. Build PPI network using Networkx
echo Building network...
python3 src/network_generator.py tests/ppi_data.json tests $queryfile
# Step 3. PPI network Analysis
echo Starting network analysis...
python3 src/network_analysis.py tests/raw_network.gml tests -stds 2 --all_metrics
# Step 4. Enrichment analysis