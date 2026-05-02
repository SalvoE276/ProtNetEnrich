#!/bin/bash

echo Starting PPI-Network Analyzer pipeline...
# Step 1. Retrieval of Protein-to-protein interaction data from STRING-DB.org
python3 src/ppi_data_retrieval.py tests/query_proteins.txt tests/ppi_data.json -il 15