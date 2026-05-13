@echo off
echo Starting PPI-Network Analyzer pipeline...

:: 1. Retrieve PPI data from STRING-DB
echo Colleting data from STRING-DB.org...
python src/ppi_data_retrieval.py query_proteins.txt output

:: 2. Build the network graph
echo Building network...
python src/network_generator.py output/ppi_data.json output query_proteins.txt

:: 3. Analyze network and identify hub proteins
echo Starting network analysis...
python -X utf8 src/network_analysis.py output/raw_network.gml output --all_metrics

:: 4. Run enrichment analysis on candidate hubs
echo Starting hubs enrichment analysis...
python src/enrichment_analysis.py output/hubs.cand.json output

:: 5. Run GSEA using the Borda ranked hubs
echo Starting hubs GSEA...
python src/GSEA.py output/borda_ranking.json output GO_26

echo "*** Process completed ***"
pause