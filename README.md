<div align="center">
    <img src="imgs/ProtNetEnrich_logo.png" width="50%"/>

**From protein lists into biological interpretation through automated network topology and enrichment analysis.**

[![Python](https://img.shields.io/badge/Python-3.13-f97316?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![STRING DB](https://img.shields.io/badge/Database-STRING.org-3776AB?style=flat-square)](https://string-db.org/)
[![ToppGene](https://img.shields.io/badge/Enrichment-ToppGene.org_&_GSEApy-8b5cf6?style=flat-square)](https://toppgene.cchmc.org/)
 </div>

---

## What is ProtNetEnrich?

**ProtNetEnrich** is a bioinformatics pipeline that transforms raw protein lists into biologically meaningful insights. It queries the [STRING database](https://string-db.org/) to build protein–protein interaction (PPI) networks, identifies topological **hub proteins**, and performs both **over-representation enrichment analysis**, using the [Toppgene](https://toppgene.cchmc.org/) API, and **Gene Set Enrichment Analysis (GSEA)** — all in a single, streamlined pipeline.

Whether investigating complex biological phenomena or specific disease states, ProtNetEnrich provides a comprehensive perspective on the interactions of your target proteins.

---

## Features

| Feature | Description |
|---|---|
| 🔗 **PPI Network Generation** | Fetches phisical protein-protein interaction data directly from STRING-DB.org |
| 🎯 **Hub Detection** | Identifies network hubs using topological metrics |
| 📊 **Enrichment Analysis** | Hubs enrichment analysis using the toppgene.org API |
| 🔎 **GSEA** | Gene Set Enrichment Analysis with ranked candidate hubs based on topological metrics |
| 📈 **Visualization** | Interactive network, diagnotic plots for hubs identification, GSEA running score curves |

---

## Table of Contents

- [Getting Started](#getting-started)
- [Pipeline Overview](#pipeline-overview)
- [Scripts](#scripts)
  - [1. PPI Data Retrieval](#1-ppi-data-retrieval-ppi_data_retrievalpy)
  - [2. Network Generation](#2-network-generation-network_generatorpy)
  - [3. Network Analysis](#3-network-analysis-network_analysispy)
  - [4. Enrichment Analysis](#4-enrichment-analysis-enrichment_analysispy)
  - [5. Gene Set Enrichment Analysis (GSEA)](#5-gene-set-enrichment-analysis-gseapy)
- [Examples/Tests](#examplestests)
- [References](#references)

---

## Getting Started

### Requirements

Install all dependencies:

```shell
pip install -r requirements.txt
```

### Quick Start

Run the full pipeline using `run_pipeline.sh `, setting the working directory and the parameter for each step.

Alternatively, run pipeline end-to-end with the following steps:

```shell
# 1. Retrieve PPI data from STRING-DB
python3 ppi_data_retrieval.py query_proteins.txt output

# 2. Build the network graph
python3 network_generator.py output/ppi_data.json output query_proteins.txt

# 3. Analyze network and identify hub proteins
python3 network_analysis.py output/raw_network.gml output --all_metrics

# 4. Run enrichment analysis on candidate hubs
python3 enrichment_analysis.py output/hubs.cand.json output

# 5. Run GSEA using the Borda ranked hubs
python3 GSEA.py output/borda_ranking.json output GO_26
```

> **Note:** `query_proteins.txt` should contain one HGNC gene symbol per line.

#### Windows only

To run this code on Windows machine use `run_pipeline_windows_only` as template.

> ⚠️ DO NOT remove `-X utf8` in the network analysis (Step 3): it is essential to generate correctly the `interactive_network.html`

---

## Pipeline Overview

```
query_proteins.txt
        │
        ▼
ppi_data_retrieval.py  ──►  ppi_data.json
        │
        ▼
network_generator.py   ──►  raw_network.gml
        │
        ▼
network_analysis.py    ──►  hubs.cand.json / borda_ranking.json / interactive_network.html
        │
        ├──► enrichment_analysis.py  ──►  enrichment_analysis/ (CSV + Excel)
        │
        └──► GSEA.py                 ──►  GSEA/ (GSEA_report.csv)
```

---

## Scripts

### 1. PPI Data Retrieval (`ppi_data_retrieval.py`)

This script provides a programmatic interface to the [STRING-DB API](https://string-db.org) to retrieve physical Protein-Protein Interaction (PPI) data. It is designed to automate the process of querying interaction partners for a specific list of proteins (using HGNC symbols) and saving the resulting interaction network in a structured JSON format.

**Syntax**

```shell
python3 ppi_data_retrieval.py <query_proteins.txt> <output_directory> [flags]
```

**Arguments**

| Argument | Description |
|---|---|
| `query_proteins.txt` | Path to the `.txt` file containing the list of proteins (HGNC symbols only). |
| `output_directory` | The folder where the results will be saved. |

**Optional Flags**

| Flag | Parameter Name | Description | Default |
|---|---|---|---|
| `-il` | `interactors_limit` | Maximum number of interaction partners to retrieve per protein. | `30` |
| `-tx` | `taxon` | NCBI taxonomy ID of the species (e.g., `9606` for Human). | `9606` |

**Output**

| File | Description |
|---|---|
| `ppi_data.json` | Protein-protein interaction data from STRING-DB. |

---

### 2. Network Generation (`network_generator.py`)

Processes the PPI data and constructs a biological network graph using [NetworkX](https://networkx.org). Query proteins are identified and highlighted within the network. The graph is exported in `.gml` format, compatible with tools like [Cytoscape](https://cytoscape.org).

**Syntax**

```shell
python3 network_generator.py <ppi_data.json> <output_folder> <query_proteins.txt>
```

**Arguments**

| Argument | Description |
|---|---|
| `ppi_data.json` | Path to the STRING PPI data file. |
| `output_folder` | The folder where the resulting network will be saved. |
| `query_proteins.txt` | The original `.txt` file of query proteins. |

**Output**

| File | Description |
|---|---|
| `raw_network.gml` | Network graph built on the interaction data. |

---

### 3. Network Analysis (`network_analysis.py`)

Analyzes the network topology and identifies candidate hub proteins  within the PPI network. The prioritization is based on statistical outlier detection with variable thresholding across multiple topological metrics (clustering coefficient, degree centrality, closeness centrality and betweenness centrality). Results are aggregated with the **Borda Count ranking** method to produce a robust, consensus-based ranking of candidate proteins. An interactive network visualization is also generated using [Pyvis](https://pyvis.readthedocs.io).

#### Statistical Prioritization of Hubs

For each selected topological metric, the script:

1. Calculates the metric value for all nodes in the network.
2. Identifies nodes exceeding a threshold defined as `mean + (n_stds × standard deviation)`, where `n_stds` is user-defined.

This approach assumes that true hub proteins follow a statistically different distribution compared to non-hubs (based on Chebyshev's inequality).

**Available metrics:**
- Clustering Coefficient (`-ccoef`)
- Degree Centrality (`-dc`)
- Closeness Centrality (`-cc`)
- Betweenness Centrality (`-bc`)

#### Borda Count Ranking

To reduce bias from any single metric, the script employs the Borda Count method:

1. Each protein is ranked by its value across individual metrics.
2. Points are assigned based on rank position.
3. The final consolidated list prioritizes proteins that consistently rank highly across multiple topological dimensions.

#### Diagnostic Plots

The script automatically generates diagnostic plots to validate network properties:

- **Degree Distribution** — visualizes the distribution of node degrees.
- **Metric Distributions** — plots topological metric score distributions to verify hub selection thresholds.
- **Interactive Network** — an HTML interface for visual network inspection. Green nodes = original query proteins; red nodes = candidate hubs. Edge width reflects the STRING interaction score.

**Syntax**

```shell
python3 network_analysis.py <raw_network.gml> <output_directory> [options]
```

**Arguments**

| Argument | Description |
|---|---|
| `raw_network.gml` | Path to the generated network. |
| `output_directory` | The folder where the results will be saved. |

**Options**

| Flag | Description | Default |
|---|---|---|
| `-stds <float>` | Number of standard deviations for the hub identification threshold. | `2` |
| `--all_metrics` | Use all available topology metrics. | `False` |
| `--remove_dead_ends` | Remove nodes with degree == 1 from the interactive visualization for a cleaner view. | `False` |
| `-ccoef` | Use clustering coefficient for hub identification. | `False` |
| `-dc` | Use degree centrality for hub identification. | `False` |
| `-cc` | Use closeness centrality for hub identification. | `False` |
| `-bc` | Use betweenness centrality for hub identification. | `False` |

**Output**

| File / Folder | Description |
|---|---|
| `diagnostic_plots/` | Folder containing all diagnostic visualizations. |
| `interactive_network.html` | Interactive HTML file for visual inspection of the network. |
| `hubs.cand.json` | Metadata of the candidate hub proteins. |
| `borda_ranking.json` | Borda ranking of hubs with scores. |
| `network.gml` | Final network with all node/edge information (excluding Borda scores). |

---

### 4. Enrichment Analysis (`enrichment_analysis.py`)

Runs an automated gene enrichment analysis pipeline. Ensembl identifiers are converted to Entrez Gene IDs via [gProfiler](https://biit.cs.ut.ee/gprofiler/), then the [ToppGene API](https://toppgene.cchmc.org) is queried across multiple biological categories (Gene Ontology, Pathways, Disease, etc.). Statistically significant results are exported for downstream analysis.

**Syntax**

```shell
python3 enrichment_analysis.py <hubs.cand.json> <output_directory>
```

**Arguments**

| Argument | Description |
|---|---|
| `hubs.cand.json` | Path to the candidate hubs JSON file. |
| `output_directory` | The folder where the results will be saved. |

**Output**

Results are saved in an `enrichment_analysis/` folder:

| File | Description |
|---|---|
| `<category>.csv` | Individual CSV file for each ToppGene category. |
| `enrichment_results.xlsx` | Excel file aggregating all categories for simplified inspection. |

---

### 5. Gene Set Enrichment Analysis (`GSEA.py`)

Performs pre-ranked Gene Set Enrichment Analysis using the [gseapy](https://gseapy.readthedocs.io) library. It takes a ranked gene list (derived from Borda Count ranking) and compares it against predefined or custom gene sets.

#### Predefined Gene Sets

| Gene Set ID | Description |
|---|---|
| `GO_26` | All gene sets in Gene Ontology 2026. |
| `KEGG_26` | Gene sets from KEGG 2026. |
| `Disease` | Expanded OMIM and Augmented Orphanet 2021 gene sets. |

**Syntax**

```shell
python3 GSEA.py <borda_ranking.json> <output_directory> <gene_set> [options]
```

**Arguments**

| Argument | Description |
|---|---|
| `borda_ranking.json` | Path to the Borda ranking JSON file. |
| `output_directory` | The folder where the results will be saved. |
| `gene_set` | Predefined gene set to use (e.g., `GO_26`, `KEGG_26`, `Disease`) or custom gene set (using `--custom_set`). |

**Options**

| Flag | Description |
|---|---|
| `--save_plots` | Save GSEA running sum plots to the output directory. |
| `--custom_set` | Use specific supported gseapy gene set |

**Output**

Results are saved in a `GSEA/` folder:

| File / Folder | Description |
|---|---|
| `GSEA_report.csv` | Filtered by statistical significance (`FDR q-value`), sorted by `NES`. |
| `plots/` | GSEA running sum plots (only if `--save_plots` is used). |
| Raw gseapy output files | Full gseapy output for further inspection. |

> **Note:** A `LookupError` during GSEA execution typically indicates that no statistically significant enrichment results were found for the provided gene set. Consider adjusting the input gene list, selecting a different gene set, or relaxing the significance threshold in `GSEA.py`.

## Examples/Tests

Two example use cases for ProtNetEnrich are provided in `tests/`. The `ATGHub/` folder contains full pipeline results from a study focused on the identification of novel candidate genes involved in the autophagic process. The `small_network/` folder contains pipeline results computed on a curated small set of immunological genes, and serves as a lightweight reference for testing and exploration.

## References
 
- **STRING-DB** — https://string-db.org
- **NetworkX** — https://networkx.org
- **Pyvis** — https://pyvis.readthedocs.io
- **Cytoscape** — https://cytoscape.org
- **gProfiler** — https://biit.cs.ut.ee/gprofiler
- **ToppGene Suite** — https://toppgene.cchmc.org
- **GSEApy** — https://gseapy.readthedocs.io