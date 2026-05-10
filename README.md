<div align="center">
    <img src="docs/ProtNetEnrich_logo.png" width="50%"/>
</div>
 
> **Transforming protein lists into biological insights through automated network topology and enrichment analysis.**
 
[![Python](https://img.shields.io/badge/Python-3.8%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-22c55e?style=flat-square)](LICENSE)
[![STRING DB](https://img.shields.io/badge/Database-STRING-f97316?style=flat-square)](https://string-db.org/)
[![ToppGene](https://img.shields.io/badge/Enrichment-ToppGene-8b5cf6?style=flat-square)](https://toppgene.cchmc.org/)
 
---
 
## What is ProtNetEnrich?
 
**ProtNetEnrich** is a bioinformatics pipeline that transforms raw protein lists into biologically meaningful insights. It queries the [STRING database](https://string-db.org/) to build protein–protein interaction (PPI) networks, identifies topological **hub proteins**, and performs both **over-representation enrichment analysis**, using the [Toppgene](https://toppgene.cchmc.org/) API, and **Gene Set Enrichment Analysis (GSEA)** — all in a single, streamlined pipeline.

Whether investigating complex biological phenomena or specific disease states, ProtNetEnrich provides a comprehensive perspective on the interactions of your target proteins.
 
---
 
## ✨ Features
 
| Module | Description |
|---|---|
| 🔗 **PPI Network Construction** | Fetches phisical protein-to-protein interaction data directly from STRING DB |
| 🎯 **Hub Detection** | Identifies network hubs using topological metrics |
| 📊 **Enrichment Analysis** | Hubs enrichment analysis using the toppgene.org API |
| 🔬 **GSEA** | Gene Set Enrichment Analysis with ranked candidate hubs based on topological metrics |
| 📈 **Visualization** | Interactive network, diagnotic plots for hubs identification, GSEA running score curves |
 
---
