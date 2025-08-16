# AR Inhibition - Computational Structural Biology Analysis

This repository contains computational analysis tools for studying protein-protein interfaces in the Androgen Receptor (AR) system. The project integrates multiple structural biology approaches to identify and analyze interface residues across different protein assemblies.

## Project Overview

This computational structural biology project focuses on:
- **Proximity-based interface analysis** using multiple distance cutoffs (4Å, 6Å, 8Å)
- **ScanNet integration** for binding site probability prediction
- **FoldX alanine scanning** for stability analysis
- **PISA interface analysis** for structural insights
- **ChimeraX visualization** scripts for structural validation

## Repository Structure

```
ARInhibition/
├── utils/                          # Analysis scripts
│   ├── generate_positions.py       # PISA XML to FoldX position converter
│   ├── position_intersection.py    # FoldX results intersection analysis
│   └── refined_interface_analysis.py # Main interface analysis pipeline
├── references/                     # PDB structure files
├── Results/                        # Analysis outputs
│   ├── [structure]_manual_script/  # Interface analysis results
│   ├── ScanNet_results_[structure]/ # ScanNet prediction outputs
│   ├── AF_binder_dimer_comp/       # specific Alphafold result analysis
│   └── AlaScan/                    # FoldX alanine scanning results
├── PISA_Output/                    # PISA interface analysis files
└── ChimeraXScripts/                # Visualization scripts
```

## Analysis Workflow

### 1. Interface Residue Identification

The `refined_interface_analysis.py` script performs proximity-based interface analysis:

- Calculates interface residues using multiple distance cutoffs (4Å, 6Å, 8Å)
- Integrates ScanNet binding site predictions
- Generates consensus residue lists across assemblies
- Creates ChimeraX visualization scripts

**Usage:**
```bash
python refined_interface_analysis.py
```

**Configuration:**
Edit the script variables to specify:
- `PDB_FILE_BASE`: Path to PDB structure files
- `ASSEMBLY_NUMBERS`: Which assemblies to analyze
- `SCAN_NET_FILE`: Path to ScanNet predictions
- `DISTANCE_CUTOFFS`: Distance thresholds for interface analysis

### 2. PISA to FoldX Position Conversion

The `generate_positions.py` script converts PISA interface analysis XML files to FoldX-compatible position strings:

**Usage:**
```bash
python generate_positions.py
```

**Input:** PISA XML files (e.g., `residue0.xml`)
**Output:** Comma-separated position strings for FoldX alanine scanning

### 3. FoldX Results Analysis

The `position_intersection.py` script analyzes FoldX alanine scanning results:

- Extracts significant positions based on ΔΔG threshold
- Finds intersecting positions across multiple assemblies
- Outputs consensus positions for further analysis

**Usage:**
```bash
python position_intersection.py
```

**Configuration:**
- `THRESHOLD`: ΔΔG threshold (default: 1.0 kcal/mol)
- `INPUT_DIR`: Directory containing FoldX output files

## Data Sources

The analysis uses the following PDB structures:
- **1t5z**: https://www.rcsb.org/structure/1T5Z
- **2am9**: https://www.rcsb.org/structure/2AM9
- **5JJM**: https://www.rcsb.org/structure/5JJM

ScanNet: http://bioinfo3d.cs.tau.ac.il/ScanNet/index_real.html

## Output Files

### Interface Analysis Results
- `assembly_X_detailed_report.txt`: Detailed interface analysis with distance cutoffs
- `surface_res_consensus.csv`: Consensus interface residues across assemblies
- `assembly_X_Y_interface_residues.txt`: Interface residues for specific assembly

### ScanNet Results
- `predictions_[PDBID].csv`: Binding site probabilities
- `annotated_[PDBID].pdb`: Structure with probabilities in B-factor field
- `annotated_[PDBID].cxc`: ChimeraX visualization script

### FoldX Results
- `PS_[structure]_scanning_output.txt`: Alanine scanning ΔΔG values
- `intersecting_positions.txt`: Consensus significant positions

## Visualization

### ChimeraX Scripts
Generated scripts for structural visualization:
- `chimeraX_assemblyX_homodimer.cxc`: Interface residue visualization
- `chimeraX_annotated_high_prob_residues_assemblyX.cxc`: High-probability residue annotation
- `chimeraX_consensus_residues_assemblyX.cxc`: Consensus residue highlighting

**Usage:**
```bash
# Open ChimeraX and load the script.
chimeraX chimeraX_assembly1_homodimer.cxc
```

## Dependencies

This project requires the following tools (assumed to be installed):
- **FoldX**: Protein stability prediction
- **ScanNet**: Binding site prediction
- **ChimeraX**: Molecular visualization
- **PISA**: Protein interface analysis
- **Python 3.9+** with BioPython, pandas, numpy

## Step-by-Step Analysis

1. **Prepare PDB structures** in the `references/` directory
2. **Run ScanNet analysis** on target structures
3. **Execute interface analysis:**
   ```bash
   python refined_interface_analysis.py
   ```
4. **Generate FoldX positions:**
   ```bash
   python generate_positions.py
   ```
5. **Run FoldX alanine scanning** using generated positions
6. **Analyze FoldX results:**
   ```bash
   python position_intersection.py
   ```
7. **Visualize results** using generated ChimeraX scripts

## File Naming Conventions

- Assembly files: `[PDBID]-assembly[number].pdb`
- ScanNet results: `predictions_[PDBID].csv`
- FoldX outputs: `PS_[structure]_scanning_output.txt`
- ChimeraX scripts: `chimeraX_[analysis_type]_assembly[number].cxc`

## Notes

- Distance cutoffs are configurable in the analysis scripts
- ScanNet probability threshold is set to 0.5 for high-probability residues
- FoldX ΔΔG threshold is set to 1.0 kcal/mol for significant mutations
- All paths in scripts are relative to the repository root

## Contributing

This project was developed as part of a Computational Structural Biology course. The analysis pipeline can be extended to other protein systems by modifying the configuration variables in the utility scripts. 