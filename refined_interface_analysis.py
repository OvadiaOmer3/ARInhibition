#!/usr/bin/env python3.9
"""
Refined Protein-Protein Interface Analysis
Focus on Assembly 2 (heterodimer) and use appropriate distance cutoffs
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
import warnings
warnings.filterwarnings('ignore')

# SCAN_NET_FILE = "ScanNet_results_2am9Test/predictions_2am9_chainA.csv"
SCAN_NET_FILE = "ScanNet_results_1t5zTest/predictions_1t5z.csv"

PDB_FILE_BASE = "5JJM-assembly"
ASSEMBLY_NUMBERS = [1, 2, 3]
RES_A_OFFSET = 0  # Residue numbering offset for PDB files
RES_B_OFFSET = 0  # Residue numbering offset for PDB files
CHIMERA_X_SCRIPT_DIR = "ChimeraXScripts/1t5z"
RESULTS_DIR = "Results/1t5z"


# PDB_FILE_BASE = "AlphaFoldRes/AF_binder_dimer_comp"
# ASSEMBLY_NUMBERS = [""]
# RES_A_OFFSET = 669  # Residue numbering offset for PDB files
# RES_B_OFFSET = 0  # Residue numbering offset for PDB files
# CHIMERA_X_SCRIPT_DIR = "ChimeraXScripts/AF_binder_dimer_comp"
# RESULTS_DIR = "Results/AF_binder_dimer_comp"



ASSEMBLY_STRACTURES_IDS = {
    1: ('1.1', '1.2'),
    2: ('1', '1'),
    3: ('1.1', '1.2'),
    "": ('1', '1')
}

# Define chain pairs for each assembly
ASSEMBLY_CHAINS = {
    1: ('A', 'A'),  # Assembly 1: two A chains
    2: ('B', 'C'),      # Assembly 2: B and C chains
    3: ('D', 'D'),   # Assembly 3: two D chains
    "": ('A', 'C')
}

DISTANCE_CUTOFFS = [4.0, 6.0, 8.0]


def calculate_refined_interface(pdb_file, chain1='B', chain2='C', distance_cutoffs=DISTANCE_CUTOFFS):
    """
    Calculate interface residues with multiple distance cutoffs
    """
    
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    
    try:
        chain_a = structure[0][chain1]
        chain_b = structure[1 if len(structure) > 1 else 0][chain2]
    except KeyError as e:
        print(f"Chain {e} not found in {pdb_file}")
        return None
    
    # Get all CA atoms
    ca_atoms_a = [(residue, residue['CA']) for residue in chain_a if 'CA' in residue]
    ca_atoms_b = [(residue, residue['CA']) for residue in chain_b if 'CA' in residue]
    
    results = {}
    # get all atoms in cuttofs range
    max_cutoff = max(distance_cutoffs)

    interface_a = dict()
    interface_b = dict()

    # Calculate distances
    for residue_a, ca_a in ca_atoms_a:
        min_distance = float('inf')
        for residue_b, ca_b in ca_atoms_b:
            for atom_a in residue_a:
                for atom_b in residue_b:
                    distance = atom_a-atom_b
                    if distance <= max_cutoff:
                        min_res_b_id = residue_b.get_id()[1] + RES_B_OFFSET
                        interface_b[min_res_b_id] = min(distance, interface_b.get(min_res_b_id, float('inf')))
                        if distance < min_distance:
                            min_distance = distance
    
        if min_distance <= max_cutoff:
            interface_a[residue_a.get_id()[1] + RES_A_OFFSET] = min_distance

    # sort by distance
    interface_a = sorted(interface_a.items(), key=lambda x: x[1])
    interface_b = sorted(interface_b.items(), key=lambda x: x[1])

    sorted_cutoffs = sorted(distance_cutoffs)

    a_idx = 0
    b_idx = 0
    for cutoff in sorted_cutoffs:
        results[f'{cutoff}A'] = {
            'chain_a': list(),
            'chain_b': list(),
            'chain_a_count': 0,
            'chain_b_count': 0
        }
        while a_idx < len(interface_a) and interface_a[a_idx][1] <= cutoff:
            results[f'{cutoff}A']['chain_a'].append(interface_a[a_idx][0])
            results[f'{cutoff}A']['chain_a_count'] += 1
            a_idx += 1

        while b_idx < len(interface_b) and interface_b[b_idx][1] <= cutoff:
            results[f'{cutoff}A']['chain_b'].append(interface_b[b_idx][0])
            results[f'{cutoff}A']['chain_b_count'] += 1
            b_idx += 1
            
    return results

def analyze_assembly_detailed(assembly_number):
    """Detailed analysis of Assembly {assembly_number} (Homodimer)"""

    report_lines = []
    report_lines.append(f"=== Detailed Analysis of Assembly {assembly_number} (Homodimer) ===\n")
    report_lines.append(f"using {SCAN_NET_FILE} as ScanNet file\n")
    # Analyze with different distance cutoffs
    interface_data = calculate_refined_interface(
        f"{PDB_FILE_BASE}{assembly_number}.pdb",
        chain1=ASSEMBLY_CHAINS[assembly_number][0],
        chain2=ASSEMBLY_CHAINS[assembly_number][1],
        distance_cutoffs=DISTANCE_CUTOFFS
    )

    if not interface_data:
        report_lines.append(f"Could not analyze Assembly {assembly_number}\n")
        # Print and save report
        report = "\n".join(report_lines)
        print(report)
        with open(f"{RESULTS_DIR}/assembly_{assembly_number}_detailed_report.txt", "w") as f:
            f.write(report)
        return

    report_lines.append("Interface residues with different distance cutoffs:\n")
    for cutoff, data in interface_data.items():
        report_lines.append(f"{cutoff} cutoff:")
        report_lines.append(f"  {ASSEMBLY_CHAINS[assembly_number][0]}: {data['chain_a_count']} residues")
        report_lines.append(f"  {ASSEMBLY_CHAINS[assembly_number][1]}: {data['chain_b_count']} residues")
        report_lines.append(f"  {ASSEMBLY_CHAINS[assembly_number][0]} residues: {data['chain_a']}")
        report_lines.append(f"  {ASSEMBLY_CHAINS[assembly_number][1]} residues: {data['chain_b']}")
        report_lines.append("")

    # Compare with ScanNet
    scan_net_df = pd.read_csv(SCAN_NET_FILE)
    high_prob_residues = set(scan_net_df[scan_net_df['Binding site probability'] >= 0.5]['Residue Index'])

    report_lines.append("=== Comparison with ScanNet ===")
    report_lines.append(f"ScanNet high-probability residues: {len(high_prob_residues)}")
    report_lines.append(f"ScanNet residues: {sorted(high_prob_residues)}\n")

    for cutoff, data in interface_data.items():
        interface_residues = set(data['chain_a'])
        overlap = high_prob_residues.intersection(interface_residues)

        report_lines.append(f"{cutoff} cutoff:")
        report_lines.append(f"  Interface residues: {len(interface_residues)}")
        report_lines.append(f"  Overlap: {len(overlap)}")
        if len(interface_residues) > 0:
            overlap_pct = len(overlap) / len(interface_residues) * 100
            report_lines.append(f"  Overlap percentage: {overlap_pct:.1f}%")
            report_lines.append(f"  Overlapping: {sorted(overlap)}")
            report_lines.append(f"  Interface-only: {sorted(interface_residues - high_prob_residues)}")
            report_lines.append(f"  ScanNet-only: {sorted(high_prob_residues - interface_residues)}")
        report_lines.append("")

    # Count the number of high-probability residues (According to ScanNet) that are in the interface
    # Over all cutoffs
    high_prob_residues_in_interface = 0
    for cutoff, data in interface_data.items():
        high_prob_residues_in_interface += len(high_prob_residues.intersection(data['chain_a']))
    report_lines.append(f"Number of high-probability residues in interface: {high_prob_residues_in_interface}")
    # Percantage of total high-probability residues that are in the interface over all cutoffs
    high_prob_residues_in_interface_pct = high_prob_residues_in_interface / len(high_prob_residues) * 100
    report_lines.append(f"Percentage of high-probability residues in interface: {high_prob_residues_in_interface_pct:.1f}%\n")

    # Print and save report
    report = "\n".join(report_lines)
    print(report)
    with open(f"{RESULTS_DIR}/assembly_{assembly_number}_detailed_report.txt", "w") as f:
        f.write(report)

    return interface_data


def create_validation_scripts(interface_data, assembly_number):
    # Create script for Assembly {assembly_number} interface
    script_content = f"""# ChimeraX script for Assembly {assembly_number} Interface Validation
# This shows the actual protein-protein interface (heterodimer)

# Open Assembly {assembly_number}
open {"../"*(CHIMERA_X_SCRIPT_DIR.count("/")+1)}{PDB_FILE_BASE}{assembly_number}.pdb

# Show surface
surface

# Color interface residues by distance cutoff: red (4Å), orange (6Å), yellow (8Å)
"""
    
    # Color interface residues by distance cutoff: red (4Å), orange (6Å), yellow (8Å)
    script_content += "\n# Chain A interface residues\n"
    for cutoff, color in zip([f"{cutoff}A" for cutoff in DISTANCE_CUTOFFS], ['red', 'orange', 'yellow']):
        residues = interface_data.get(cutoff, {}).get('chain_a', [])
        for residue in residues:
            script_content += f"color #{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}/{ASSEMBLY_CHAINS[assembly_number][0]}:{residue} {color}\n"

    script_content += "\n# Chain B interface residues\n"
    for cutoff, color in zip([f"{cutoff}A" for cutoff in DISTANCE_CUTOFFS], ['red', 'orange', 'yellow']):
        residues = interface_data.get(cutoff, {}).get('chain_b', [])
        for residue in residues:
            script_content += f"color #{ASSEMBLY_STRACTURES_IDS[assembly_number][1]}/{ASSEMBLY_CHAINS[assembly_number][1]}:{residue} {color}\n"

    script_content += f"""
# Set surface transparency
surface transparency 50

# Print information
echo "Assembly {assembly_number} Interface Validation - Homodimer"
echo "#{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}/{ASSEMBLY_CHAINS[assembly_number][0]} interface residues"
echo "#{ASSEMBLY_STRACTURES_IDS[assembly_number][1]}/{ASSEMBLY_CHAINS[assembly_number][1]} interface residues"
echo "Distance cutoff: {DISTANCE_CUTOFFS[-1]} Å"
echo ""
"""
    
    with open(f'{CHIMERA_X_SCRIPT_DIR}/chimeraX_assembly{assembly_number}_homodimer.cxc', 'w') as f:
        f.write(script_content)
    
    print(f"\nGenerated validation script: {CHIMERA_X_SCRIPT_DIR}/chimeraX_assembly{assembly_number}_homodimer.cxc")


def build_consensus_summary(interfaces) -> pd.DataFrame:
    # search for the most common interface residues
    surface_res_consensus = dict()
    scan_net_df = pd.read_csv(SCAN_NET_FILE)
    for assembly_number in ASSEMBLY_NUMBERS:
        for cutoff, data in interfaces[assembly_number].items():
            for residue in data['chain_a']:
                res_dict = surface_res_consensus.get(residue, dict([("Residue", residue)]))
                res_dict['Count'] = res_dict.get('Count', 0) + 1
                # Add ScanNet Probability
                if any(scan_net_df['Residue Index'] == residue):
                    res_dict['ScanNet Probability'] = scan_net_df[scan_net_df['Residue Index'] == residue]['Binding site probability'].iloc[0]
                else:
                    res_dict['ScanNet Probability'] = 0
                res_dict[f'Assembly {assembly_number}'] = cutoff
                surface_res_consensus[residue] = res_dict
    
    df =  pd.DataFrame(surface_res_consensus.values())
    for assembly_number in ASSEMBLY_NUMBERS:
        df[f"Assembly {assembly_number}"] = df[f"Assembly {assembly_number}"].fillna(np.inf)
    # Sort by Count, ScanNet Probability, Residue
    df = df.sort_values(by=['Count', 'ScanNet Probability', 'Residue'] + [f"Assembly {assembly_number}" for assembly_number in ASSEMBLY_NUMBERS], ascending=[False, False, True] + [True] * len(ASSEMBLY_NUMBERS))
    return df


def create_chimeraX_script_for_consensus_residues(assembly_number, df):
    # Create a ChimeraX script that colors the consensus residues
    script_content = f"""# ChimeraX script for Consensus Residues
    # Open Assembly {assembly_number}
    open {"../"*(CHIMERA_X_SCRIPT_DIR.count("/")+1)}{PDB_FILE_BASE}{assembly_number}.pdb

    # Show surface
    surface

    color gray

    hide #{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}
    hide #{ASSEMBLY_STRACTURES_IDS[assembly_number][1]}
    surface hide #{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}
    surface hide #{ASSEMBLY_STRACTURES_IDS[assembly_number][1]}

    surface show #{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}/{ASSEMBLY_CHAINS[assembly_number][0]}
    surface show #{ASSEMBLY_STRACTURES_IDS[assembly_number][1]}/{ASSEMBLY_CHAINS[assembly_number][1]}

    transparency #{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}/{ASSEMBLY_CHAINS[assembly_number][0]} 0 surface 
    transparency #{ASSEMBLY_STRACTURES_IDS[assembly_number][1]}/{ASSEMBLY_CHAINS[assembly_number][1]} 70 surface 
    
    """

    # color consensus residues
    for residue in df['Residue'][df['Count'] == 3]:
        script_content += f"color #{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}/{ASSEMBLY_CHAINS[assembly_number][0]}:{residue} red\n"

    # annotate high prob residues
    for residue in df['Residue'][(df['ScanNet Probability'] >= 0.5) & (df['Count'] == 3)]:
        label = f"ScanNet Probability: {df['ScanNet Probability'][df['Residue'] == residue].iloc[0]:.2f}"
        script_content += f"label #{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}/{ASSEMBLY_CHAINS[assembly_number][0]}:{residue} text '{label}'\n"
        # color high prob residues
        script_content += f"color #{ASSEMBLY_STRACTURES_IDS[assembly_number][0]}/{ASSEMBLY_CHAINS[assembly_number][0]}:{residue} yellow\n"

    with open(f'{CHIMERA_X_SCRIPT_DIR}/chimeraX_annotated_high_prob_residues_assembly{assembly_number}.cxc', 'w') as f:
        f.write(script_content)

def main():
    """Main analysis function"""

    print("=== Refined Interface Analysis ===")

    interfaces = dict()
    
    for assembly_number in ASSEMBLY_NUMBERS:
        interface_data = analyze_assembly_detailed(assembly_number)
        if interface_data:
            interfaces[assembly_number] = interface_data
            # Create validation scripts
            create_validation_scripts(interface_data, assembly_number)
    
    # Build consensus summary DataFrame
    surface_res_consensus = build_consensus_summary(interfaces)
    for assembly_number in ASSEMBLY_NUMBERS:
        create_chimeraX_script_for_consensus_residues(assembly_number, surface_res_consensus)

    # save surface_res_consensus to a csv file
    surface_res_consensus.to_csv(f'{RESULTS_DIR}/surface_res_consensus.csv', index=False)
    print(f"Saved surface_res_consensus to {RESULTS_DIR}/surface_res_consensus.csv")


if __name__ == "__main__":
    main() 