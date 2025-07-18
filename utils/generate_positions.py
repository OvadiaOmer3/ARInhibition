import xml.etree.ElementTree as ET

def generate_foldx_positions(xml_file):
    """
    Parses a pdb-pisa interface residue XML file to generate a position string
    for an alanine scan.

    Args:
        xml_file (str): The path to the input XML file.

    Returns:
        str: A comma-separated string of positions for the FoldX command.
    """
    # Standard 3-letter to 1-letter amino acid code mapping
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    positions = []
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except FileNotFoundError:
        return f"Error: The file '{xml_file}' was not found."
    except ET.ParseError:
        return f"Error: The file '{xml_file}' is not a valid XML file."


    # Loop through both RESIDUE1 and RESIDUE2 sections
    for residue_group in root:
        for residue in residue_group.findall('RESIDUE'):
            # Find the score and check if it's greater than 0
            score_element = residue.find('BURIEDSURFACEAREASCORE')
            if score_element is not None and int(score_element.text) > 0:
                
                # Extract structure string, e.g., " C:PRO 671    "
                structure_str = residue.find('STRUCTURE').text.strip()
                
                # Split into parts, e.g., "C:PRO" and "671"
                parts = structure_str.split()
                chain_res_part = parts[0]
                res_num = parts[1]

                # Split chain/residue part, e.g., "C" and "PRO"
                chain_id, res_name = chain_res_part.split(':')
                
                # Convert 3-letter code to 1-letter code
                one_letter_code = three_to_one.get(res_name)

                if one_letter_code:
                    # Format the string for FoldX and add to our list
                    # e.g., P + A + 671 + a -> PA671a
                    foldx_position = f"{one_letter_code}{chain_id}{res_num}a"
                    positions.append(foldx_position)

    # Join all the individual position strings with a comma
    return ",".join(positions)

# --- Main execution ---
if __name__ == "__main__":
    # Replace 'residue0.xml' with the actual path to your file if different
    filename = 'residue0.xml'
    position_string = generate_foldx_positions(filename)
    print(position_string)