import os
import glob

"""This script processes FoldX mutation data from multiple text files,
extracts significant positions based on a threshold,
and finds the intersection of these positions across all files.
It outputs the intersecting positions to a text file.
the format of the output is ready for RFdiffusion input."""

# Parameters
THRESHOLD = 1.0  # kcal/mol
INPUT_DIR = "results"  # folder with .txt files
EXT = "*.txt"

def extract_position(mutation_str):
    """
    Parses FoldX mutation strings like PROA752A into PRO752 (ignores chain).
    """
    aa3 = "A"
    residue_number = ''.join(c for c in mutation_str[4:-1] if c.isdigit())
    return f"{aa3}{residue_number}"

def extract_significant_positions(file_path, threshold):
    positions = set()
    with open(file_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith("POS"):
                continue
            try:
                mutation_str, ddg = line.strip().split()
                ddg = float(ddg)
                if ddg >= threshold:
                    pos = extract_position(mutation_str)
                    positions.add(pos)
            except ValueError:
                continue
    return positions

def main():
    all_files = glob.glob(os.path.join(INPUT_DIR, EXT))
    all_significant = []

    for file in all_files:
        positions = extract_significant_positions(file, THRESHOLD)
        print(f"{os.path.basename(file)}: {len(positions)} significant positions")
        all_significant.append(positions)

    if len(all_significant) < 2:
        print("Need at least 2 files to find intersections.")
        return

    # Find intersection
    intersect = set.intersection(*all_significant)
    print("\nIntersecting significant positions across all files:")
    for pos in sorted(intersect):
        print(pos)

    # Save to a one-line CSV
    with open("intersecting_positions.txt", "w") as f:
        f.write(",".join(intersect))

if __name__ == "__main__":
    main()