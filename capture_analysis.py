
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv

# Read peptide sequences from the text file
with open("peptides.txt", "r") as file:
    sequences = [
        line.strip().upper()
        for line in file
        if line.strip() and not line.startswith(">")
    ]

valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
results = []

for i, seq in enumerate(sequences, start=1):

    if not set(seq).issubset(valid_amino_acids):
        print(f"Skipping invalid sequence {i}: {seq}")
        continue

    peptide = ProteinAnalysis(seq)

    result = {
        "Peptide_ID": f"P{i:02d}",
        "Sequence": seq,
        "Length": len(seq),
        "Molecular_Weight_Da": round(
            peptide.molecular_weight(), 2
        ),
        "Isoelectric_Point": round(
            peptide.isoelectric_point(), 2
        ),
        "Net_Charge_pH7": round(
            peptide.charge_at_pH(7), 2
        ),
        "GRAVY": round(peptide.gravy(), 3),
        "Hydrophobic_Percent": round(
            sum(aa in "AVILMFW" for aa in seq)
            / len(seq) * 100, 2
        ),
        "Lysine_Count": seq.count("K"),
        "Arginine_Count": seq.count("R"),
        "Aspartate_Count": seq.count("D"),
        "Glutamate_Count": seq.count("E"),
        "Cysteine_Count": seq.count("C")
    }

    results.append(result)

if results:
    with open("capture_analysis.csv", "w", newline="") as file:
        writer = csv.DictWriter(
            file, fieldnames=results[0].keys()
        )
        writer.writeheader()
        writer.writerows(results)

    print("\nPEPTIDE ANALYSIS COMPLETE!")
    print(f"Successfully analyzed {len(results)} peptides.")
    print("Results saved to capture_analysis.csv")

    for result in results:
        print(
            result["Peptide_ID"],
            result["Sequence"],
            "Charge:", result["Net_Charge_pH7"],
            "GRAVY:", result["GRAVY"]
        )
else:
    print("No valid peptide sequences found.")
