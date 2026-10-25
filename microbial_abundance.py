import csv
from pathlib import Path
from collections import defaultdict

# -----------------------------------------
# 1. INPUT AND OUTPUT DIRECTORIES
# -----------------------------------------

base = Path("results/bracken")

out = Path("results/abundance/top5_rebuilt")
out.mkdir(parents=True, exist_ok=True)

# -----------------------------------------
# 2. SAMPLE IDENTIFICATION
# -----------------------------------------

sample_map = {
    "SRR7961333": "SKIN_1",
    "SRR7961317": "SKIN_2",
    "SRR7961304": "SKIN_3",
    "SRR7961281": "SKIN_4",
    "SRR33915754": "SKIN_5",

    "SRR7961222": "LUNG_1",
    "SRR7961255": "LUNG_2",
    "SRR7961291": "LUNG_3",
    "SRR7961223": "LUNG_4",
    "SRR7961221": "LUNG_5",

    "SRR7961292": "COLON_1",
    "SRR7961262": "COLON_2",
    "SRR7961278": "COLON_3",
    "SRR7961254": "COLON_4",

    "SRR2012208_1": "COLON_5A",
    "SRR2012208_2": "COLON_5B"
}

# -----------------------------------------
# 3. STORAGE
# -----------------------------------------

reports = {}

# Dictionary preserves insertion order.
# It prevents alphabetical rearrangement.
selected = {}

top5_rows = []

# -----------------------------------------
# 4. PROCESS ALL 16 BRACKEN REPORTS
# -----------------------------------------

for accession, sample in sample_map.items():

    file = base / f"{accession}.bracken"

    if not file.exists():
        raise FileNotFoundError(
            f"Missing Bracken report: {file}"
        )

    with file.open(newline="") as f:
        rows = list(
            csv.DictReader(f, delimiter="\t")
        )

    # Retain species-level assignments
    species = [
        r for r in rows
        if r["taxonomy_lvl"].strip() == "S"
    ]

    # Sort by fraction_total_reads
    # Homo sapiens remains included here.
    species.sort(
        key=lambda r: float(
            r["fraction_total_reads"]
        ),
        reverse=True
    )

    # Select top five
    top5 = species[:5]

    # Record top-five taxa
    for rank, r in enumerate(top5, 1):

        name = r["name"].strip()

        top5_rows.append({
            "Sample": sample,
            "Rank": rank,
            "Species": name,
            "new_est_reads": r["new_est_reads"],
            "fraction_total_reads":
                r["fraction_total_reads"]
        })

        # Exclude Homo sapiens only
        # AFTER selecting the top five.
        if name != "Homo sapiens":

            # Add species only once.
            # Preserve first appearance.
            if name not in selected:
                selected[name] = None

    # -------------------------------------
    # 5. RETRIEVE FULL BRACKEN ABUNDANCES
    # -------------------------------------

    full_counts = defaultdict(int)

    # IMPORTANT:
    # Use every species in the report,
    # not only the top five.
    for r in species:

        name = r["name"].strip()

        count = int(
            r["new_est_reads"]
        )

        full_counts[name] += count

    reports[sample] = full_counts

# Preserve original selection order
selected_species = list(selected.keys())

# -----------------------------------------
# 6. SAVE TOP-FIVE SELECTION
# -----------------------------------------

with (
    out / "top5_per_report.csv"
).open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Sample",
            "Rank",
            "Species",
            "new_est_reads",
            "fraction_total_reads"
        ]
    )

    writer.writeheader()
    writer.writerows(top5_rows)

# -----------------------------------------
# 7. COMBINE COLON_5A AND COLON_5B
# -----------------------------------------

reports["COLON_5"] = defaultdict(int)

all_colon5_species = (
    set(reports["COLON_5A"])
    | set(reports["COLON_5B"])
)

for name in all_colon5_species:

    forward = reports["COLON_5A"].get(
        name, 0
    )

    reverse = reports["COLON_5B"].get(
        name, 0
    )

    reports["COLON_5"][name] = (
        forward + reverse
    )

# -----------------------------------------
# 8. FINAL 15 BIOLOGICAL SAMPLES
# -----------------------------------------

samples = (
    [f"SKIN_{i}" for i in range(1, 6)]
    + [f"LUNG_{i}" for i in range(1, 6)]
    + [f"COLON_{i}" for i in range(1, 6)]
)

# -----------------------------------------
# 9. CREATE FINAL ABUNDANCE TABLE
# -----------------------------------------

with (
    out / "abundance_rebuilt.csv"
).open("w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow(
        ["Taxon"] + samples
    )

    for name in selected_species:

        abundances = [
            reports[sample].get(name, 0)
            for sample in samples
        ]

        writer.writerow(
            [name] + abundances
        )

# -----------------------------------------
# 10. VERIFY COLON_5 ADDITION
# -----------------------------------------

with (
    out / "COLON_5_verification.csv"
).open("w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "Species",
        "COLON_5A",
        "COLON_5B",
        "COLON_5_combined"
    ])

    for name in selected_species:

        forward = reports["COLON_5A"].get(
            name, 0
        )

        reverse = reports["COLON_5B"].get(
            name, 0
        )

        combined = forward + reverse

        writer.writerow([
            name,
            forward,
            reverse,
            combined
        ])

# -----------------------------------------
# 11. FINAL CHECKS
# -----------------------------------------

assert len(sample_map) == 16
assert len(samples) == 15

for name in selected_species:

    assert reports["COLON_5"][name] == (
        reports["COLON_5A"].get(name, 0)
        + reports["COLON_5B"].get(name, 0)
    )

print()
print("===================================")
print("TOP-FIVE ABUNDANCE REBUILD COMPLETE")
print("===================================")
print()

print("Bracken reports processed:", len(sample_map))
print("Unique non-human species:", len(selected_species))
print("Final biological samples:", len(samples))

print()
print("Species order: FIRST APPEARANCE")
print("Alphabetical sorting: DISABLED")

print()
print("COLON_5A + COLON_5B: VERIFIED")

print()
print("Files saved to:", out)

print()
print("1. top5_per_report.csv")
print("2. abundance_rebuilt.csv")
print("3. COLON_5_verification.csv")
