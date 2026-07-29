# =========================================================
# LNC-SCOPE FULL PIPELINE
# PART 1: Custom lncRNA Folding + Accessibility
# PART 2: miRNA Binding + Functional Prediction
# =========================================================

import json
import time
import csv
import os
import math
from collections import defaultdict
from turner_folding import fold_turner, can_pair, mutation_sensitivity_fast
# =========================================================
# PART 1 : CUSTOM RNA FOLDING
# =========================================================

# -------------------------------
# TURNER FOLDING WRAPPER
# -------------------------------

def fold_rna(seq):
    structure, mfe = fold_turner(seq)
    fold_rna.last_mfe = mfe
    return structure


def fold_rna_with_energy(seq):
    return fold_turner(seq)

# -------------------------------
# ACCESSIBLE REGIONS (6-mer seed)
# -------------------------------

def find_accessible_regions(structure, seed_length=6, min_open=4):
    regions = []

    for i in range(len(structure) - seed_length + 1):
        window = structure[i:i+seed_length]

        if window.count('.') >= min_open:
            regions.append({
                "start": i,
                "end": i + seed_length,
                "length": seed_length,
                "open_count": window.count('.')
            })

    return regions


def merge_accessible_regions(regions):
    if not regions:
        return []

    regions = sorted(regions, key=lambda r: r["start"])
    merged = []
    current = {"start": regions[0]["start"], "end": regions[0]["end"]}

    for r in regions[1:]:
        if r["start"] <= current["end"]:
            current["end"] = max(current["end"], r["end"])
        else:
            merged.append(current)
            current = {"start": r["start"], "end": r["end"]}
    merged.append(current)

    for m in merged:
        m["length"] = m["end"] - m["start"]

    return merged


# -------------------------------
# FUNCTIONAL ZONES
# -------------------------------

def classify_functional_zones(regions, seed_length=7):
    zones = []

    for region in regions:
        l = region["length"]

        if l >= seed_length * 3:
            zone_type = "Sponge Zone"
        elif l < seed_length:
            zone_type = "Scaffold Candidate"
        else:
            zone_type = "Hybrid Zone"

        zones.append({
            "start": region["start"],
            "end": region["end"],
            "length": l,
            "type": zone_type
        })

    return zones


# -------------------------------
# MUTATION SENSITIVITY
# -------------------------------

def mutate_base(base):
    return {'A': 'U', 'U': 'G', 'G': 'C'}.get(base, 'A')


# =========================================================
# PART 2 : miRNA MODULE
# =========================================================

# -------------------------------
# Tunable thresholds
# -------------------------------

SEED_START = 1
SEED_END = 8

MIN_PAIRS = 6              # minimum WC/wobble pairs in the (7nt) seed window
DELTAG_THRESHOLD = -3.0    # kcal/mol cutoff to call a candidate "high confidence"
SPONGE_DISTANCE = 30       # nt: sites within this distance are considered clustered
SPONGE_MIN_CLUSTER = 2     # minimum clustered sites to call a region a sponge cluster

CONF_DELTAG_SCALE = 300.0      # kcal/mol of cumulative |deltaG| at which the deltaG component reaches ~63% of its max
CONF_DIVERSITY_SCALE = 20.0    # distinct miRNA species at which the diversity component reaches ~63% of its max
CONF_DELTAG_WEIGHT = 0.6
CONF_DIVERSITY_WEIGHT = 0.4


def load_mirna_csv(file_path):
    mirna_dict = {}
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['name']
            seq = row['sequence'].upper().replace("T", "U")
            if len(seq) >= 7:
                mirna_dict[name] = seq
    return mirna_dict


def is_pair(a, b):
    return (
        (a == 'A' and b == 'U') or
        (a == 'U' and b == 'A') or
        (a == 'G' and b == 'C') or
        (a == 'C' and b == 'G') or
        (a == 'G' and b == 'U') or
        (a == 'U' and b == 'G')
    )


def reverse_complement(seq):
    comp = {'A':'U','U':'A','G':'C','C':'G'}
    return ''.join(comp[b] for b in seq[::-1])


def build_kmer_set(sequence, k=7):
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}


def prefilter_mirna(sequence, miRNA_data):
    seed_len = SEED_END - SEED_START
    kmer_set = build_kmer_set(sequence, seed_len)
    filtered = {}
    for name, seq in miRNA_data.items():
        seed = seq[SEED_START:SEED_END]
        seed_rc = reverse_complement(seed)
        if seed_rc in kmer_set:
            filtered[name] = seq
    return filtered


def get_pair_energy(a, b):
    pair = a + b
    energy_table = {
        "AU": -1.1, "UA": -1.1,
        "GC": -2.3, "CG": -2.3,
        "GU": -0.9, "UG": -0.9
    }
    return energy_table.get(pair, 0)


def evaluate_duplex(window, seed):
    rev_window = window[::-1]
    pair_count = 0
    deltaG = 0.0
    for a, b in zip(rev_window, seed):
        if is_pair(a, b):
            pair_count += 1
        deltaG += get_pair_energy(a, b)
    return pair_count, deltaG


def find_and_score_matches(sequence, regions, miRNA_data):
    scored = []
    n_mirnas = len(miRNA_data)
    n_regions = len(regions)

    for idx, (name, mir_seq) in enumerate(miRNA_data.items(), 1):
        seed = mir_seq[SEED_START:SEED_END]
        seed_len = len(seed)

        for region in regions:
            sub_seq = sequence[region["start"]:region["end"]]

            if len(sub_seq) < seed_len:
                continue

            for i in range(len(sub_seq) - seed_len + 1):
                window = sub_seq[i:i + seed_len]
                pair_count, deltaG = evaluate_duplex(window, seed)

                if pair_count >= MIN_PAIRS and deltaG < DELTAG_THRESHOLD:
                    scored.append({
                        "miRNA": name,
                        "pos": region["start"] + i,
                        "pair_count": pair_count,
                        "deltaG": deltaG,
                        "site": window
                    })

    return scored


def merge_overlapping_hits(scored, seed_len=SEED_END - SEED_START):    
    by_mirna = defaultdict(list)
    for s in scored:
        by_mirna[s["miRNA"]].append(s)

    merged = []
    for mirna, hits in by_mirna.items():
        hits.sort(key=lambda h: h["pos"])
        cluster = [hits[0]]
        for h in hits[1:]:
            if h["pos"] - cluster[-1]["pos"] <= seed_len:
                cluster.append(h)
            else:
                merged.append(min(cluster, key=lambda c: c["deltaG"]))
                cluster = [h]
        merged.append(min(cluster, key=lambda c: c["deltaG"]))

    return merged


def classify_sponge_distance(interactions, max_gap=SPONGE_DISTANCE, min_cluster=SPONGE_MIN_CLUSTER):
    
    if not interactions:
        return "non-sponge", []

    positions = sorted(i["pos"] for i in interactions)

    clusters = []
    current_cluster = [positions[0]]

    for pos in positions[1:]:
        if pos - current_cluster[-1] <= max_gap:
            current_cluster.append(pos)
        else:
            clusters.append(current_cluster)
            current_cluster = [pos]
    clusters.append(current_cluster)

    sponge_clusters = [c for c in clusters if len(c) >= min_cluster]

    label = "miRNA sponge" if sponge_clusters else "non-sponge"
    return label, sponge_clusters


def compute_confidence(scored):
    if not scored:
        return 0.0

    total_abs_deltaG = sum(abs(s["deltaG"]) for s in scored)
    distinct_mirnas = len(set(s["miRNA"] for s in scored))

    deltaG_component = 1 - math.exp(-total_abs_deltaG / CONF_DELTAG_SCALE)
    diversity_component = 1 - math.exp(-distinct_mirnas / CONF_DIVERSITY_SCALE)

    confidence = (CONF_DELTAG_WEIGHT * deltaG_component +
                  CONF_DIVERSITY_WEIGHT * diversity_component)

    return confidence


def reconcile_classification(scored, sponge_label, sponge_clusters, functional_zones):
    
    def zone_type_at(pos):
        for z in functional_zones:
            if z["start"] <= pos < z["end"]:
                return z["type"]
        return None

    structural_hit = any(zone_type_at(s["pos"]) == "Sponge Zone" for s in scored)
    functional_hit = sponge_label == "miRNA sponge"

    if structural_hit and functional_hit:
        verdict = "Validated miRNA Sponge"
    elif functional_hit and not structural_hit:
        verdict = "Putative miRNA Sponge (functional evidence only, no structural Sponge Zone overlap)"
    elif structural_hit and not functional_hit:
        verdict = "Structural Sponge Candidate (unconfirmed : no clustered binding evidence)"
    elif scored:
        zone_types = [zone_type_at(s["pos"]) for s in scored if zone_type_at(s["pos"])]
        if zone_types:
            dominant = max(set(zone_types), key=zone_types.count)
            verdict = dominant
        else:
            verdict = "Hybrid Zone (unclassified overlap)"
    else:
        verdict = "non-sponge"

    return {
        "verdict": verdict,
        "structural_sponge_zone_hit": structural_hit,
        "functional_cluster_hit": functional_hit,
        "sponge_clusters": sponge_clusters
    }


# =========================================================
# MAIN PROGRAM
# =========================================================

pipeline_start = time.time()

print("=== LNC-SCOPE PIPELINE STARTED ===")

# -------------------------------
# INPUT
# -------------------------------

choice = input("Do you have dot-bracket structure? (yes/no): ").strip().lower()

sequence = ""
structure = ""

if choice == "yes":
      structure = input("Paste dot-bracket structure: ")
      mfe = float(input("Enter MFE from RNAfold (kcal/mol): "))
      sequence = input("Paste the corresponding RNA sequence:\n").strip().upper()
      sequence = sequence.replace("T", "U")
      if len(sequence) != len(structure):
          print(f"WARNING: sequence length ({len(sequence)}) != structure length ({len(structure)}) — indices will misalign.")
else:
    file_choice = input("Do you have FASTA file input? (yes/no): ").strip().lower()

    if file_choice == "yes":
        path = input("Enter FASTA file path:\n").strip().strip('"')

        with open(path, "r") as f:
            lines = f.readlines()

        sequence = ''.join([l.strip() for l in lines if not l.startswith(">")])
        sequence = sequence.upper().replace("T", "U")

        print("\nSequence preview:")
        print(sequence[:60] + "...")

    else:
        sequence = input("Paste RNA sequence:\n").strip().upper()
        sequence = sequence.replace("T", "U")

# -------------------------------
# FOLD (only if structure wasn't already provided)
# -------------------------------

if choice != "yes":
    start = time.time()
    structure, mfe = fold_rna_with_energy(sequence)
    single_fold = time.time() - start

    print("\nSingle fold time:", round(single_fold, 2), "sec")

    print("\nPredicted structure (Turner 2004 model)")
    print("Predicted MFE:", round(mfe, 2), "kcal/mol")

# -------------------------------
# ANALYSIS
# -------------------------------

print(f"\nSequence length: {len(sequence)} nt")

accessible_regions = find_accessible_regions(structure, seed_length=7, min_open=5)
print(f"Accessible regions found: {len(accessible_regions)}")
merged_regions = merge_accessible_regions(accessible_regions)
functional_zones = classify_functional_zones(merged_regions, seed_length=7)

sensitive_sites = mutation_sensitivity_fast(sequence, window=80, min_delta=0.5)

# Save structural stage
with open("structure_output.json", "w") as f:
    json.dump({
        "sequence": sequence,
        "structure": structure,
        "mfe_kcal_mol": mfe,
        "accessible_regions": accessible_regions,
        "merged_regions": merged_regions,
        "functional_zones": functional_zones,
        "mutation_sensitive_positions": sensitive_sites
    }, f, indent=4)

print("\nStructure stage complete.")

# =========================================================
# miRNA STAGE
# =========================================================

miRNA_data = load_mirna_csv("mirna_sequences.csv")
print("Total miRNAs:", len(miRNA_data))

miRNA_data = prefilter_mirna(sequence, miRNA_data)
print("After filtering:", len(miRNA_data))

raw_scored = find_and_score_matches(sequence, accessible_regions, miRNA_data)
scored = merge_overlapping_hits(raw_scored)
sponge_label, sponge_clusters = classify_sponge_distance(scored)
confidence = compute_confidence(scored)
if scored:
    _total_abs_deltaG = sum(abs(s["deltaG"]) for s in scored)
    _distinct_mirnas = len(set(s["miRNA"] for s in scored))
    print(f" confidence inputs: total|deltaG|={_total_abs_deltaG:.1f} kcal/mol, "
          f"distinct miRNAs={_distinct_mirnas}, raw confidence={confidence!r}")
combined = reconcile_classification(scored, sponge_label, sponge_clusters, functional_zones)
hotspots = [s["pos"] for s in scored]

final_output = {
    "lncRNA": "custom_input",
    "mfe_kcal_mol": mfe,
    "function": combined["verdict"],
    "structural_sponge_zone_hit": combined["structural_sponge_zone_hit"],
    "functional_cluster_hit": combined["functional_cluster_hit"],
    "miRNA": list(set(i["miRNA"] for i in scored)),
    "binding_sites": hotspots,
    "sponge_clusters": sponge_clusters,
    "confidence": confidence,
    "accessible_regions": accessible_regions,
    "merged_regions": merged_regions,
    "functional_zones": functional_zones
}

# -------------------------------
# FINAL SAVE
# -------------------------------

with open("final_lnc_scope_output.json", "w") as f:
    json.dump(final_output, f, indent=4)

# -------------------------------
# RESULTS
# -------------------------------

print("\n=== PIPELINE COMPLETE ===")
print("Function prediction:", final_output["function"])
print("Confidence:", round(final_output["confidence"], 4))
print("Results saved to final_lnc_scope_output.json")

print("\nRaw window hits before merging:", len(raw_scored))
print("Distinct binding events after merging:", len(scored))

print("\nTotal pipeline time:",
      round(time.time() - pipeline_start, 2),
      "seconds")
