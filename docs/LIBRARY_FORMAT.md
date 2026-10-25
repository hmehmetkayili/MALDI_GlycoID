# N-Glycan Library Format & Schema — MALDI-GlycoID

## 1. Default Human N-Glycan Library
The default library is stored at:
- `data/libraries/human_nglycan_native_v1.json`
- `data/libraries/human_nglycan_native_v1.csv`

It contains **333 biosynthetically constrained human N-glycan compositions** generated from the conserved trimannosyl chitobiose core (`Hex3HexNAc2`) across `paucimannose`, `high_mannose`, `hybrid`, and `complex` classes.

## 2. Required Fields per Library Entry
| Field | Type | Description |
|---|---|---|
| `library_id` | `string` | Unique identifier (e.g., `NG_0001`) |
| `composition` | `string` | Canonical composition string (e.g., `Hex5HexNAc4dHex1Neu5Ac2`) |
| `residue_counts` | `object` | Mapping of monosaccharide residue counts |
| `glycan_class` | `string` | `paucimannose`, `high_mannose`, `hybrid`, or `complex` |
| `antennary_hypothesis` | `string` | Composition-inferred branching/isomer group description |
| `theoretical_neutral_mass` | `float` | Exact monoisotopic neutral mass (Da) verified against Stage 1 |
| `hill_formula` | `string` | Elemental formula in Hill notation |
| `adduct_mzs` | `object` | Precomputed `m/z` values for standard adducts (`[M+Na]+`, `[M+H]+`, `[M+K]+`, etc.) |
| `structural_annotation_status` | `string` | `COMPOSITION_CANDIDATE`, `TOPOLOGY_HYPOTHESIS`, `LINKAGE_SPECIFIC_DERIVATIZED`, or `STRUCTURALLY_VALIDATED_REFERENCE` |
| `source_provenance` | `string` | Provenance record |
| `license_info` | `string` | License metadata (`CC0-1.0 / Theoretical Calculation`) |

## 3. Regenerating & Validating the Library
```bash
python scripts/validate_library.py
```
