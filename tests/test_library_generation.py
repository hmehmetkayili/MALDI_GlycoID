"""Unit tests for Stage 3: Automatic human N-glycan library generation, search, validation, and I/O."""

from __future__ import annotations

import unittest

from app.core.mass_error import MassToleranceUnit
from app.libraries import (
    CompositionConstraints,
    GlycanClass,
    GlycanLibrary,
    GlycanLibraryEntry,
    LibrarySearchFilter,
    StructuralAnnotationStatus,
    detect_library_duplicates,
    export_library_to_delimited,
    export_library_to_json,
    generate_human_n_glycan_library,
    import_library_from_delimited,
    import_library_from_json,
    query_glycan_library,
    validate_glycan_library,
)


class TestHumanNGlycanLibrary(unittest.TestCase):
    """Comprehensive verification of Stage 3 library generation, constraints, search, and I/O."""

    def setUp(self) -> None:
        self.default_constraints = CompositionConstraints(
            min_hex=3,
            max_hex=9,
            min_hexnac=2,
            max_hexnac=6,
            min_dhex=0,
            max_dhex=2,
            min_neu5ac=0,
            max_neu5ac=3,
            include_neu5gc=False,
            labeling="NONE",
        )
        self.lib = generate_human_n_glycan_library(
            constraints=self.default_constraints,
            library_name="Test Human N-Glycan Library",
            library_version="1.0.0",
        )

    def test_biosynthetic_composition_generation_and_classes(self) -> None:
        """Verify generated library contains valid high-mannose, hybrid, and complex compositions."""
        self.assertGreater(self.lib.entry_count, 50)

        # Ensure all entries are marked as COMPOSITION_CANDIDATE (not claiming topological resolution)
        for entry in self.lib.entries:
            self.assertEqual(entry.structural_annotation_status, StructuralAnnotationStatus.COMPOSITION_CANDIDATE)
            self.assertFalse(entry.contains_neu5gc)

        classes_present = {e.glycan_class for e in self.lib.entries}
        self.assertIn(GlycanClass.HIGH_MANNOSE, classes_present)
        self.assertIn(GlycanClass.HYBRID, classes_present)
        self.assertIn(GlycanClass.COMPLEX, classes_present)

        # Verify nonsensical compositions (e.g. Hex5HexNAc2Neu5Ac2: oligomannose with sialic acids) are excluded
        compositions_set = {e.composition for e in self.lib.entries}
        self.assertIn("Hex5HexNAc2", compositions_set)  # Man5
        self.assertIn("Hex5HexNAc4dHex1", compositions_set)  # G2F
        self.assertIn("Hex5HexNAc4Neu5Ac2", compositions_set)  # A2G2S2
        self.assertNotIn("Hex5HexNAc2Neu5Ac2", compositions_set)  # Biosynthetically invalid

    def test_optional_neu5gc_and_linkage_specific_expansion(self) -> None:
        """Verify optional Neu5Gc generation and alpha2,3/alpha2,6 sialic linkage expansion."""
        c_gc = CompositionConstraints(
            min_hex=5,
            max_hex=5,
            min_hexnac=4,
            max_hexnac=4,
            min_dhex=0,
            max_dhex=0,
            min_neu5ac=0,
            max_neu5ac=1,
            include_neu5gc=True,
            max_neu5gc=1,
            labeling="EEC_SIALIC",
            sialic_linkage_expansion=True,
        )
        lib_eec = generate_human_n_glycan_library(constraints=c_gc)
        comps = {e.composition for e in lib_eec.entries}
        # Hex5HexNAc4Neu5Ac1 should expand into [a2,3:0/a2,6:1] and [a2,3:1/a2,6:0]
        self.assertIn("Hex5HexNAc4Neu5Ac1[a2,3:0/a2,6:1]", comps)
        self.assertIn("Hex5HexNAc4Neu5Ac1[a2,3:1/a2,6:0]", comps)
        # Check Neu5Gc warning is recorded
        neu5gc_entries = [e for e in lib_eec.entries if e.contains_neu5gc]
        self.assertGreater(len(neu5gc_entries), 0)
        self.assertTrue(any("Neu5Gc" in w for w in neu5gc_entries[0].warnings))

    def test_library_validation_and_duplicate_isomer_detection(self) -> None:
        """Test library mass consistency validation and distinction between isomers and duplicates."""
        val_res = validate_glycan_library(self.lib)
        self.assertTrue(val_res.is_valid)
        self.assertEqual(len(val_res.errors), 0)
        self.assertEqual(len(val_res.duplicate_groups), 0)

        # Create two structural isomers sharing Hex5HexNAc4dHex1 (e.g., core fucosylated G2F vs antennary Lewis-x)
        base_entry = next(e for e in self.lib.entries if e.composition == "Hex5HexNAc4dHex1")
        isomer_a = GlycanLibraryEntry(
            **{
                **base_entry.__dict__,
                "library_id": "ISO_001",
                "reference_id": "G12345AA",
                "structural_annotation_status": StructuralAnnotationStatus.STRUCTURALLY_VALIDATED_REFERENCE,
            }
        )
        isomer_b = GlycanLibraryEntry(
            **{
                **base_entry.__dict__,
                "library_id": "ISO_002",
                "reference_id": "G67890BB",
                "structural_annotation_status": StructuralAnnotationStatus.STRUCTURALLY_VALIDATED_REFERENCE,
            }
        )
        dups = detect_library_duplicates([isomer_a, isomer_b])
        self.assertEqual(len(dups), 1)
        self.assertTrue(dups[0].are_structural_isomers)
        self.assertIn("Isomeric composition group", dups[0].description)

    def test_offline_library_search_modes(self) -> None:
        """Verify composition search, mass-range search, exact m/z search, and class filters."""
        # 1. Exact composition query (supporting alias 'Fuc')
        q1 = query_glycan_library(
            self.lib,
            LibrarySearchFilter(composition_query="Hex5HexNAc4Fuc1"),
        )
        self.assertEqual(len(q1), 1)
        self.assertEqual(q1[0].entry.composition, "Hex5HexNAc4dHex1")

        # 2. Exact m/z search for G2F [M+Na]+ (1809.6393 Da) within 10 ppm
        q2 = query_glycan_library(
            self.lib,
            LibrarySearchFilter(
                target_mass=1809.6393,
                target_mass_is_mz=True,
                adduct_code="[M+Na]+",
                tolerance=10.0,
                tolerance_unit=MassToleranceUnit.PPM,
            ),
        )
        self.assertEqual(len(q2), 1)
        self.assertEqual(q2[0].entry.composition, "Hex5HexNAc4dHex1")
        assert q2[0].abs_error_ppm is not None
        self.assertLess(q2[0].abs_error_ppm, 1.0)

        # 3. Class & fucosylation filter
        q3 = query_glycan_library(
            self.lib,
            LibrarySearchFilter(
                glycan_classes=[GlycanClass.HIGH_MANNOSE],
                fucosylated_only=False,
            ),
        )
        self.assertTrue(all(m.entry.glycan_class == GlycanClass.HIGH_MANNOSE and not m.entry.is_fucosylated for m in q3))

    def test_json_csv_tsv_roundtrip_consistency(self) -> None:
        """Verify exporting and re-importing a library in JSON, CSV, and TSV preserves all masses and metadata."""
        json_str = export_library_to_json(self.lib)
        re_json = import_library_from_json(json_str)
        self.assertEqual(re_json.entry_count, self.lib.entry_count)
        self.assertAlmostEqual(
            re_json.entries[0].theoretical_neutral_mass,
            self.lib.entries[0].theoretical_neutral_mass,
            places=7,
        )

        csv_str = export_library_to_delimited(self.lib, delimiter=",")
        re_csv = import_library_from_delimited(csv_str, delimiter=",")
        self.assertEqual(re_csv.entry_count, self.lib.entry_count)
        self.assertTrue(validate_glycan_library(re_csv).is_valid)

        tsv_str = export_library_to_delimited(self.lib, delimiter="\t")
        re_tsv = import_library_from_delimited(tsv_str, delimiter="\t")
        self.assertEqual(re_tsv.entry_count, self.lib.entry_count)
        self.assertAlmostEqual(
            re_tsv.entries[-1].adduct_mzs["[M+Na]+"],
            self.lib.entries[-1].adduct_mzs["[M+Na]+"],
            places=6,
        )


if __name__ == "__main__":
    unittest.main()
