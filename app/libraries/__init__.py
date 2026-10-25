"""Public interface for Stage 3: Human N-glycan library generation, search, validation, and I/O."""

from __future__ import annotations

from app.libraries.composition_generator import (
    CompositionConstraints,
    GlycanClass,
    GlycanLibrary,
    GlycanLibraryEntry,
    StructuralAnnotationStatus,
    classify_n_glycan_composition,
    generate_human_n_glycan_library,
)
from app.libraries.library_io import (
    export_library_to_delimited,
    export_library_to_json,
    import_library_from_delimited,
    import_library_from_json,
)
from app.libraries.library_search import (
    LibraryQueryMatch,
    LibrarySearchFilter,
    get_entry_adduct_mz,
    query_glycan_library,
)
from app.libraries.library_validation import (
    DuplicateGroup,
    LibraryValidationResult,
    detect_library_duplicates,
    validate_glycan_library,
)

__all__ = [
    "CompositionConstraints",
    "DuplicateGroup",
    "GlycanClass",
    "GlycanLibrary",
    "GlycanLibraryEntry",
    "LibraryQueryMatch",
    "LibrarySearchFilter",
    "LibraryValidationResult",
    "StructuralAnnotationStatus",
    "classify_n_glycan_composition",
    "detect_library_duplicates",
    "export_library_to_delimited",
    "export_library_to_json",
    "generate_human_n_glycan_library",
    "get_entry_adduct_mz",
    "import_library_from_delimited",
    "import_library_from_json",
    "query_glycan_library",
    "validate_glycan_library",
]
