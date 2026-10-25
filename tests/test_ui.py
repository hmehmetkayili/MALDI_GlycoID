"""Unit tests for Stage 11 User Interface and Workflow Integration (`app/templates/index.html`, `app/static/css/styles.css`, `app/static/js/app.js`)."""

from __future__ import annotations

from pathlib import Path
import unittest

from app.exports import REQUIRED_12_EXPORT_CATEGORIES


class TestUserInterfaceAndWorkflow(unittest.TestCase):
    """Verify the 12-page responsive English UI, 12 export cards, and workflow controller bindings."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.html_path = Path("app/templates/index.html")
        cls.css_path = Path("app/static/css/styles.css")
        cls.js_path = Path("app/static/js/app.js")
        cls.html_text = cls.html_path.read_text(encoding="utf-8")
        cls.css_text = cls.css_path.read_text(encoding="utf-8")
        cls.js_text = cls.js_path.read_text(encoding="utf-8")

    def test_all_12_navigation_views_present(self) -> None:
        """Verify all 12 required navigation buttons and view sections exist in index.html."""
        required_views = [
            "dashboard",
            "import",
            "library",
            "search",
            "glycopeptides",
            "batch",
            "quantification",
            "qc",
            "statistics",
            "exports",
            "settings",
            "help",
        ]
        for v in required_views:
            self.assertIn(f'data-view="{v}"', self.html_text)
            self.assertIn(f'id="view-{v}"', self.html_text)

    def test_all_12_export_category_cards_present_in_ui(self) -> None:
        """Verify all 12 standardized export categories have download cards in the UI."""
        for cat in REQUIRED_12_EXPORT_CATEGORIES:
            self.assertIn(f'data-category="{cat}"', self.html_text)
            self.assertIn(f'data-cat="{cat}"', self.html_text)

    def test_scientific_safeguards_and_workflow_bindings_in_ui(self) -> None:
        """Verify key scientific disclaimers, residue tables, and JS API endpoints exist."""
        self.assertIn("162.0528234315", self.html_text)
        self.assertIn("203.0793725330", self.html_text)
        self.assertIn("146.0579088090", self.html_text)
        self.assertIn("291.0954165270", self.html_text)
        self.assertIn("307.0903311490", self.html_text)

        # Verify JS endpoints
        required_endpoints = [
            "/api/spectra/upload",
            "/api/spectra/example",
            "/api/libraries/default",
            "/api/libraries/generate",
            "/api/libraries/query",
            "/api/search/ms1",
            "/api/search/quantify",
            "/api/glycopeptides/screen",
            "/api/batch/run",
            "/api/qc/calibrate",
            "/api/qc/evaluate",
            "/api/exports/download",
        ]
        for ep in required_endpoints:
            self.assertIn(ep, self.js_text)

        # Verify CSS responsive layout rules
        self.assertIn(".app-shell", self.css_text)
        self.assertIn(".dropzone", self.css_text)
        self.assertIn("@media (max-width: 900px)", self.css_text)


if __name__ == "__main__":
    unittest.main()
