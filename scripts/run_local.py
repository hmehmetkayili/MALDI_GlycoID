"""CLI launcher script to run MALDI-GlycoID locally on Windows, macOS, or Linux (Stage 12)."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import threading
import webbrowser

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from app import __version__ as SOFTWARE_VERSION  # noqa: E402
from app.api.compat import run_http_server  # noqa: E402
from app.config import settings  # noqa: E402
from app.main import app  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Run MALDI-GlycoID local web application server.")
    parser.add_argument("--host", default=settings.host, help="Bind host address (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=settings.port, help="Bind TCP port (default: 8000)")
    parser.add_argument(
        "--open-browser",
        action="store_true",
        help="Automatically open the default web browser once the server starts.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Verify application initialization and default library loading, then exit 0.",
    )
    args = parser.parse_args()

    if args.check_only:
        from app.api.compat import TestClient

        client = TestClient(app)
        resp = client.get("/api/health")
        if resp.status_code != 200:
            print(f"[ERROR] Health check failed with status {resp.status_code}", file=sys.stderr)
            return 1
        data = resp.json()
        print(
            f"[OK] {data['app_name']} v{data['version']} initialized successfully "
            f"(Default library: {data['active_library_name']} with {data['active_library_entries_count']} entries)."
        )
        return 0

    browser_host = "127.0.0.1" if args.host in ("0.0.0.0", "::") else args.host
    url = f"http://{browser_host}:{args.port}"
    print(f"MALDI-GlycoID v{SOFTWARE_VERSION} listening on {url}")
    if args.open_browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    run_http_server(app, host=args.host, port=args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
