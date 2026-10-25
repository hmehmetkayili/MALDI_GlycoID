# Installation Guide — MALDI-GlycoID

## 1. System Requirements
- **Operating System:** Windows 10/11, Linux (Ubuntu 20.04+), or macOS 12+
- **Python Version:** Python 3.11 or newer
- **Memory & Disk:** Minimum 2 GB RAM, 200 MB free disk space

## 2. Standard Python Virtual Environment Installation

```bash
# Clone or extract repository
cd MALDI-GLYCOID2

# Create virtual environment
python -m venv .venv

# Activate virtual environment (Windows PowerShell)
.\.venv\Scripts\Activate.ps1

# Activate virtual environment (Linux / macOS)
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

## 3. Offline / Scientific Python Distribution Execution
MALDI-GlycoID is engineered so that all core scientific modules, multi-sheet `.xlsx` ECMA-376 OpenXML export, vector `.pdf` report generation, and the local HTTP web server (`scripts/run_local.py`) run out-of-the-box in standard scientific Python 3.11 environments (`numpy`, `scipy`, `pandas`, `scikit-learn`, `matplotlib`, `lxml`, `PyPDF2`) even without internet access.

## 4. Verifying Your Installation

```bash
python scripts/validate_library.py
python scripts/run_local.py --check-only
python -m unittest discover -s tests -p "test_*.py" -v
```
