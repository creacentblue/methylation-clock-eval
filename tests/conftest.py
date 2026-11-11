"""Make ``scripts/`` importable; tests run without any downloaded data."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
