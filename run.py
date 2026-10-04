"""Run this workflow from any working directory; no arguments shows help."""
import runpy
import sys
from pathlib import Path

if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    sys.path.insert(0, str(directory))
    target = directory / "amorphous_md.py"
    sys.argv = [str(target), *(sys.argv[1:] or ["--help"])]
    runpy.run_path(str(target), run_name="__main__")
