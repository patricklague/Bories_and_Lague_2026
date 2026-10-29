"""Cylinder version: reads distribution_data/total_cyl, writes total_cyl/ (same pooled-count method)."""
import importlib.util
import os

_spec = importlib.util.spec_from_file_location(
    "pmf_counting", os.path.join(os.path.dirname(os.path.abspath(__file__)), "pmf-from-distribution.py"))
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)

if __name__ == "__main__":
    _module.run("total_cyl")
    print("Done.")
