"""Miscellaneous utilities for PyRosettaCluster unit tests."""

__author__ = "Jason C. Klima"


import platform


ROSETTACOMMONS_CONDA_CHANNEL = "https://conda.rosettacommons.org" # "https://conda.graylab.jhu.edu"


def detect_platform():
    """Detect system platform string used by GitHub Actions."""
    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == "linux":
        plat = "linux-64" if "64" in machine else "linux-32"
    elif system == "darwin":
        plat = "osx-arm64" if "arm" in machine else "osx-64"
    elif system == "windows":
        plat = "win-64"
    else:
        raise RuntimeError(f"Unsupported platform: {system} ({machine})")

    return plat
