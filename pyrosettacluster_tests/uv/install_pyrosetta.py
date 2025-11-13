__author__ = "Jason C. Klima"


import pyrosetta_installer


try:
    pyrosetta_installer.install_pyrosetta(
        distributed=False,
        serialization=True,
        skip_if_installed=False,
        mirror=0
    )
except Exception as ex:
    print(f"PyRosetta installation with 'mirror=0' failed: {ex}. Retrying with 'mirror=1'.")
    pyrosetta_installer.install_pyrosetta(
        distributed=False,
        serialization=True,
        skip_if_installed=False,
        mirror=1
    )
