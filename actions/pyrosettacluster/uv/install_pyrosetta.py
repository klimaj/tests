__author__ = "Jason C. Klima"


import pyrosetta_installer

kwargs = dict(
    distributed=False,
    serialization=True,
    mirror=0,
    type="Release",
    use_setup_py_package=False,
    skip_if_installed=False,
    extras="",
    silent=False,
)

try:
    pyrosetta_installer.install_pyrosetta(**kwargs)
except Exception as ex:
    print(f"PyRosetta installation with 'mirror=0' failed: {ex}. Retrying with 'mirror=1'.")
    kwargs.update(mirror=1)
    pyrosetta_installer.install_pyrosetta(**kwargs)
