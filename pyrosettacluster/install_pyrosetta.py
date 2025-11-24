__author__ = "Jason C. Klima"


import pyrosetta_installer


# Order of mirrors to prioritize
MIRROR_ORDER = (1,) # (1, 0)


kwargs = dict(
    distributed=False,
    serialization=True,
    type="Release",
    use_setup_py_package=False,
    skip_if_installed=False,
    extras="",
    silent=False,
)

try:
    kwargs.update(mirror=MIRROR_ORDER[0])
    pyrosetta_installer.install_pyrosetta(**kwargs)
except Exception as ex:
    _err_msg = f"PyRosetta installation with 'mirror={MIRROR_ORDER[0]}' failed: {ex}"
    if len(MIRROR_ORDER) == 2:
        print(f"{_err_msg}. Retrying with 'mirror={MIRROR_ORDER[1]}'.")
        kwargs.update(mirror=MIRROR_ORDER[1])
        pyrosetta_installer.install_pyrosetta(**kwargs)
    else:
        raise NotImplementedError(f"{_err_msg} Skipping installation attempts with mirrors.")
