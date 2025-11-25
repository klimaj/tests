__author__ = "Jason C. Klima"


import pyrosetta_installer


# Order of mirrors to prioritize
MIRROR_ORDER = (0, 1)

# Base `pyrosetta_installer.install_pyrosetta` kwargs
BASE_KWARGS = dict(
    distributed=False,
    serialization=True,
    type="Release",
    use_setup_py_package=False,
    skip_if_installed=False,
    extras="",
    silent=False,
)


def install_pyrosetta_with_mirrors(mirrors, base_kwargs):
    """
    Attempt to install PyRosetta using the `pyrosetta-installer` package with a
    prioritized list of mirrors. Returns cleanly on first success, or raises if
    all attempts fail.
    """
    prev_exception = None
    for mirror in mirrors:
        try:
            print(f"Attempting PyRosetta installation using `mirror={mirror}`...")
            kwargs = {**base_kwargs, "mirror": mirror}
            pyrosetta_installer.install_pyrosetta(**kwargs)
            print(f"PyRosetta installation succeeded using `mirror={mirror}`!")
            break
        except Exception as ex:
            prev_exception = ex
            print(f"{type(ex).__name__}: PyRosetta installation failed with `mirror={mirror}`: {ex}")
            continue
    else:
        raise RuntimeError(
            f"All PyRosetta installation attempts failed. "
            f"Mirrors tried: {mirrors}. Most recent exception: {prev_exception}"
        ) from prev_exception


if __name__ == "__main__":
    install_pyrosetta_with_mirrors(
        mirrors=MIRROR_ORDER,
        base_kwargs=BASE_KWARGS,
    )
