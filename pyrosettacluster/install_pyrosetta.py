"""
Perform PyRosetta installation with the PyPI `pyrosetta-installer` package.
"""

__author__ = "Jason C. Klima"


import argparse
import pyrosetta_installer


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
            f"Mirrors tried: {mirrors}. Last exception: {prev_exception}"
        ) from prev_exception


def parse_args():
    parser = argparse.ArgumentParser(
        description="Install PyRosetta using the PyPI 'pyrosetta-installer' package with automatic mirror fallback."
    )
    parser.add_argument(
        "--mirror_order",
        nargs="+",
        type=int,
        default=[0, 1],
        help=(
            "Optionally specify the PyRosetta installer mirror order to try, e.g. `--mirror_order 0 1`. "
            "See the PyPI 'pyrosetta-installer' package website for details."
        ),
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    mirror_order = tuple(args.mirror_order)
    install_pyrosetta_with_mirrors(
        mirrors=mirror_order,
        base_kwargs=BASE_KWARGS,
    )
