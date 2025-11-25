"""
This script extracts and writes a PyRosettaCluster environment file
based on metadata from a PyRosettaCluster result. It requires PyRosetta
to be installed for the PyRosettaCluster environment file extraction.
"""

__author__ = "Jason C. Klima"


import argparse
import pyrosetta
import pyrosetta.distributed
import os

from datetime import datetime
from pyrosetta.distributed.cluster import get_instance_kwargs
from typing import Optional


def main(
    input_file: Optional[str],
    scorefile: Optional[str],
    decoy_name: Optional[str],
    env_dir: Optional[str],
    pyrosetta_init_flags: Optional[str],
    pandas_secure: bool,
) -> None:
    """
    Dump the PyRosettaCluster environment file(s) based on metadata from a PyRosettaCluster result.
    """
    if (
        isinstance(input_file, str)
        and input_file.endswith((".pkl_pose", ".pkl_pose.bz2", ".b64_pose", ".b64_pose.bz2"))
        and pyrosetta_init_flags
    ):
        pyrosetta.distributed.init(pyrosetta_init_flags)

    if (
        isinstance(scorefile, str)
        and not scorefile.endswith(".json")
    ):
        if pandas_secure:
            pyrosetta.secure_unpickle.add_secure_package("pandas")
        else:
            raise RuntimeError(
                "Please also pass the `--pandas_secure` flag to retrieve data from a `pandas.DataFrame`. "
                "Note that the input pickled `pandas.DataFrame` scorefile will be deserialized using the "
                "`pickle` module, so please only run this with data sources you trust."
            )

    instance_kwargs, metadata_kwargs = get_instance_kwargs(
        input_file=input_file,
        scorefile=scorefile,
        decoy_name=decoy_name,
        skip_corrections=False,
        with_metadata_kwargs=True,
    )

    environment = instance_kwargs.get("environment")
    if not environment:
        raise RuntimeError(
            "The PyRosettaCluster result contains an empty environment file string; "
            "the environment cannot be recreated."
        )

    toml = metadata_kwargs.get("toml")
    toml_format = metadata_kwargs.get("toml_format")
    env_manager = metadata_kwargs.get("environment_manager")  # may be `None` in legacy cases

    sha1 = instance_kwargs.get("sha1")
    print(f"[INFO] GitHub SHA1: {sha1}")

    # Determine output files based on manager
    if env_manager == "pixi":
        env_file = os.path.join(env_dir, "pixi.lock")
        write_file(env_file, environment)
        if toml:
            toml_file = os.path.join(env_dir, toml_format)
            write_file(toml_file, toml)
        else:
            print(
                f"[WARNING] The PyRosettaCluster result contains an empty pixi manifest file string. "
                "Please generate a pixi manifest file for the written `pixi.lock` file, or ensure "
                "the original pixi manifest file is put in the output directory before recreating "
                "the pixi project."
            )

    elif env_manager == "uv":
        env_file = os.path.join(env_dir, "requirements.txt")
        write_file(env_file, environment)
        if toml:
            toml_file = os.path.join(env_dir, toml_format)
            write_file(toml_file, toml)
        else:
            print(
                f"[WARNING] The PyRosettaCluster result contains an empty uv TOML file string. "
                "Please generate a uv TOML file for the written `requirements.txt` file, or ensure "
                "the original uv TOML file is put in the output directory before recreating "
                "the uv project."
            )

    elif env_manager in ("conda", "mamba"):
        env_file = os.path.join(env_dir, "environment.yml")
        write_file(env_file, environment)

    else:
        # Legacy fallback if 'environment_manager' key is missing assumes conda
        env_file = os.path.join(env_dir, "environment.yml")
        print(
            f"[WARNING] Environment manager not detected from PyRosettaCluster result.",
            f"Writing '{os.path.basename(env_file)}' file.",
        )
        write_file(env_file, environment)


def write_file(path: str, content: str) -> None:
    """Utility function for writing a file."""
    with open(path, "w") as f:
        f.write(content)
    print(f"[INFO] Wrote: '{path}'")


def ensure_env_dir(path: Optional[str]) -> str:
    """
    Normalize the directory path, create it if it does not exist, and ensure it is writable.
    """
    if path is None:
        # Create a timestamped directory name
        dirname = "PyRosettaCluster_" + datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f")
        path = os.path.abspath(dirname)
    else:
        if not isinstance(path, str):
            raise argparse.ArgumentTypeError(
                f"The 'env_dir' parameter must be of type `str`. Received: {type(path)}"
            )
        path = os.path.abspath(os.path.expanduser(path))

    # Create directory if missing
    if not os.path.isdir(path):
        try:
            os.makedirs(path, exist_ok=True)
        except Exception as ex:
            raise argparse.ArgumentTypeError(f"Could not create directory '{path}': {ex}")

    # Ensure writable
    if not os.access(path, os.W_OK):
        raise argparse.ArgumentTypeError(f"The directory '{path}' is not writable.")

    print(f"[INFO] Environment directory: '{path}'")

    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Dump a PyRosettaCluster environment file using metadata contained in a "
            "PyRosettaCluster output decoy file or PyRosettaCluster output scorefile."
        )
    )

    # -------------------------------------------------------------------------
    # Inputs for `pyrosetta.distributed.cluster.get_instance_kwargs`
    # -------------------------------------------------------------------------
    parser.add_argument(
        "--input_file",
        type=str,
        required=False,
        default=None,
        help=(
            "Path to a PyRosettaCluster output decoy file (a '.pdb', '.pdb.bz2', '.pkl_pose', "
            "'.pkl_pose.bz2', '.b64_pose', or '.b64_pose.bz2' file) or an output PyRosetta "
            "initialization file ('.init' or '.init.bz2' file). If provided, the `--scorefile` "
            "and `--decoy_name` flags are ignored."
        ),
    )

    parser.add_argument(
        "--scorefile",
        type=str,
        required=False,
        default=None,
        help=(
            "Path to a PyRosettaCluser output scorefile ('.json' file or a pickled `pandas.DataFrame` "
            "scorefile containing full simulation records. Required if the `--decoy_name` flag is used "
            "and the `--input_file` flag is not provided."
        ),
    )

    parser.add_argument(
        "--decoy_name",
        type=str,
        required=False,
        default=None,
        help=(
            "The PyRosettaCluster output decoy name in the provided scorefile for which to extract "
            "instance kwargs. Requires the `--scorefile` flag, and is ignored if the `--input_file` "
            "flag is provided."
        ),
    )

    # -------------------------------------------------------------------------
    # Environment directory
    # -------------------------------------------------------------------------
    parser.add_argument(
        "--env_dir",
        type=ensure_env_dir,
        required=False,
        default=None,
        help=(
            "Directory where environment files will be written. If not provided, "
            "a timestamped directory named 'PyRosettaCluster_<timestamp>' will be "
            "created in the current working directory."
        ),
    )

    # -------------------------------------------------------------------------
    # PyRosetta initialization options
    # -------------------------------------------------------------------------
    parser.add_argument(
        "--pyrosetta_init_flags",
        type=str,
        required=False,
        default=None,
        help=(
            "Optional PyRosetta initialization flags to use when loading a pickled or base64-encoded "
            "PyRosettaCluster output decoy file (i.e., '.pkl_pose', '.pkl_pose.bz2', '.b64_pose', "
            "or '.b64_pose.bz2' file type extensions) with the `--input_file` flag, specifying any "
            "extra Rosetta topology `.params` files or patch files required to load the `Pose` object "
            "into memory. Ignored for all other input file types or scorefile inputs."
        ),
    )

    # -------------------------------------------------------------------------
    # Security option for loading pickled `pandas.DataFrames` objects
    # -------------------------------------------------------------------------
    parser.add_argument(
        "--pandas_secure",
        action="store_true",
        default=False,
        help=(
            "Allow loading a pickled `pandas.DataFrame` scorefile. "
            "*Warning*: This will deserialize the scorefile using Python's `pickle` module, "
            "which can execute arbitrary code. Please only enable this if the scorefile "
            "comes from a trusted source."
        ),
    )

    args = parser.parse_args()

    main(
        input_file=args.input_file,
        scorefile=args.scorefile,
        decoy_name=args.decoy_name,
        env_dir=args.env_dir,
        pyrosetta_init_flags=args.pyrosetta_init_flags,
        pandas_secure=args.pandas_secure,
    )
