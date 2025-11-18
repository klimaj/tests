__author__ = "Jason C. Klima"


import argparse
import pyrosetta
import os

from datetime import datetime
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Uncomment after PR #567 is merged
# from pyrosetta.distributed.cluster import get_instance_kwargs
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
from typing import Optional

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Temporary bootstrap - delete after PR #567 is merged
import json
import warnings

from functools import singledispatch
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.distributed.cluster import get_scores_dict
from pyrosetta.distributed.cluster.io import secure_read_pickle
from pyrosetta.distributed.cluster.exceptions import InputFileError

from typing import Any, Dict, NoReturn, Optional, Tuple, Union


@singledispatch
def parse_input_file_to_instance_kwargs(obj: Any) -> NoReturn:
    raise InputFileError(obj)

@singledispatch
def parse_input_file_to_instance_metadata_kwargs(obj: Any) -> NoReturn:
    raise InputFileError(obj)

@parse_input_file_to_instance_metadata_kwargs.register(PackedPose)
@parse_input_file_to_instance_metadata_kwargs.register(Pose)
@parse_input_file_to_instance_metadata_kwargs.register(str)
def _parse_str(obj: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    scores_dict = get_scores_dict(obj)
    return scores_dict["instance"], scores_dict["metadata"]

@singledispatch
def parse_scorefile(obj: Any) -> NoReturn:
    raise TypeError(
        "The `scorefile` argument parameter must be of type `str`, "
        + "not of type {0}.".format(type(obj))
    )

@parse_scorefile.register(str)
def _parse_str(obj: str) -> Union[str, NoReturn]:
    if not os.path.exists(obj):
        raise ValueError(
            "The `scorefile` argument parameter must exist! Received {0}".format(obj)
        )
    return obj

@singledispatch
def parse_decoy_name(obj: Any) -> NoReturn:
    raise TypeError(
        "The `decoy_name` argument parameter must be of type `str`, "
        + "not of type {0}.".format(type(obj))
    )

@parse_decoy_name.register(str)
def _from_str(obj: str) -> str:
    return obj

def get_instance_kwargs(
    input_file: Optional[Union[str, Pose, PackedPose]] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
    skip_corrections: Optional[bool] = None,
    with_metadata_kwargs: Optional[bool] = None,
) -> Union[Dict[str, Any], Tuple[Dict[str, Any], Dict[str, Any]], NoReturn]:
    """
    Given an input file that was written by PyRosettaCluster, or a scorefile
    and a decoy name that was written by PyRosettaCluster, return the PyRosettaCluster
    instance kwargs needed to reproduce the decoy using PyRosettaCluster.

    Args:
        input_file: A `str` object specifying the path to the '.pdb', '.pdb.bz2', '.pkl_pose',
            '.pkl_pose.bz2', '.b64_pose', '.b64_pose.bz2', '.init', or '.init.bz2' file, or a
            `Pose` or `PackedPose` object, from which to extract PyRosettaCluster instance kwargs.
            If 'input_file' is provided, then ignore the 'scorefile' and 'decoy_name' keyword
            argument parameters.
            Default: None
        scorefile: A `str` object specifying the path to the JSON-formatted scorefile
            (or pickled `pandas.DataFrame` scorefile) from a PyRosettaCluster simulation
            from which to extract PyRosettaCluster instance kwargs. If 'scorefile'
            is provided, 'decoy_name' must also be provided. In order to use a scorefile,
            it must contain full simulation records from the original production
            run; i.e., the attribute 'simulation_records_in_scorefile' was set to True.
            Default: None
        decoy_name: A `str` object specifying the decoy name for which to extract
            PyRosettaCluster instance kwargs. If 'decoy_name' is provided, 'scorefile'
            must also be provided.
            Default: None
        skip_corrections: A `bool` object specifying whether or not to skip any ScoreFunction
            corrections specified in the PyRosettaCluster task initialization options
            (extracted from the 'input_file' or 'scorefile' keyword argument parameter).
            Default: None
        with_metadata_kwargs: A `bool` object specifying whether or not to return a `tuple`
            object with the instance kwargs as the first element and the metadata kwargs as
            the second element.
            Default: None

    Returns:
        A `dict` object of PyRosettaCluster instance kwargs, or a `tuple` object of `dict`
        objects with the PyRosettaCluster instance kwargs as the first element and the
        PyRosettaCluster metadata kwargs as the second element when `with_metadata_kwargs=True`.
    """
    _simulation_records_in_scorefile_msg = (
        "The 'scorefile' argument parameter does not contain the full simulation records. "
        + "In order to reproduce a decoy using a 'scorefile', the PyRosettaCluster "
        + "attribute 'simulation_records_in_scorefile' must have been set to `True` in "
        + "the original simulation. Please provide an 'input_file' generated by PyRosettaCluster, "
        + "or a 'scorefile' with full simulation records generated by PyRosettaCluster, "
        + "in order to reproduce."
    )
    if input_file:
        if scorefile or decoy_name:
            warnings.warn(
                "Received 'input_file' and either 'scorefile' or 'decoy_name' keyword argument parameters. "
                + "Ignoring 'scorefile' and 'decoy_name' and using 'input_file' keyword argument parameter!",
                UserWarning,
                stacklevel=2,
            )
        if with_metadata_kwargs:
            instance_kwargs, metadata_kwargs = parse_input_file_to_instance_metadata_kwargs(input_file)
        else:
            instance_kwargs = parse_input_file_to_instance_kwargs(input_file)
    elif scorefile and decoy_name:
        scorefile = parse_scorefile(scorefile)
        decoy_name = parse_decoy_name(decoy_name)
        instance_kwargs = None
        if scorefile.endswith(".json"):
            with open(scorefile, "r") as f:
                lines = f.readlines()
                for line in lines:
                    try:
                        scorefile_entry = json.loads(line)
                    except:
                        raise TypeError(
                            "`get_instance_kwargs()` received `scorefile` which does not appear to be JSON-formatted."
                        )
                    if all(k in scorefile_entry for k in ("metadata", "instance")):
                        if "decoy_name" in scorefile_entry["metadata"]:
                            if scorefile_entry["metadata"]["decoy_name"] == decoy_name:
                                instance_kwargs = scorefile_entry["instance"]
                                if with_metadata_kwargs:
                                    metadata_kwargs = scorefile_entry["metadata"]
                                break
                    else:
                        raise NotImplementedError(_simulation_records_in_scorefile_msg)
        else:
            try:
                df = secure_read_pickle(scorefile, compression="infer")
            except:
                raise TypeError(
                    "`get_instance_kwargs()` received `scorefile` which does not appear to be "
                    + "readable by `pyrosetta.distributed.cluster.io.secure_read_pickle(compression='infer')`."
                )
            if all(k in df.columns for k in ("metadata", "instance")):
                for instance, metadata in df[["instance", "metadata"]].values:
                    if "decoy_name" in metadata:
                        if metadata["decoy_name"] == decoy_name:
                            instance_kwargs = dict(instance)
                            if with_metadata_kwargs:
                                metadata_kwargs = dict(metadata)
                            break
            else:
                raise NotImplementedError(_simulation_records_in_scorefile_msg)
        if instance_kwargs is None:
            raise KeyError(
                "Error in `get_instance_kwargs()`! The provided `decoy_name` is not in the provided `scorefile`."
            )
    else:
        raise NotImplementedError(
            "`get_instance_kwargs()` requires either `input_file` (or `scorefile` and `decoy_name`) argument parameter inputs."
        )
    assert isinstance(
        instance_kwargs, dict
    ), "Returned instance keyword arguments are not of type `dict`."
    if with_metadata_kwargs:
        assert isinstance(
            metadata_kwargs, dict
        ), "Returned metadata keyword arguments are not of type `dict`."

    if skip_corrections:
        assert isinstance(
            instance_kwargs["tasks"], dict
        ), "PyRosettaCluster 'tasks' keyword argument parameter must be an instance of `dict`."
        for option in ("extra_options", "options"):
            if option in instance_kwargs["tasks"]:
                if isinstance(instance_kwargs["tasks"][option], dict):
                    instance_kwargs["tasks"][option] = toolz.dicttoolz.keyfilter(
                        lambda k: not k.startswith(("corrections:", "-corrections:")),
                        instance_kwargs["tasks"][option],
                    )
                elif isinstance(instance_kwargs["tasks"][option], str):
                    if "corrections:" in instance_kwargs["tasks"][option]:
                        raise NotImplementedError(
                            "Cannot skip ScoreFunction corrections because the original simulation did not output "
                            + "PyRosettaCluster results with normalized PyRosetta initialization options or configure "
                            + "the task's PyRosetta initialization options as an instance of `dict`. Please disable the "
                            + "'skip_corrections' keyword argument to continue with the reproduction."
                        )
                else:
                    raise TypeError(
                        f"PyRosettaCluster task key '{option}' must have a value of type `dict` or `str`. "
                        + f"Received: {type(instance_kwargs['tasks'][option])}"
                    )

    if with_metadata_kwargs:
        return instance_kwargs, metadata_kwargs
    else:
        return instance_kwargs
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def main(
    input_file: Optional[str],
    scorefile: Optional[str],
    decoy_name: Optional[str],
    env_dir: Optional[str],
    pandas_secure: bool,
) -> None:
    """
    Dump the PyRosettaCluster environment file(s) based on metadata from a PyRosettaCluster result.
    """
    if scorefile and (not scorefile.endswith(".json")):
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

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Temporary bootstrap - delete after PR #567 is merged
    import sys
    from pathlib import Path

    # sys.path.append(str(Path(__file__).resolve().parent.parent))
    from actions.pyrosettacluster.utils import (
        ROSETTACOMMONS_CONDA_CHANNEL,
        detect_platform,
    )

    if env_manager == "pixi":
        if not toml:
            # Detect Python version
            # py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
            # py_feature = f"py{sys.version_info.major}{sys.version_info.minor}"
            # # Detect platform
            # plat = detect_platform()
            # toml_format = "pixi.toml"
            # template_toml_file = Path(__file__).resolve().parent.parent / "actions" / "pyrosettacluster" / "pixi" / toml_format
            # with open(template_toml_file, "r") as f:
            #     toml = f.read().format(
            #         rosettacommons_conda_channel=ROSETTACOMMONS_CONDA_CHANNEL,
            #         name=os.path.basename(env_dir),
            #         plat=plat,
            #         py_version=py_version,
            #         py_feature=py_feature,
            #     )
            original_env_dir = os.path.join(os.path.dirname(env_dir), os.path.basename(env_dir).split("_reproduce")[0])
            toml_format = "pixi.toml"
            original_toml_file = os.path.join(original_env_dir, toml_format)
            with open(original_toml_file, "r") as f:
                toml = f.read()

    elif env_manager == "uv":
        if not toml:
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            toml_format = "pyproject.toml"
            template_toml_file = Path(__file__).resolve().parent.parent / "actions" / "pyrosettacluster" / "uv" / toml_format
            with open(template_toml_file, "r") as f:
                toml = f.read().format(
                    name=os.path.basename(env_dir),
                    py_version=py_version,
                )
            # original_env_dir = os.path.join(os.path.dirname(env_dir), os.path.basename(env_dir).split("_reproduce")[0])
            # toml_format = "pyproject.toml"
            # original_toml_file = os.path.join(original_env_dir, toml_format)
            # with open(original_toml_file, "r") as f:
            #     toml = f.read()
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
        pandas_secure=args.pandas_secure,
    )
