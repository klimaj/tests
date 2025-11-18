"""
*Warning*: This function runs a subprocess with one of the following commands:
    - `conda env create ...`: when 'conda' is an executable
    - `mamba env create ...`: when 'mamba' is an executable
    - `uv pip install ...`: when 'uv' is an executable
    - `pixi install ...`: when 'pixi' is an executable
Installing certain packages may not be secure, so please only run with input files you trust.
Learn more about PyPI security `here <https://pypi.org/security>`_ and conda security `here <https://www.anaconda.com/docs/reference/security>`_.

Given an input file that was written by PyRosettaCluster, or a scorefile
and a decoy name that was written by PyRosettaCluster, recreate the
environment that was used to generate the decoy with a new environment name.

The environment manager used (i.e., either 'conda', 'mamba', 'uv', or 'pixi') is
automatically determined from the operating system environment variable
'PYROSETTACLUSTER_ENVIRONMENT_MANAGER' if exported, or otherwise it must be
provided using the `--env_manager` flag.
"""

__author__ = "Jason C. Klima"


# import argparse
# import os
# import shutil

# from datetime import datetime


# # -----------------------------------------
# import argparse
# import os
# import subprocess

# from pathlib import Path

# # from pyrosetta.distributed.cluster import recreate_environment

# import json
# import logging
# import os
# import shutil
# import subprocess
# import sys
# import tempfile
# import warnings

# from datetime import datetime
# from functools import lru_cache
# # from pyrosetta.distributed.packed_pose.core import PackedPose
# # from pyrosetta.rosetta.core.pose import Pose
# from typing import (
#     Dict,
#     Generic,
#     NoReturn,
#     Optional,
#     Tuple,
#     TypeVar,
#     Union,
# )

# # from pyrosetta.distributed.cluster.tools import get_instance_kwargs

# sys.path.append(str(Path(__file__).resolve().parent.parent.parent))
# from pyrosettacluster_tests.utils import (
#     ROSETTACOMMONS_CONDA_CHANNEL,
#     detect_platform,
# )

# G = TypeVar("G")

# def get_instance_kwargs(scorefile, decoy_name):
#     instance_kwargs = None
#     if scorefile.endswith(".json"):
#         with open(scorefile, "r") as f:
#             lines = f.readlines()
#             for line in lines:
#                 try:
#                     scorefile_entry = json.loads(line)
#                 except:
#                     raise TypeError(
#                         "`get_instance_kwargs()` received `scorefile` which does not appear to be JSON-formatted."
#                     )
#                 if all(k in scorefile_entry for k in ("metadata", "instance")):
#                     if "decoy_name" in scorefile_entry["metadata"]:
#                         if scorefile_entry["metadata"]["decoy_name"] == decoy_name:
#                             instance_kwargs = scorefile_entry["instance"]
#                             break
#                 else:
#                     raise NotImplementedError("Could not locate decoy.")
#     assert instance_kwargs is not None, instance_kwargs

#     return instance_kwargs


# class EnvironmentConfig(Generic[G]):
#     _ENV_VAR: str = "PYROSETTACLUSTER_ENVIRONMENT_MANAGER"
#     _ENV_MANAGERS: Tuple[str, ...] = ("pixi", "uv", "mamba", "conda")

#     def env_create_cmd(
#         self, environment_name: str, raw_spec: str, tmp_dir: str, base_dir: str
#     ) -> Union[str, NoReturn]:
#         # Create a project directory for uv/pixi, or prefix directory for conda/mamba
#         project_dir = os.path.join(base_dir, environment_name)
#         # Raise exception if the project directory exists
#         if os.path.isdir(project_dir):
#             if self.environment_manager in ("conda", "mamba"):
#                 _err_msg = f"The {self.environment_manager} environment prefix directory already exists: '{project_dir}'"
#             elif self.environment_manager in ("uv", "pixi"):
#                 _err_msg = f"The {self.environment_manager} project directory already exists: '{project_dir}'"
#             else:
#                 raise RuntimeError(f"Unsupported environment manager: '{self.environment_manager}'")
#             raise IsADirectoryError(_err_msg)
#         os.makedirs(project_dir, exist_ok=False)

#         if self.environment_manager in ("conda", "mamba"):
#             yml_file = os.path.join(tmp_dir, f"{environment_name}.yml")
#             with open(yml_file, "w") as f:
#                 f.write(raw_spec)

#             if self.environment_manager == "conda":
#                 return f"conda env create -f '{yml_file}' -p '{project_dir}'"

#             elif self.environment_manager == "mamba":
#                 return f"mamba env create -f '{yml_file}' -p '{project_dir}'"

#         elif self.environment_manager == "uv":
#             # Write the requirements.txt file
#             req_file = os.path.join(tmp_dir, f"{environment_name}.txt")
#             with open(req_file, "w") as f:
#                 f.write(raw_spec)
#             # Determine the current Python version
#             py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
#             return (
#                 f"uv init '{project_dir}' --python {py_version} && "
#                 f"uv add --requirements '{req_file}' --project '{project_dir}'"
#             )

#         elif self.environment_manager == "pixi":
#             py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
#             py_feature = f"py{sys.version_info.major}{sys.version_info.minor}"

#             # Detect platform
#             plat = detect_platform()
#             template_toml_file = os.path.join(os.path.dirname(__file__), os.pardir, "pixi", "pixi.toml")
#             with open(template_toml_file, "r") as f:
#                 pixi_toml = f.read().format(
#                     rosettacommons_conda_channel=ROSETTACOMMONS_CONDA_CHANNEL,
#                     name=environment_name,
#                     plat=plat,
#                     py_version=py_version,
#                     py_feature=py_feature,
#                 )
#             toml_file = os.path.join(project_dir, "pixi.toml")
#             with open(toml_file, "w") as f:
#                 f.write(pixi_toml)
#             lock_file = os.path.join(project_dir, "pixi.lock")
#             with open(lock_file, "w") as f:
#                 f.write(raw_spec)

#             return "pixi install --frozen"

#         raise RuntimeError(f"Unsupported environment manager: '{self.environment_manager}'")


# @lru_cache(maxsize=1)
# def get_environment_config() -> EnvironmentConfig:
#     """Return an instance of the `EnvironmentConfig` class on the host process."""
#     return EnvironmentConfig()


# def get_environment_manager() -> str:
#     """Get the configured environment manager."""
#     return get_environment_config().environment_manager


# def get_environment_cmd() -> str:
#     """Get the configured environment export command."""
#     return get_environment_config().env_export_cmd


# def get_environment_var() -> str:
#     """Get the PyRosettaCluster operating system environment variable name."""
#     return EnvironmentConfig._ENV_VAR


# def recreate_environment(
#     environment_name: Optional[str] = None,
#     input_file: Optional[str] = None,
#     scorefile: Optional[str] = None,
#     decoy_name: Optional[str] = None,
#     timeout: Optional[int] = None,
#     base_dir: Optional[str] = None,
# ) -> Optional[NoReturn]:
#     """
#     *Warning*: This function runs a subprocess with one of the following commands:
#       - `conda env create ...`: when 'conda' is an executable
#       - `mamba env create ...`: when 'mamba' is an executable
#       - `uv pip ...`: when 'uv' is an executable
#       - `pixi install ...`: when 'pixi' is an executable
#     Installing certain packages may not be secure, so please only run with input files you trust.
#     Learn more about PyPI security `here <https://pypi.org/security>`_ and conda security `here <https://www.anaconda.com/docs/reference/security>`_.

#     Given an input file that was written by PyRosettaCluster, or a scorefile
#     and a decoy name that was written by PyRosettaCluster, recreate the
#     environment that was used to generate the decoy with a new environment name.

#     The environment manager used (i.e., either 'conda', 'mamba', 'uv', or 'pixi') is
#     automatically determined from the operating system environment variable
#     'PYROSETTACLUSTER_ENVIRONMENT_MANAGER' if exported, or otherwise it is automatically
#     selected based on which manager is available on the system. If no executables
#     are explicitly found, then 'conda' is used by default.

#     Args:
#         environment_name: A `str` object specifying the new name of the environment
#             to recreate. If using 'conda' and 'mamba', this is the prefix directory that will
#             be created in the 'base_dir' directory. If using 'uv' or 'pixi', this is the
#             local project directory name that will be created in the 'base_dir' directory.
#             Default: 'PyRosettaCluster_' + datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f")
#         input_file: A `str` object specifying the path to the '.pdb', '.pdb.bz2', '.pkl_pose',
#             '.pkl_pose.bz2', '.b64_pose', or '.b64_pose.bz2' file, or a `Pose` or `PackedPose`
#             object, from which to extract PyRosettaCluster instance kwargs. If 'input_file' is
#             provided, then ignore the 'scorefile' and 'decoy_name' keyword argument parameters.
#             Default: None
#         scorefile: A `str` object specifying the path to the JSON-formatted scorefile
#             (or pickled `pandas.DataFrame` scorefile) from a PyRosettaCluster simulation
#             from which to extract PyRosettaCluster instance kwargs. If 'scorefile'
#             is provided, 'decoy_name' must also be provided. The scorefile must
#             contain full simulation records from the original production run
#             (i.e., 'simulation_records_in_scorefile' was set to `True`).
#             Default: None
#         decoy_name: A `str` object specifying the decoy name for which to extract
#             PyRosettaCluster instance kwargs. Must be provided if 'scorefile'
#             is used.
#             Default: None
#         timeout: An `int` object specifying the timeout in seconds before any
#             subprocesses are terminated.
#             Default: None
#         base_dir: A `str` object specifying the base directory in which to create
#             the environment.
#             Default: `.`

#     Returns:
#         None
#     """
#     def _run_subprocess(cmd, cwd=None, env=None) -> str:
#         try:
#             return subprocess.check_output(
#                 cmd,
#                 shell=True,
#                 stderr=subprocess.STDOUT,
#                 timeout=timeout,
#                 text=True,
#                 cwd=cwd,
#                 env=env,
#                 executable="/bin/bash", # Ensure `&&` works properly
#             )
#         except subprocess.CalledProcessError as ex:
#             print(
#                 f"Command failed: `{cmd}`\n",
#                 f"Return code: {ex.returncode}\n",
#                 f"Output:\n{ex.output}",
#             )
#             raise RuntimeError(cmd)

#     if not environment_name:
#         environment_name = "PyRosettaCluster_" + datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f")
#     elif not isinstance(environment_name, str):
#         raise TypeError(
#             f"The 'environment_name' keyword argument parameter must be of type `str`. Received: {type(environment_name)}"
#         )

#     if not base_dir:
#         base_dir = os.path.abspath(os.curdir)
#     elif not isinstance(base_dir, str):
#         raise TypeError(f"The 'base_dir' keyword argument parameter must be of type `str`. Received: {type(base_dir)}")
#     else:
#         base_dir = os.path.abspath(os.path.expanduser(base_dir))
#     if not os.path.isdir(base_dir):
#         raise NotADirectoryError(
#             f"The 'base_dir' keyword argument parameter must be an existing directory. Received: '{base_dir}'"
#         )

#     _env_config = get_environment_config()
#     environment_manager = _env_config.environment_manager
#     environment_var = get_environment_var()

#     # Extract environment spec from record
#     _instance_kwargs = get_instance_kwargs(scorefile, decoy_name)
#     if "environment" not in _instance_kwargs:
#         raise RuntimeError(
#             "PyRosettaCluster 'environment' instance attribute does not exist. "
#             + "`recreate_environment()` cannot create environment!"
#         )
#     raw_spec = _instance_kwargs.get("environment", None)
#     if not raw_spec:
#         raise RuntimeError(
#             "PyRosettaCluster 'environment' instance attribute is empty. "
#             + "`recreate_environment()` cannot create environment!"
#         )
#     # # Get original environment manager
#     # for line in raw_spec.splitlines():
#     #     if line.startswith(f"# {environment_var}"):
#     #         original_environment_manager = line.split("=")[-1].strip()
#     #         break
#     # else:
#     #     original_environment_manager = "conda" # For legacy PyRosettaCluster simulations with a non-empty raw spec

#     # # Test that current environment manager can recreate original environment
#     # if original_environment_manager == "uv" and environment_manager != "uv":
#     #     raise RuntimeError(
#     #         f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
#     #         + f"an environment manager, but the current environment manager is configured to use "
#     #         + f"'{environment_manager}' as an environment manager! The environment specification "
#     #         + f"is saved in a 'requirements.txt' format, and therefore requires '{original_environment_manager}' "
#     #         + f"to recreate the environment. Please ensure '{original_environment_manager}' is installed "
#     #         + f"and run `export {environment_var}={original_environment_manager}` to properly configure "
#     #         + f"the '{original_environment_manager}' environment manager, then try again. For installation "
#     #         + "instructions, please visit:\n"
#     #         + "https://docs.astral.sh/uv/guides/install-python\n"
#     #     )
#     # elif original_environment_manager != "uv" and environment_manager == "uv":
#     #     raise RuntimeError(
#     #         f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
#     #         + "an environment manager, but the current environment manager is configured to use "
#     #         + f"'{environment_manager}' as an environment manager! The environment specification "
#     #         + "is saved in a YAML file format, and therefore requires 'conda', 'mamba', or 'pixi' to "
#     #         + "recreate the environment. Please ensure that 'conda', 'mamba', or 'pixi' is installed, then "
#     #         + "configure the environment manager by running one of the following commands, then try again:\n"
#     #         + f"To configure 'conda', run: `export {environment_var}=conda`\n"
#     #         + f"To configure 'mamba', run: `export {environment_var}=mamba`\n"
#     #         + f"To configure 'pixi', run:  `export {environment_var}=pixi`\n"
#     #         + "For installation instructions, please visit:\n"
#     #         + "https://docs.anaconda.com/anaconda/install\n"
#     #         + "https://github.com/conda-forge/miniforge\n"
#     #         + "https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html\n"
#     #         + "https://pixi.sh/latest/installation\n"
#     #     )
#     # elif original_environment_manager != environment_manager:
#     #     warnings.warn(
#     #         f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
#     #         + "an environment manager, but the current environment manager is configured to use "
#     #         + f"'{environment_manager}' as an environment manager! Now attempting to recreate the "
#     #         + "environment with cross-manager compatibility using the universal YAML file string...",
#     #         RuntimeWarning,
#     #         stacklevel=2,
#     #     )

#     # Recreate the environment
#     with tempfile.TemporaryDirectory() as tmp_dir:
#         env_create_cmd = _env_config.env_create_cmd(environment_name, raw_spec, tmp_dir, base_dir)
#         env_dir = os.path.join(base_dir, environment_name)
#         print(f"Running environment create command: `{env_create_cmd}`")
#         output = _run_subprocess(env_create_cmd, cwd=env_dir, env=None)
#         print(
#             f"\nEnvironment successfully created using {environment_manager}: '{environment_name}'\nOutput:\n{output}\n",
#             flush=True,
#         )

# # ---------------------------------------

import argparse
import os
import shutil
import subprocess

from typing import Any, Dict, Optional


def run_subprocess(
    cmd: str,
    cwd: Optional[str] = None,
    env: Optional[Dict[str, Any]] = None,
    timeout: float = 3600,
) -> str:
    """Run a shell command and return its stdout, raising RuntimeError on failure."""
    print(f"[INFO] Running command: `{cmd}`")
    try:
        output = subprocess.check_output(
            cmd,
            shell=True,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            text=True,
            cwd=cwd,
            env=env,
            executable="/bin/bash",
        )
        print("[INFO] Subprocess stdout:\n" + output)
        return output
    except subprocess.CalledProcessError as ex:
        print(
            f"Command failed: `{cmd}`\n"
            f"Return code: {ex.returncode}\n"
            f"Output:\n{ex.output}"
        )
        raise RuntimeError(cmd) from ex


def recreate_environment(env_dir: str, env_manager: str, timeout: float):
    """
    Recreate an environment using pixi, uv, conda, or mamba inside `env_dir`.
    The directory must already exist.
    """

    if env_manager == "pixi":
        lock_file = os.path.join(env_dir, "pixi.lock")
        if not os.path.isfile(lock_file):
            raise FileNotFoundError(
                "Please ensure that the pixi 'pixi.lock' file "
                "is in the pixi project directory, then try again."
            )

        for filename in ("pixi.toml", "pyproject.toml"):
            if os.path.isfile(os.path.join(env_dir, filename)):
                break
        else:
            raise FileNotFoundError(
                "Please ensure that the pixi manifest 'pixi.toml' or 'pyproject.toml' file "
                "is in the pixi project directory, then try again."
            )

        env_create_cmd = "pixi install --frozen"

    elif env_manager == "uv":
        req_file = os.path.join(env_dir, "requirements.txt")
        if not os.path.isfile(req_file):
            raise FileNotFoundError(
                "Please ensure that the uv project 'requirements.txt' file "
                "is in the uv project directory, then try again."
            )

        # Install packages strictly from requirements.txt
        env_create_cmd = f"uv pip install -r '{req_file}'"
        # env_create_cmd = f"uv add --requirements '{req_file}'"

    elif env_manager in ("conda", "mamba"):
        yml_file = os.path.join(env_dir, "environment.yml")
        if not os.path.isfile(yml_file):
            raise FileNotFoundError(
                f"Please ensure that the {env_manager} environment 'environment.yml' file "
                f"is in the {env_manager} environment prefix directory, then try again."
            )

        env_create_cmd = f"{env_manager} env create -f '{yml_file}' -p '{env_dir}'"

    else:
        raise ValueError(f"Unsupported environment manager: {env_manager}")

    run_subprocess(env_create_cmd, cwd=env_dir, timeout=timeout)

    print(
        f"[INFO] Environment successfully created using {env_manager} in directory: '{env_dir}'",
        flush=True,
    )


def parse_env_dir(path: str | None) -> str:
    """Validate and normalize the environment directory path."""
    if path is None:
        path = os.path.abspath(os.curdir)
    else:
        if not isinstance(path, str):
            raise argparse.ArgumentTypeError(
                f"The 'env_dir' parameter must be of type `str`. Received: {type(path)}"
            )
        path = os.path.abspath(os.path.expanduser(path))

    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(
            f"The 'env_dir' parameter must be an existing directory. Received: '{path}'"
        )

    if not os.access(path, os.W_OK):
        raise argparse.ArgumentTypeError(
            f"The directory '{path}' is not writable."
        )

    return path


def validate_env_manager(manager: str) -> str:
    """Validate that the environment manager exists on PATH."""
    allowed = ["pixi", "uv", "conda", "mamba"]
    manager = manager.lower()

    if manager not in allowed:
        raise argparse.ArgumentTypeError(
            f"Invalid environment manager '{manager}'. Must be one of: {', '.join(allowed)}"
        )

    if shutil.which(manager) is None:
        raise argparse.ArgumentTypeError(
            f"The environment manager executable '{manager}' was not found on PATH. "
            "Please ensure it is installed."
        )

    return manager


if __name__ == "__main__":
    env_manager_default = os.getenv("PYROSETTACLUSTER_ENVIRONMENT_MANAGER")

    parser = argparse.ArgumentParser(
        description=(
            "Recreate a PyRosettaCluster environment using one of the supported "
            "environment managers ('pixi', 'uv', 'conda', 'mamba')."
        )
    )

    parser.add_argument(
        "--env_dir",
        type=parse_env_dir,
        required=False,
        default=None,
        help=(
            "Directory in which the environment will be created. Must exist and be "
            "writable. Defaults to the current working directory."
        ),
    )

    parser.add_argument(
        "--env_manager",
        type=str,
        required=False,
        default=env_manager_default,
        help=(
            "Environment manager to use: pixi, uv, conda, or mamba.\n"
            "If omitted, the environment variable "
            "`PYROSETTACLUSTER_ENVIRONMENT_MANAGER` will be used if set."
        ),
    )

    parser.add_argument(
        "--timeout",
        type=float,
        required=False,
        default=1800.0,
        help=(
            "Timeout specifying the amount of time in seconds "
            "before any subprocesses are terminated."
        ),
    )

    args = parser.parse_args()

    if args.env_manager is None:
        raise SystemExit(
            "Error: No environment manager was provided.\n"
            "Provide the `--env_manager` flag, or otherwise set the "
            "environment variable 'PYROSETTACLUSTER_ENVIRONMENT_MANAGER'."
        )

    args.env_manager = validate_env_manager(args.env_manager)

    recreate_environment(
        env_dir=args.env_dir,
        env_manager=args.env_manager,
        timeout=args.timeout,
    )
