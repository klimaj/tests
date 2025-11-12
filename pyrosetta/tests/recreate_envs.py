__author__ = "Jason C. Klima"


import argparse
import os
import subprocess
import textwrap

# from pyrosetta.distributed.cluster import recreate_environment

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

from datetime import datetime
from functools import lru_cache
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.core.pose import Pose
from typing import (
    Dict,
    Generic,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
)

from pyrosetta.distributed.cluster.tools import get_instance_kwargs

G = TypeVar("G")

class EnvironmentConfig(Generic[G]):
    _ENV_VAR: str = "PYROSETTACLUSTER_ENVIRONMENT_MANAGER"
    _ENV_MANAGERS: Tuple[str, ...] = ("pixi", "uv", "mamba", "conda")
    _ENV_EXPORT_CMDS: Dict[str, str] = {
        "pixi": "pixi workspace export conda-environment",
        "uv": "uv export --format requirements-txt --frozen",
        "mamba": f"mamba env export --prefix '{sys.prefix}'",
        "conda": f"conda env export --prefix '{sys.prefix}'",
    }

    def __init__(self) -> None:
        _env_var_manager = os.environ.get(EnvironmentConfig._ENV_VAR, None)
        if _env_var_manager:
            self.environment_manager = _env_var_manager
            logging.debug(
                "Configuring environment manager for PyRosettaCluster from operating system "
                + f"environment variable: {EnvironmentConfig._ENV_VAR}={self.environment_manager}"
            )
            if self.environment_manager not in EnvironmentConfig._ENV_MANAGERS:
                raise ValueError(
                    "The '{0}' environment variable must be in: '{1}'. Received: '{2}'.".format(
                        EnvironmentConfig._ENV_VAR,
                        EnvironmentConfig._ENV_MANAGERS,
                        self.environment_manager,
                    )
                )
        else:
            for _manager in EnvironmentConfig._ENV_MANAGERS:
                if shutil.which(_manager):
                    self.environment_manager = _manager
                    logging.debug(f"Configuring environment manager for PyRosettaCluster: '{_manager}'")
                    break
            else:
                self.environment_manager = "conda"
                warnings.warn(
                    f"Warning: could not configure an environment manager for PyRosettaCluster. "
                    + "Please ensure that either of 'pixi', 'uv', 'mamba', or 'conda' is installed. "
                    + "Using 'conda' as the default environment manager.",
                    UserWarning,
                    stacklevel=7,
                )

    @property
    def env_export_cmd(self) -> str:
        return self._ENV_EXPORT_CMDS[self.environment_manager]

    def env_create_cmd(
        self, environment_name: str, raw_spec: str, tmp_dir: str, base_dir: str
    ) -> Union[str, NoReturn]:
        # Create a project directory for uv/pixi, or prefix directory for conda/mamba
        project_dir = os.path.join(base_dir, environment_name)
        # Raise exception if the project directory exists
        if os.path.isdir(project_dir):
            if self.environment_manager in ("conda", "mamba"):
                _err_msg = f"The {self.environment_manager} environment prefix directory already exists: '{project_dir}'"
            elif self.environment_manager in ("uv", "pixi"):
                _err_msg = f"The {self.environment_manager} project directory already exists: '{project_dir}'"
            else:
                raise RuntimeError(f"Unsupported environment manager: '{self.environment_manager}'")
            raise IsADirectoryError(_err_msg)
        os.makedirs(project_dir, exist_ok=False)

        if self.environment_manager in ("conda", "mamba"):
            yml_file = os.path.join(tmp_dir, f"{environment_name}.yml")
            with open(yml_file, "w") as f:
                f.write(raw_spec)

            if self.environment_manager == "conda":
                return f"conda env create -f '{yml_file}' -p '{project_dir}'"

            elif self.environment_manager == "mamba":
                return f"mamba env create -f '{yml_file}' -p '{project_dir}'"
            # Updated

        elif self.environment_manager == "uv":
            # Write the requirements.txt file
            req_file = os.path.join(tmp_dir, f"{environment_name}.txt")
            with open(req_file, "w") as f:
                f.write(raw_spec)
            # Determine the current Python version
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
            return (
                f"uv init '{project_dir}' --python {py_version} && "
                f"uv add --project '{project_dir}' pip && "
                f"uv pip sync '{req_file}' --project '{project_dir}'"
            ) # Updated

        elif self.environment_manager == "pixi": # Updated
            subprocess.run(f"pixi init '{project_dir}'", shell=True, stderr=subprocess.STDOUT)
            lock_file = os.path.join(project_dir, "pixi.lock")
            with open(lock_file, "w") as f:
                f.write(raw_spec)
            # return f"pixi init '{project_dir}' && pixi install --locked --manifest-path '{project_dir}'" # Updated
            return f"pixi install --locked --manifest-path '{project_dir}'" # Updated

        raise RuntimeError(f"Unsupported environment manager: '{self.environment_manager}'")


@lru_cache(maxsize=1)
def get_environment_config() -> EnvironmentConfig:
    """Return an instance of the `EnvironmentConfig` class on the host process."""
    return EnvironmentConfig()


def get_environment_manager() -> str:
    """Get the configured environment manager."""
    return get_environment_config().environment_manager


def get_environment_cmd() -> str:
    """Get the configured environment export command."""
    return get_environment_config().env_export_cmd


def get_environment_var() -> str:
    """Get the PyRosettaCluster operating system environment variable name."""
    return EnvironmentConfig._ENV_VAR


def recreate_environment(
    environment_name: Optional[str] = None,
    input_file: Optional[Union[str, Pose, PackedPose]] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
    timeout: Optional[int] = None,
    base_dir: Optional[str] = None,
) -> Optional[NoReturn]:
    """
    *Warning*: This function runs a subprocess with one of the following commands:
      - `conda env create ...`: when 'conda' is an executable
      - `mamba env create ...`: when 'mamba' is an executable
      - `uv pip ...`: when 'uv' is an executable
      - `pixi install ...`: when 'pixi' is an executable
    Installing certain packages may not be secure, so please only run with input files you trust.
    Learn more about PyPI security `here <https://pypi.org/security>`_ and conda security `here <https://www.anaconda.com/docs/reference/security>`_.

    Given an input file that was written by PyRosettaCluster, or a scorefile
    and a decoy name that was written by PyRosettaCluster, recreate the
    environment that was used to generate the decoy with a new environment name.

    The environment manager used (i.e., either 'conda', 'mamba', 'uv', or 'pixi') is
    automatically determined from the operating system environment variable
    'PYROSETTACLUSTER_ENVIRONMENT_MANAGER' if exported, or otherwise it is automatically
    selected based on which manager is available on the system. If no executables
    are explicitly found, then 'conda' is used by default.

    Args:
        environment_name: A `str` object specifying the new name of the environment
            to recreate. If using 'conda' and 'mamba', this is the prefix directory that will
            be created in the 'base_dir' directory. If using 'uv' or 'pixi', this is the
            local project directory name that will be created in the 'base_dir' directory.
            Default: 'PyRosettaCluster_' + datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f")
        input_file: A `str` object specifying the path to the '.pdb', '.pdb.bz2', '.pkl_pose',
            '.pkl_pose.bz2', '.b64_pose', or '.b64_pose.bz2' file, or a `Pose` or `PackedPose`
            object, from which to extract PyRosettaCluster instance kwargs. If 'input_file' is
            provided, then ignore the 'scorefile' and 'decoy_name' keyword argument parameters.
            Default: None
        scorefile: A `str` object specifying the path to the JSON-formatted scorefile
            (or pickled `pandas.DataFrame` scorefile) from a PyRosettaCluster simulation
            from which to extract PyRosettaCluster instance kwargs. If 'scorefile'
            is provided, 'decoy_name' must also be provided. The scorefile must
            contain full simulation records from the original production run
            (i.e., 'simulation_records_in_scorefile' was set to `True`).
            Default: None
        decoy_name: A `str` object specifying the decoy name for which to extract
            PyRosettaCluster instance kwargs. Must be provided if 'scorefile'
            is used.
            Default: None
        timeout: An `int` object specifying the timeout in seconds before any
            subprocesses are terminated.
            Default: None
        base_dir: A `str` object specifying the base directory in which to create
            the environment.
            Default: `.`

    Returns:
        None
    """
    def _run_subprocess(cmd) -> str:
        try:
            return subprocess.check_output(
                cmd,
                shell=True,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                text=True,
                executable="/bin/bash", # Ensure `&&` works properly
            )
        except subprocess.CalledProcessError as ex:
            raise RuntimeError(
                f"Command failed: `{cmd}`\n"
                f"Return code: {ex.returncode}\n"
                f"Output:\n{ex.output}"
            ) from ex

    if not environment_name:
        environment_name = "PyRosettaCluster_" + datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f")
    elif not isinstance(environment_name, str):
        raise TypeError(
            f"The 'environment_name' keyword argument parameter must be of type `str`. Received: {type(environment_name)}"
        )

    if not base_dir:
        base_dir = os.path.abspath(os.curdir)
    elif not isinstance(base_dir, str):
        raise TypeError(f"The 'base_dir' keyword argument parameter must be of type `str`. Received: {type(base_dir)}")
    else:
        base_dir = os.path.abspath(os.path.expanduser(base_dir))
    if not os.path.isdir(base_dir):
        raise NotADirectoryError(
            f"The 'base_dir' keyword argument parameter must be an existing directory. Received: '{base_dir}'"
        )

    _env_config = get_environment_config()
    environment_manager = _env_config.environment_manager
    environment_var = get_environment_var()

    # Extract environment spec from record
    _instance_kwargs = get_instance_kwargs(
        input_file=input_file,
        scorefile=scorefile,
        decoy_name=decoy_name,
        skip_corrections=None,
    )
    if "environment" not in _instance_kwargs:
        raise RuntimeError(
            "PyRosettaCluster 'environment' instance attribute does not exist. "
            + "`recreate_environment()` cannot create environment!"
        )
    raw_spec = _instance_kwargs.get("environment", None)
    if not raw_spec:
        raise RuntimeError(
            "PyRosettaCluster 'environment' instance attribute is empty. "
            + "`recreate_environment()` cannot create environment!"
        )
    # # Get original environment manager
    # for line in raw_spec.splitlines():
    #     if line.startswith(f"# {environment_var}"):
    #         original_environment_manager = line.split("=")[-1].strip()
    #         break
    # else:
    #     original_environment_manager = "conda" # For legacy PyRosettaCluster simulations with a non-empty raw spec

    # # Test that current environment manager can recreate original environment
    # if original_environment_manager == "uv" and environment_manager != "uv":
    #     raise RuntimeError(
    #         f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
    #         + f"an environment manager, but the current environment manager is configured to use "
    #         + f"'{environment_manager}' as an environment manager! The environment specification "
    #         + f"is saved in a 'requirements.txt' format, and therefore requires '{original_environment_manager}' "
    #         + f"to recreate the environment. Please ensure '{original_environment_manager}' is installed "
    #         + f"and run `export {environment_var}={original_environment_manager}` to properly configure "
    #         + f"the '{original_environment_manager}' environment manager, then try again. For installation "
    #         + "instructions, please visit:\n"
    #         + "https://docs.astral.sh/uv/guides/install-python\n"
    #     )
    # elif original_environment_manager != "uv" and environment_manager == "uv":
    #     raise RuntimeError(
    #         f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
    #         + "an environment manager, but the current environment manager is configured to use "
    #         + f"'{environment_manager}' as an environment manager! The environment specification "
    #         + "is saved in a YAML file format, and therefore requires 'conda', 'mamba', or 'pixi' to "
    #         + "recreate the environment. Please ensure that 'conda', 'mamba', or 'pixi' is installed, then "
    #         + "configure the environment manager by running one of the following commands, then try again:\n"
    #         + f"To configure 'conda', run: `export {environment_var}=conda`\n"
    #         + f"To configure 'mamba', run: `export {environment_var}=mamba`\n"
    #         + f"To configure 'pixi', run:  `export {environment_var}=pixi`\n"
    #         + "For installation instructions, please visit:\n"
    #         + "https://docs.anaconda.com/anaconda/install\n"
    #         + "https://github.com/conda-forge/miniforge\n"
    #         + "https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html\n"
    #         + "https://pixi.sh/latest/installation\n"
    #     )
    # elif original_environment_manager != environment_manager:
    #     warnings.warn(
    #         f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
    #         + "an environment manager, but the current environment manager is configured to use "
    #         + f"'{environment_manager}' as an environment manager! Now attempting to recreate the "
    #         + "environment with cross-manager compatibility using the universal YAML file string...",
    #         RuntimeWarning,
    #         stacklevel=2,
    #     )

    # Recreate the environment
    with tempfile.TemporaryDirectory() as tmp_dir:
        env_create_cmd = _env_config.env_create_cmd(environment_name, raw_spec, tmp_dir, base_dir)
        print(f"Running environment create command: `{env_create_cmd}`")
        output = _run_subprocess(env_create_cmd)
        print(
            f"\nEnvironment successfully created using {environment_manager}: '{environment_name}'\nOutput:\n{output}\n",
            flush=True,
        )
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

ROSETTACOMMONS_CONDA_CHANNEL = "https://conda.rosettacommons.org"


def run_recreate_environment(
    env_manager,
    reproduce_env_dir,
    original_scorefile_path,
    original_decoy_name,
):
    if env_manager in ("conda", "mamba"):
        # Write .condarc file dynamically
        os.environ["CONDARC"] = os.path.join(os.getcwd(), ".condarc")
        with open(os.environ["CONDARC"], "w") as f:
            f.write("channels:\n")
            f.write(f"  - {ROSETTACOMMONS_CONDA_CHANNEL}\n")
    # Set environment manager
    os.environ["PYROSETTACLUSTER_ENVIRONMENT_MANAGER"] = env_manager
    # Setup parameters
    environment_name = os.path.basename(reproduce_env_dir)
    base_dir = os.path.dirname(reproduce_env_dir)
    # Recreate environment
    recreate_environment(
        environment_name=environment_name,
        input_file=None,
        scorefile=original_scorefile_path,
        decoy_name=original_decoy_name,
        timeout=999,
        base_dir=base_dir,
    )
    if env_manager == "uv":
        # The recreated uv environment uses the PyPI 'pyrosetta-installer' package, which does not allow specifying PyRosetta version.
        # Therefore, installing the correct PyRosetta version in the recreated uv environment depends fortuitously on a prompt
        # uv environment recreation after the original uv environment creation.
        print("Running PyRosetta installer in recreated uv environment...")
        # Run PyRosetta installer with mirror fallback
        install_script = textwrap.dedent("""
            import pyrosetta_installer
            try:
                pyrosetta_installer.install_pyrosetta(
                    distributed=False,
                    serialization=True,
                    skip_if_installed=True,
                    mirror=0
                )
            except Exception as e:
                print(f"Recreated PyRosetta installation with 'mirror=0' failed: {e}. Retrying with 'mirror=1'.")
                pyrosetta_installer.install_pyrosetta(
                    distributed=False,
                    serialization=True,
                    skip_if_installed=True,
                    mirror=1
                )
        """)
        subprocess.run(
            ["uv", "run", "--project", str(reproduce_env_dir), "--active", "python", "-c", install_script],
            check=True,
        )

if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--env_manager', type=str)
    parser.add_argument('--reproduce_env_dir', type=str)
    parser.add_argument('--original_scorefile_path', type=str)
    parser.add_argument('--original_decoy_name', type=str)
    args = parser.parse_args()
    run_recreate_environment(
        args.env_manager,
        args.reproduce_env_dir,
        args.original_scorefile_path,
        args.original_decoy_name,
    )
