#!/usr/bin/env python3
"""
*Warning*: This script runs a subprocess with one of the following commands:
    - `conda env create ...`: when 'conda' is an executable
    - `mamba env create ...`: when 'mamba' is an executable
    - `uv pip sync ...`: when 'uv' is an executable
    - `pixi install ...`: when 'pixi' is an executable
Installing certain packages may not be secure, so please only run with input files you trust.
Learn more about PyPI security at <https://pypi.org/security> and conda security
at <https://www.anaconda.com/docs/reference/security>.

Given an environment directory with dumped files that were written by PyRosettaCluster,
recreate the environment that was used to generate the decoy with a new environment name.
The environment manager used (i.e., either 'conda', 'mamba', 'uv', or 'pixi') is
automatically determined from the operating system environment variable
'PYROSETTACLUSTER_ENVIRONMENT_MANAGER' if exported, or otherwise it must be
provided using the `--env_manager` flag. Run `./recreate_env.py --help` for more details.
"""

__author__ = "Jason C. Klima"


import argparse
import os
import shutil
import subprocess
import tempfile

from pathlib import Path
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
        env_create_cmd = f"uv venv --project '{env_dir}' && uv pip sync --project '{env_dir}' '{req_file}'"

    elif env_manager == "conda":
        yml_file = os.path.join(env_dir, "environment.yml")
        if not os.path.isfile(yml_file):
            raise FileNotFoundError(
                "Please ensure that the conda environment 'environment.yml' file is in the environment directory."
            )
        env_create_cmd = f"conda env create -f '{yml_file}' -p '{env_dir}'"

    elif env_manager == "mamba":
        yml_file = os.path.join(env_dir, "environment.yml")
        if not os.path.isfile(yml_file):
            raise FileNotFoundError(
                "Please ensure that the mamba environment 'environment.yml' file is in the environment directory."
            )
        # Transactional Mamba logic
        with tempfile.TemporaryDirectory(prefix="env_backup_") as temp_dir:
            # Move all files in `env_dir` into `temp_dir`
            for f in os.listdir(env_dir):
                shutil.move(os.path.join(env_dir, f), temp_dir)
            temp_yml = os.path.join(temp_dir, "environment.yml")
            env_create_cmd = f"mamba env create -f '{temp_yml}' -p '{env_dir}'"
            try:
                run_subprocess(env_create_cmd, cwd=env_dir, timeout=timeout)
            except Exception:
                # On failure, restore all original files from `temp_dir`
                for f in os.listdir(temp_dir):
                    shutil.move(os.path.join(temp_dir, f), env_dir)
                raise
            else:
                # On success, restore original files that don't conflict
                for f in os.listdir(temp_dir):
                    target = os.path.join(env_dir, f)
                    if not os.path.exists(target):
                        shutil.move(os.path.join(temp_dir, f), env_dir)
                    else:
                        print(f"[WARNING] Existing file is being overwritten in the created environment directory: '{target}'")

    else:
        raise ValueError(f"Unsupported environment manager: {env_manager}")

    # For all managers except mamba
    if env_manager != "mamba":
        run_subprocess(env_create_cmd, cwd=env_dir, timeout=timeout)
    
    if env_manager == "uv":
        # The recreated uv environment uses the PyPI 'pyrosetta-installer' package, which does not allow specifying PyRosetta version.
        # Therefore, installing the correct PyRosetta version in the recreated uv environment depends fortuitously on a prompt
        # uv environment recreation after the original uv environment creation.
        print("[INFO] Running PyRosetta installer in uv environment...")
        install_pyrosetta_file = Path(__file__).resolve().parent / "install_pyrosetta.py"
        run_subprocess(
            f"uv run --project '{env_dir}' python '{install_pyrosetta_file}'",
            cwd=env_dir,
            timeout=timeout,
        )

    print(
        f"[INFO] Environment successfully created using {env_manager} in directory: '{env_dir}'",
        flush=True,
    )


def parse_env_dir(path: Optional[str]) -> str:
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
            "If omitted, the environment variable 'PYROSETTACLUSTER_ENVIRONMENT_MANAGER' "
            "will be used if set."
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
            "No environment manager was provided. Please provide "
            "the `--env_manager` flag, or otherwise set the "
            "environment variable 'PYROSETTACLUSTER_ENVIRONMENT_MANAGER'."
        )

    args.env_manager = validate_env_manager(args.env_manager)

    recreate_environment(
        env_dir=args.env_dir,
        env_manager=args.env_manager,
        timeout=args.timeout,
    )
