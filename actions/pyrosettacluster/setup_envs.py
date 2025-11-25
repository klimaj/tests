__author__ = "Jason C. Klima"


import argparse
import os
import subprocess
import sys
import tempfile

from pathlib import Path

from actions.pyrosettacluster.utils import (
    ROSETTACOMMONS_CONDA_CHANNEL,
    detect_platform,
)


def setup_pixi_environment(env_dir, timeout):
    """
    Create a fresh pixi environment containing 'pyrosetta' and 'pyrosetta-distributed' packages.

    Note: this requires that `pixi` is an executable installed and on `${PATH}`. This function:
    - detects the current Python version
    - detects the current platform (linux/mac/windows)
    - writes a compatible 'pixi.toml' file
    - runs `pixi install` to build a new pixi environment
    """
    # Detect current Python version
    # Pixi/Conda only publishes minor versions (e.g., `3.9.x`), not every patch release (e.g., `3.9.25`)
    # Therefore, use a wild-card to fetch the latest micro version available
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}.*"
    py_feature = f"py{sys.version_info.major}{sys.version_info.minor}"

    # Detect platform
    plat = detect_platform()

    template_toml_file = Path(__file__).resolve().parent / "pixi" / "pixi.toml"
    with open(template_toml_file, "r") as f:
        pixi_toml = f.read().format(
            rosettacommons_conda_channel=ROSETTACOMMONS_CONDA_CHANNEL,
            name=os.path.basename(env_dir),
            plat=plat,
            py_version=py_version,
            py_feature=py_feature,
        )

    # Create environment directory
    env_path = Path(env_dir)
    env_path.mkdir(parents=True, exist_ok=False)
    toml_path = env_path / "pixi.toml"
    toml_path.write_text(pixi_toml.strip() + "\n")
    print(f"Created '{toml_path}' file for platform '{plat}' and Python-{py_version}")

    # Install pixi environment
    print("Running `pixi install`...")
    subprocess.run(["pixi", "install"], cwd=env_path, check=True, timeout=timeout)

    print(f"Pixi environment setup complete in directory: '{env_path}'.")


def setup_uv_environment(env_dir, timeout):
    """
    Create a fresh uv environment using the 'pyrosetta-installer' package.

    Note: this requires that `uv` is an executable installed and on `${PATH}`. This function:
    - detects the current Python version
    - adds 'pyrosetta-installer' via `uv add ...`
    - runs the PyRosetta installer using `uv run python -c ...`
    """
    env_path = Path(env_dir)
    if env_path.exists():
        raise FileExistsError(f"The specified uv environment path already exists: '{env_path}'.")

    # Create uv environment using the current Python
    print(f"Creating uv environment at '{env_path}'...")
    os.makedirs(env_dir, exist_ok=True)
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    template_toml_file = Path(__file__).resolve().parent / "uv" / "pyproject.toml"
    with open(template_toml_file, "r") as f:
        toml_data = f.read().format(
            name=os.path.basename(env_dir),
            py_version=py_version,
        )
    toml_file = env_path / "pyproject.toml"
    with open(toml_file, "w") as f:
        f.write(toml_data)

    # Install pyrosetta-installer
    print("Adding 'pyrosetta-installer', 'pip', and `pyrosetta.distributed` depedencies to uv environment...")
    requirements_txt_file = Path(__file__).resolve().parent / "uv" / "requirements.txt"
    subprocess.run(
        [
            "uv",
            "add",
            "--project", str(env_path),
            "--requirements", str(requirements_txt_file),
        ],
        check=True,
        cwd=str(env_path),
        timeout=timeout,
    )

    # Run PyRosetta installer with mirror fallback
    print("Running PyRosetta installer in uv environment...")
    install_pyrosetta_file = Path(__file__).resolve().parent.parent.parent / "pyrosettacluster" / "install_pyrosetta.py"
    subprocess.run(
        ["uv", "run", "--project", str(env_path), "python", install_pyrosetta_file],
        check=True,
        timeout=timeout,
    )

    print(f"Uv environment setup complete in directory: '{env_path}'.")


def setup_conda_environment(env_dir, timeout, env_manager="conda"):
    """
    Create a fresh conda/mamba environment containing 'pyrosetta' and 'pyrosetta-distributed' packages.

    Note: this requires that `conda` or `mamba` is an executable installed and on `${PATH}`. This function:
    - detects the current Python version
    - detects the current platform (linux/mac/windows)
    - writes a temporary 'environment.yml' file
    - runs `{env_manager} env create -f environment.yml ...` to build a new environment
    """

    # Validate environment manager
    if env_manager not in ("conda", "mamba"):
        raise ValueError("The 'env_manager' keyword argument parameter must be either 'conda' or 'mamba'.")

    # Detect Python version
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}"

    # Detect platform
    plat = detect_platform()

    # Build the conda environment file dynamically
    name = os.path.basename(env_dir)
    env_yaml_file = Path(__file__).resolve().parent / "conda" / "environment.yml"
    env_yaml = env_yaml_file.read_text().format(
        name=name,
        rosettacommons_conda_channel=ROSETTACOMMONS_CONDA_CHANNEL,
        py_version=py_version,
    )

    # Create environment directory
    env_path = Path(env_dir)
    with tempfile.TemporaryDirectory() as tmp_dir:
        env_file = os.path.join(tmp_dir, "environment.yml")
        with open(env_file, "w") as f:
            f.write(env_yaml.strip() + "\n")
        print(f"Created '{env_file}' for platform '{plat}' and Python-{py_version}")

        # Run conda/mamba to build the environment
        print(f"Running: {env_manager} env create -f {env_file} -p {env_path}")
        subprocess.run(
            [env_manager, "env", "create", "-f", env_file, "-p", str(env_path)],
            check=True,
            timeout=timeout,
        )

    print(f"{env_manager.title()} environment setup complete in directory: '{env_path}'.")


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--env_manager', type=str)
    parser.add_argument('--env_dir', type=str)
    parser.add_argument('--timeout', type=float)
    args = parser.parse_args()
    if args.env_manager == "pixi":
        setup_pixi_environment(args.env_dir, args.timeout)
    elif args.env_manager == "uv":
        setup_uv_environment(args.env_dir, args.timeout)
    elif args.env_manager in ("conda", "mamba"):
        setup_conda_environment(args.env_dir, args.timeout, env_manager=args.env_manager)
    else:
        raise NotImplementedError(args.env_manager)
