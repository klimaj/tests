__author__ = "Jason C. Klima"


import argparse
import os
import subprocess
import sys
import tempfile
import textwrap

from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)))
from pyrosettacluster_tests.utils import (
    ROSETTACOMMONS_CONDA_CHANNEL,
    detect_platform,
)


def setup_pixi_environment(env_dir):
    """
    Create a fresh pixi environment containing 'pyrosetta' and 'pyrosetta-distributed' packages.

    Note: this requires that `pixi` is an executable installed and on `${PATH}`. This function:
    - detects the current Python version
    - detects the current platform (linux/mac/windows)
    - writes a compatible 'pixi.toml' file
    - runs `pixi install` to build a new pixi environment
    """
    # Detect Python version
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    py_feature = f"py{sys.version_info.major}{sys.version_info.minor}"

    # Detect platform
    plat = detect_platform()

    # Build 'pixi.toml' file dynamically
    # pixi_toml = textwrap.dedent(f"""
    # [workspace]
    # channels = ["{ROSETTACOMMONS_CONDA_CHANNEL}", "conda-forge"]
    # name = "pyrosetta-pixi"
    # platforms = ["{plat}"]
    # version = "1.0.0"

    # [dependencies]
    # pyrosetta = "*"

    # [pypi-dependencies]
    # pyrosetta-distributed = "*"

    # [feature.{py_feature}.dependencies]
    # python = "{py_version}.*"

    # [environments]
    # {py_feature} = ["{py_feature}"]
    # """)

    template_toml_file = os.path.join(os.path.dirname(__file__), os.pardir, "pixi", "pixi.toml")
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
    subprocess.run(["pixi", "install"], cwd=env_path, check=True)

    print(f"Pixi environment setup complete in directory: '{env_path}'.")


def setup_uv_environment(env_dir):
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

    # Detect Python version
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}"

    # Create uv environment using the current Python
    print(f"Creating uv environment at '{env_path}'...")
    subprocess.run(
        ["uv", "init", str(env_path), "--python", py_version],
        check=True,
    )

    # Install pyrosetta-installer
    print("Adding 'pyrosetta-installer', 'pip', and `pyrosetta.distributed` depedencies to uv environment...")
    # subprocess.run(
    #     [
    #         "uv",
    #         "add",
    #         "--project", str(env_path),
    #         "pyrosetta-installer>=0.1.2",
    #         "pip>=25.3",
    #         "attrs>=19.3.0",
    #         "billiard>=3.6.3.0",
    #         "blosc>=1.8.3",
    #         "cloudpickle>=1.5.0",
    #         "cryptography>=2.8",
    #         "dask>=2.16.0",
    #         "dask-jobqueue>=0.7.0",
    #         "distributed>=2.16.0",
    #         "gitpython>=3.1.1",
    #         "jupyter>=1.0.0",
    #         "numpy>=1.17.3",
    #         "pandas>=0.25.2",
    #         "py3Dmol>=0.8.0",
    #         "python-xz>=0.4.0",
    #         "scipy>=1.4.1",
    #         "traitlets>=4.3.3",
    #     ],
    #     check=True,
    #     cwd=str(env_path),
    # )
    requirements_txt_file = Path(__file__).resolve().parent.parent / "uv" / "requirements.txt"
    subprocess.run(
        [
            "uv",
            "add",
            "--project", str(env_path),
            "--requirements", str(requirements_txt_file),
        ],
        check=True,
        cwd=str(env_path),
    )

    # Run PyRosetta installer with mirror fallback
    print("Running PyRosetta installer in uv environment...")
    # install_script = textwrap.dedent("""
    #     import pyrosetta_installer
    #     try:
    #         pyrosetta_installer.install_pyrosetta(
    #             distributed=False,
    #             serialization=True,
    #             skip_if_installed=False,
    #             mirror=0
    #         )
    #     except Exception as e:
    #         print(f"PyRosetta installation with 'mirror=0' failed: {e}. Retrying with 'mirror=1'.")
    #         pyrosetta_installer.install_pyrosetta(
    #             distributed=False,
    #             serialization=True,
    #             skip_if_installed=False,
    #             mirror=1
    #         )
    # """)
    install_pyrosetta_file = Path(__file__).resolve().parent.parent / "uv" / "install_pyrosetta.py"
    install_pyrosetta_script = install_pyrosetta_file.read_text()
    subprocess.run(
        ["uv", "run", "--project", str(env_path), "python", "-c", install_pyrosetta_script],
        check=True,
    )

    print(f"Uv environment setup complete in directory: '{env_path}'.")


def setup_conda_environment(env_dir, env_manager="conda"):
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
    env_yaml = textwrap.dedent(f"""
    name: {name}
    channels:
      - {ROSETTACOMMONS_CONDA_CHANNEL}
      - conda-forge
    dependencies:
      - python={py_version}.*
      - pyrosetta
      - pip
      - pip:
        - pyrosetta-distributed
    """)

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
        )

    print(f"{env_manager.title()} environment setup complete in directory: '{env_path}'.")


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--env_manager', type=str)
    parser.add_argument('--env_dir', type=str)
    args = parser.parse_args()
    if args.env_manager == "pixi":
        setup_pixi_environment(args.env_dir)
    elif args.env_manager == "uv":
        setup_uv_environment(args.env_dir)
    elif args.env_manager in ("conda", "mamba"):
        setup_conda_environment(args.env_dir, env_manager=args.env_manager)
    else:
        raise NotImplementedError(args.env_manager)
