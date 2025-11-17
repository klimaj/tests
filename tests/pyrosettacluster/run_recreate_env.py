__author__ = "Jason C. Klima"


import argparse
import os
import subprocess

from pathlib import Path

# from pyrosetta.distributed.cluster import recreate_environment

sys.path.append(str(Path(__file__).resolve().parent.parent.parent))
from pyrosettacluster_tests.utils import (
    ROSETTACOMMONS_CONDA_CHANNEL,
    detect_platform,
)


def run_recreate_environment(
    env_manager,
    reproduce_env_dir,
    original_scorefile_path,
    original_decoy_name,
):
    # TODO: remove custom .condarc implementation
    if env_manager in ("conda", "mamba"):
        # Write .condarc file dynamically
        template_condarc_file = Path(__file__).resolve().parent.parent / "conda" / ".condarc"
        condarc_file = Path.cwd() / ".condarc"
        with open(template_condarc_file, "r") as f1, open(condarc_file, "w") as f2:
            f2.write(
                f1.read().format(rosettacommons_conda_channel=ROSETTACOMMONS_CONDA_CHANNEL)
            )
        os.environ["CONDARC"] = str(condarc_file)
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
        install_pyrosetta_file = Path(__file__).resolve().parent.parent / "uv" / "install_pyrosetta.py"
        install_pyrosetta_script = install_pyrosetta_file.read_text()
        subprocess.run(
            ["uv", "run", "--project", str(reproduce_env_dir), "python", "-c", install_pyrosetta_script],
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
