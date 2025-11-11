__author__ = "Jason C. Klima"


import argparse
import os
import subprocess
import textwrap

from pyrosetta.distributed.cluster import recreate_environment


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
            ["uv", "run", "-p", str(reproduce_env_dir), "python", "-c", install_script],
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
