__author__ = "Jason C. Klima"

import os
import sys

# print("### Debug")
# print("PYTHONPATH:", os.environ.get("PYTHONPATH"))
# print("sys.executable:", sys.executable)
# print("sys.path:")
# for p in sys.path:
#     print(" ", p)
# print("CWD:", os.getcwd())
# print("### Debug")

import argparse
import os
import tempfile

from pyrosetta.distributed.cluster import reproduce, run

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
import pyrosetta
import subprocess
import sys

from pyrosetta.bindings.utility import bind_method
from pyrosetta.distributed.cluster.config import (
    get_environment_manager,
    # get_environment_var,
    # source_domains,
)


def get_env_export_cmd(env_manager: str) -> str:
    """
    Return the appropriate environment export command for the given environment manager.
    Automatically adjusts for Pixi and uv when a project/manifest path is set via
    environment variables or custom run paths.

    Args:
        env_manager: A `str` object of either "pixi", "uv", "mamba", "conda".

    Returns:
        A shell command string for subprocess execution.
    """

    # Default commands
    base_cmds = {
        "pixi": "pixi lock --check || pixi lock --no-install",
        "uv": "uv export --format requirements-txt --frozen",
        "mamba": f"mamba env export --prefix '{sys.prefix}'",
        "conda": f"conda env export --prefix '{sys.prefix}'",
    }

    # Pixi
    if env_manager == "pixi":
        # https://pixi.sh/dev/reference/environment_variables/#environment-variables-set-by-pixi
        manifest_path = os.environ.get("PIXI_PROJECT_MANIFEST")
        if manifest_path:
            # Append --manifest-path flag to both commands in the OR clause
            return (
                f"pixi lock --check --manifest-path '{manifest_path}' || "
                f"pixi lock --no-install --manifest-path '{manifest_path}'"
            )

        return base_cmds["pixi"]

    # Uv
    elif env_manager == "uv":
        project_dir = (
            os.environ.get("UV_PROJECT")
            or os.environ.get("UV_PROJECT_ENVIRONMENT")
        )
        if project_dir:
            # Append --project flag
            return (
                f"uv export --format requirements-txt --frozen --project '{project_dir}'"
            )

        return base_cmds["uv"]

    # Conda / Mamba
    else:
        return base_cmds.get(env_manager, "")


# @bind_method(pyrosetta.distributed.cluster.toolkit)
# @bind_method(pyrosetta.distributed.cluster.converters)
# @bind_method(pyrosetta.distributed.cluster.converter_tasks)
# def get_yml() -> str:
#     """
#     Run environment export command to return a YML file string with the current virtual
#     environment, excluding certain source domains.
#     """
#     def remove_comments(text: str) -> str:
#         """Remove lines starting with '#' from a requirements.txt file string."""
#         return "\n".join(
#             line for line in text.splitlines()
#             if not line.strip().startswith("#")
#         )

#     # _ENV_EXPORT_CMDS = {
#     #     "pixi": "pixi lock --check || pixi lock --no-install", # Updated
#     #     "uv": "uv export --format requirements-txt --frozen",
#     #     "mamba": f"mamba env export --prefix '{sys.prefix}'",
#     #     "conda": f"conda env export --prefix '{sys.prefix}'",
#     # }
#     env_manager = get_environment_manager()
#     # environment_cmd = _ENV_EXPORT_CMDS[env_manager]
#     environment_cmd = get_env_export_cmd(env_manager)
#     print(f"Running environment command: `{environment_cmd}`")
#     if env_manager == "pixi": # Updated
#         subprocess.run(
#             environment_cmd,
#             shell=True,
#             check=True,
#             stderr=subprocess.DEVNULL,
#         )
#         with open("pixi.lock") as f:
#             yml = f.read()
#     else:
#         try:
#             raw_yml = subprocess.check_output(
#                 environment_cmd,
#                 shell=True,
#                 stderr=subprocess.DEVNULL,
#             ).decode()
#         except subprocess.CalledProcessError:
#             raw_yml = ""

#         yml = (
#             (
#                 os.linesep.join(
#                     # [f"# {get_environment_var()}={get_environment_manager()}"]
#                     # + [ # Updated
#                     [
#                         line
#                         for line in raw_yml.split(os.linesep)
#                         if all(
#                             source_domain not in line for source_domain in source_domains
#                         )
#                         and all(not line.startswith(s) for s in ["name:", "prefix:"])
#                         and line
#                     ]
#                 )
#                 + os.linesep
#             )
#             if raw_yml
#             else raw_yml
#         )

#     if env_manager == "uv":
#         yml = remove_comments(yml)

#     return yml


@bind_method(pyrosetta.distributed.cluster.toolkit)
@bind_method(pyrosetta.distributed.cluster.converters)
@bind_method(pyrosetta.distributed.cluster.converter_tasks)
def get_yml() -> str:
    """
    Export the current environment to a YAML string, excluding metadata lines
    depending on the environment manager.
    """

    def remove_comments(text: str) -> str:
        """Remove lines starting with '#'."""
        return "\n".join(
            line for line in text.splitlines() if not line.strip().startswith("#")
        )

    def remove_metadata(text: str) -> str:
        """Remove 'name:' and 'prefix:' lines."""
        filtered_lines = [
            line
            for line in text.splitlines()
            if not line.startswith(("name:", "prefix:")) and line.strip()
        ]
        return "\n".join(filtered_lines) + "\n"


    env_manager = get_environment_manager()
    environment_cmd = get_env_export_cmd(env_manager)

    print(f"Running environment command: `{environment_cmd}`")

    # Handle Pixi separately since it writes a pixi.lock file
    if env_manager == "pixi":
        try:
            subprocess.run(
                environment_cmd,
                shell=True,
                check=True,
                stderr=subprocess.DEVNULL,
            )
            # https://pixi.sh/dev/reference/environment_variables/#environment-variables-set-by-pixi
            manifest_path = os.environ.get("PIXI_PROJECT_MANIFEST")
            lock_path = os.path.join(
                os.path.dirname(manifest_path) if manifest_path else os.getcwd(),
                "pixi.lock",
            )
            with open(lock_path, encoding="utf-8") as f:
                return f.read()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return ""

    # For all other environment managers, run the export command and process the output
    try:
        result = subprocess.run(
            environment_cmd,
            shell=True,
            check=True,
            stderr=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError:
        return ""

    raw_yml = result.stdout.strip()
    if not raw_yml:
        return ""

    if env_manager == "uv":
        return remove_comments(raw_yml)
    elif env_manager in ("conda", "mamba"):
        return remove_metadata(raw_yml)

    raise NotImplementedError(env_manager)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

def create_tasks():
    yield {
        "options": "-ex1 0 -ex1aro 0 -ex1aro_exposed 0 -ex2 0 -ex2aro 0 -ex2aro_exposed 0 -ex3 0 -ex4 0 -lazy_ig 1",
        "extra_options": "-out:level 300 -multithreading:total_threads 1 -ignore_unrecognized_res 1 -load_PDB_components 0",
        "set_logging_handler": "logging",
        "silent": False,
        "seq": "NEW/ENV",
    }


def my_protocol(packed_pose, **kwargs):
    import math
    import pyrosetta
    import pyrosetta.distributed.io as io
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

    assert packed_pose is None, f"The input `PackedPose` object must be `None`. Received: {packed_pose}"
    pose = pyrosetta.pose_from_sequence(kwargs["seq"])

    pose.cache["SEQUENCE"] = kwargs["seq"]
    pose.cache["VALUE"] = math.pi # JSON-compatible

    scorefxn = pyrosetta.create_score_function("ref2015.wts")
    pack_rotamers = PackRotamersMover(
        scorefxn=scorefxn,
        task=pyrosetta.standard_packer_task(pose),
        nloop=3,
    )
    pack_rotamers.apply(pose)
    scorefxn(pose)

    return pose


def print_logs(prc_log, protocol_log):
    for log_file in (prc_log, protocol_log):
        if os.path.isfile(log_file):
            with open(log_file, "r") as f:
                print(f"Output: '{log_file}':", f.read(), sep="\n")
        else:
            print(f"Warning: missing PyRosettaCluster log file: '{log_file}'")


def run_original_simulation(
    env_manager,
    output_path,
    scorefile_name,
    verbose=True,
):
    # Set environment manager
    os.environ["PYROSETTACLUSTER_ENVIRONMENT_MANAGER"] = env_manager
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Setup simulation
        scratch_dir = os.path.join(tmp_dir, "scratch")
        tasks = list(create_tasks())
        decoy_dir_name = "test_decoys"
        logs_dir_name = "logs"
        project_name = "Original_Environment_Simulation"
        simulation_name = env_manager
        protocols = [my_protocol,]
        instance_kwargs = dict(
            tasks=tasks,
            input_packed_pose=None,
            seeds=None,
            decoy_ids=None,
            client=None,
            scheduler=None,
            scratch_dir=scratch_dir,
            cores=None,
            processes=None,
            memory=None,
            min_workers=1,
            max_workers=1,
            nstruct=1,
            dashboard_address=None,
            compressed=False,
            compression=True,
            logging_level="INFO",
            scorefile_name=scorefile_name,
            project_name=project_name,
            simulation_name=simulation_name,
            environment=None,
            output_path=output_path,
            simulation_records_in_scorefile=True,
            decoy_dir_name=decoy_dir_name,
            logs_dir_name=logs_dir_name,
            ignore_errors=False,
            timeout=0.1,
            max_delay_time=1.0,
            sha1=None,
            dry_run=False,
            save_all=False,
            system_info=None,
            pyrosetta_build=None,
            output_decoy_types=[".pdb", ".pkl_pose", ".b64_pose"],
            output_scorefile_types=[".json",],
            norm_task_options=None,
            output_init_file=None,
            protocols=protocols,
            clients_indices=None,
            resources=None,
            author="Alice",
            email=None,
            license="LICENSE.PyRosetta.md"
        )
        # Run simulation
        run(**instance_kwargs)
        # Maybe print log files
        if verbose:
            protocol_log = os.path.join(output_path, logs_dir_name, f"{project_name}_{simulation_name}.log")
            prc_log = os.path.join(output_path, logs_dir_name, "PyRosettaCluster.log")
            print_logs(prc_log, protocol_log)


def run_reproduce_simulation(
    env_manager,
    output_path,
    scorefile_name,
    original_scorefile,
    original_decoy_name,
    verbose=True,
):
    # Set environment manager
    os.environ["PYROSETTACLUSTER_ENVIRONMENT_MANAGER"] = env_manager

    # print("*" * 50)
    # print("Reproduced environment environment file:")
    # print(get_yml())
    # print("*" * 50)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Setup simulation
        scratch_dir = os.path.join(tmp_dir, "scratch")
        logs_dir_name = "logs"
        project_name = "Recreated_Environment_Simulation"
        simulation_name = env_manager
        # Run simulation
        reproduce(
            input_file=None,
            scorefile=original_scorefile,
            decoy_name=original_decoy_name,
            protocols=None, # Automatically detect protocols in frame from scorefile 
            input_packed_pose=None,
            client=None,
            clients=None,
            resources=None,
            instance_kwargs={
                "output_path": output_path,
                "scratch_dir": scratch_dir,
                "sha1": None,
                "scorefile_name": scorefile_name,
                "logs_dir_name": logs_dir_name,
                "project_name": project_name,
                "simulation_name": simulation_name,
                "output_decoy_types": [".pdb", ".pkl_pose", ".b64_pose"],
                "output_scorefile_types": [".json",],
                "author": "Bob",
                "email": None,
                "license": "LICENSE.PyRosetta.md"
            },
        )
        # Maybe print log files
        if verbose:
            protocol_log = os.path.join(output_path, logs_dir_name, f"{project_name}_{simulation_name}.log")
            prc_log = os.path.join(output_path, logs_dir_name, "PyRosettaCluster.log")
            print_logs(prc_log, protocol_log)


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--env_manager', type=str)
    parser.add_argument('--output_path', type=str)
    parser.add_argument('--scorefile_name', type=str)
    parser.add_argument('--original_scorefile', type=str, default=None)
    parser.add_argument('--original_decoy_name', type=str, default=None)
    parser.add_argument('--reproduce', dest='reproduce', action='store_true')
    parser.set_defaults(reproduce=False)
    args = parser.parse_args()
    if not args.reproduce:
        run_original_simulation(
            args.env_manager,
            args.output_path,
            args.scorefile_name,
        )
    else:
        run_reproduce_simulation(
            args.env_manager,
            args.output_path,
            args.scorefile_name,
            args.original_scorefile,
            args.original_decoy_name,
        )
