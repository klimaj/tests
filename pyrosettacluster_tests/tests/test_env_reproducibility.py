"""
PyRosettaCluster environment reproducibility unit test suite using the `unittest` framework.
"""

__author__ = "Jason C. Klima"


import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
import uuid

from pathlib import Path


class TestEnvironmentReproducibility(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workdir = tempfile.TemporaryDirectory()
        cls.run_tag = uuid.uuid4().hex[:12]

    @classmethod
    def tearDownClass(cls):
        try:
            cls.workdir.cleanup()
        except Exception as ex:
            print(f"Warning: failed to cleanup temporary directory: {ex}... Continuing.")
        os.environ.pop("PYROSETTACLUSTER_ENVIRONMENT_MANAGER", None)

    @staticmethod
    def run_subprocess(cmd, module_dir=None, cwd=None, live_output=False):
        print("Running command:", cmd)
        if module_dir:
            env = os.environ.copy()
            pythonpath = os.environ.get("PYTHONPATH")
            env["PYTHONPATH"] = f"{module_dir}{os.pathsep + pythonpath if pythonpath else ''}"
        else:
            env = None
        cmd_list = shlex.split(cmd)
        try:
            if live_output:
                # Use live output streaming for GitHub Actions visibility
                process = subprocess.Popen(
                    cmd_list,
                    cwd=cwd,
                    env=env,
                    shell=False,
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                    text=True,
                )
                returncode = process.wait()
            else:
                # Capture output for debugging
                result = subprocess.run(
                    cmd_list,
                    cwd=cwd,
                    env=env,
                    shell=False,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                returncode = result.returncode

                print("------", flush=True)
                print("\nStdout:\n", flush=True)
                print(result.stdout, flush=True)
                print("\nStderr:\n", flush=True)
                print(result.stderr, flush=True)
                print("------", flush=True)

            if returncode != 0:
                raise subprocess.CalledProcessError(returncode, cmd)
        except subprocess.CalledProcessError as ex:
            print(f"Subprocess command failed (return code: {ex.returncode}): {cmd}", flush=True)
            raise
        except Exception as ex:
            print(f"Unexpected error in subprocess: {ex}", flush=True)
            raise
        else:
            print(f"Return code: {returncode}", flush=True)

        return returncode

    def recreate_environment_test(self, environment_manager="conda"):
        """Test for PyRosettaCluster decoy reproducibility in a recreated virtual environment."""
        self.assertIn(environment_manager, ("conda", "mamba", "uv", "pixi"))

        test_script = os.path.join(os.path.dirname(__file__), "recreate_environment_test_runs.py")

        # Create new environment
        original_env_name = f"{environment_manager}_env_{self.run_tag}"
        original_env_dir = os.path.join(self.workdir.name, original_env_name)
        setup_env_script = os.path.join(os.path.dirname(__file__), "setup_envs.py")
        module = os.path.splitext(os.path.basename(setup_env_script))[0]
        cmd = "{0} -u -m {1} --env_manager '{2}' --env_dir '{3}'".format(
            sys.executable,
            module,
            environment_manager,
            original_env_dir,
        )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=os.path.dirname(setup_env_script),
            cwd=None,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Run original simulation inside new environment
        original_output_path = os.path.join(original_env_dir, f"{environment_manager}_original_outputs")
        original_scorefile_name = "test_scores.json"
        module = os.path.splitext(os.path.basename(test_script))[0]
        if environment_manager == "pixi":
            cmd = (
                f"pixi run python -u -m {module} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{original_output_path}' "
                f"--scorefile_name '{original_scorefile_name}'"
            )
        elif environment_manager == "uv":
            cmd = (
                f"uv run --project {original_env_dir} python -u -m {module} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{original_output_path}' "
                f"--scorefile_name '{original_scorefile_name}'"
            )
        elif environment_manager in ("conda", "mamba"):
            cmd = (
                f"conda run -p {original_env_dir} python -u -m {module} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{original_output_path}' "
                f"--scorefile_name '{original_scorefile_name}'"
            )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=os.path.dirname(test_script),
            # For pixi, activate the original pixi environment context
            # For conda/mamba/uv, run from environment directory for consistency with pixi workflow
            cwd=original_env_dir,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Recreate new environment from output scorefile
        original_scorefile_path = os.path.join(original_output_path, original_scorefile_name)
        self.assertTrue(os.path.isfile(original_scorefile_path), msg=f"Missing original output scorefile: {original_scorefile_path}")
        with open(original_scorefile_path, "r") as f:
            original_data = [json.loads(line) for line in f]
        self.assertEqual(len(original_data), 1)
        original_record = original_data[0]
        self.assertIn("environment_manager", original_record["metadata"])
        self.assertEqual(original_record["metadata"]["environment_manager"], environment_manager)
        self.assertIn("decoy_name", original_record["metadata"])
        original_decoy_name = original_record["metadata"]["decoy_name"]

        # print("*" * 50)
        # print("Cached environment file:")
        # print(original_record["instance"]["environment"])
        # print("*" * 50)

        # Recreate environment
        reproduce_env_name = f"{original_env_name}_reproduce"
        reproduce_env_dir = os.path.join(self.workdir.name, reproduce_env_name)
        recreate_env_script = os.path.join(os.path.dirname(__file__), "recreate_envs.py")
        if environment_manager == "pixi":
            cmd = (
                # f"pixi run python -u {recreate_env_script} "
                f"{sys.executable} -u {recreate_env_script} "
                f"--env_manager '{environment_manager}' "
                f"--reproduce_env_dir '{reproduce_env_dir}' "
                f"--original_scorefile_path '{original_scorefile_path}' "
                f"--original_decoy_name {original_decoy_name}"
            )
        elif environment_manager == "uv":
            cmd = (
                # f"uv run --project {original_env_dir} python -u {recreate_env_script} "
                f"{sys.executable} -u {recreate_env_script} "
                f"--env_manager '{environment_manager}' "
                f"--reproduce_env_dir '{reproduce_env_dir}' "
                f"--original_scorefile_path '{original_scorefile_path}' "
                f"--original_decoy_name {original_decoy_name}"
            )
        elif environment_manager in ("conda", "mamba"):
            cmd = (
                f"conda run -p {original_env_dir} python -u {recreate_env_script} "
                f"--env_manager '{environment_manager}' "
                f"--reproduce_env_dir '{reproduce_env_dir}' "
                f"--original_scorefile_path '{original_scorefile_path}' "
                f"--original_decoy_name {original_decoy_name}"
            )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=None,
            # For pixi, activate the original pixi environment context
            # For conda/mamba/uv, run from environment directory for consistency with pixi workflow
            cwd=original_env_dir,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")
        self.assertTrue(
            os.path.isdir(reproduce_env_dir),
            f"Reproduced '{environment_manager}' environment directory was not created: '{reproduce_env_dir}'",
        )

        # Run reproduction simulation inside recreated environment
        reproduce_output_path = os.path.join(reproduce_env_dir, f"{environment_manager}_reproduce_outputs")
        reproduce_scorefile_name = "test_scores.json"
        # module = os.path.splitext(os.path.basename(test_script))[0]
        module = "pyrosettacluster_tests.tests.recreate_environment_test_runs"

        def find_test_root(start_path):
            path = Path(start_path).resolve()
            for parent in path.parents:
                if parent.name == "pyrosettacluster_tests":
                    return parent
            raise FileNotFoundError(f"Could not find 'pyrosetta' root directory from: '{start_path}'")

        # print("Reproduced environment directory files:", os.listdir(reproduce_env_dir))
        # if environment_manager == "uv":
        #     venv_dir = os.path.join(reproduce_env_dir, ".venv")
        #     if os.path.isdir(venv_dir):
        #         print("Reproduced environment .venv files:", os.listdir(venv_dir))
        #     bin_dir = os.path.join(reproduce_env_dir, ".venv", "bin")
        #     if os.path.isdir(bin_dir):
        #         print("Reproduced environment .venv/bin files:", os.listdir(bin_dir))
        #     pyproject_toml = os.path.join(reproduce_env_dir, "pyproject.toml")
        #     if os.path.isfile(pyproject_toml):
        #         with open(pyproject_toml, "r") as f:
        #             print("Reproduced environment pyproject_toml.toml file:", f.read())
        # elif environment_manager == "pixi":
        #     pixi_toml = os.path.join(reproduce_env_dir, "pixi.toml")
        #     if os.path.isfile(pixi_toml):
        #         with open(pixi_toml, "r") as f:
        #             print("Reproduced environment pixi.toml file:", f.read())
        #     pixi_dir = os.path.join(reproduce_env_dir, ".pixi")
        #     if os.path.isdir(pixi_dir):
        #         print("Reproduced environment .pixi files:", os.listdir(pixi_dir))
        #     pixi_envs_dir = os.path.join(reproduce_env_dir, ".pixi", "envs")
        #     if os.path.isdir(pixi_envs_dir):
        #         print("Reproduced environment .pixi/envs files:", os.listdir(pixi_envs_dir))
        #     pixi_envs_default_dir = os.path.join(reproduce_env_dir, ".pixi", "envs", "default")
        #     if os.path.isdir(pixi_envs_default_dir):
        #         print("Reproduced environment .pixi/envs/default files:", os.listdir(pixi_envs_default_dir))
        #     pixi_envs_default_bin_dir = os.path.join(reproduce_env_dir, ".pixi", "envs", "default", "bin")
        #     if os.path.isdir(pixi_envs_default_bin_dir):
        #         print("Reproduced environment .pixi/envs/default/bin files:", os.listdir(pixi_envs_default_bin_dir))

        #     def pixi_which_python(project_dir):
        #         cmd = ["pixi", "run", "which", "python"]
        #         result = subprocess.run(cmd, cwd=project_dir, capture_output=True, text=True, check=True)
        #         return result.stdout.strip()

        #     def pixi_python_sys(project_dir):
        #         code = "import sys; print(sys.executable); print(sys.path)"
        #         cmd = ["pixi", "run", "python", "-c", code]
        #         result = subprocess.run(cmd, cwd=project_dir, capture_output=True, text=True, check=True)
        #         lines = result.stdout.strip().splitlines()
        #         executable = lines[0]
        #         # Evaluate the printed list safely
        #         path_list = eval(lines[1])
        #         return executable, path_list

        #     python_executable = pixi_which_python(reproduce_env_dir)
        #     print("Reproduced environment Pixi Python executable:", python_executable)

        #     executable, sys_path = pixi_python_sys(reproduce_env_dir)
        #     print("Reproduced environment sys.executable:", executable)
        #     print("Reproduced environment sys.path:", sys_path)


        src_test_root = find_test_root(__file__)
        dst_test_root = os.path.join(reproduce_env_dir, "pyrosettacluster_tests")
        shutil.copytree(src_test_root, dst_test_root, dirs_exist_ok=True)

        if environment_manager == "pixi":
            cmd = (
                f"pixi run python -u -m {module} "
                # f"{reproduce_env_dir}/.pixi/envs/default/bin/python -u -m {module} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
            returncode = TestEnvironmentReproducibility.run_subprocess(
                cmd,
                module_dir=reproduce_env_dir,
                # For pixi, activate the recreated pixi environment context
                cwd=reproduce_env_dir,
            )
        elif environment_manager == "uv":
            cmd = (
                f"uv run --project {reproduce_env_dir} python -u -m {module} "
                # f"{reproduce_env_dir}/.venv/bin/python -u -m {module} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
            returncode = TestEnvironmentReproducibility.run_subprocess(
                cmd,
                module_dir=reproduce_env_dir,
                # For uv, activate the recreated uv environment context
                cwd=reproduce_env_dir,
            )
        elif environment_manager in ("conda", "mamba"):
            cmd = (
                f"conda run -p {reproduce_env_dir} python -u -m {module} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
            returncode = TestEnvironmentReproducibility.run_subprocess(
                cmd,
                module_dir=os.path.dirname(test_script),
                # For conda/mamba, run from recreated environment directory for consistency with pixi/uv workflows
                cwd=reproduce_env_dir,
            )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Validate reproduced decoy is identical to original decoy
        reproduce_scorefile_path = os.path.join(reproduce_output_path, reproduce_scorefile_name)
        self.assertTrue(os.path.isfile(reproduce_scorefile_path), msg=f"Missing reproduced output scorefile: {reproduce_scorefile_path}")
        with open(reproduce_scorefile_path, "r") as f:
            reproduce_data = [json.loads(line) for line in f]
        self.assertEqual(len(reproduce_data), 1)
        reproduce_record = reproduce_data[0]
        self.assertIn("environment_manager", reproduce_record["metadata"])
        self.assertEqual(reproduce_record["metadata"]["environment_manager"], environment_manager)

        self.assertEqual(
            original_record["scores"]["SEQUENCE"],
            reproduce_record["scores"]["SEQUENCE"],
        )
        self.assertEqual(
            original_record["scores"]["VALUE"],
            reproduce_record["scores"]["VALUE"],
        )
        self.assertEqual(
            original_record["scores"]["total_score"],
            reproduce_record["scores"]["total_score"],
        )
        self.assertListEqual(
            original_record["instance"]["seeds"],
            reproduce_record["instance"]["seeds"],
        )
        self.assertListEqual(
            original_record["instance"]["decoy_ids"],
            reproduce_record["instance"]["decoy_ids"],
        )
        self.assertNotEqual(
            original_record["instance"]["author"],
            reproduce_record["instance"]["author"],
        )
        self.assertNotEqual(
            original_record["metadata"]["decoy_name"],
            reproduce_record["metadata"]["decoy_name"],
        )

        original_output_file = original_record["metadata"]["output_file"]
        reproduce_output_file = reproduce_record["metadata"]["output_file"]
        assert_coordinates_script = os.path.join(os.path.dirname(__file__), "assert_coordinates.py")
        module = os.path.splitext(os.path.basename(assert_coordinates_script))[0]
        if environment_manager == "pixi":
            cmd = (
                f"pixi run python -u -m {module} "
                # f"{reproduce_env_dir}/.pixi/envs/default/bin/python -u -m {module} "
                f"--original_output_file '{original_output_file}' "
                f"--reproduce_output_file '{reproduce_output_file}' "
            )
        elif environment_manager == "uv":
            cmd = (
                f"uv run --project {reproduce_env_dir} python -u -m {module} "
                # f"{reproduce_env_dir}/.venv/bin/python -u -m {module} "
                f"--original_output_file '{original_output_file}' "
                f"--reproduce_output_file '{reproduce_output_file}' "
            )
        elif environment_manager in ("conda", "mamba"):
            cmd = (
                f"conda run -p {reproduce_env_dir} python -u -m {module} "
                f"--original_output_file '{original_output_file}' "
                f"--reproduce_output_file '{reproduce_output_file}' "
            )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=os.path.dirname(assert_coordinates_script),
            cwd=reproduce_env_dir,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")


    @unittest.skipIf(shutil.which("conda") is None, "The executable 'conda' is not available.")
    def test_recreate_environment_conda(self):
        return self.recreate_environment_test(environment_manager="conda")

    @unittest.skipIf(shutil.which("mamba") is None, "The executable 'mamba' is not available.")
    def test_recreate_environment_mamba(self):
        return self.recreate_environment_test(environment_manager="mamba")

    @unittest.skipIf(shutil.which("uv") is None, "The executable 'uv' is not available.")
    def test_recreate_environment_uv(self):
        return self.recreate_environment_test(environment_manager="uv")

    @unittest.skipIf(shutil.which("pixi") is None, "The executable 'pixi' is not available.")
    def test_recreate_environment_pixi(self):
        return self.recreate_environment_test(environment_manager="pixi")
