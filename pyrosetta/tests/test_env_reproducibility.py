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
    def run_subprocess(cmd, module_dir=None, cwd=None):
        print("Running command:", cmd)
        if module_dir:
            env = os.environ.copy()
            pythonpath = os.environ.get("PYTHONPATH")
            env["PYTHONPATH"] = f"{module_dir}{os.pathsep + pythonpath if pythonpath else ''}"
        else:
            env = None
        try:
            # Use live output streaming for GitHub Actions visibility
            process = subprocess.Popen(
                shlex.split(cmd),
                cwd=cwd,
                env=env,
                shell=False,
                stdout=sys.stdout,
                stderr=sys.stderr,
                text=True,
            )
            returncode = process.wait()
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
        cmd = "{0} -m {1} --env_manager '{2}' --env_dir '{3}'".format(
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

        print("*" * 50)
        print("Cached environment file:")
        print(original_record["instance"]["environment"])
        print("*" * 50)

        # Recreate environment
        reproduce_env_name = f"{original_env_name}_reproduce"
        reproduce_env_dir = os.path.join(self.workdir.name, reproduce_env_name)
        recreate_env_script = os.path.join(os.path.dirname(__file__), "recreate_envs.py")
        if environment_manager == "pixi":
            cmd = (
                f"pixi run python -u {recreate_env_script} "
                f"--env_manager '{environment_manager}' "
                f"--reproduce_env_dir '{reproduce_env_dir}' "
                f"--original_scorefile_path '{original_scorefile_path}' "
                f"--original_decoy_name {original_decoy_name}"
            )
        elif environment_manager == "uv":
            cmd = (
                f"uv run --project {original_env_dir} python -u {recreate_env_script} "
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
        module = os.path.splitext(os.path.basename(test_script))[0]
        if environment_manager == "pixi":
            cmd = (
                f"pixi run python -u -m {module} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
        elif environment_manager == "uv":
            cmd = (
                f"uv run --project {reproduce_env_dir} --active python -u -m {module} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
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
            # For pixi, activate the recreated pixi environment context
            # For conda/mamba/uv, run from recreated environment directory for consistency with pixi workflow
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
                f"--original_output_file '{original_output_file}' "
                f"--reproduce_output_file '{reproduce_output_file}' "
            )
        elif environment_manager == "uv":
            cmd = (
                f"uv run --project {reproduce_env_dir} python -u -m {module} "
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
