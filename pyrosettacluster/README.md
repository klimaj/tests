# 👋 PyRosettaCluster
PyRosettaCluster is a python framework for reproducible, high-throughput job distribution of user-defined PyRosetta protocols efficiently parallelized on the user's local computer, high-performance computing (HPC) cluster, or elastic cloud computing infrastructure with available compute resources.

# 🏠 Creating Environments for PyRosettaCluster
PyRosettaCluster supports running reproducible PyRosetta simulations from reproducible virtual environments created with [Conda](https://docs.conda.io/), [Mamba](https://mamba.readthedocs.io/), [uv](https://docs.astral.sh/uv/), and [Pixi](https://pixi.sh/) environment managers. Please install [PyRosetta](https://www.pyrosetta.org/downloads) (with `cxx11thread.serialization` support) and the following packages to get started!
- `attrs`
- `billiard`
- `blosc`
- `cloudpickle`
- `cryptography`
- `dask`
- `dask-jobqueue`
- `distributed`
- `gitpython`
- `numpy`
- `pandas`
- `python-xz`
- `scipy`
- `traitlets`

[Official Full List of Packages](https://github.com/RosettaCommons/rosetta/blob/main/tests/benchmark/tests/__init__.py#L69-L84)

> [!IMPORTANT]
> It is recommended to install the required packages individually when using the `pixi`/`uv` environment managers.
>
> _Explanation:_ If using `pixi`/`uv` environment managers, it is highly recommended to avoid using the [PyPI pyrosetta-distributed](https://pypi.org/project/pyrosetta-distributed/) package to install the required `pyrosetta.distributed` framework packages (which installs subpackages using `pip install ...`, so the subpackages will _not_ get registered as installed in the exported `pixi`/`uv` environment file). Instead, please install the required `pyrosetta.distributed` framework packages individually to register them properly in the environment.

> [!IMPORTANT]
> It is recommended to use either the `pixi`, `conda` or `mamba` environment manager.
> 
> _Explanation:_ If using the `uv` environment manager, currently the only way to install PyRosetta is using the [PyPI pyrosetta-installer](https://pypi.org/project/pyrosetta-installer/) package, which installs the `pyrosetta` package using `pip install ...` (instead of `uv pip install ...`), so `pyrosetta` will _not_ get registered as installed in the exported environment file. Furthermore, when recreating a `uv` environment (see below) using the [PyPI pyrosetta-installer](https://pypi.org/project/pyrosetta-installer/) package, the syntax does not allow specifying an exact PyRosetta version (automatically defaulting to the latest published PyRosetta version), and therefore installing the correct PyRosetta version in the recreated `uv` environment depends fortuitously on a prompt `uv` environment recreation after the original `uv` environment creation (typically within the same week). For this reason, `pixi`, `conda` and `mamba` are highly preferred environment managers instead of `uv` for reproducible environments used for reproducible PyRosettaCluster simulations. Note that the PyRosetta developers are working on supporting a bona fide `uv pip install pyrosetta==X.Y` syntax in the future — stay tuned!

# 💨 Running PyRosettaCluster simulations
Please see the [official PyRosettaCluster documentation](https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.distributed.cluster.html#pyrosetta.distributed.cluster.PyRosettaCluster) for all of the configuration options.
> [!CAUTION]
> PyRosettaCluster uses the [cloudpickle module](https://github.com/cloudpipe/cloudpickle), which can lead to arbitrary code execution resulting in a critical security vulnerability. **Please only run (1) with user-provided PyRosetta protocols from trusted sources, and (2) over a secure network (unless running on a local network).**
> 
> A primary feature of PyRosettaCluster is that arbitrary user-provided PyRosetta protocols are pickled, sent over a network, and unpickled, which allows the user to run customized macromolecular design and modeling workflows. Unless running solely on a local network, it is _highly_ recommended to operate PyRosettaCluster behind a trusted private network segment (i.e., a firewall) or setup TLS/SSL communication between network endpoints for authenticated and encrypted transmission of data. In PyRosettaCluster, there are two easy ways to setup TLS/SSL communication:
>
> **(1)** Use the `PyRosettaCluster(security=True)` option to invoke the `cryptography` package to automatically generate a `Security.temporary()` object on-the-fly through the `dask` or `dask-jobqueue` package.
>
> **(2)** Pre-generate a `dask.distributed.Security()` object using _OpenSSL_ (instead of the `cryptography` package) and pass it to PyRosettaCluster. See the `pyrosetta.distributed.cluster.generate_dask_tls_security()` function docstring for more information:
> ```
> security = pyrosetta.distributed.cluster.generate_dask_tls_security()
> PyRosettaCluster(security=security)
> ```

> [!IMPORTANT]
> PyRosettaCluster is a tool for reproducible computational macromolecular modeling and design. It is up to the user to define their PyRosetta protocols with reproducibility in mind – meaning user-provided PyRosetta protocols ought to be _deterministic_:
> 
> **(1)** Set seeds for any _impure functions_ (i.e., non-deterministic functions) implemented.
>  - Pseudo-random seeds can either _(i)_ be hard-coded into PyRosetta protocols, _(ii)_ distributed as `kwargs` key values with each submitted task, or _(iii)_ be dynamically set based on each PyRosetta protocol's automatically-assigned random seed accessible through each PyRosetta protocol's `kwargs["PyRosettaCluster_seed"]` key value.
> 
> **(2)** If impure functions cannot be made pure through controlling the underlying randomness, please do _not_ rely on them to update the `Pose` in meaningful ways.
>  - i.e., a randomly-named score key might be alright, but not randomly selecting the number of `Pose` conformational updates.
>
> **(3)** If implementing third-party software in PyRosetta protocols that does not support determinism, please ask the developers to support determinism in their software.
>
> In general, the determinism can be (and ought to be) strictly controlled when developing PyRosetta protocols. Note that PyRosettaCluster can still be used as a job distributor, even if PyRosetta protocol determinism is impossible for a specific application.

# 🏘️ *Re*creating Environments for PyRosettaCluster

## 1️⃣ ✂ Extract environment file ️
The virtual environment configuration used for the original simulation is cached in the PyRosettaCluster output decoy file and in the _full-record_ output scorefile. 
Please refer to the following table to select _one_ environment file extraction method based on the file type being used to recreate the original virtual environment:

| File type extension | Output type | Extraction method #1<br>(_without_ PyRosetta) | Extraction method #2<br>(_requires_ PyRosetta) |
| --- | --- | --- | --- |
| `.pdb` | Decoy | Read file → Copy → Paste into new file | Run `dump_env_file.py` helper |
| `.pdb.bz2` | Decoy | Unzip with `bzip2` → Read file → Copy → Paste into new file | Run `dump_env_file.py` helper |
| `.pkl_pose`, `.pkl_pose.bz2`, `.b64_pose`, `.b64_pose.bz2` | Decoy | | Run `dump_env_file.py` helper |
| `.json` | Full-record scorefile | Read file → Copy → Paste into new file | Run `dump_env_file.py` helper |
| Pickled `pandas.DataFrame`<br>(`.gz`, `.xz`, `.tar`, etc.) | Full-record scorefile | | Run `dump_env_file.py` helper |
| `.init`, `.init.bz2` | PyRosetta initialization file | | Run `dump_env_file.py` helper |

> [!NOTE]  
> **Extraction method #1:** If copy/pasting into a new file, the environment file string is located in the `record["instance"]["environment"]` nested key value of the PyRosettaCluster full record. Please paste it into one of the following file names (as expected in the next step) in a new folder, depending on the environment manager you're using to recreate the environment:
> | Environment manager | New file name |
> | --- | --- |
> | `pixi` | `pixi.lock` |
> | `uv` | `requirements.txt` |
> | `conda` | `environment.yml` |
> | `mamba` | `environment.yml` |
> 
> If using `pixi`/`uv` environment managers, please also extract the manifest file string (`pixi`) or project file string (`uv`) located in the `record["metadata"]["toml"]` nested key value of the PyRosettaCluster full record. The `record["metadata"]["toml_format"]` nested key value also specifies the TOML file format. Please paste it into one of the following file names (as expected in the next step) in the same new folder, depending on the environment manager you're using to recreate the environment:
> | Environment manager | New file name |
> | --- | --- |
> | `pixi` | `pixi.toml` / `pypyroject.toml` |
> | `uv` | `pyproject.toml` |
>
> Also note the `record["instance"]["sha1"]` nested key value holding the GitHub SHA1 required to reproduce the PyRosettaCluster simulation!

> [!NOTE]  
> **Extraction method #2:** If running `dump_env_file.py`, the `pyrosetta` package must be installed in any existing virtual environment, and that virtual environment's python interpreter used to run the script.
>
> Also note the printed GitHub SHA1 required to reproduce the PyRosettaCluster simulation!

> [!TIP]
> **Extraction method #2:** See `python dump_env_file.py --help` for details.

## 2️⃣ 🛠️ Recreate environment
Run `python recreate_env.py` to recreate the virtual environment.

> [!CAUTION]
> This script runs a subprocess with one of the following commands:<br>
> - `conda env create ...`: when using the `conda` environment manager<br>
> - `mamba env create ...`: when using the `mamba` environment manager<br>
> - `uv pip sync ...`: when using the `uv` environment manager<br>
> - `pixi install ...`: when using the `pixi` environment manager<br>
> Installing certain packages may not be secure, so please only run with an input environment file you trust!<br>
> 👉 Learn more about [PyPI security](https://pypi.org/security) and [conda security](https://www.anaconda.com/docs/reference/security).

> [!IMPORTANT]
> If using `pixi`/`uv` environment managers, please use the _system python interpreter_, since the script creates a new `pixi`/`uv` environment and cannot be run from an existing virtual environment.
> If using `conda`/`mamba`, any python interpreter may be used.

> [!NOTE]
> If using the `uv` environment manager, the PyRosetta installation step may be subsequently required if using the `pyrosetta-installer` package installation method.  Note that installing the identical PyRosetta version of the original `uv` environment in the recreated `uv` environment depends fortuitously on a prompt `uv` environment recreation after the original `uv` environment creation (typically within the same week). See the [PyPI pyrosetta-installer](https://pypi.org/project/pyrosetta-installer/) for details.

> [!TIP]
> See `python recreate_env.py --help` for details.

# 🚀 Reproducing PyRosettaCluster simulations
In order to re-run the same user-provided PyRosetta protocols, clone the original GitHub repository and checkout the original GitHub SHA1 used by the original PyRosettaCluster simulation. You'll need to know the owner and repository name (and if not, don't worry, there are ways to search GitHub by SHA1):
```
git clone --no-checkout https://github.com/<owner>/<repo>.git
cd <repo>
git fetch origin <SHA1>
git checkout <SHA1>
```
Then, use the python interpreter of the recreated environment to run your PyRosettaCluster simulation reproduction script. Here's a template script to get started!
```
from pyrosetta.distributed.cluster import reproduce

# Import (or copy/paste) the original user-provided PyRosetta protocols
from my_protocols import protocol_1, protocol_2  # Change depending on the original GitHub repository structure

def main():
    reproduce(
        # Input either a PyRosettaCluster output decoy file or output PyRosetta initialization file
        input_file=...,

        # Or input a PyRosettaCluster scorefile and decoy name
        scorefile=...,
        decoy_name=...,

        # Optional configurations:
        protocols=[protocol_1, protocol_2], # Can be `None` for auto-detection
        clients=...,
        input_packed_pose=...,
        instance_kwargs={
            "output_path": ...,
            "scratch_dir": ...,
            "project_name": ...,
            "simulation_name": ...,
        },
        clients_indices=...,
        resources=...,
        skip_corrections=...,
        init_from_file_kwargs=...,
    )

if __name__ == "__main__":
    main()
```
✅ Save your PyRosettaCluster simulation reproduction script, and run it with the _recreated environment's python interpreter_ (with the local repository `HEAD` at that same commit SHA1 for PyRosettaCluster SHA1 validation). The PyRosetta build string and the environment file string will also be validated against the original record at this step.

🎉 Congrats! You have now recreated a virtual environment and used it to successfully reproduce a distributed PyRosettaCluster simulation.




