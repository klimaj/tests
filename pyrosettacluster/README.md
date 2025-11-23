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
> _Explanation:_ If using the `uv` environment manager, currently the only way to install PyRosetta is using the [PyPI pyrosetta-installer](https://pypi.org/project/pyrosetta-installer/) package, which installs the `pyrosetta` package using `pip install ...` (instead of `uv pip install ...`), so `pyrosetta` will _not_ get registered as installed in the exported environment file. Furthermore, when recreating a `uv` environment (see below) using the [PyPI pyrosetta-installer](https://pypi.org/project/pyrosetta-installer/) package, the syntax does not allow specifying an exact PyRosetta version (automatically defaulting to the latest published PyRosetta version), and therefore installing the correct PyRosetta version in the recreated `uv` environment depends fortuitously on a prompt `uv` environment recreation after the original `uv` environment creation (typically within the same week). For this reason, `pixi`, `conda` and `mamba` are highly preferred environment managers instead of `uv` for reproducible environments used for reproducible PyRosettaCluster simulations. Note that the PyRosetta developers are working on supporting a bona fide `uv pip install pyrosetta=X.Y` syntax in the future — stay tuned!

# 🏘️ *Re*creating Environments for PyRosettaCluster

## ✂ 1️⃣ Extract environment file ️
The virtual environment configuration used for the original simulation is cached in the PyRosettaCluster output decoy file and in the _full-record_ output scorefile. 
Please refer to the following table to select _one_ environment file extraction method based on the file type being used to recreate the original virtual environment:

| File type extension | Output type | Extraction method #1<br>(_without_ PyRosetta) | Extraction method #2<br>(_requires_ PyRosetta) |
| --- | --- | --- | --- |
| `.pdb` | Decoy | Read file → Copy → Paste into new file | Run `dump_env_file.py` helper |
| `.pdb.bz2` | Decoy | Unzip with `bzip2` → Read file → Copy → Paste into new file | Run `dump_env_file.py` helper |
| `.pkl_pose`, `.pkl_pose.bz2`, `.b64_pose`, `.b64_pose.bz2` | Decoy | | Run `dump_env_file.py` helper |
| `.json` | Full-record scorefile | Read file → Copy → Paste into new file | Run `dump_env_file.py` helper |
| Pickled `pandas.DataFrame`<br>(`.gz`, `.xz`, `.tar`, etc.) | Full-record scorefile | | Run `dump_env_file.py` helper |

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

> [!NOTE]  
> **Extraction method #2:** If running `dump_env_file.py`, the `pyrosetta` package must be installed in any existing virtual environment, and that virtual environment's python interpreter used to run the script.

> [!TIP]
> **Extraction method #2:** See `python dump_env_file.py --help` for details.

## 🛠️ 2️⃣ Recreate environment
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

## 🚀 3️⃣ Reproduce PyRosettaCluster simulation!
Use the python interpreter of the recreated environment to run your PyRosettaCluster simulation reproduction script.










