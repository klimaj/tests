# 🏘️ Recreating Environments for PyRosettaCluster

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
> If running `dump_env_file.py`, the `pyrosetta` package must be installed in any existing virtual environment, and that virtual environment's python interpreter used to run the script.

> [!TIP]
> See `python dump_env_file.py --help` for details.

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

> [!NOTE]
> If using `pixi`/`uv` environment managers, please use the _system python interpreter_, since the script creates a new `pixi`/`uv` environment and cannot be run from an existing virtual environment.
> If using `conda`/`mamba`, any python interpreter may be used.

> [!NOTE]
> If using the `uv` environment manager, the PyRosetta installation step may be subsequently required if using the `pyrosetta-installer` package installation method. See the [PyPI pyrosetta-installer](https://pypi.org/project/pyrosetta-installer/) for details.

> [!TIP]
> See `python recreate_env.py --help` for details.

## 🚀 3️⃣ Reproduce PyRosettaCluster simulation!
Use the python interpreter of the recreated environment to run your PyRosettaCluster simulation reproduction script.





