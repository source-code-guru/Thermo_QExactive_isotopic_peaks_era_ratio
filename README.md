# FTMS RAW Isotope Ratio Extractor (Q Exactive Orbitrap)

A Python script to extract **MS1 peak intensities** and compute the **mono / first isotopic peak integral ratio** from **Thermo Fisher Q Exactive (Orbitrap) LC–MS(/MS) `.raw` files** using the **ThermoFisher CommonCore RawFileReader (.NET)** libraries.

---

## What this does

Given a folder of Thermo `.raw` files and a CSV list of target ions (m/z and optional RT), the script:

- Loads Thermo `.raw` files via **ThermoFisher.CommonCore.RawFileReader** (through **pythonnet / CLR**)
- Extracts centroid **MS1** spectra across the run
- For each target ion:
  - Finds the best matching **monoisotopic peak** within a ppm window and (optionally) within a retention time window
  - Searches for the **first isotope** (mass + 1.003355)
  - Computes a **mono / first isotope integral ratio** using interpolation + trapezoidal integration
- Writes results to CSV output files:
  - `FTMS_search_output_trapz.csv`
  - `FTMS_int_search_output_trapz.csv`

---

## Tested environments

- **Ubuntu 24.04 LTS** with **conda 24.11.3** and **Python 3.x**
- **Windows (10/11)** with **conda 24.11.3** and **Python 3.x**

> The script requires the **Thermo Fisher CommonCore / RawFileReader DLLs** to be available on the machine.

---

## Requirements

### Core requirements
- Python **3.x**
- Conda (recommended) or pip
- **pythonnet** (>= 2.5)
- **.NET SDK** (>= 8.x)
- Thermo Fisher **RawFileReader / CommonCore DLLs**:
  - `ThermoFisher.CommonCore.Data`
  - `ThermoFisher.CommonCore.RawFileReader`
  - `ThermoFisher.CommonCore.BackgroundSubtraction`
  - `ThermoFisher.CommonCore.MassPrecisionEstimator`

### Python dependencies
- numpy
- pandas
- scipy

---

## Installation

### 1) Clone the repository
```bash
git clone Thermo_QExactive_isotopic_peaks_era_ratio.git
cd Thermo_QExactive_isotopic_peaks_era_ratio
```

### 2) Create and activate a conda environment
**Ubuntu / Linux / macOS (bash/zsh)**
```bash
conda create -n ftms-raw python=3.11 -y
conda activate ftms-raw
```

**Windows (Anaconda Prompt / PowerShell)**
```bash
conda create -n ftms-raw python=3.11 -y
conda activate ftms-raw
```

### 3) Install Python dependencies
- Recommended (conda-forge)
```bash
conda install -c conda-forge numpy pandas scipy pythonnet -y
```

Or using pip
```bash
python -m pip install numpy pandas scipy pythonnet
```

- Install .NET SDK (>= 8.x)
**Ubuntu 24.04 LTS**
```bash
sudo apt update
sudo apt install -y dotnet-sdk-8.0
dotnet --list-sdks
```

**Windows**

Install a .NET SDK 8+ from Microsoft, then verify:

```bash
dotnet --list-sdks
```

**Thermo Fisher CommonCore / RawFileReader DLLs**

You must provide the Thermo Fisher RawFileReader / CommonCore assemblies (DLLs). If clr.AddReference(...) fails, it usually means the DLLs cannot be found.

**Common approaches:**

Copy the DLLs next to the script (same folder as the .py file), or make sure the DLL directory is discoverable by the runtime / pythonnet.
If you encounter errors such as:

Could not load file or assembly 'ThermoFisher.CommonCore.RawFileReader'
then verify:
- the DLLs are present on the machine
- the process can locate them at runtime

## Input
### Target list CSV format

The script expects a CSV file with ; as separator, containing at least:
- Mass (float)
- optionally RT (float, minutes)
- Example:

```csv
Mass;RT
445.1203;5.20
300.2001;8.75
```

## Configuration

Edit these variables at the top of the script:

```python
delta_mass_ppm = 10.0
detection_threshold = 0
delta_rt = 10.0

csv_file_mass = "..."
working_dir = "..."
```

### Path examples
**Windows**
```python
csv_file_mass = r"C:\Users\Marc\Documents\Test\list_mass_test.csv"
working_dir   = r"C:\Users\Marc\Documents\Test\"
```

**Linux (Ubuntu)** 
```python
csv_file_mass = "/home/marc/Documents/Test/list_mass_test.csv"
working_dir   = "/home/marc/Documents/Test/"
```

### Running the script 
**Ubuntu / Linux** 
```bash
conda activate ftms-raw
python Thermo_QExactive_isotopic_peaks_era_ratio.py
```

**Windows (Anaconda Prompt / PowerShell)**
```bash
conda activate ftms-raw
python Thermo_QExactive_isotopic_peaks_era_ratio.py
```

## Output

Two CSV files are written in working_dir:

### 1) Per-ion / per-file results

File name : FTMS_search_output_trapz.csv

Includes:
- file name
- mass
- charge
- scan number
- retention time
- monoisotopic intensity
- first isotope intensity
- first isotope charge
- first/mono isotopic integral ratio
- basepeak

### 2) Intensity matrix across RAW files

File name : FTMS_int_search_output_trapz.csv

Columns:
- Mass
- one column per RAW file (intensity values)

## Troubleshooting
- ModuleNotFoundError: No module named 'clr'

Install pythonnet:

```bash
conda install -c conda-forge pythonnet -y
```
or
```bash
python -m pip install pythonnet
```

- No such file or directory: 'dotnet'

Install .NET SDK 8+ and verify:

```bash
dotnet --list-sdks
```

- ThermoFisher DLL load errors

If you see errors like:
Could not load file or assembly 'ThermoFisher.CommonCore.RawFileReader'
then:
- ensure the Thermo Fisher DLLs are present
- ensure they are in a location discoverable by the runtime

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).
You are free to:
- use
- study
- modify
- redistribute

this software under the terms of the GPL-3.0 license.

Any derivative work must also be distributed under GPL-3.0.

See the LICENSE file for the full license text, or
https://www.gnu.org/licenses/gpl-3.0.html

## Contributing

Issues and pull requests are welcome.
For bug reports, please include:
- operating system and version
- conda environment details (conda list)
- Python version
- full traceback / log output
