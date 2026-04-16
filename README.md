# sndMuonFlux

Codebase for the SND@LHC muon flux measurement analysis. This repository contains the code for track counting, tracking efficiency estimation, and luminosity extraction, combining C++ extensions with Python analysis scripts.

## Project Structure

- `muonFluxUtils/`: C++ utilities for muon flux calculations.
- `nTracks/`: C++ implementation for counting tracks in reconstructed data.
- `trkeff/`: C++ and Python implementations for tracking efficiency estimation.
- `pythonHelpers/`: Shared Python modules for data processing, ROOT interaction, and library loading.
- `notebooks/`: Jupyter notebooks for exploratory data analysis.
- `scripts/`: Main execution scripts for the analysis pipeline.

## Prerequisites

- **SNDSW**: The SND@LHC software envronment.

Required environment variables:

- `SNDSW_ROOT`
- `ROOTSYS`

## Installation

1. **Clone the repository**:

    ```bash
    git clone https://github.com/idionisov/sndMuonFlux.git
    cd sndMuonFlux
    ```

2. **Environment Setup**:
   `alienv enter sndsw/latest`

3. **Build the C++ extensions**:

    ```bash
    mkdir build && cd build
    cmake ..
    make -j$(nproc)
    cd ..
    ```

## Usage

The script `getMuonFlux.py` containts the entire pipeline for muon flux calculation.

Independent control of track counting and efficiency estimation is provided by the scripts `getNtracks.py` and `getTrkeff.py`. Precalculated values for track counts and efficiency can be given to `getMuonFlux.py` with the corresponding flags.

#### Notes:

- Assumed implemented methods `sndRecoTrack::getPlaneAtZ(float Z)` and `sndRecoTrack::getDoca(const MuFilterHit* mfHit) const` [https://github.com/idionisov/sndsw/blob/354bee7d9d7f8490a6b4f0124d3b5a95cac430e8/shipLHC/sndRecoTrack.cxx#L270-L304]
- Bash expands the wildcard "\*". Input is read by TChain, so make sure to escape it with "\\\*".
- Output is either appended to an existing file, or a new file is created.
    - Only implemented for data and not Monte Carlo samples.
    - MC results printed to stdout instead in all cases.

### 1. Track Counting (minimal options)

```bash
python getNtracks.py -i /path/to/recoData/\*.root
```

See `python getNtracks.py --help` for all options.

### 2. Tracking Efficiency (minimal options)

```bash
python getTrkeff.py -i /path/to/recoData/\*.root -g /path/to/geofile.root
```

For both MC-Truth and tagging track methods.
See `python getTrkeff.py --help` for all options.

### 3. Muon Flux Calculation (minimal options)

```bash
python getMuonFlux.py -i /path/to/recoData/\*.root -g /path/to/geofile.root
```

### 4. Bunch structure extraction

```bash
python extractBunchStruct.py -i /path/to/recoData/\*.root -o /path/to/output.root
```

### 5. Luminosity extraction (1/nb)

```bash
python getLumi.py /path/to/inputData
```
