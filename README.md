# IE_prediction

This julia package enables the prediction of ionization efficiency for compounds with known and unknown structures.

This allows semi-quantifation in electrospray ionization mass spectrometry (ESI-MS) for identified and unidentified compounds without relying on reference standards.

### Installation
```julia

# Install dependencies
using PyCall
using Conda
Conda.add("joblib")         # Package to load Python models
Conda.add("padelpy")        # Package to calculate molecular fingerprints
Conda.add("pubchempy")      # Package to calculate canonical SMILES from InCHiKey
Conda.add("catboost")       # Package to run CatBoost regression models for IE IE_prediction

using Pkg
Pkg.build("PyCall")         # Build PyCall with the installed packages

# --Restart julia--

using Pkg
Pkg.add(url="https://github.com/pockos56/IE_prediction.jl")

```

### Documentation

#### Ionization efficiency prediction for compounds with known structures


```julia
# e.g. for isopropyl propyl ether measured in positive ionization mode and at pH 7

# Using the canonical SMILES, the ionization mode and the pH of the mobile phase
using IE_prediction
logIE_1 = logIE_from_SMILES("CCCOC(C)C","positive", 7)

# Using the InChIKey instead
logIE_2 = logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N","positive", 7)
logIE_1 == logIE_2      # True

```

#### Ionization efficiency prediction for compounds with unknown structures

```julia
# e.g. for an unknown compound with a m/z of 240.25 measured in negative ionization mode and at pH 9

# In case of a multiple fragments with m/z 222.24, 210.24 and 179.9
using IE_prediction
logIE_from_CNLs([222.24, 210.24, 179.9], 240.25, "negative", 5)

# In case of a single fragment with m/z 210.24
logIE_from_CNLs([210.24], 240.25, "negative", 5)

```