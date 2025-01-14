# IE_prediction

This julia package enables the prediction of ionization efficiency for compounds with known and unknown structures.

This allows semi-quantifation in electrospray ionization mass spectrometry (ESI-MS) for identified and unidentified compounds without relying on reference standards.

### Installation
```julia

# Install dependencies
using Pkg
Pkg.add("PyCall")
Pkg.add("Conda")

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
using IE_prediction
# Using the canonical SMILES and the pH of the mobile phase
julia> logIE_from_structure("CC(=O)OC1=CC=CC=C1C(=O)O", 2.7)
1.2200884862831083

# Using the InChIKey and pH
julia> logIE_from_structure("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", 2.7)
1.2200884862831083

# Predicting the minimum or maximum value of ionization efficiency
julia> logIE_from_structure("CCCOC(C)C", 7, "min")
1.5592564226940018

# For multiple compounds, using the batch mode
julia> logIE_from_structure(["CCCOC(C)C", "CC(=O)OC1=CC=CC=C1C(=O)O"], [10, 10], "min", FP_calculation="batch")
2-element Vector{Float64}:
 1.7958649829457431       
 1.2213394324250775
```

#### Ionization efficiency prediction for compounds with unknown structures

```julia
# e.g. for an unknown compound with a m/z of 240.25 measured at pH 9

# In case of multiple fragments with m/z 222.24, 210.24 and 179.9
using IE_prediction
logIE_from_CNLs([222.24, 210.24, 179.9], 240.25, 9, "mean")

# In case of a single fragment with m/z 210.24
logIE_from_CNLs([210.24], 240.25, 9, "max")

```