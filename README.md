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

### Documentation (Basics)

#### Ionization efficiency prediction for compounds with known structures


```julia
using IE_prediction
# Using canonical SMILES and the pH of the mobile phase
logIE_from_structure("CC(=O)OC1=CC=CC=C1C(=O)O", 2.7)
> 1.2200884862831083

# Using InChIKey and pH
logIE_from_structure("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", 2.7)
> 1.2200884862831083

# Predicting the minimum or maximum value of ionization efficiency
logIE_from_structure("CCCOC(C)C", 7, "min")
> 1.5592564226940018

# For multiple compounds, using the batch mode
logIE_from_structure(["CCCOC(C)C", "CC(=O)OC1=CC=CC=C1C(=O)O"], [10, 10], "min", FP_calculation="batch")
> 2-element Vector{Float64}:
 1.7958649829457431       
 1.2213394324250775
```

#### Ionization efficiency prediction for compounds with unknown structures

```julia
# e.g. for an unknown compound with a m/z of 240.25 measured at pH 9

# In case of multiple fragments with m/z 222.24, 210.24 and 179.9
using IE_prediction
logIE_from_CNLs([222.24, 210.24, 179.9], 240.25, 9, "mean")
> 3.524817834899175

# In case of a single fragment with m/z 210.24
logIE_from_CNLs([210.24], 240.25, 9, "max")
> 3.977743567907343

```

### Documentation (Advanced)

``` julia
logIE_from_structure(identifier::Union{String, Vector{String}, DataFrame}, pH, data_mode::String="mean"; FP_calculation::String="default", identifier_type::String="auto")
```

#### Parameters
1. data_mode
- Defines the model used for the IE prediction.
- Possible values: "min", "mean", "max"
- Default value: "mean"

2. FP_calculation
- Defines the mode in which the input data is generated.  
*default* - Suitable for a single identifier (Expected type for identifier variable: String)  
*batch* - Suitable for multiple identifiers (Expected type for identifier variable: Vector{String})  
*pre-calculated* - Suitable in cases when the PubChem fingerprints have been already calculated for the compounds of interest. This bypasses completely the calculation of fingerprints and the prediction model is applied directly. The **PubChem fingerprints** can be calculated with [padelpy](https://github.com/Hamada-Noreldeen/PaDELPy) (Expected type for identifier variable: DataFrame)  
- Possible values: "default", "batch", "pre-calculated"
- Default value: "default"

3. identifier_type
- Defines which type of identifier is provided. The function automatically detects the type of fingerprint by default, but there might be cases when manual selection is useful.
- It doesn't affect the process when in pre-calculated fingerprints.
- Possible values: "auto", "SMILES", "INCHIKEY"
- Default value: "auto"