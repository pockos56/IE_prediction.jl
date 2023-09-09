# IE_prediction

[![Build Status](https://github.com/pockos56/IE_prediction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pockos56/IE_prediction.jl/actions/workflows/CI.yml?query=branch%3Amain)


### Documentation

#### Ionization efficiency prediction for compounds with known structures


```julia

# For Isopropyl propyl ether measured in negative ionization mode and at pH 7

# Using the canonical SMILES, the ionization mode and the pH of the mobile phase
logIE_1 = logIE_from_SMILES("CCCOC(C)C","positive", 7)

# Using the InChIKey instead
logIE_2 = logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N","positive", 7)
logIE_1 == logIE_2      # True

```

#### Ionization efficiency prediction for compounds with unknown structures

```julia
# For an unknown compound with a m/z of 240.25 measured in positive ionization mode and at pH 9

# In case of a multiple fragments with m/z 222.24, 210.24 and 179.9
logIE_from_CNLs([222.24, 210.24, 179.9], 240.25, "positive", 9)

# In case of a single fragment with m/z 210.24
logIE_from_CNLs([210.24], 240.25, "positive", 9)

```