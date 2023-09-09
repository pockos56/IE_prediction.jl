# IE_prediction

This julia package enables the prediction of ionization efficiency for compounds with known and unknown structures.

This allows semi-quantifation in electrospray ionization mass spectrometry (ESI-MS) for identified and unidentified compounds without relying on reference standards.


### Documentation

#### Ionization efficiency prediction for compounds with known structures


```julia
# e.g. for isopropyl propyl ether measured in positive ionization mode and at pH 7

# Using the canonical SMILES, the ionization mode and the pH of the mobile phase
logIE_1 = logIE_from_SMILES("CCCOC(C)C","positive", 7)

# Using the InChIKey instead
logIE_2 = logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N","positive", 7)
logIE_1 == logIE_2      # True

```

#### Ionization efficiency prediction for compounds with unknown structures

```julia
# e.g. for an unknown compound with a m/z of 240.25 measured in negative ionization mode and at pH 9

# In case of a multiple fragments with m/z 222.24, 210.24 and 179.9
logIE_from_CNLs([222.24, 210.24, 179.9], 240.25, "negative", 5)

# In case of a single fragment with m/z 210.24
logIE_from_CNLs([210.24], 240.25, "negative", 5)

```