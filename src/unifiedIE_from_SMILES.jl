## import packages ##
using ScikitLearn, CSV, Statistics, DataFrames, DataStructures, PyCall, Conda

"""
Predict the ionization efficiency (in log units) using canonical SMILES and the pH of the mobile phase.

# Examples
```julia-repl
julia> logIE_from_SMILES("CCCOC(C)C", 7)
1.43249846779072
```
"""
function ulogIE_from_SMILES(SMILES::String, pH)
    
    jblb = pyimport("joblib")
    pd = pyimport("padelpy")
    cd(@__DIR__)

    # Loading models
    reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg.joblib"))

    if pH > 14 || pH < 0 
        error("Set pH to a valid value between 0 and 14")
    end

    fingerprint = DataFrame()
    fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(SMILES,fingerprints=true, descriptors=false)))))
    if size(fingerprint,2) != 881
        error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set PubChem to true")
    end

    input = Matrix(hcat(pH, fingerprint))[1,:]
    IE_pred = predict(reg, input)
    return IE_pred        
end