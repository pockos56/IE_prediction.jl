## import packages ##
using ScikitLearn
using CSV
using Statistics
using DataFrames
using PyCall
using Conda
using DataStructures

"""
Predict the ionization efficiency (in log units) using canonical SMILES, the ionization mode of the ESI source, and the pH of the mobile phase.

# Examples
```julia-repl
julia> logIE_from_SMILES("CCCOC(C)C","positive", 7)
1.1706905569667065
```
"""
function logIE_from_SMILES(SMILES::String, ESI_mode::String, pH)
    
    jblb = pyimport("joblib")
    pd = pyimport("padelpy")
    cd(@__DIR__)

    if pH > 14 || pH < 0 
        error("Set pH to a valid value between 0 and 14")
    end

    # Loading models
    if ESI_mode == "negative" || ESI_mode == "neg"
        reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg_neg.joblib"))
    elseif ESI_mode == "positive" || ESI_mode == "pos"
        if mode == None
            reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg_pos.joblib"))
        else
            reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg_$mode.joblib"))
        end
    else error("ESI_mode should be set to positive or negative")
    end
    
    fingerprint = DataFrame()
    fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(SMILES,fingerprints=true, descriptors=false)))))
    if size(fingerprint,2) != 780
        error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set to 2DAPC")
    end

    input = Matrix(hcat(pH, fingerprint))[1,:]
    IE_pred = predict(reg, input)
    return IE_pred        
end