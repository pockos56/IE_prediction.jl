## import packages ##
using ScikitLearn, CSV, Statistics, DataFrames, DataStructures, PyCall, Conda

"""
Predict the ionization efficiency (in log units) using canonical SMILES and the pH of the mobile phase.

# Examples
```julia-repl
julia> logIE_from_SMILES("CCCOC(C)C", 7)
1.56133081771261
```
"""
function logIE_from_SMILES(SMILES::Union{String, Vector{String}}, pH, data_mode::String)
    # Packages
    jblb = pyimport("joblib")
    pd = pyimport("padelpy")
    cd(@__DIR__)

    # Loading models
    if data_mode != "min" && data_mode != "mean" && data_mode != "max"
        error("Set data_mode to min, mean, or max")
    else
        reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg_$data_mode.joblib"))
    end

    # Test pH
    if pH > 14 || pH < 0 
        error("Set pH to a valid value between 0 and 14")
    end

    if typeof(SMILES) == String
        fingerprint = DataFrame()
        fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(SMILES,fingerprints=true, descriptors=false)))))
        if size(fingerprint,2) != 881
            error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set PubChem to true")
        end

        input = Matrix(hcat(pH, fingerprint))[1,:]
        IE_pred = predict(reg, input)
        return IE_pred        
    end
    if  typeof(SMILES) == Vector{String}
        error("Functionality to be implemented. Currently only single string values are accepted")
        fingerprint = DataFrame()
        fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(SMILES,fingerprints=true, descriptors=false)))))
        if size(fingerprint,2) != 881
            error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set PubChem to true")
        end

        input = Matrix(hcat(pH, fingerprint))[1,:]
        IE_pred = predict(reg, input)
        println("This function is deprecated. Please use the function logIE_from_structure instead for new features, including batch calculations and logIE predictions for pre-calculated PubChem fingerprints")
        return IE_pred        
    end 
end