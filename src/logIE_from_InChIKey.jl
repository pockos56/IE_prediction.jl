## import packages ##
using ScikitLearn, CSV, Statistics, DataFrames, DataStructures, PyCall, Conda

"""
Predict the ionization efficiency (in log units) using InChIKey and the pH of the mobile phase.

# Examples
```julia-repl
julia> logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 7, "mean")
1.56133081771261
```
"""
function logIE_from_InChIKey(INCHIKEY::Union{String, Vector{String}}, pH::Float64, data_mode::String)
    # Packages
    jblb = pyimport("joblib")
    pcp = pyimport("pubchempy")
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
    if typeof(INCHIKEY) == Vector{String}
        error("Functionality to be implemented. Currently only single string values are accepted")
    end

    fingerprint = DataFrame()
    cids = pcp.get_compounds(INCHIKEY, "inchikey")
    if isempty(cids)
        error("CID not found: The INCHIKEY ($INCHIKEY) could not be associated with a PubChem compound.")
    end
    shortest_cid = cids[argmin([cids[y].cid for y in 1:length(cids)])]
    fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(shortest_cid.isomeric_smiles,fingerprints=true, descriptors=false)))))
    if size(fingerprint,2) != 881
        error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set PubChem to true")
    end

    input = Matrix(hcat(pH, fingerprint))[1,:]
    IE_pred = predict(reg, input)
    return IE_pred        
end