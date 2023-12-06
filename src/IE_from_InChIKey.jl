## import packages ##
using ScikitLearn
using CSV
using Statistics
using DataFrames
using PyCall
using Conda
using DataStructures

"""
Predict the ionization efficiency (in log units) using InChIKey, the ionization mode of the ESI source, and the pH of the mobile phase.

# Examples
```julia-repl
julia> logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N","positive", 7)
1.1706905569667065
```
"""
function logIE_from_InChIKey(INCHIKEY::String, ESI_mode::String, pH; mode::String=None)

    jblb = pyimport("joblib")
    pcp = pyimport("pubchempy")
    pd = pyimport("padelpy")
    cd(@__DIR__)

    if pH > 14 || pH < 0 
        error("Set pH to a valid value between 0 and 14")
    end

    # Loading models
    if ESI_mode == "negative" || ESI_mode == "neg"
        reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg_neg.joblib"))
    elseif ESI_mode == "positive" || ESI_mode == "pos"
        reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg_pos.joblib"))
    else error("ESI_mode should be set to positive or negative")
    end

    fingerprint = DataFrame()
    cids = pcp.get_compounds(INCHIKEY, "inchikey")
    if isempty(cids)
        error("CID not found: The INCHIKEY ($INCHIKEY) could not be associated with a PubChem compound.")
    end
    shortest_cid = cids[argmin([cids[y].cid for y in 1:length(cids)])]
    fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(shortest_cid.isomeric_smiles,fingerprints=true, descriptors=false)))))
    if size(fingerprint,2) != 780
        error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set to 2DAPC")
    end

    input = Matrix(hcat(pH, fingerprint))[1,:]
    IE_pred = predict(reg, input)
    return IE_pred        
end