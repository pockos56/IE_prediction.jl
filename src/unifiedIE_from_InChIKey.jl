## import packages ##
using ScikitLearn, CSV, Statistics, DataFrames, DataStructures, PyCall, Conda

"""
Predict the ionization efficiency (in log units) using InChIKey and the pH of the mobile phase.

# Examples
```julia-repl
julia> logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 7)
1.43249846779072
```
"""
INCHIKEY = "JIEJJGMNDWIGBJ-UHFFFAOYSA-N"
pH = 7
function ulogIE_from_InChIKey(INCHIKEY::String, pH)

    jblb = pyimport("joblib")
    pcp = pyimport("pubchempy")
    pd = pyimport("padelpy")
    cd(@__DIR__)

    # Loading models
    reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg.joblib"))

    if pH > 14 || pH < 0 
        error("Set pH to a valid value between 0 and 14")
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