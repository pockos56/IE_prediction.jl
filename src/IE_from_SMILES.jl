## import packages ##
using ScikitLearn
using CSV
using Statistics
using DataFrames
using PyCall
using Conda
using DataStructures
jblb = pyimport("joblib")
pcp = pyimport("pubchempy")
pd = pyimport("padelpy")


function logIE_from_SMILES(SMILES::String, ESI_mode::String, pH)

    # Loading models
    FP_reg_neg = jblb.load(joinpath(@__DIR__, "data", "FP_reg_neg.joblib"))
    FP_reg_pos = jblb.load(joinpath(@__DIR__, "data", "FP_reg_pos.joblib"))

    if pH > 14 || pH < 0 
        error("Set pH to a valid value between 0 and 14")
    end

    if ESI_mode == "negative" || ESI_mode == "neg"
        reg = FP_reg_neg
    elseif ESI_mode == "positive" || ESI_mode == "pos"
        reg = FP_reg_pos
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