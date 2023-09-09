module IE_prediction

    # Loading packages
    using ScikitLearn
    using CSV
    using Statistics
    using DataFrames
    using PyCall
    #using Conda
    using DataStructures
    jblb = pyimport("joblib")
    #using JSON2
    #model = JSON2.read("r",x, "src\\data\\CNL_reg_neg.json")
    #using BSON
    # Saving the models (BSON)
    #BSON.@load "C:\\Users\\alex_\\Documents\\GitHub\\IE_prediction-project\\Models\\FP_reg_neg.bson" reg_neg
    #BSON.@load "C:\\Users\\alex_\\Documents\\GitHub\\IE_prediction-project\\Models\\FP_reg_pos.bson" reg_pos 
    pcp = pyimport("pubchempy")
    pd = pyimport("padelpy")
    #using Pkg
    #Pkg.build("PyCall")

    include("IE_from_SMILES.jl")
    include("IE_from_InChIKey.jl")
    include("IE_from_CNLs.jl")

    export logIE_from_SMILES, logIE_from_InChIKey, logIE_from_CNLs

end