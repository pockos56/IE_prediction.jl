module IE_prediction

# Loading packages
    #using Conda
    using ScikitLearn
    using CSV
    using Statistics
    using DataFrames
    using PyCall
    using DataStructures
    jblb = pyimport("joblib")
    pcp = pyimport("pubchempy")
    pd = pyimport("padelpy")

    include("IE_from_SMILES.jl")
    include("IE_from_InChIKey.jl")
    include("unifiedIE_from_SMILES.jl")
    include("unifiedIE_from_InChIKey.jl")
    include("IE_from_CNLs.jl")

    export logIE_from_SMILES, logIE_from_InChIKey, ulogIE_from_SMILES, ulogIE_from_InChIKey, logIE_from_CNLs

end