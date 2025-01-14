module IE_prediction

# Loading packages
    #using Conda
    using ScikitLearn
    using CSV
    using Statistics
    using DataFrames
    using PyCall
    using DataStructures
    using ProgressBars
    
    jblb = pyimport("joblib")
    pcp = pyimport("pubchempy")
    pd = pyimport("padelpy")

    include("logIE_from_SMILES.jl")
    include("logIE_from_InChIKey.jl")
    include("logIE_from_CNLs.jl")
    include("logIE_from_structure.jl")

    export logIE_from_SMILES, logIE_from_InChIKey, logIE_from_CNLs, logIE_from_structure
end