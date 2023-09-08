module IE_prediction

    # Loading packages
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
    using Pkg
    Pkg.build("PyCall")

    include("IE_from_SMILES.jl")
    include("IE_from_InChIKey.jl")
    include("IE_from_CNLs.jl")

    export logIE_from_SMILES, logIE_from_InChIKey, logIE_from_CNLs

end