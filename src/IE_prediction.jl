module IE_prediction

    export logIE_from_SMILES
    include("IE_from_SMILES.jl")

    export logIE_from_InChIKey
    include("IE_from_InChIKey.jl")

    export logIE_from_CNLs
    include("IE_from_CNLs.jl")

end