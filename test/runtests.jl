using IE_prediction
using Test

@testset "IE_prediction.jl" begin
    #@test IE_prediction.logIE_from_SMILES("CCCOC(C)C", 7) == 1.56133081771261
    #@test IE_prediction.logIE_from_SMILES("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 13) == 1.9962543539670365 
    #@test IE_prediction.logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 7) == 1.56133081771261
    #@test IE_prediction.logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 2.7) == 1.7215733476286195
    #@test IE_prediction.logIE_from_CNLs([30.01,30,22.19,16.21], 40.2, 7) == 1.7351528403607068
    #@test IE_prediction.logIE_from_CNLs([22.19], 40.2, 7) == 1.78702178814918
    #@test IE_prediction.logIE_from_CNLs([300.01,30,29,2.13], 400.52, 5) == 3.8657452530743783
    @test IE_prediction.logIE_from_structure("CC(=O)OC1=CC=CC=C1C(=O)O", 2.7) == 1.2200884862831083
end