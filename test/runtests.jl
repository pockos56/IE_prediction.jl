using IE_prediction
using Test

@testset "IE_prediction.jl" begin
    @test IE_prediction.logIE_from_SMILES("CCCOC(C)C", 7, "mean") == 1.56133081771261
    @test IE_prediction.logIE_from_SMILES("CCCOC(C)C", 7, "max") == 1.6364964476168147
    @test IE_prediction.logIE_from_SMILES("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 13, "mean") == 1.9962543539670365 
    @test IE_prediction.logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 7, "mean") == 1.56133081771261
    @test IE_prediction.logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 7, "min") == 1.59729888233392
    @test IE_prediction.logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 2.7, "mean") == 1.7215733476286195
    @test IE_prediction.logIE_from_CNLs([30.01,30,22.19,16.21], 40.2, 7, "mean") == 1.7351528403607068
    @test IE_prediction.logIE_from_CNLs([22.19], 40.2, 7, "max") == 1.78702178814918
    @test IE_prediction.logIE_from_CNLs([300.01,30,29,2.13], 400.52, 5, "max") == 3.8657452530743783
end