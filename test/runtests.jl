using IE_prediction
using Test

import IE_prediction
@testset "IE_prediction.jl" begin
    @test IE_prediction.logIE_from_SMILES("CCCOC(C)C","positive", 7) == 1.1706905569667065
    @test IE_prediction.logIE_from_SMILES("CCCOC(C)C","negative", 7) == -1.3006570357939462
    @test IE_prediction.logIE_from_SMILES("CN1C=NC2=C1C(=O)N(C(=O)N2C)C","positive", 13) == 2.0745916429350912 
    @test IE_prediction.logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N","positive", 7) == 1.1706905569667065
    @test IE_prediction.logIE_from_InChIKey("JIEJJGMNDWIGBJ-UHFFFAOYSA-N","negative", 7) == -1.3006570357939462
    @test IE_prediction.logIE_from_InChIKey("RYYVLZVUVIJVGH-UHFFFAOYSA-N","positive", 13) == 2.0745916429350912 
    @test IE_prediction.logIE_from_CNLs([30.01,30,29,2.1], 40.2,"positive", 7) == 1.3916559822536396
    @test IE_prediction.logIE_from_CNLs([2.1], 40.2,"negative", 7) == -1.3750886819949022
    @test IE_prediction.logIE_from_CNLs([300.01,30,29,2.13], 400.52,"positive", 3) == 3.799384180522007
end


