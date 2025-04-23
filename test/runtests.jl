using IE_prediction
using Test

@testset "IE_prediction.jl" begin
    @test IE_prediction.logIE_from_structure("CCCOC(C)C", 7) == 1.4138901778271273
    @test IE_prediction.logIE_from_structure("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 13) == 2.1577971679270567 
    @test IE_prediction.logIE_from_structure("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 7) == 1.4138901778271273
    @test IE_prediction.logIE_from_structure("JIEJJGMNDWIGBJ-UHFFFAOYSA-N", 2.7) == 1.2989089872507973
    @test IE_prediction.logIE_from_CNLs([30.01,30,22.19,16.21], 40.2, 7) == 1.1460003005991275
    @test IE_prediction.logIE_from_CNLs([22.19], 40.2, 2.7) == 2.2589066230991515
    @test IE_prediction.logIE_from_CNLs([300.01,30,29,2.13], 400.52, 5) == 3.579897331838181
    @test IE_prediction.logIE_from_structure("CC(=O)OC1=CC=CC=C1C(=O)O", 2.7) == 1.2200884862831083
end