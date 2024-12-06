using ToolBox
using Test

@testset verbose=true "setup" begin
    @test verbose(21,12) == 4
end