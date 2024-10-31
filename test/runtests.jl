using Test


using ToolBox
using FHist


@testset verbose=true "ToolBox.jl" begin
    @testset verbose=true "mathutils.jl" begin
        @testset "cosbetween" begin
        end
    end

    @testset verbose=true "statisticsutils.jl" begin
        @testset "efficiency" begin
            sample = [
                Hist1D([1,1,4,5,4.1,3.1,7,8];binedges=0.0:10.0), # signal
                Hist1D([1,1,5,5,5,8];binedges=0.0:10.0) # background
            ]
            x1,x2 = 3,6
            signal_index = 1
            @test efficiency(sample,signal_index,x1,x2)
        end
    end
end