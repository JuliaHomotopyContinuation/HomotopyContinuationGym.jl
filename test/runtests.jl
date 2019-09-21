using HomotopyContinuationGym, Test

@testset "HomotopyContinuationGym" begin

    @testset "Systems" begin
        @test Steiner() isa TestSystem
        @test Cyclooctane() isa TestSystem
        @test FourBar() isa TestSystem
    end
end
