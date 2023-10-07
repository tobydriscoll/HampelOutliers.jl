using HampelOutliers
using Test,Random

@testset "Identification" begin
    rng = MersenneTwister(123)
    x = rand(rng,1000)
    k = [1,2,100,500,700,900,950]
    x[k] .= 4
    @test all(findall(identify(x)).==k)
    x[k] .= -1
    @test all(findall(identify(x)).==k)
end

@testset "Filtering" begin
    t = range(0,10,length=200)
    x = @. cos(t) + 0.1sin(4t)
    k = [1,5,100,101,102,103,199]
    x[k] .= 2
    
    y = filter(x,7)
    @test findall(y.!=x)==k
    z = filter(x,fill(1,15))
    @test all(y.==z)

    u = filter(x,[1,1,1,3,1,1,1])
end


