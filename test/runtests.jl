using HampelOutliers
using Test, Random

# Some results are validated against python's hampel. But that package uses unnormalized
# MAD and doesn't offer the repeated-value option. It also handles boundaries differently.

# These tests aren't truly adequate, but they're better than nothing.

@testset "Identification" begin
    rng = MersenneTwister(123)
    x = rand(rng,1000)
    k = [1,2,100,500,700,900,950]
    x[k] .= 4
    @test all(findall(Hampel.identify(x)).==k)
    x[k] .= -1
    @test all(findall(Hampel.identify(x)).==k)
end

@testset "Filtering" begin
    t = range(0, 10, length=200)
    x = @. cos(t) + 0.1sin(4t)
    k = [1,5,100,101,102,103,199]
    x[k] .= 2

    y = Hampel.filter(x, 3)
    @test findall(y.!=x) == [1, 5, 199]
    y = Hampel.filter(x, 5)
    @test findall(y.!=x) == k
    z = Hampel.filter(x, fill(1, 11))
    @test all(y.==z)
    y = Hampel.filter(x, 5, threshold=3)
    @test findall(y.!=x) == [1, 5, 199]

    u = Hampel.filter(x,[1,1,1,3,1,1,1])
    v = Hampel.filter(x[[94,95,96,97,97,97,98,99,100]], 4)
    @test u[97] ≈ v[5]
end
