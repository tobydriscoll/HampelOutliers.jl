module HampelOutliers

using StatsBase: mad,median

import Base: filter,filter!
export identify,identify!,screen,filter,filter!

# default function to measure spread of data
spread = x -> mad(x,normalize=true)

"""
    identify(x; spread=mad(x), threshold=3)

Identify outliers using the Hampel criterion.

Given vector `x`, identify elements xₖ such that

```math
|xₖ - m| > t S,
```

where ``m`` is the median of the elements, the dispersion scale ``S`` is provided by the 
function `spread`, and the parameter `t`` is given by `threshold`. The return value
is a `Bool` vector.

By default, `spread` is `StatsBase.mad` and `threshold` is 3.

# See also
[`screen`](@ref)
"""
identify(x::AbstractVector;kwargs...) = identify!(falses(size(x)),x;kwargs...)

"""
    identify(y,x; spread, threshold)

In-place version of `identify`.

# See also 
[`identify`](@ref)
"""
function identify!(y::AbstractVector,x::AbstractVector;spread=spread,threshold=3)
    S = spread(x)
    m = median(x)
    for k in eachindex(x)
        y[k] = abs(x[k]-m) > threshold*S
    end
    return y
end


"""
    screen(x; spread=mad(x), threshold=3)

Remove outliers using the Hampel criterion.

Returns a new copy of `x`, with elements flagged by `identify` removed.

# See also
[`identify`](@ref)
"""
function screen(x;kwargs...)
    idx = .!identify(x;kwargs...)
    return x[idx],findall[idx]
end

"""
    filter(x,K; spread=mad, threshold=2)
    filter(x,weights; ...)

Apply a weighted Hampel filter to a time series.

Given vector `x` and half-width `K`, apply a Hampel criterion within a 
sliding window of width 2K+1. The median ``m`` of the window replaces the element
``xₖ`` at the center of the window if it satisfies 

```math
|xₖ - m| > t S,
```

where the dispersion scale ``S`` is provided by the function `spread` and the parameter 
``t`` is given by `threshold`. The window shortens near the beginning and end of the vector 
to avoid referencing fictitious elements. Larger values of ``t`` make the filter less agressive,
while ``t=0`` is the standard median filter.

For recursive filtering, see `filter!` 

Given vector `x` and a vector `weights` of positive intgers, before computing the criterion 
each element in the window is repeated by the number of times given by its corresponding 
weight. This is typically used to make the central element more influential than the others.

# See also
[`filter!`](@ref), [`identify`](@ref)
"""
filter(x,K::Integer;kwargs...) = filter(x,fill(1,2K+1);kwargs...)
filter(x::AbstractVector,weights::Vector{<:Integer};kwargs...) = filter!(similar(x),x,weights;kwargs...)

"""
    filter!(y,x,K; spread=mad, threshold=2)
    filter(y,x,weights; ...)

Apply a weighted Hampel filter in-place.

The idiom `filter!(x,x,...)` will make the filter *recursive*, i.e., vector elements are replaced as 
they are found, possibly affecting future results. 
"""
filter!(y::AbstractVector,x::AbstractVector,K::Integer;kwargs...) = filter!(y,x,fill(1,2K+1);kwargs...)
function filter!(y::AbstractVector,x::AbstractVector,weights::Vector{<:Integer};
                       spread=spread,threshold=2)

    KK1 = length(weights)
    @assert isodd(KK1) "Must specify an odd number of weights"
    K = (KK1-1)÷2  # filter half-width
    N = length(x)
    # Set up duplications as specified by the weights.
    idx = vcat( [fill(i,m) for (i,m) in zip(-K:K,weights)]... )
    for n in 1:N
        @. idx += 1
        valid = @. 1 <= idx <= N
        window = view(x,idx[valid])
        wmed = median(window)
        S = spread(window)
        out = abs(x[n] - wmed) > threshold*S
        y[n] =  out ? wmed : x[n]
    end

    return y
end

end 

