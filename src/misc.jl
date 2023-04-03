ind(v::AbstractVector{T}, a::U) where {T<:Number,U<:Number} = argmin(abs.(v .- a))

meshgrid(x::AbstractVector,y::AbstractVector) =
    repeat(x', outer=(length(y),1)),repeat(y, outer=(1,length(x)))

export ind, meshgrid
