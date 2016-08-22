ind{T<:Number,U<:Number}(v::AbstractVector{T}, a::U) = indmin(abs(v-a))

export ind
