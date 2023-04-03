using RecipesBase

function replace_attr!(r::RecipeData, src::Symbol, dst::Symbol, trf=identity)
    if src ∈ keys(r.plotattributes)
        r.plotattributes[dst] = trf(r.plotattributes[src])
        delete!(r.plotattributes, src)
    end
end

function rplot(v, args...; kwargs...)
    map(RecipesBase.apply_recipe(Dict{Symbol,Any}(), v, args...)) do r
        replace_attr!(r, :markershape, :marker,
                      s -> Dict(:circle => "o")[s])
        plot(r.args...; r.plotattributes..., kwargs...)
    end
end

export rplot
