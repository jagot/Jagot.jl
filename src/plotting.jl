# * Setup
module plotting
using PyPlot
using PyCall

const COL = PyNULL()
const masked_array = PyNULL()
const pgf_backend = PyNULL()
const pyplot_collections = PyNULL()
const RegularPolyCollection = PyNULL()

function __init__()
    copy!(COL, pyimport_conda("matplotlib.colors", "matplotlib"))
    copy!(masked_array, pyimport_conda("numpy.ma", "numpy"))
    copy!(pgf_backend, pyimport_conda("matplotlib.backends.backend_pgf", "matplotlib"))
    copy!(pyplot_collections, pyimport_conda("matplotlib.collections", "matplotlib"))
    copy!(RegularPolyCollection, pyplot_collections.RegularPolyCollection)
end

using UnicodeFun

import Jagot: meshgrid

# using GSL
using Printf

plot_style(style::String) = matplotlib.style.use(style)

# * Figure wrappers
import PyPlot: figure, subplot

function figure(fun::Function, args...; kwargs...)
    old_fig = gcf()
    figure(args...; kwargs...)
    fun()
    figure(old_fig.number)
end

function subplot(fun::Function, args...; kwargs...)
    old_ax = gca()
    subplot(args...; kwargs...)
    fun()
    sca(old_ax)
end

# * Colormaps

include("colormaps.jl")

function colorbar_hack(v::AbstractVector,cmap)
    Z = [0 0
         0 0]
    CS3 = contourf(Z, v, cmap=cmap)
    clf()
    CS3
end

# * Map plots

function filter_kwargs!(kwargs, sym, default=nothing)
    v = get(kwargs, sym, default)
    filter!(kv -> kv[1] != sym, kwargs)
    v
end

function Base.size(pyobj::PyCall.PyObject)
    :shape ∈ keys(pyobj) ||
        throw(ArgumentError("`:shape` not present among `PyObject`s `keys`."))
    pyobj.shape
end
Base.size(pyobj::PyCall.PyObject, i) = size(pyobj)[i]

function plot_map(args...; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)
    set_default!(key, val) = (kwargs[key] = get(kwargs, key, val))
    kwargs[:cmap] = get_cmap(get(kwargs, :cmap, "viridis"))
    set_default!(:rasterized, true)
    get(kwargs, :norm, :lin) == :log && (kwargs[:norm] = COL.LogNorm())
    aw = filter_kwargs!(kwargs, :align_ticks, false)
    xtl = filter_kwargs!(kwargs, :xticklabels, nothing)
    ytl = filter_kwargs!(kwargs, :yticklabels, nothing)
    tickrotation = filter_kwargs!(kwargs, :tickrotation, 0)
    xtickrotation = filter_kwargs!(kwargs, :xtickrotation, tickrotation)
    ytickrotation = filter_kwargs!(kwargs, :ytickrotation, tickrotation)
    xticker,yticker = if aw != false
        x,y = if length(args) == 3
            x,y,z = args
            x,y
        elseif length(args) == 1
            z = args[1]
            x = 1:size(z,2)
            y = 1:size(z,1)
            x,y
        else
            @warn "Don't know how to shift x/y axes in case of $(length(args)) `args`"
            nothing,nothing
        end
        xt = !isnothing(x) && aw != :y ? () -> xticks(x .- (x[2]-x[1])/2,
                                                     isnothing(xtl) ? x : xtl,
                                                     rotation=xtickrotation) : () -> ()
        yt = !isnothing(y) && aw != :x ? () -> yticks(y .- (y[2]-y[1])/2,
                                                     isnothing(ytl) ? x : xtl,
                                                     rotation=ytickrotation) : () -> ()
        xt,yt
    else
        () -> (), () -> ()
    end
    p = pcolormesh(args...; kwargs...)
    margins(0,0)
    xticker()
    yticker()
    p
end

function plot_polar_map(r::AbstractVector, ϕ::AbstractVector, V::AbstractMatrix,
                        args...; kwargs...)
    nϕ = length(ϕ)
    R = repeat(r,1,nϕ)
    Φ = repeat(ϕ',length(r),1)
    X = R.*cos.(Φ)
    Y = R.*sin.(Φ)
    plot_map(X,Y,V, args...; kwargs...)
end

plot_polar_map(r::AbstractVector, V::AbstractMatrix, nϕ::Integer, args...; kwargs...) =
    plot_polar_map(r, range(0,stop=2π,length=nϕ), V, args...; kwargs...)

plot_polar_map(r::AbstractVector, v::AbstractVector, nϕ::Integer, args...; kwargs...) =
    plot_polar_map(r, repeat(v,1,nϕ), nϕ, args...; kwargs...)

# function spherical_harmonic_plot(fun::Function,
#                                  r::AbstractVector,
#                                  J::AbstractVector{Int},
#                                  v::AbstractMatrix,
#                                  nθ::Int,
#                                  args...;
#                                  wfn = false,
#                                  J_max = Inf,
#                                  m_min = 0, m_max = Inf,
#                                  kwargs...)
#     nr = length(r)
#     (nr,length(J)) == size(v) || (wfn && (nr,length(J)^2) == size(v)) || error("Dimension mismatch!")

#     θ = range(0,stop=2π,length=nθ)
#     cosθ = cos.(θ)

#     R = repmat(r,1,nθ)
#     Θ = repmat(θ',nr,1)

#     X = R.*cos.(Θ)
#     Y = R.*sin.(Θ)

#     V = zeros(eltype(v), nr, nθ)

#     leg_fun = (wfn ? sf_legendre_sphPlm : sf_legendre_Plm)

#     if size(v,2) == length(J)
#         for j in eachindex(J)
#             J[j] > J_max && break
#             PJ = map(x -> leg_fun(J[j], 0, x), cosθ)'
#             V += broadcast(*, v[:,j], PJ)
#         end
#     elseif !wfn
#         error("3d plot only implemented for wavefunctions")
#     else
#         j = 0
#         for ell in J
#             ell > J_max && break
#             for m in -ell:ell
#                 j += 1
#                 abs(m) > m_max && continue
#                 abs(m) < m_min && continue
#                 Pell_m = (m<0 ? (-1)^abs(m) : 1
#                           )*map(x -> leg_fun(ell, abs(m), x), cosθ)'
#                 V += broadcast(*, v[:,j], Pell_m)
#             end
#         end
#     end

#     plot_map(X, Y, fun.(V), args...; kwargs...)
# end
# spherical_harmonic_plot(r::AbstractVector,
#                         J::AbstractVector{Int},
#                         v::AbstractMatrix,
#                         nθ::Int,
#                         args...; kwargs...) = spherical_harmonic_plot(identity,
#                                                                       r::AbstractVector,
#                                                                       J::AbstractVector{Int},
#                                                                       v::AbstractMatrix,
#                                                                       nθ,
#                                                                       args...; kwargs...)
# * Matrix plots

function plot_matrix(a, args...; align_ticks=true, kwargs...)
    aa = pycall(masked_array.masked_equal, Any, Matrix(a), 0)
    plot_map(aa, args...; align_ticks=align_ticks, kwargs...)
    gca().invert_yaxis()
    square_axs()
end

# # Inspired by/stolen from
# # https://github.com/tonysyu/mpltools/blob/master/mpltools/special/hinton.py

# @pydef mutable struct SquareCollection <: pyplot_collections.RegularPolyCollection
#     __init__(self; kwargs...) = begin
#         pyplot_collections.RegularPolyCollection.__init__(
#             self, 4, rotation=pi/4; kwargs...
#         )
#     end
#     get_transform(self) = begin
#         """Return transform scaling circle areas to data space."""
#         ax = self.axes
#         pts2pixels = 72.0 / ax.figure.dpi
#         scale_x = pts2pixels * ax.bbox.width / ax.viewLim.width
#         scale_y = pts2pixels * ax.bbox.height / ax.viewLim.height
#         matplotlib.transforms.Affine2D().scale(scale_x, scale_y)
#     end
# end

# function hinton_plot_matrix(a::AbstractMatrix;
#                             max_weight=nothing,
#                             pos_color=nothing,
#                             neg_color=nothing,
#                             bg_color="none",
#                             traditional=false)
#     if traditional
#         pos_color = "white"
#         neg_color = "black"
#         bg_color = "gray"
#     else
#         pos_color == nothing && (pos_color="C0")
#         neg_color == nothing && (neg_color="C1")
#     end

#     ax = gca()
#     ax.set_facecolor(bg_color)

#     height, width = size(a)
#     if max_weight == nothing
#         max_weight = 2^ceil(log2(maximum(abs, a)))
#     end

#     vals = clamp.(a/max_weight, -1, 1)
#     neg = vals .< 0
#     pos = vals .> 0

#     cols,rows = meshgrid(1:width,1:height)

#     for (sel,color) in zip([neg,pos], [neg_color,pos_color])
#         if any(sel)
#             xy = collect(zip(cols[sel], rows[sel]))
#             circle_areas = π/2 * abs.(vals[sel])
#             squares = SquareCollection(sizes=circle_areas,
#                                        offsets=xy, transOffset=ax.transData,
#                                        facecolor=color, edgecolor=color)
#             ax.add_collection(squares, autolim=true)
#         end
#     end
#     square_axs()
#     grid(false)
#     ax.set_xlim(0.5, width+0.5)
#     ax.set_ylim(0.5, height+0.5)
#     gca().invert_yaxis()
#     ax
# end
# hinton_plot_matrix(a;kwargs...) = hinton_plot_matrix(full(a);kwargs...)

# * LaTeX/font setup

function set_pgf_to_pdf(preamble=[]; texsystem = "xelatex")
    matplotlib.rcdefaults()
    ion()
    matplotlib.backend_bases.register_backend("pdf", pgf_backend.FigureCanvasPgf)

    set_latex_serif()
    rc("pgf";
       texsystem = texsystem,
       preamble = preamble)
end

function set_font(; kwargs...)
    if :serif in keys(kwargs)
        kwargs[:family] = "serif"
    elseif :sans_serif in keys(kwargs)
        kwargs[:family] = "sans-serif"
    end

    kwargs = Dict(Symbol(replace(string(k), "_" => "-")) => v
                  for (k,v) in kwargs)

    rc("font"; kwargs...)
end

function set_times_new_roman(; kwargs...)
    set_font(serif = "Times New Roman"; kwargs...)
    rc("mathtext";
       rm = "serif",
       it = "serif:italic",
       bf = "serif:bold",
       fontset = "stix")
end

function set_latex_serif(; kwargs...)
    rc("text";
       usetex = true)
    rc("font";
       family = "serif")
    rc("pgf";
       rcfonts = false)
    rc("mathtext";
       rm = "serif",
       it = "serif:italic",
       bf = "serif:bold",
       fontset = "cm")
    set_font(; kwargs...)
end

function latex(o)
    io = IOBuffer()
    show(io, MIME"text/latex"(), o)
    takebuf_string(io)
end

# * Scientific notation

function latex_base10(v)
    if v==0
        return "0"
    end
    ex = floor(log10(abs(v)))
    v /= 10^ex
    if ex == 0
        @sprintf("%0.7f",v)
    elseif ex == 1
        @sprintf("%0.7f\\cdot10",v)
    else
        @sprintf("%0.7f\\cdot10^{%i}",v,ex)
    end
end
latex_base10(v::Vector) = map(latex_base10, v)

function base10(v)
    if v==0
        return "0"
    end
    ex = floor(log10(abs(v)))
    v /= 10^ex
    ex = floor(Int, ex)
    if v==1
        @sprintf("10%s",to_superscript(ex))
    else
        @sprintf("%0.3f×10%s",v,to_superscript(ex))
    end
end
base10(v::Vector) = map(base10, v)

# * Axes/tick(label)s

function axis_add_ticks(ticks, labels, ax = :x; kwargs...)
    a = gca()
    a2 = a[ax == :x ? :twiny : :twinx]()
    a2[ax == :x ? :set_xlim : :set_ylim](a[ax == :x ? :get_xlim : :get_ylim]())
    a2[ax == :x ? :set_xscale : :set_yscale](a[ax == :x ? :get_xscale : :get_yscale]())
    a2[ax == :x ? :set_xticks : :set_yticks](vec(collect(ticks)))
    a2[ax == :x ? :set_xticklabels : :set_yticklabels](vec(collect(labels)); kwargs...)
    sca(a)
end

function set_ticklabel_props(ax=:x; kwargs...)
    labels = gca()[ax == :x ? :get_xticklabels : :get_yticklabels]()
    setp(labels; kwargs...)
end

function π_frac_string(i; pi_sym = "\\pi")
    den_str = v -> denominator(v) != 1 ? "/$(denominator(v))" : ""
    if i == 0
        L"$0$"
    elseif numerator(i) == 1
        latexstring("\$$(pi_sym)$(den_str(i))\$")
    elseif numerator(i) == -1
        latexstring("\$-$(pi_sym)$(den_str(i))\$")
    else
        latexstring("\$$(numerator(i))$(pi_sym)$(den_str(i))\$")
    end
end

function π_labels(ax = :x, max_count = 10; pi_sym = "\\pi")
    lims = gca()[ax == :x ? :get_xlim : :get_ylim]()
    f = 2
    mi,ma = map(l -> trunc(Int, l/(π/2)), lims)
    if ma-mi > max_count
        f = 1
        mi,ma = map(l -> trunc(Int, l/π), lims)
    end
    d = round(Int, max((ma-mi)/max_count,1))
    r = (mi:d:ma)//f
    gca()[ax == :x ? :set_xticks : :set_yticks](collect(r)*π)

    tick_labels = map(r) do i
        π_frac_string(i, pi_sym = pi_sym)
    end
    gca()[ax == :x ? :set_xticklabels : :set_yticklabels](tick_labels)
end

function frac_ticks(ts::AbstractVector{Rational{Int}}, axis = :x; sfrac = false)
    if axis == :x
        xticks(ts)
    else
        yticks(ts)
    end
    tls = map(ts) do t
        if t == 0
            L"0"
        elseif denominator(t)==1
            latexstring("$(numerator(t))")
        else
            if sfrac
                latexstring("\\sfrac{$(numerator(t))}{$(denominator(t))}")
            else
                latexstring("$(numerator(t))/$(denominator(t))")
            end
        end
    end
    gca()[axis == :x ? :set_xticklabels : :set_yticklabels](tls)
end

function sci_ticks(ax, min=0, max=0)
    ticklabel_format(style="sci", axis="$ax", scilimits=(min,max))
end

function colorbar_sci_ticks(cb, min=0, max=0)
    cb.formatter.set_scientific(true)
    cb.formatter.set_powerlimits((min, max))
    cb.update_ticks()
end

function square_axs()
    gca().set_aspect("equal")
    gca().autoscale(tight=true)
end

function axes_labels_opposite(axis = :y, ax = gca())
    if axis == :y
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    else
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
    end
end

function no_tick_labels(axis = :x, ax = gca(), ticks = false)
    a,b,la,lb = if axis == :x
        :bottom,:top,:labelbottom,:labeltop
    else
        :left,:right,:labelleft,:labelright
    end
    ax.tick_params(; which="both",
                   Dict(a => ticks, b => ticks,
                        la => false, lb => false)...)
end

# * Misc

pyslice(args...) = pycall(pybuiltin("slice"), PyObject, args...)

function savefig_f(filename, args...; kwargs...)
    savefig(filename, args...; kwargs...)
    filename
end

reltext(x,y,string,args...;kwargs...) = text(x,y,string,transform=gca().transAxes,args...;kwargs...)

function disp()
    display(gcf())
    println()
end

GridSpec = matplotlib.gridspec.GridSpec

# * ICC support
include("save_pgf_with_icc.jl")

# * Default settings
hide_toolbar() = (matplotlib.rcParams["toolbar"] = nothing)

# * Exports

export plot_style, figure, subplot,
    colormaps, colorbar_hack,
    plot_map, plot_polar_map, spherical_harmonic_plot, plot_matrix, hinton_plot_matrix,
    set_pgf_to_pdf, set_font, set_times_new_roman, set_latex_serif,
    latex, latex_base10, base10,
    axis_add_ticks, set_ticklabel_props, π_frac_string, π_labels, frac_ticks, sci_ticks, colorbar_sci_ticks,
    square_axs, axes_labels_opposite, no_tick_labels,
    pyslice, savefig_f, reltext, disp, hide_toolbar,
    GridSpec

end
