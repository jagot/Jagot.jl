# New matplotlib colormaps by Nathaniel J. Smith, Stefan van der Walt,
# and (in the case of viridis) Eric Firing.
#
# This file and the colormaps in it are released under the CC0 license /
# public domain dedication. We would appreciate credit if you use or
# redistribute these colormaps, but do not impose any legal restrictions.
#
# To the extent possible under law, the persons who associated CC0 with
# mpl-colormaps have waived all copyright and related or neighboring rights
# to mpl-colormaps.
#
# You should have received a copy of the CC0 legalcode along with this
# work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
#
# Adapted for Julia by Jagot

module colormaps
using PyCall
using PyPlot
@pyimport matplotlib.colors as COL

include("colormaps_data.jl")

macro define_colormap(name)
    name_str = string(name)
    data_sym = symbol("_$(name)_data")
    eval(quote
        const $name = COL.ListedColormap($data_sym, name=$name_str)
        export $name
    end)
end

@define_colormap(magma)
@define_colormap(inferno)
@define_colormap(plasma)
@define_colormap(viridis)

end
