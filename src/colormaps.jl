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

function get_cmap(args...)
    plt[:cm][:get_cmap](args...)
end

import Base.call
function call(cmap::PyPlot.ColorMap, i::Int)
    if :colors in keys(cmap)
        vec(cmap[:colors][i,:])
    else
        pycall(cmap, PyAny, i)
    end
end

function lerp(a,b,t)
    (1-t)*a + t*b
end

function lerp(a::Tuple, b::Tuple, t)
    tuple([lerp(a[i],b[i],t) for i = 1:length(a)]...)
end

function call(cmap::PyPlot.ColorMap, f::Real)
    i = clamp(1+f*(cmap[:N]-1),1,cmap[:N])
    fl,ce = floor(Integer,i), ceil(Integer, i)
    lerp(cmap(fl),cmap(ce),f)
end

call(cmap::PyPlot.ColorMap, i::Integer, l::Integer) =
    l != 0 ? cmap((i-1)/(l-1)) : cmap(0.5)

export get_cmap, call

end
