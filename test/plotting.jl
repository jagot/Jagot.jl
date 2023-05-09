using PythonPlot
using Jagot.plotting

@testset "Plotting" begin
    plot_style("ggplot")

    @testset "Simple plots" begin
        fig = cfigure("my figure") do
            csubplot(211, nox=true) do
                plot(rand(5))
            end
            csubplot(2,2,(2,1)) do
                plot(rand(5))
            end
            csubplot(2,2,(2,2)) do
                plot(rand(5))
            end
        end
        @test fig isa Figure

        cfigure("my second figure") do
            m,n = subplot_ratio(16,1)
            @test (m,n) == (4,4)
            gr = grid_spec(gcf(), m, n)
            axs = gr_sub_plots(gr)
            @test length(axs) == 16
            ci = CartesianIndices((m,n))
            for I in ci
                csubplot(axs[I]) do
                    plot(rand(10))
                    frac_ticks([1//3, 18//5])
                    π_labels(:y)
                end
            end
        end
    end

    A = Matrix(reshape(1:16, 4, 4))
    A[:,4] .= 0

    @testset "Map plots" begin
        @testset "Map plot" begin
            cfigure("my map plot") do
                plot_map(A)
            end
        end

        @testset "Polar map plots" begin
            r = range(0.3, stop=1.0, length=101)
            θ = range(0, stop=π, length=101)
            z = r' .* sin.(θ)

            cfigure("my polar map plot") do
                plot_polar_map(r, θ, z)
                axis("equal")
            end
        end
    end

    @testset "Matrix plots" begin
        cfigure("my matrix plot") do
            plot_matrix(A)
        end
    end

    @testset "Colormaps" begin
        cmap = get_cmap("viridis")
        @test cmap(0.0) ≈ [0.267004, 0.004874, 0.329415]
        @test cmap(0.5) ≈ [0.1281485, 0.565107, 0.5508925]
        @test cmap(1.0) ≈ [0.993248, 0.906157, 0.143936]

        @test cmap(000) ≈ [0.267004, 0.004874, 0.329415]
        @test cmap(255) ≈ [0.993248, 0.906157, 0.143936]

        @test cmap(1, 7) ≈ [0.267004, 0.004874, 0.329415]
        @test cmap(4, 7) ≈ [0.1281485, 0.565107, 0.5508925]
        @test cmap(7, 7) ≈ [0.993248, 0.906157, 0.143936]
    end
end
