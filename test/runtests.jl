using Test
using TethysPlots
using Tethys
using Plots
using Graphs
using DataFrames
using SparseArrays
using LinearAlgebra

@testset "TethysPlots Test Suite" begin
    # 1. Synthetic Geometry & Tessellation Setup
    poly1 = [(-64.0, 44.0), (-63.0, 44.0), (-63.0, 45.0), (-64.0, 45.0)]
    poly2 = [(-63.0, 44.0), (-62.0, 44.0), (-62.0, 45.0), (-63.0, 45.0)]
    polys = [poly1, poly2]
    cents = [(-63.5, 44.5), (-62.5, 44.5)]
    g = SimpleGraph(2)
    add_edge!(g, 1, 2)
    W = spzeros(Float64, 2, 2)
    W[1, 2] = 1.0
    W[2, 1] = 1.0

    au = (
        polygons = polys,
        centroids = cents,
        graph = g,
        W = W,
        hull_coords = [(-64.0, 44.0), (-62.0, 44.0), (-62.0, 45.0), (-64.0, 45.0)]
    )

    @testset "Coordinate Transformations" begin
        # Zone 20 UTM coordinate conversion test
        lon, lat = utm_to_lonlat(500.0, 4900.0; zone = 20, is_km = true)
        @test -66.0 <= lon <= -60.0
        @test 43.0 <= lat <= 46.0

        x_km, y_km = lonlat_to_utm(lon, lat; zone = 20, is_km = true)
        @test isapprox(x_km, 500.0, atol = 1e-2)
        @test isapprox(y_km, 4900.0, atol = 1e-2)
    end

    @testset "Static Choropleth & Hexagonal Field" begin
        vals = [12.5, 18.2]
        p_choro = choropleth(polys, vals; title = "Test Choropleth")
        @test p_choro isa Plots.Plot

        p_hex = plot_hexagonal_field(vals, au; title = "Test Hex Field")
        @test p_hex isa Plots.Plot

        p_hsi = plot_hsi_choropleth([0.3, 0.8], au; show_centroids = true)
        @test p_hsi isa Plots.Plot

        p_diff = plot_diffusion_map([0.5, 1.2], au)
        @test p_diff isa Plots.Plot

        Gamma = [0.8 0.2; 0.3 0.7]
        p_res = plot_residence_time_map(Gamma, au)
        @test p_res isa Plots.Plot
    end

    @testset "Spatial Partitioning & Graph Plot" begin
        p_graph = spatial_graph_plot(cents, g; polygons = polys)
        @test p_graph isa Plots.Plot

        p_part = plot_partition(au; title = "Areal Units")
        @test p_part isa Plots.Plot
    end

    @testset "Continuous Spatial Surface" begin
        coords = [-63.5 44.5; -62.5 44.5; -63.0 44.8]
        z_vals = [2.0, 5.0, 3.5]
        p_surf = plot_spatial_surface(coords, z_vals; grid_res = 30)
        @test p_surf isa Plots.Plot

        df = DataFrame(lon = [-63.5, -62.5, -63.0], lat = [44.5, 44.5, 44.8])
        p_surf_df = plot_spatial_surface(df, z_vals; grid_res = 30)
        @test p_surf_df isa Plots.Plot
    end

    @testset "Advection & Velocity Fields" begin
        p_adv = plot_advection_arrows(au; hsi = [0.2, 0.9], velocity = 1.5)
        @test p_adv isa Plots.Plot

        p_vel = plot_velocity_field(au; hsi = [0.2, 0.9])
        @test p_vel isa Plots.Plot
    end

    @testset "Hydrodynamic Diagnostics" begin
        hydro = (
            N2 = [1e-4 2e-4; 1.5e-4 2.5e-4],
            salinity = [31.5 32.8; 31.8 33.1],
            temperature = [14.0 6.0; 13.5 5.5],
            kappa_v = [1e-3 1e-5; 2e-3 1.2e-5],
            nu_v = [1e-2 1e-4; 1.5e-2 1.1e-4],
            depth_levels = [5.0, 50.0],
            depths = [80.0, 110.0]
        )

        p_strat = plot_hydrodynamic_stratification(hydro, au; depth_level = 1)
        @test p_strat isa Plots.Plot

        p_diff = plot_hydrodynamic_diffusion(hydro, au; depth_level = 1)
        @test p_diff isa Plots.Plot

        p_sec = plot_hydrodynamic_section(hydro, au; var = :temperature)
        @test p_sec isa Plots.Plot
    end

    @testset "Timeseries CI & Export" begin
        x = 2010:2020
        y_mean = sin.(1:11)
        y_low = y_mean .- 0.5
        y_up = y_mean .+ 0.5

        p_ci = timeseries_ci(x, y_mean, y_low, y_up; title = "Trend")
        @test p_ci isa Plots.Plot

        theme = create_theme(; fontsize = 11)
        @test theme.legendfontsize == 11

        tmp_file = joinpath(tempdir(), "test_tethys_plot.png")
        saved = save_plot(p_ci, tmp_file)
        @test isfile(saved)
        rm(saved; force = true)
    end

    @testset "Interactive Leaflet HTML Maps" begin
        vals = [12.5, 18.2]
        m_choro = choropleth(polys, vals; mode = :leaflet)
        @test m_choro isa LeafletMap
        @test occursin("<html", m_choro.html)

        m_graph = spatial_graph_plot(au = au, mode = :leaflet)
        @test m_graph isa LeafletMap

        m_hsi = plot_hsi_choropleth([0.3, 0.8], au; mode = :leaflet)
        @test m_hsi isa LeafletMap

        m_diff = plot_diffusion_map([0.5, 1.2], au; mode = :leaflet)
        @test m_diff isa LeafletMap

        Gamma = [0.8 0.2; 0.3 0.7]
        m_res = plot_residence_time_map(Gamma, au; mode = :leaflet)
        @test m_res isa LeafletMap

        m_adv = plot_advection_arrows(au; hsi = [0.2, 0.9], mode = :leaflet)
        @test m_adv isa LeafletMap

        tmp_html = joinpath(tempdir(), "test_leaflet_map.html")
        saved_html = save_html(m_choro, tmp_html)
        @test isfile(saved_html)
        rm(saved_html; force = true)
    end
end
