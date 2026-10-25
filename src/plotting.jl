"""
    plotting.jl

Unified static and interactive plotting library for the Tethys modeling ecosystem.
Provides publication-grade spatial maps, discrete tessellation choropleths, continuous
random field surfaces, vector drift fields, hydrodynamic water-column diagnostics, and
time series credible interval ribbons.
"""

# =============================================================================
# Section 1: Geometric & Theme Helpers
# =============================================================================

# Canonical polygon decomposition and centroid extraction routines are provided by Tethys
const _poly_xy = poly_xy
const _centroids_xy = centroids_xy
const _extract_polygons = extract_polygons
const _extract_centroids = extract_centroids

"""
    _get_cscheme(cmap)

Resolves a `ColorSchemes.ColorScheme` object from a `Symbol`, `ColorScheme`, or
fallback.
"""
function _get_cscheme(cmap)
    if cmap isa ColorSchemes.ColorScheme
        return cmap
    elseif cmap isa Symbol
        try
            return get(ColorSchemes, cmap, ColorSchemes.viridis)
        catch
            return ColorSchemes.viridis
        end
    else
        return ColorSchemes.viridis
    end
end

"""
    _map_color_from_range(cscheme, v, vmin, vmax)

Linearly maps a scalar value `v` in `[vmin, vmax]` to an RGB color from `cscheme`.
Returns light gray for `NaN` values without using clamp truncation.
"""
function _map_color_from_range(cscheme, v, vmin, vmax)
    if isnan(v)
        return RGB(0.9, 0.9, 0.9)
    end
    denom = vmax - vmin
    t = if abs(denom) <= 1e-12
        0.5
    elseif v < vmin
        0.0
    elseif v > vmax
        1.0
    else
        (v - vmin) / denom
    end
    return get(cscheme, t)
end

"""
    create_theme(; fontsize=10, titlefontsize=12, font="DejaVu Sans", size=(900,600))

Returns a `NamedTuple` of standard plotting attributes for pass-through to Plots functions.
"""
function create_theme(;
    fontsize::Int = 10,
    titlefontsize::Int = 12,
    font::String = "DejaVu Sans",
    size::Tuple{Int, Int} = (900, 600)
)
    return (;
        legendfontsize = fontsize,
        guidefontsize = fontsize,
        tickfontsize = max(8, fontsize - 1),
        titlefontsize = titlefontsize,
        fontfamily = font,
        size = size
    )
end

# =============================================================================
# Section 2: Spatial Primitives (Choropleth, Mesh, Surface)
# =============================================================================

"""
    _choropleth_impl(polygons, values; cmap=:viridis, vmin=nothing, vmax=nothing,
                     clims=nothing, center_zero=false, show_colorbar=true,
                     title="", border_color=:black, border_alpha=0.6,
                     lw=0.4, closed=true, theme_kwargs=NamedTuple(),
                     colorbar_label=nothing, kwargs...)

Internal rendering function generating a static polygon choropleth figure using Plots.jl.
"""
function _choropleth_impl(
    polygons,
    values;
    cmap = :viridis,
    vmin = nothing,
    vmax = nothing,
    clims = nothing,
    center_zero::Bool = false,
    show_colorbar::Bool = true,
    title::String = "",
    border_color = :black,
    border_alpha::Real = 0.6,
    lw::Real = 0.4,
    closed::Bool = true,
    theme_kwargs = NamedTuple(),
    colorbar_label::Union{Nothing, String} = nothing,
    kwargs...
)
    n = length(polygons)
    @assert length(values) == n "length(values) ($(length(values))) must equal length(polygons) ($n)"

    vals = collect(Float64, values)
    mask = .!isnan.(vals)

    if !any(mask)
        p = Plots.plot(
            aspect_ratio = :equal,
            title = title,
            legend = false;
            theme_kwargs...,
            kwargs...
        )
        for poly in polygons
            if length(poly) < 3
                continue
            end
            xs, ys = _poly_xy(poly)
            if closed && !isempty(xs) && (xs[1], ys[1]) != (xs[end], ys[end])
                push!(xs, xs[1])
                push!(ys, ys[1])
            end
            if !isempty(xs)
                Plots.plot!(
                    p, xs, ys,
                    seriestype = :shape,
                    fillcolor = :white,
                    linecolor = border_color,
                    lw = lw,
                    alpha = border_alpha,
                    label = nothing
                )
            end
        end
        return p
    end

    cs_sym = center_zero ? :RdBu : (cmap isa Symbol ? cmap : :viridis)
    cs = _get_cscheme(cs_sym)
    vvals = vals[mask]
    if clims !== nothing
        vmin = clims[1]
        vmax = clims[2]
    end
    if vmin === nothing
        vmin = length(vvals) > 1 ? quantile(vvals, 0.02) : minimum(vvals)
    end
    if vmax === nothing
        vmax = length(vvals) > 1 ? quantile(vvals, 0.98) : maximum(vvals)
    end

    if center_zero
        limit = max(abs(vmin), abs(vmax))
        vmin, vmax = -limit, limit
    end

    p = Plots.plot(
        aspect_ratio = :equal,
        title = title,
        legend = false;
        theme_kwargs...,
        kwargs...
    )

    for (i, poly) in enumerate(polygons)
        if length(poly) < 3
            continue
        end
        xs, ys = _poly_xy(poly)
        if closed && !isempty(xs) && (xs[1], ys[1]) != (xs[end], ys[end])
            push!(xs, xs[1])
            push!(ys, ys[1])
        end
        if !isempty(xs)
            fillcolor = mask[i] ?
                _map_color_from_range(cs, vals[i], vmin, vmax) :
                RGB(0.9, 0.9, 0.9)
            Plots.plot!(
                p, xs, ys,
                seriestype = :shape,
                fillcolor = fillcolor,
                linecolor = border_color,
                lw = lw,
                alpha = 0.95,
                label = nothing
            )
        end
    end

    if show_colorbar && abs(vmax - vmin) > 1e-12
        Plots.scatter!(
            p, [NaN], [NaN],
            zcolor = [vmin],
            clims = (vmin, vmax),
            markersize = 0,
            markeralpha = 0,
            c = cs_sym,
            colorbar = true,
            label = nothing
        )
        if colorbar_label !== nothing
            Plots.plot!(p, colorbar_title = colorbar_label)
        end
    end

    return p
end

"""
    choropleth(polygons, values; mode=:plots, kwargs...)
    choropleth(values, polygons; mode=:plots, kwargs...)

Creates a high-quality choropleth map from spatial `polygons` and aligned numeric `values`.
Supports automated percentile-based scaling, divergent center-zero palettes, custom themes,
and interactive Leaflet HTML dashboard generation via `mode=:leaflet`.

# Inputs
- `polygons`: Collection of polygon vertex sequences.
- `values`: Numerical vector aligned with `polygons`.

# Keyword Arguments
- `mode::Symbol`: Backend mode (`:plots` for static vector graphic, `:leaflet` for HTML).
- `cmap::Symbol`: Colormap symbol (e.g. `:viridis`, `:plasma`, `:RdBu`).
- `clims::Tuple`: Explicit value limits `(vmin, vmax)`.
- `center_zero::Bool`: When `true`, symmetrically centers colormap about 0.
- `show_colorbar::Bool`: Whether to render a value colorbar.
- `title::String`: Plot title.

# Outputs
- `Plots.Plot` or `LeafletMap`.
"""
function choropleth(
    polygons::AbstractVector{<:AbstractVector},
    values::AbstractVector{<:Real};
    mode::Symbol = :plots,
    kwargs...
)
    if mode == :leaflet || mode == :html
        return leaflet_choropleth(polygons, values; kwargs...)
    end
    return _choropleth_impl(polygons, values; kwargs...)
end

function choropleth(
    values::AbstractVector{<:Real},
    polygons::AbstractVector{<:AbstractVector};
    mode::Symbol = :plots,
    kwargs...
)
    if mode == :leaflet || mode == :html
        return leaflet_choropleth(polygons, values; kwargs...)
    end
    return _choropleth_impl(polygons, values; kwargs...)
end

# Canonical spatial domain grid and continuous surface interpolation routines are provided by Tethys
const _build_spatial_domain_grid = build_spatial_domain_grid
const _reconstruct_spatial_surface = reconstruct_spatial_surface

"""
    plot_spatial_surface(coords, values; grid_res=150, ext_factor=0.05,
                         cmap=:viridis, show_points=true, show_contours=false,
                         title="Continuous Spatial Random Field",
                         colorbar_title=nothing, polygons=nothing,
                         point_color=:white, point_size=2.5, point_alpha=0.6,
                         clims=nothing, theme_kwargs=NamedTuple(), kwargs...)

Renders a smooth continuous random surface over the spatial domain as a high-resolution
heatmap, with optional contour lines, boundary polygon overlays, and observation markers.
"""
function plot_spatial_surface(
    coords::AbstractMatrix{<:Real},
    values::AbstractVector{<:Real};
    grid_res::Int = 150,
    ext_factor::Real = 0.05,
    cmap = :viridis,
    show_points::Bool = true,
    show_contours::Bool = false,
    title::String = "Continuous Spatial Random Field",
    colorbar_title::Union{String, Nothing} = nothing,
    polygons = nothing,
    point_color = :white,
    point_size::Real = 2.5,
    point_alpha::Real = 0.6,
    clims = nothing,
    theme_kwargs = NamedTuple(),
    kwargs...
)
    c_mat = Matrix{Float64}(coords)
    v_vec = vec(collect(Float64, values))

    recon = _reconstruct_spatial_surface(
        c_mat, v_vec; grid_res = grid_res, ext_factor = ext_factor
    )
    grid_x, grid_y, Z = recon.grid_x, recon.grid_y, recon.surface

    c_limits = if !isnothing(clims)
        clims
    else
        valid_z = filter(!isnan, Z)
        !isempty(valid_z) ? (minimum(valid_z), maximum(valid_z)) : (0.0, 1.0)
    end

    cs_sym = cmap isa Symbol ? cmap : :viridis

    p = Plots.heatmap(
        grid_x, grid_y, Z';
        c = cs_sym,
        clims = c_limits,
        title = title,
        xlabel = "Coordinate X",
        ylabel = "Coordinate Y",
        aspect_ratio = :equal,
        colorbar = true,
        colorbar_title = colorbar_title,
        theme_kwargs...,
        kwargs...
    )

    if show_contours
        Plots.contour!(p, grid_x, grid_y, Z', color = :black, alpha = 0.35, lw = 0.6)
    end

    if !isnothing(polygons)
        for poly in polygons
            if length(poly) >= 3
                xs, ys = _poly_xy(poly)
                if !isempty(xs)
                    Plots.plot!(
                        p, xs, ys,
                        linecolor = :black,
                        lw = 0.8,
                        alpha = 0.7,
                        label = nothing
                    )
                end
            end
        end
    end

    if show_points && size(c_mat, 1) > 0
        Plots.scatter!(
            p, c_mat[:, 1], c_mat[:, 2],
            markersize = point_size,
            markercolor = point_color,
            markeralpha = point_alpha,
            markerstrokewidth = 0.4,
            markerstrokecolor = :black,
            label = nothing
        )
    end

    return p
end

function plot_spatial_surface(
    coords::AbstractVector,
    values::AbstractVector{<:Real};
    kwargs...
)
    xs, ys = _centroids_xy(coords)
    return plot_spatial_surface([xs ys], values; kwargs...)
end

function plot_spatial_surface(
    df::DataFrame,
    values::AbstractVector{<:Real};
    x::Union{Symbol, AbstractString, Nothing} = nothing,
    y::Union{Symbol, AbstractString, Nothing} = nothing,
    kwargs...
)
    col_x, col_y = _resolve_df_coords(df; x = x, y = y)
    coords_mat = Matrix{Float64}(df[!, [col_x, col_y]])
    return plot_spatial_surface(coords_mat, values; kwargs...)
end

"""
    _resolve_df_coords(df; x=nothing, y=nothing)

Internal coordinate resolver checking DataFrame column names.
"""
function _resolve_df_coords(df::DataFrame; x = nothing, y = nothing)
    if !isnothing(x) && !isnothing(y)
        return Symbol(x), Symbol(y)
    end
    for (cx, cy) in [
        (:s_x, :s_y), (:plon, :plat), (:lon, :lat), (:longitude, :latitude),
        (:easting, :northing), (:x, :y)
    ]
        if hasproperty(df, cx) && hasproperty(df, cy)
            return cx, cy
        end
    end
    error("Could not automatically detect 2D spatial coordinate columns in DataFrame.")
end

# =============================================================================
# Section 3: Spatial Partitioning & Network Graph Visualizations
# =============================================================================

"""
    spatial_graph_plot(centroids, g; polygons=nothing, hull_coords=nothing, pts=nothing,
                       node_size=3, node_color=:black, edge_color=:red, edge_alpha=0.6,
                       title="Spatial Partitioning", theme_kwargs=NamedTuple(), mode=:plots,
                       kwargs...)

Plots spatial partitioning structures: polygons, adjacency graph edges, centroids,
boundary hulls, and raw data points. Supports `mode=:plots` (Plots.jl) or
`mode=:leaflet` (interactive HTML).
"""
function spatial_graph_plot(au::NamedTuple, g; kwargs...)
    polys = _extract_polygons(au)
    cents = _extract_centroids(au)
    return spatial_graph_plot(cents, g; polygons = polys, kwargs...)
end

function spatial_graph_plot(
    centroids,
    g;
    polygons = nothing,
    hull_coords = nothing,
    pts = nothing,
    node_size = 3,
    node_color = :black,
    edge_color = :red,
    edge_alpha = 0.6,
    title::String = "Spatial Partitioning",
    theme_kwargs = NamedTuple(),
    mode::Symbol = :plots,
    kwargs...
)
    # If first argument is vector of polygons and centroids is in kwargs
    if !isempty(centroids) && (centroids[1] isa AbstractVector) && length(centroids[1]) > 2 && haskey(kwargs, :centroids)
        actual_polys = centroids
        actual_cents = kwargs[:centroids]
        filtered_kwargs = pairs(Base.structdiff((; kwargs...), (; centroids = nothing)))
        return spatial_graph_plot(
            actual_cents, g;
            polygons = actual_polys,
            hull_coords = hull_coords,
            pts = pts,
            node_size = node_size,
            node_color = node_color,
            edge_color = edge_color,
            edge_alpha = edge_alpha,
            title = title,
            theme_kwargs = theme_kwargs,
            mode = mode,
            filtered_kwargs...
        )
    end

    if mode == :leaflet || mode == :html
        return leaflet_spatial_graph(
            centroids, g;
            polygons = polygons,
            hull_coords = hull_coords,
            pts = pts,
            title = title,
            kwargs...
        )
    end
    p = Plots.plot(
        aspect_ratio = :equal,
        title = title,
        legend = false;
        theme_kwargs...,
        kwargs...
    )

    if !isnothing(polygons)
        for poly in polygons
            if length(poly) > 2
                xs, ys = _poly_xy(poly)
                if !isempty(xs) && (xs[1], ys[1]) != (xs[end], ys[end])
                    push!(xs, xs[1])
                    push!(ys, ys[1])
                end
                if !isempty(xs)
                    Plots.plot!(
                        p, xs, ys,
                        seriestype = :shape,
                        fillalpha = 0.1,
                        linecolor = :black,
                        lw = 0.5,
                        label = nothing
                    )
                end
            end
        end
    end

    if !isnothing(g)
        edges_list, _ = _extract_graph_edges_and_degrees(g, length(centroids))
        for (u, v) in edges_list
            if u <= length(centroids) && v <= length(centroids)
                p1, p2 = centroids[u], centroids[v]
                Plots.plot!(
                    p, [p1[1], p2[1]], [p1[2], p2[2]],
                    color = edge_color,
                    lw = 1.5,
                    alpha = edge_alpha,
                    label = nothing
                )
            end
        end
    end

    if !isnothing(pts)
        Plots.scatter!(
            p, [pt[1] for pt in pts], [pt[2] for pt in pts],
            markersize = 1,
            color = :gray,
            alpha = 0.3,
            label = "Points"
        )
    end

    if !isnothing(centroids) && !isempty(centroids)
        xs, ys = _centroids_xy(centroids)
        Plots.scatter!(
            p, xs, ys,
            markersize = node_size,
            color = node_color,
            markerstrokecolor = :white,
            label = "Centroids"
        )
    end

    if !isnothing(hull_coords) && length(hull_coords) > 2
        bx, by = _poly_xy(hull_coords)
        if !isempty(bx)
            Plots.plot!(p, bx, by, color = :black, lw = 2, ls = :dash, label = nothing)
        end
    end

    return p
end

function spatial_graph_plot(;
    au = nothing,
    pts = nothing,
    plot_title = "Spatial Partitioning",
    title = plot_title,
    mode::Symbol = :plots,
    kwargs...
)
    if !isnothing(au)
        if mode == :leaflet || mode == :html
            return leaflet_spatial_graph(au; pts = pts, title = title, kwargs...)
        end
        polygons = hasproperty(au, :polygons) ? au.polygons :
            (haskey(au, :polygons) ? au[:polygons] : nothing)
        centroids = hasproperty(au, :centroids) ? au.centroids :
            (haskey(au, :centroids) ? au[:centroids] : nothing)
        g = hasproperty(au, :graph) ? au.graph :
            (haskey(au, :graph) ? au[:graph] : nothing)
        hull = hasproperty(au, :hull_coords) ? au.hull_coords :
            (haskey(au, :hull_coords) ? au[:hull_coords] : nothing)
        return spatial_graph_plot(
            centroids, g;
            polygons = polygons,
            hull_coords = hull,
            pts = pts,
            title = title,
            mode = mode,
            kwargs...
        )
    else
        error("spatial_graph_plot requires either (centroids, graph) or `au=(...)`.")
    end
end

"""
    plot_partition(au; title="Spatial Partitioning", mode=:plots, kwargs...)

Canonical visualization of a spatial tessellation structure (areal units, hexagonal
mesh, boundary hulls, and adjacency graph).
"""
function plot_partition(
    au;
    title::String = "Spatial Partitioning",
    mode::Symbol = :plots,
    kwargs...
)
    return spatial_graph_plot(; au = au, title = title, mode = mode, kwargs...)
end

"""
    plot_kde_simple(s_coord_tuple; grid_res=600, title="Spatial Intensity (KDE)")
    plot_kde_simple(df::DataFrame; x=nothing, y=nothing, grid_res=600, title="Spatial Intensity (KDE)")

Generates a 2D intensity heatmap of spatial observation locations using Gaussian kernel
density evaluation.
"""
function plot_kde_simple(
    s_coord_tuple;
    grid_res::Int = 150,
    title::String = "Spatial Intensity (KDE)"
)
    xs = [Float64(p[1]) for p in s_coord_tuple]
    ys = [Float64(p[2]) for p in s_coord_tuple]
    n_pts = length(xs)

    dx = max(maximum(xs) - minimum(xs), 1e-6)
    dy = max(maximum(ys) - minimum(ys), 1e-6)
    h_x = 1.06 * std(xs) * n_pts^(-0.2)
    h_y = 1.06 * std(ys) * n_pts^(-0.2)
    h_x = max(h_x, 0.05 * dx)
    h_y = max(h_y, 0.05 * dy)

    gx = range(minimum(xs) - 0.1 * dx, maximum(xs) + 0.1 * dx, length = grid_res)
    gy = range(minimum(ys) - 0.1 * dy, maximum(ys) + 0.1 * dy, length = grid_res)

    dens = zeros(Float64, length(gx), length(gy))
    inv_2h2_x = 1.0 / (2.0 * h_x^2)
    inv_2h2_y = 1.0 / (2.0 * h_y^2)

    for (j, y_val) in enumerate(gy)
        for (i, x_val) in enumerate(gx)
            s = 0.0
            for k in 1:n_pts
                s += exp(-(x_val - xs[k])^2 * inv_2h2_x - (y_val - ys[k])^2 * inv_2h2_y)
            end
            dens[i, j] = s / (2.0 * pi * h_x * h_y * n_pts)
        end
    end

    p = Plots.heatmap(
        gx, gy, dens';
        c = :viridis,
        title = title,
        aspect_ratio = :equal,
        xlabel = "X",
        ylabel = "Y"
    )
    Plots.scatter!(
        p, xs, ys,
        markersize = 2,
        markercolor = :white,
        markeralpha = 0.5,
        label = "Observations"
    )
    return p
end

function plot_kde_simple(
    df::DataFrame;
    x = nothing,
    y = nothing,
    grid_res::Int = 150,
    title::String = "Spatial Intensity (KDE)"
)
    col_x, col_y = _resolve_df_coords(df; x = x, y = y)
    s_coord_tuple = tuple.(df[!, col_x], df[!, col_y])
    return plot_kde_simple(s_coord_tuple; grid_res = grid_res, title = title)
end

"""
    render_paths!(p::Plots.Plot, paths; au=nothing, centroids=nothing, polygons=nothing,
                  labels=nothing, color=:black, lw=1.0, markersize=2.0)

Overlays trajectory paths across spatial units or coordinates onto an existing plot.
"""
function render_paths!(
    p::Plots.Plot,
    paths;
    au = nothing,
    centroids = nothing,
    polygons = nothing,
    labels = nothing,
    color = :black,
    lw::Real = 1.0,
    markersize::Real = 2.0
)
    poly_list = if !isnothing(au) && hasproperty(au, :polygons)
        au.polygons
    elseif !isnothing(polygons)
        polygons
    else
        nothing
    end

    if !isnothing(poly_list)
        for poly in poly_list
            if length(poly) > 2
                px, py = _poly_xy(poly)
                Plots.plot!(
                    p, px, py,
                    seriestype = :shape,
                    fillalpha = 0.03,
                    linecolor = :gray,
                    linewidth = 0.3,
                    label = nothing
                )
            end
        end
    end

    cents = if !isnothing(au) && hasproperty(au, :centroids)
        au.centroids
    elseif !isnothing(centroids)
        centroids
    else
        nothing
    end

    if isa(paths, AbstractMatrix)
        if isnothing(cents)
            error("When `paths` is an integer unit matrix, `au` or `centroids` must be provided.")
        end
        n_indiv, n_steps = size(paths)
        for i in 1:n_indiv
            xs = Float64[cents[paths[i, t]][1] for t in 1:n_steps]
            ys = Float64[cents[paths[i, t]][2] for t in 1:n_steps]
            lab = (labels === nothing || i > length(labels)) ? nothing : labels[i]
            Plots.plot!(
                p, xs, ys,
                marker = :circle,
                markersize = markersize,
                lw = lw,
                label = lab,
                color = color
            )
        end
        return p
    end

    if isa(paths, AbstractVector) && !isempty(paths) && !isa(paths[1], Number) && isa(paths[1], AbstractVector)
        for (i, path) in enumerate(paths)
            xs = [Float64(pt[1]) for pt in path]
            ys = [Float64(pt[2]) for pt in path]
            lab = (labels === nothing || i > length(labels)) ? nothing : labels[i]
            Plots.plot!(
                p, xs, ys,
                marker = :circle,
                markersize = markersize,
                lw = lw,
                label = lab,
                color = color
            )
        end
        return p
    end
    error("render_paths! expects paths as Matrix of unit indices or Vector of Vector of (x,y) points.")
end

"""
    map_point_occupancy(polygons, centroids, id::Integer, step::Integer=1;
                        highlight_color=:red, title="")

Highlights the specific spatial polygon unit occupied by index `id` on a map.
"""
function map_point_occupancy(
    polygons,
    centroids,
    id::Integer,
    step::Integer = 1;
    highlight_color = :red,
    title::String = ""
)
    title_str = isempty(title) ? "Occupancy (Unit $id, Step $step)" : title
    p = Plots.plot(aspect_ratio = :equal, legend = false, title = title_str)

    if !isnothing(polygons)
        for poly in polygons
            if length(poly) > 2
                xs, ys = _poly_xy(poly)
                if !isempty(xs)
                    Plots.plot!(
                        p, xs, ys,
                        seriestype = :shape,
                        fillcolor = :white,
                        linecolor = :gray,
                        lw = 0.3,
                        alpha = 0.6,
                        label = nothing
                    )
                end
            end
        end
    end

    if !isnothing(polygons) && 1 <= id <= length(polygons)
        poly = polygons[id]
        if length(poly) > 2
            xs, ys = _poly_xy(poly)
            if !isempty(xs)
                Plots.plot!(
                    p, xs, ys,
                    seriestype = :shape,
                    fillcolor = highlight_color,
                    linecolor = :black,
                    lw = 0.8,
                    alpha = 0.5,
                    label = nothing
                )
            end
        end
    end

    if !isnothing(centroids) && !isempty(centroids)
        xs, ys = _centroids_xy(centroids)
        Plots.scatter!(p, xs, ys, markersize = 3, color = :black, label = nothing)
    end

    return p
end

function map_individual_occupancy(tracts, id::Integer, step::Integer, au_context)
    polys = get(au_context, :polygons, nothing)
    cents = get(au_context, :centroids, nothing)
    unit_idx = tracts[id, step]
    return map_point_occupancy(
        polys, cents, unit_idx, step;
        title = "Occupancy of Mark $id at Step $step"
    )
end

# =============================================================================
# Section 4: Specialized Spatial & Movement Fields
# =============================================================================

"""
    plot_hexagonal_field(values, au; mode=:plots, cmap=:viridis, title="Hexagonal Field",
                         colorbar_label=nothing, clims=nothing, center_zero=false, kwargs...)

Renders a spatial scalar field across polygonal areal units, specifically designed
for regular hexagonal lattices (e.g. constructed via `build_hex_mesh_planar`).
"""
function plot_hexagonal_field(
    values::AbstractVector{<:Real},
    au;
    mode::Symbol = :plots,
    cmap::Symbol = :viridis,
    title::String = "Hexagonal Spatial Field",
    colorbar_label::Union{Nothing, String} = nothing,
    clims = nothing,
    center_zero::Bool = false,
    kwargs...
)
    polys = _extract_polygons(au)
    n_p = length(polys)
    n_v = length(values)
    if n_p != n_v
        throw(DimensionMismatch(
            "Dimension mismatch: length(values) ($n_v) != length(polygons) ($n_p)."
        ))
    end

    if mode == :leaflet || mode == :html
        return leaflet_choropleth(
            polys, values;
            cmap = cmap,
            title = title,
            clims = clims,
            kwargs...
        )
    else
        return _choropleth_impl(
            polys, values;
            cmap = cmap,
            title = title,
            colorbar_label = colorbar_label,
            clims = clims,
            center_zero = center_zero,
            kwargs...
        )
    end
end

"""
    plot_hsi_choropleth(hsi, au; title="Habitat Suitability Index (HSI)",
                        cmap=:viridis, show_centroids=false, show_hull=true,
                        show_colorbar=true, border_color=:gray, lw=0.4, kwargs...)

Renders a spatial polygon choropleth map of the Habitat Suitability Index (HSI)
over the tessellated spatial domain.
"""
function plot_hsi_choropleth(
    hsi::AbstractVector{<:Real},
    au;
    mode::Symbol = :plots,
    title::String = "Habitat Suitability Index (HSI)",
    cmap::Symbol = :viridis,
    show_centroids::Bool = false,
    show_hull::Bool = true,
    show_colorbar::Bool = true,
    border_color = :gray,
    lw::Real = 0.4,
    kwargs...
)
    if mode == :leaflet || mode == :html
        return leaflet_hsi_map(
            hsi, au;
            title = title,
            cmap = cmap,
            show_centroids = show_centroids,
            show_hull = show_hull,
            kwargs...
        )
    end
    polys = _extract_polygons(au)
    p = _choropleth_impl(
        polys, hsi;
        title = title,
        cmap = cmap,
        show_colorbar = show_colorbar,
        colorbar_label = "HSI",
        border_color = border_color,
        lw = lw,
        kwargs...
    )

    if show_centroids
        cents = _extract_centroids(au)
        xs = [c[1] for c in cents]
        ys = [c[2] for c in cents]
        Plots.scatter!(
            p, xs, ys,
            markersize = 2.5,
            color = :black,
            markerstrokecolor = :white,
            markerstrokewidth = 0.5,
            label = nothing
        )
    end

    if show_hull && hasproperty(au, :hull_coords) && !isnothing(au.hull_coords)
        bx, by = _poly_xy(au.hull_coords)
        if !isempty(bx)
            Plots.plot!(p, bx, by, color = :black, lw = 1.5, ls = :dash, label = nothing)
        end
    end

    return p
end

"""
    plot_diffusion_map(diffusion, au; title="Spatial Diffusion / Dispersal Field",
                       cmap=:plasma, show_colorbar=true, border_color=:gray, lw=0.4, kwargs...)

Renders a polygon choropleth map representing the spatial distribution of
dispersal or diffusion rates across spatial units.
"""
function plot_diffusion_map(
    diffusion::Union{Real, AbstractVector{<:Real}},
    au;
    mode::Symbol = :plots,
    title::String = "Spatial Diffusion / Dispersal Field",
    cmap::Symbol = :plasma,
    show_colorbar::Bool = true,
    border_color = :gray,
    lw::Real = 0.4,
    kwargs...
)
    if mode == :leaflet || mode == :html
        return leaflet_diffusion_map(
            diffusion, au;
            title = title,
            cmap = cmap,
            kwargs...
        )
    end
    polys = _extract_polygons(au)
    S = length(polys)
    diff_vec = diffusion isa Real ? fill(Float64(diffusion), S) : Float64.(diffusion)

    return _choropleth_impl(
        polys, diff_vec;
        title = title,
        cmap = cmap,
        show_colorbar = show_colorbar,
        colorbar_label = "Diffusion (D)",
        border_color = border_color,
        lw = lw,
        kwargs...
    )
end

"""
    plot_residence_time_map(Gamma, au; title="Stationary Residence Distribution (π)",
                            cmap=:inferno, show_colorbar=true, kwargs...)

Computes and maps the stationary distribution ``\\pi`` satisfying ``\\pi \\Gamma = \\pi``,
reflecting long-term spatial occupancy and residence probabilities under a Markov transition kernel.
"""
function plot_residence_time_map(
    Gamma::AbstractMatrix{<:Real},
    au;
    mode::Symbol = :plots,
    title::String = "Stationary Residence Distribution (π)",
    cmap::Symbol = :inferno,
    show_colorbar::Bool = true,
    kwargs...
)
    if mode == :leaflet || mode == :html
        return leaflet_residence_time_map(
            Gamma, au;
            title = title,
            cmap = cmap,
            kwargs...
        )
    end
    S = size(Gamma, 1)
    pi_vec = fill(1.0 / S, S)
    for _ in 1:250
        pi_next = vec(pi_vec' * Gamma)
        s_sum = sum(pi_next)
        s_sum > 1e-12 && (pi_next ./= s_sum)
        pi_vec = pi_next
    end

    polys = _extract_polygons(au)
    return _choropleth_impl(
        polys, pi_vec;
        title = title,
        cmap = cmap,
        show_colorbar = show_colorbar,
        colorbar_label = "Stationary π",
        kwargs...
    )
end

"""
    plot_advection_arrows(au; hsi=nothing, velocity=1.0, Gamma=nothing,
                          relationship=:exponential, arrow_scale=1.0,
                          arrow_color=:white, background=:hsi, cmap=:viridis,
                          title="Advection Drift & Velocity Field",
                          min_speed_quantile=0.05, kwargs...)

Generates a spatial vector field (quiver plot) depicting the magnitude and
direction of directed movement or advective drift across spatial units.
"""
function plot_advection_arrows(
    au;
    mode::Symbol = :plots,
    hsi::Union{Nothing, AbstractVector{<:Real}} = nothing,
    velocity::Real = 1.0,
    Gamma::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
    relationship::Symbol = :exponential,
    arrow_scale::Real = 1.0,
    arrow_color = :white,
    background::Symbol = :hsi,
    cmap::Symbol = :viridis,
    title::String = "Advection Drift & Velocity Field",
    min_speed_quantile::Real = 0.05,
    kwargs...
)
    if mode == :leaflet || mode == :html
        col_str = arrow_color isa Symbol ? string(arrow_color) : string(arrow_color)
        return leaflet_advection_arrows(
            au;
            hsi = hsi,
            velocity = velocity,
            Gamma = Gamma,
            relationship = relationship,
            arrow_scale = arrow_scale,
            arrow_color = col_str,
            background = background,
            cmap = cmap,
            title = title,
            min_speed_quantile = min_speed_quantile,
            kwargs...
        )
    end

    cents = _extract_centroids(au)
    S = length(cents)
    vx = zeros(Float64, S)
    vy = zeros(Float64, S)

    if !isnothing(hsi) && hasproperty(au, :W) && !isnothing(au.W)
        W = au.W
        rows = SparseArrays.rowvals(W)
        vals = SparseArrays.nonzeros(W)
        hsi_vec = Float64.(hsi)

        for i in 1:S
            ci = cents[i]
            for j_idx in SparseArrays.nzrange(W, i)
                j = rows[j_idx]
                i == j && continue
                cj = cents[j]
                dh = hsi_vec[j] - hsi_vec[i]
                if dh > 0.0
                    dx = cj[1] - ci[1]
                    dy = cj[2] - ci[2]
                    dist = sqrt(dx^2 + dy^2) + 1e-9
                    weight = if relationship == :exponential
                        vals[j_idx] * exp(dh)
                    elseif relationship == :logistic
                        vals[j_idx] / (1.0 + exp(-4.0 * dh))
                    else
                        vals[j_idx] * dh
                    end
                    vx[i] += weight * (dx / dist)
                    vy[i] += weight * (dy / dist)
                end
            end
        end
        vx .*= Float64(velocity)
        vy .*= Float64(velocity)
    elseif !isnothing(Gamma)
        for i in 1:S
            ci = cents[i]
            for j in 1:S
                i == j && continue
                cj = cents[j]
                p_ij = Gamma[i, j]
                vx[i] += p_ij * (cj[1] - ci[1])
                vy[i] += p_ij * (cj[2] - ci[2])
            end
        end
    end

    p = if background == :hsi && !isnothing(hsi)
        plot_hsi_choropleth(
            hsi, au;
            title = title,
            cmap = cmap,
            border_color = :gray40,
            lw = 0.3,
            show_centroids = false,
            kwargs...
        )
    elseif background == :polygons
        polys = _extract_polygons(au)
        p_base = Plots.plot(aspect_ratio = :equal, title = title, legend = false; kwargs...)
        for poly in polys
            if length(poly) > 2
                xs, ys = _poly_xy(poly)
                Plots.plot!(
                    p_base, xs, ys,
                    seriestype = :shape,
                    fillcolor = :white,
                    linecolor = :gray60,
                    lw = 0.4,
                    label = nothing
                )
            end
        end
        p_base
    else
        Plots.plot(aspect_ratio = :equal, title = title, legend = false; kwargs...)
    end

    speeds = sqrt.(vx.^2 .+ vy.^2)
    max_sp = maximum(speeds)

    if max_sp > 1e-12
        dists = Float64[]
        for i in 1:min(S, 20)
            ci = cents[i]
            nn_d = minimum([
                sqrt((cj[1] - ci[1])^2 + (cj[2] - ci[2])^2)
                for (j, cj) in enumerate(cents) if j != i
            ])
            push!(dists, nn_d)
        end
        avg_spacing = !isempty(dists) ? mean(dists) : 10.0
        target_max_arrow = avg_spacing * 0.75 * arrow_scale
        scale_fac = target_max_arrow / max_sp

        if min_speed_quantile < 0.0 || min_speed_quantile > 1.0
            throw(ArgumentError(
                "min_speed_quantile must reside within [0.0, 1.0], got $(min_speed_quantile)"
            ))
        end
        q_val = min(Float64(min_speed_quantile), 0.5)
        speed_cutoff = quantile(speeds, q_val)

        xs_arr = Float64[]
        ys_arr = Float64[]
        qx = Float64[]
        qy = Float64[]

        for i in 1:S
            if speeds[i] >= speed_cutoff
                push!(xs_arr, cents[i][1])
                push!(ys_arr, cents[i][2])
                push!(qx, vx[i] * scale_fac)
                push!(qy, vy[i] * scale_fac)
            end
        end

        if !isempty(xs_arr)
            Plots.quiver!(
                p, xs_arr, ys_arr,
                quiver = (qx, qy),
                color = arrow_color,
                lw = 1.5,
                label = nothing
            )
            Plots.scatter!(
                p, xs_arr, ys_arr,
                markersize = 1.8,
                color = arrow_color,
                label = nothing
            )
        end
    end

    return p
end

"""
    plot_advection_arrows(au, vx::AbstractVector{<:Real}, vy::AbstractVector{<:Real};
                          arrow_scale=1.0, arrow_color=:navy, title="Velocity Vectors",
                          background=:polygons, kwargs...)

Renders directed velocity arrows `(vx, vy)` located at unit centroids over spatial units `au`.
"""
function plot_advection_arrows(
    au,
    vx::AbstractVector{<:Real},
    vy::AbstractVector{<:Real};
    arrow_scale::Real = 1.0,
    arrow_color = :navy,
    title::String = "Velocity Vectors",
    background::Symbol = :polygons,
    kwargs...
)
    polys = _extract_polygons(au)
    cents = haskey(kwargs, :centroids) ? kwargs[:centroids] : _extract_centroids(au)
    S = length(cents)
    @assert length(vx) == S && length(vy) == S "length(vx) and length(vy) must equal length(centroids) ($S)"

    p_base = if background == :polygons
        p_b = Plots.plot(aspect_ratio = :equal, title = title, legend = false; kwargs...)
        for poly in polys
            if length(poly) > 2
                xs, ys = _poly_xy(poly)
                Plots.plot!(
                    p_b, xs, ys,
                    seriestype = :shape,
                    fillcolor = :white,
                    linecolor = :gray60,
                    lw = 0.4,
                    label = nothing
                )
            end
        end
        p_b
    else
        Plots.plot(aspect_ratio = :equal, title = title, legend = false; kwargs...)
    end

    speeds = sqrt.(vx.^2 .+ vy.^2)
    max_sp = maximum(speeds)
    if max_sp > 1e-12
        dists = Float64[]
        for i in 1:min(S, 20)
            ci = cents[i]
            nn_d = minimum([
                sqrt((cj[1] - ci[1])^2 + (cj[2] - ci[2])^2)
                for (j, cj) in enumerate(cents) if j != i
            ])
            push!(dists, nn_d)
        end
        avg_spacing = !isempty(dists) ? mean(dists) : 10.0
        target_max_arrow = avg_spacing * 0.75 * arrow_scale
        scale_fac = target_max_arrow / max_sp

        xs_arr = Float64[]
        ys_arr = Float64[]
        us_arr = Float64[]
        vs_arr = Float64[]
        for i in 1:S
            ci = cents[i]
            push!(xs_arr, float(ci[1]))
            push!(ys_arr, float(ci[2]))
            push!(us_arr, float(vx[i]) * scale_fac)
            push!(vs_arr, float(vy[i]) * scale_fac)
        end

        Plots.quiver!(
            p_base, xs_arr, ys_arr,
            quiver = (us_arr, vs_arr),
            color = arrow_color,
            lw = 1.4
        )
    end
    return p_base
end

const plot_velocity_field = plot_advection_arrows

# =============================================================================
# Section 5: Hydrodynamic Diagnostics
# =============================================================================

"""
    plot_hydrodynamic_stratification(hydro_data, au; depth_level=1, mode=:plots,
                                     title="Ocean Stratification Diagnostics", kwargs...)

Generates multi-panel hydrographic stratification diagnostics, visualizing the
squared Brunt-Väisälä buoyancy frequency (\$N^2\$), practical salinity (\$S\$), and
representative vertical water-column profiles across depth levels.
"""
function plot_hydrodynamic_stratification(
    hydro_data::NamedTuple,
    au;
    depth_level::Integer = 1,
    mode::Symbol = :plots,
    title::String = "Ocean Stratification Diagnostics",
    kwargs...
)
    if mode == :leaflet || mode == :html
        return leaflet_hydrodynamic_dashboard(hydro_data, au; title = title, kwargs...)
    end

    z_coords = hasproperty(hydro_data, :depth_levels) ?
        collect(Float64, hydro_data.depth_levels) : [0.0]
    n_z = length(z_coords)

    if depth_level < 1 || depth_level > n_z
        throw(BoundsError(z_coords, depth_level))
    end
    z_val = z_coords[depth_level]

    n2_mat = hydro_data.N2
    sal_mat = hydro_data.salinity

    n2_slice = ndims(n2_mat) >= 2 ? vec(n2_mat[:, depth_level]) : vec(n2_mat)
    sal_slice = ndims(sal_mat) >= 2 ? vec(sal_mat[:, depth_level]) : vec(sal_mat)

    p_a = plot_hexagonal_field(
        n2_slice, au;
        title = "(A) Buoyancy Freq N² (z=$(round(z_val, digits=1))m)",
        cmap = :viridis,
        colorbar_label = "N² (s⁻²)"
    )

    p_b = plot_hexagonal_field(
        sal_slice, au;
        title = "(B) Salinity S (z=$(round(z_val, digits=1))m)",
        cmap = :haline,
        colorbar_label = "Salinity (PSU)"
    )

    p_c = Plots.plot(
        title = "(C) Mean Water-Column Profile",
        xlabel = "Value",
        ylabel = "Depth (m)",
        yflip = true,
        legend = :bottomright
    )

    if ndims(n2_mat) >= 2 && n_z > 1
        mean_n2 = [Statistics.mean(filter(!isnan, n2_mat[:, k])) for k in 1:n_z]
        mean_sal = [Statistics.mean(filter(!isnan, sal_mat[:, k])) for k in 1:n_z]

        Plots.plot!(
            p_c, mean_n2 .* 1e4, z_coords,
            lw = 2.4, color = :darkblue, label = "N² × 10⁴ (s⁻²)"
        )
        Plots.plot!(
            p_c, mean_sal, z_coords,
            lw = 2.2, color = :darkgreen, linestyle = :dash, label = "Salinity (PSU)"
        )
    else
        Plots.scatter!(
            p_c, [Statistics.mean(filter(!isnan, n2_slice)) * 1e4], [z_val],
            color = :darkblue, markersize = 6, label = "N² × 10⁴"
        )
    end

    layout_spec = Plots.@layout [a{0.38w} b{0.38w} c{0.24w}]
    return Plots.plot(
        p_a, p_b, p_c;
        layout = layout_spec,
        size = (1250, 440),
        plot_title = title,
        kwargs...
    )
end

"""
    plot_hydrodynamic_diffusion(hydro_data, au; depth_level=1, mode=:plots,
                                title="Turbulent Eddy Diffusivity & Viscosity", kwargs...)

Visualizes turbulent vertical eddy diffusivity \$\\kappa_v\$ (\$m^2/s\$) and eddy viscosity
\$\\nu_v\$ (\$m^2/s\$) across hexagonal units and vertical water-column profiles.
"""
function plot_hydrodynamic_diffusion(
    hydro_data::NamedTuple,
    au;
    depth_level::Integer = 1,
    mode::Symbol = :plots,
    title::String = "Turbulent Eddy Diffusivity & Viscosity",
    kwargs...
)
    if mode == :leaflet || mode == :html
        return leaflet_hydrodynamic_dashboard(hydro_data, au; title = title, kwargs...)
    end

    z_coords = hasproperty(hydro_data, :depth_levels) ?
        collect(Float64, hydro_data.depth_levels) : [0.0]
    n_z = length(z_coords)

    if depth_level < 1 || depth_level > n_z
        throw(BoundsError(z_coords, depth_level))
    end
    z_val = z_coords[depth_level]

    kv_mat = hydro_data.kappa_v
    nv_mat = hydro_data.nu_v

    kv_slice = ndims(kv_mat) >= 2 ? vec(kv_mat[:, depth_level]) : vec(kv_mat)
    nv_slice = ndims(nv_mat) >= 2 ? vec(nv_mat[:, depth_level]) : vec(nv_mat)

    p_a = plot_hexagonal_field(
        kv_slice, au;
        title = "(A) Eddy Diffusivity κ_v (z=$(round(z_val, digits=1))m)",
        cmap = :plasma,
        colorbar_label = "κ_v (m²/s)"
    )

    p_b = plot_hexagonal_field(
        nv_slice, au;
        title = "(B) Eddy Viscosity ν_v (z=$(round(z_val, digits=1))m)",
        cmap = :cividis,
        colorbar_label = "ν_v (m²/s)"
    )

    p_c = Plots.plot(
        title = "(C) Vertical Diffusion Profile",
        xlabel = "log₁₀(Coefficient) [m²/s]",
        ylabel = "Depth (m)",
        yflip = true,
        legend = :bottomright
    )

    if ndims(kv_mat) >= 2 && n_z > 1
        mean_kv = [
            Statistics.mean(filter(x -> !isnan(x) && x > 0, kv_mat[:, k]))
            for k in 1:n_z
        ]
        mean_nv = [
            Statistics.mean(filter(x -> !isnan(x) && x > 0, nv_mat[:, k]))
            for k in 1:n_z
        ]

        Plots.plot!(
            p_c, log10.(max.(mean_kv, 1e-8)), z_coords,
            lw = 2.4, color = :purple, label = "log₁₀(κ_v)"
        )
        Plots.plot!(
            p_c, log10.(max.(mean_nv, 1e-8)), z_coords,
            lw = 2.2, color = :teal, linestyle = :dash, label = "log₁₀(ν_v)"
        )
    else
        val_kv = Statistics.mean(filter(x -> !isnan(x) && x > 0, kv_slice))
        Plots.scatter!(
            p_c, [log10(max(val_kv, 1e-8))], [z_val],
            color = :purple, markersize = 6, label = "log₁₀(κ_v)"
        )
    end

    layout_spec = Plots.@layout [a{0.38w} b{0.38w} c{0.24w}]
    return Plots.plot(
        p_a, p_b, p_c;
        layout = layout_spec,
        size = (1250, 440),
        plot_title = title,
        kwargs...
    )
end

"""
    plot_hydrodynamic_section(hydro_data, au; var=:temperature, coordinate=:lon,
                              slice_val=nothing, tolerance=nothing, cmap=nothing,
                              title=nothing, kwargs...)

Renders a 2D vertical hydrographic cross-section (depth \$z\$ versus horizontal position)
along a transect corridor, overlaid with seafloor bathymetry bedrock masking.
"""
function plot_hydrodynamic_section(
    hydro_data::NamedTuple,
    au;
    var::Symbol = :temperature,
    coordinate::Symbol = :lon,
    along::Union{Nothing, Symbol} = nothing,
    z_coords = nothing,
    slice_val::Union{Nothing, Real} = nothing,
    tolerance::Union{Nothing, Real} = nothing,
    cmap = nothing,
    title = nothing,
    kwargs...
)
    cents = _extract_centroids(au)
    lons = [Float64(c[1]) for c in cents]
    lats = [Float64(c[2]) for c in cents]

    eff_coord = !isnothing(along) ? ((along == :lon || along == :longitude) ? :lat : :lon) : coordinate
    fixed = eff_coord == :lon ? lons : lats
    varying = eff_coord == :lon ? lats : lons
    along_label = eff_coord == :lon ? "Latitude (°N)" : "Longitude (°E)"

    target_val = !isnothing(slice_val) ? Float64(slice_val) : Statistics.median(fixed)

    tol = if !isnothing(tolerance)
        Float64(tolerance)
    else
        u_fixed = sort(unique(fixed))
        spacing = length(u_fixed) > 1 ? Statistics.median(diff(u_fixed)) : 0.2
        max(0.04, 1.25 * spacing)
    end

    matched_idx = findall(abs.(fixed .- target_val) .<= tol)
    if isempty(matched_idx)
        matched_idx = [argmin(abs.(fixed .- target_val))]
    end

    perm = sortperm(varying[matched_idx])
    matched_idx = matched_idx[perm]
    x_coords = varying[matched_idx]

    field_mat = getproperty(hydro_data, var)

    z_coords_vec = if !isnothing(z_coords)
        collect(Float64, z_coords)
    elseif hasproperty(hydro_data, :depth_levels)
        collect(Float64, hydro_data.depth_levels)
    elseif hasproperty(hydro_data, :depths_vec)
        collect(Float64, hydro_data.depths_vec)
    elseif ndims(field_mat) >= 2
        collect(Float64, range(0.0, 500.0, length = size(field_mat, 2)))
    else
        [0.0, 50.0, 100.0, 200.0]
    end
    n_z = length(z_coords_vec)
    n_stn = length(matched_idx)

    sec_data = Matrix{Float64}(undef, n_stn, n_z)

    for (si, s_idx) in enumerate(matched_idx)
        for (zi, _) in enumerate(z_coords_vec)
            val = ndims(field_mat) >= 2 ? field_mat[s_idx, zi] : field_mat[s_idx]
            sec_data[si, zi] = val
        end
    end

    cs = if !isnothing(cmap)
        cmap
    elseif var in [:temp, :temperature, :bottom_temp]
        :thermal
    elseif var in [:sal, :salinity]
        :haline
    elseif var in [:N2, :buoyancy]
        :viridis
    else
        :plasma
    end

    fig_title = !isnothing(title) ? title :
        "Vertical Transect: $var along $(eff_coord)=$(round(target_val, digits=2))"

    p = Plots.heatmap(
        x_coords, z_coords_vec, sec_data';
        c = cs,
        title = fig_title,
        xlabel = along_label,
        ylabel = "Depth (m)",
        yflip = true,
        colorbar = true,
        colorbar_title = string(var),
        size = (900, 480),
        kwargs...
    )

    seafloor_vals = if hasproperty(hydro_data, :seafloor)
        Float64.(hydro_data.seafloor[matched_idx])
    elseif hasproperty(hydro_data, :depths)
        Float64.(hydro_data.depths[matched_idx])
    else
        nothing
    end

    if !isnothing(seafloor_vals)
        Plots.plot!(
            p, x_coords, seafloor_vals,
            color = :black,
            lw = 2.5,
            label = "Seafloor Bedrock"
        )
    end

    return p
end

# =============================================================================
# Section 6: Timeseries & Export Primitives
# =============================================================================

"""
    timeseries_ci(x, mean, lower, upper; color=:royalblue, title="", xlabel="", ylabel="",
                  theme_kwargs=NamedTuple(), kwargs...)

Plots a timeseries or 1D effect with a shaded credible/confidence interval ribbon.
If `x` is empty or nothing, defaults to integer indices `1:length(mean)`.
"""
function timeseries_ci(
    x,
    mean,
    lower,
    upper;
    color = :royalblue,
    title::String = "",
    xlabel::String = "",
    ylabel::String = "",
    theme_kwargs = NamedTuple(),
    kwargs...
)
    meanv = collect(Float64, mean)
    lowv = collect(Float64, lower)
    upv = collect(Float64, upper)

    xvals = (isnothing(x) || isempty(x)) ? collect(1:length(meanv)) : x

    return Plots.plot(
        xvals, meanv,
        ribbon = (meanv .- lowv, upv .- meanv),
        color = color,
        lw = 2,
        fillalpha = 0.2,
        legend = false,
        title = title,
        xlabel = xlabel,
        ylabel = ylabel;
        kwargs...,
        theme_kwargs...
    )
end

"""
    save_plot(p::Union{Plots.Plot, LeafletMap}, path::AbstractString; fmt=nothing, dpi::Integer=150)

Saves plot `p` (either a `Plots.Plot` figure or an interactive `LeafletMap` HTML document)
to file path `path`. Creates parent directories automatically if needed.
"""
function save_plot(
    p::Union{Plots.Plot, LeafletMap},
    path::AbstractString;
    fmt = nothing,
    dpi::Integer = 150
)
    if p isa LeafletMap
        return save_html(p, path)
    end
    dir = dirname(path)
    if !isempty(dir)
        try
            isdir(dir) || mkpath(dir)
        catch
        end
    end
    out_path = path
    if !isnothing(fmt)
        ext = "." * lowercase(string(lstrip(string(fmt), '.')))
        if !endswith(lowercase(out_path), ext)
            out_path = string(out_path, ext)
        end
    end
    Plots.savefig(p, out_path)
    return out_path
end

"""
    save_plots(plots_obj, save_dir::AbstractString; prefix::String="", fmt::String="png",
               dpi::Integer=150)

Saves all plots in a `NamedTuple`, `Dict`, or container to `save_dir`.
Automatically saves `Plots.Plot` objects as `fmt` and `LeafletMap` objects as `.html`.
"""
function save_plots(
    plots_obj,
    save_dir::AbstractString;
    prefix::String = "",
    fmt::String = "png",
    dpi::Integer = 150
)
    if !isdir(save_dir)
        mkpath(save_dir)
    end

    plots_to_save = hasproperty(plots_obj, :plots) ? plots_obj.plots : plots_obj
    saved_paths = String[]
    clean_fmt = lowercase(string(lstrip(fmt, '.')))

    for (k, p) in pairs(plots_to_save)
        if p isa LeafletMap
            fn = isempty(prefix) ? "$(k).html" : "$(prefix)_$(k).html"
            out_path = joinpath(save_dir, fn)
            save_html(p, out_path)
            push!(saved_paths, out_path)
        elseif p isa Plots.Plot
            fn = isempty(prefix) ? "$(k).$(clean_fmt)" : "$(prefix)_$(k).$(clean_fmt)"
            out_path = joinpath(save_dir, fn)
            save_plot(p, out_path; fmt = nothing, dpi = dpi)
            push!(saved_paths, out_path)
        elseif p isa NamedTuple || p isa Dict
            for (sub_k, sub_p) in pairs(p)
                if sub_p isa LeafletMap
                    fn = isempty(prefix) ? "$(k)_$(sub_k).html" : "$(prefix)_$(k)_$(sub_k).html"
                    out_path = joinpath(save_dir, fn)
                    save_html(sub_p, out_path)
                    push!(saved_paths, out_path)
                elseif sub_p isa Plots.Plot
                    fn = isempty(prefix) ? "$(k)_$(sub_k).$(clean_fmt)" : "$(prefix)_$(k)_$(sub_k).$(clean_fmt)"
                    out_path = joinpath(save_dir, fn)
                    save_plot(sub_p, out_path; fmt = nothing, dpi = dpi)
                    push!(saved_paths, out_path)
                end
            end
        end
    end
    return saved_paths
end
