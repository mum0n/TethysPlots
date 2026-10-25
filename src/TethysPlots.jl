"""
    TethysPlots

Unified visualization, mapping, and geospatial diagnostics library for the Tethys
oceanographic, ecological, and spatial partitioning ecosystem.
Provides both publication-grade static graphics (via Plots.jl) and standalone
interactive GIS dashboards (via Leaflet.js).
"""
module TethysPlots

using CoordRefSystems
using Unitful
using DataFrames
using Dates
using Graphs
using LinearAlgebra
using Printf
using SparseArrays
using Statistics
using StatsBase
using ColorSchemes
using Plots
using Tethys

# Include modules
include("leaflet.jl")
include("plotting.jl")

# Static plotting exports
export choropleth,
       spatial_graph_plot,
       plot_partition,
       plot_hexagonal_field,
       plot_spatial_surface,
       plot_advection_arrows,
       plot_velocity_field,
       plot_hsi_choropleth,
       plot_diffusion_map,
       plot_residence_time_map,
       plot_hydrodynamic_stratification,
       plot_hydrodynamic_diffusion,
       plot_hydrodynamic_section,
       timeseries_ci,
       plot_kde_simple,
       render_paths!,
       map_point_occupancy,
       map_individual_occupancy,
       create_theme,
       save_plot,
       save_plots,
       poly_xy,
       centroids_xy,
       extract_polygons,
       extract_centroids,
       build_spatial_domain_grid,
       reconstruct_spatial_surface,
       extract_graph_edges_and_degrees,
       _poly_xy,
       _centroids_xy,
       _extract_polygons,
       _extract_centroids,
       _build_spatial_domain_grid,
       _reconstruct_spatial_surface,
       _extract_graph_edges_and_degrees

# Interactive Leaflet exports
export LeafletMap,
       save_html,
       leaflet_choropleth,
       leaflet_spatial_graph,
       leaflet_hsi_map,
       leaflet_diffusion_map,
       leaflet_residence_time_map,
       leaflet_advection_arrows,
       leaflet_spacetime_map,
       leaflet_current_density_map,
       leaflet_hydrodynamic_dashboard,
       leaflet_particle_tracks,
       utm_to_lonlat,
       lonlat_to_utm,
       is_geographic_coordinates,
       CoordinateTransformer,
       build_coordinate_transformer,
       transform_point,
       transform_polygon,
       extract_wkt,
       is_wgs84_wkt

end # module TethysPlots
