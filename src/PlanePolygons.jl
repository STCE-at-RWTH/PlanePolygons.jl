module PlanePolygons
using StaticArrays
using LinearAlgebra

export Point, is_in_neighborhood
export Vec, vectors_parallel
export Line,
    point_on,
    direction_of,
    right_normal,
    left_normal,
    point_in_right_half_plane,
    point_in_left_half_plane,
    is_other_point_on_line,
    lines_coincident,
    line_intersect,
    lines_parallel

export ClockwiseOrientedPolygon, SizedClockwiseOrientedPolygon
export MClosedPolygon, SClosedPolygon, ClosedPolygon
export num_vertices,
    edge_starts,
    edge_ends,
    edge_directions,
    edge_lines,
    edge_tangents,
    inward_edge_normals,
    outward_edge_normals,
    point_inside,
    poly_area,
    poly_line_intersections,
    poly_intersection,
    are_polygons_intersecting,
    cut_poly_with_line,
    polygons_equal

export unpack_polygon_tangent

const _HOW_CLOSE_IS_TOO_CLOSE = 1.0e-12

include("dtypes.jl")

_POINT_DOES_NOT_EXIST(T) = Point{T}(T(NaN), T(NaN))

include("interface.jl")



end
