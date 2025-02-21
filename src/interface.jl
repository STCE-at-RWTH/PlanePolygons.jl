using LinearAlgebra
using StaticArrays

"""

Inserted by extensions if AD is available.
"""
function unpack_polygon_tangent end

function is_in_neighborhood(p0::Point, p1::Point)
    return isapprox(norm(p1 - p0), 0; atol = _HOW_CLOSE_IS_TOO_CLOSE)
end

"""
    vectors_parallel(u::Vec, v::Vec)

Test if vectors `u` and `v` are parallel`.
"""
function vectors_parallel(u::Vec, v::Vec)
    scalarprod2 = (u ⋅ v)^2
    norm2prod = (u[1]^2 + u[2]^2) * (v[1]^2 + v[2]^2)
    return isapprox(scalarprod2, norm2prod; atol = _HOW_CLOSE_IS_TOO_CLOSE) ||
           isapprox(scalarprod2, -1 * norm2prod; atol = _HOW_CLOSE_IS_TOO_CLOSE)
end

"""
    point_on_line(ℓ, point; atol=1.0e-12)

Test if the point `pt` is on the line `ℓ`.
"""
function is_other_point_on_line(ℓ, pt)
    return vectors_parallel(pt - point_on(ℓ), direction_of(ℓ))
end

"""
    projected_component(ℓ, pt)

What is the projected length from the base point of `ℓ` to the point `pt`?
"""
function projected_component(ℓ, pt)
    v = pt - point_on(ℓ)
    return (v ⋅ direction_of(ℓ)) / (direction_of(ℓ) ⋅ direction_of(ℓ))
end

"""
    project_point_onto(ℓ, pt)

Project point `pt` onto line `ℓ`.
"""
function project_point_onto(ℓ, pt)
    return point_on(ℓ) + direction_of(ℓ) * projected_component(ℓ, pt)
end

"""
    lines_parallel(ℓ1, ℓ2; atol)

Test if two lines are parallel.
"""
function lines_parallel(ℓ1, ℓ2;)
    return vectors_parallel(ℓ1.dir, ℓ2.dir)
end

"""
    lines_coincident(ℓ1, ℓ2)

Test if two lines are coincendent.
"""
function lines_coincident(ℓ1, ℓ2)
    return is_other_point_on_line(ℓ1, point_on(ℓ2)) &&
           vectors_parallel(direction_of(ℓ1), direction_of(ℓ2))
end

"""
    right_normal(ℓ)

Get a vector normal to `ℓ` that points into its right half-plane.
"""
function right_normal(ℓ)
    dir = direction_of(ℓ)
    return Vec{eltype(dir)}(dir[2], -dir[1])
end

"""
    left_normal(ℓ)

Get a vector normal to `ℓ` that points into its left half-plane.
"""
function left_normal(ℓ)
    dir = direction_of(ℓ)
    return Vec{eltype(dir)}(-dir[2], dir[1])
end

"""
    line_intersect(ℓ_1, ℓ_2)

Computes the point of intersection between the lines ``\\ell_1`` and ``\\ell_2``.

Returns `PlanePolygons._POINT_DOES_NOT_EXIST` if there is no intersection. May throw a SingularException.
"""
function line_intersect(ℓ1, ℓ2)
    if vectors_parallel(direction_of(ℓ1), direction_of(ℓ2))
        return _POINT_DOES_NOT_EXIST(eltype(point_on(ℓ1)))
    end
    d = point_on(ℓ2) - point_on(ℓ1)
    A = hcat(direction_of(ℓ1), -1 * direction_of(ℓ2))
    (_, s) = A \ d
    return point_on(ℓ2) + s * direction_of(ℓ2)
end

"""
    line_intersect(p, q, u, v; atol=1.0e-12); atol=atol

end
Computes the point of intersection between the lines
``
    \\vec p + t\\vec q
``
and
``
    \\vec u + s\\vec v
``.
"""
function line_intersect(p::Point, q::Vec, u::Point, v::Vec;)
    return line_intersect(Line(p, q), Line(u, v))
end

"""
Helper function. Computes the LHS & RHS of `E⋅n = P⋅n`
"""
function _right_half_plane(ℓ, pt)
    v1 = right_normal(ℓ) ⋅ pt
    v2 = right_normal(ℓ) ⋅ point_on(ℓ)
    return (v1, v2)
end

"""
    point_in_right_half_plane(ℓ, pt)

Test if the point `point` is to the right of the hyperplane defined by the line `ℓ`.
"""
function point_in_right_half_plane(ℓ, pt)
    (v1, v2) = _right_half_plane(ℓ, pt)
    return v1 > v2 || isapprox(v1 - v2, 0; atol = _HOW_CLOSE_IS_TOO_CLOSE)
end

function point_in_right_half_plane_strict(ℓ, pt)
    (v1, v2) = _right_half_plane(ℓ, pt)
    return v1 > v2 && !isapprox(v1 - v2, 0; atol = _HOW_CLOSE_IS_TOO_CLOSE)
end

"""
    point_in_left_half_plane(ℓ, pt)

Test if the point `point` is to the left of the hyperplane defined by the line `ℓ`.
"""
function point_in_left_half_plane(ℓ, pt;)
    return (is_other_point_on_line(ℓ, pt) || !point_in_right_half_plane(ℓ, pt))
end

function point_in_left_half_plane_strict(ℓ, pt;)
    return !point_in_right_half_plane(ℓ, pt)
end

"""
    edge_directions(poly)

Returns an iterator of `Vec`s that are parellel to the corresponding edges of `poly`.
"""
function edge_directions(poly)
    return Iterators.map(zip(edge_starts(poly), edge_ends(poly))) do (p1, p2)
        return p2 - p1
    end
end

function edge_directions(poly::SVector{TWONV,T}) where {TWONV,T}
    return edge_ends(poly) - edge_starts(poly)
end

function edge_directions(poly::SClosedPolygon{NV,T}) where {NV,T}
    return edge_ends(poly) - edge_starts(poly)
end

"""
    edge_lines(poly)

Returns an iterator of `Line`s that coincide with the edges of `poly`. 
The interior of `poly` is to the right of every `Line` in this iterator.
"""
function edge_lines(poly::ClockwiseOrientedPolygon)
    return Iterators.map(edge_starts(poly), edge_ends(poly)) do a, b
        return Line(a, b - a)
    end
end

# default method 
function edge_lines(poly)
    return Iterators.map(edge_starts(poly), edge_ends(poly)) do a, b
        return vcat(a, b - a)
    end
end

function edge_lines(poly::SClosedPolygon{NV,T}) where {NV,T}
    return mapreduce(vcat, edge_starts(poly), edge_ends(poly)) do a, b
        return SVector(Line(a, b - a))
    end
end

function edge_lines(poly::SVector{TWONV,T}) where {TWONV,T}
    return map(edge_starts(poly), edge_ends(poly)) do a, b
        return SVector(a..., (b - a)...)
    end
end

function edge_tangents(poly)
    return Iterators.map(normalize, edge_directions(poly))
end

function edge_tangents(poly::SVector{TWONV,T}) where {TWONV,T}
    return map(normalize, edge_directions(poly))
end

function edge_tangents(poly::SClosedPolygon{NV,T}) where {NV,T}
    return map(normalize, edge_directions(poly))
end

function inward_edge_normals(poly)
    return Iterators.map(edge_tangents(poly)) do t̂
        return Point(t̂[2], -t̂[1])
    end
end

function inward_edge_normals(poly::SVector{TWONV,T}) where {TWONV,T}
    return map(edge_tangents(poly)) do t̂
        return Point(t̂[2], -t̂[1])
    end
end

function inward_edge_normals(poly::SClosedPolygon{NV,T}) where {NV,T}
    return map(edge_tangents(poly)) do t̂
        return Point(t̂[2], -t̂[1])
    end
end

function outward_edge_normals(poly)
    return Iterators.map(v -> -1 * v, inward_edge_normals(poly))
end

function outward_edge_normals(poly::SVector{TWONV,T}) where {TWONV,T}
    return map(n -> -1 * n, inward_edge_normals(poly))
end

function outward_edge_normals(poly::SClosedPolygon{NV,T}) where {NV,T}
    return map(n -> -1 * n, inward_edge_normals(poly))
end

function point_inside(poly, pt)
    return all(Base.Fix2(point_in_right_half_plane, pt), edge_lines(poly))
end

function point_inside_strict(poly, pt)
    return all(Base.Fix2(point_in_right_half_plane_strict, pt), edge_lines(poly))
end

"""
    poly_area(poly)

Computes the area of clockwise-oriented polygon `poly`.
"""
function poly_area(poly)
    twoA = zero(_numeric_dtype(poly))
    dne = _POINT_DOES_NOT_EXIST(_numeric_dtype(poly))
    if any(==(dne), edge_starts(poly))
        return twoA
    end
    for (p1, p2) ∈ zip(edge_starts(poly), edge_ends(poly))
        twoA -= p1[1] * p2[2] - p1[2] * p2[1]
    end
    return twoA / 2
end

"""
    poly_line_intersections(poly, ℓ)

Computes the intersections of each edge of polygon `poly` with the line ``\\ell``.
"""
function poly_line_intersections(poly, ℓ)
    isections = Iterators.map(edge_lines(poly)) do ℓ1
        line_intersect(ℓ, ℓ1)
    end
    return isections
end

function poly_line_intersections(poly::SVector{TWONV,T}, ℓ) where {TWONV,T}
    return map(edge_lines(poly)) do ℓ1
        line_intersect(ℓ, ℓ1)
    end
end

function poly_line_intersections(poly::SClosedPolygon{NV,T}, ℓ) where {NV,T}
    return map(edge_lines(poly)) do ℓ1
        line_intersect(ℓ, ℓ1)
    end
end

function _cut_yields_new_poly(poly, isect_pts)
    return any(isect_pts) do pt
        return pt != _POINT_DOES_NOT_EXIST(eltype(pt)) && point_inside(poly, pt)
    end
end

"""
    cut_poly_with_line(poly, ℓ)
    
Cuts polygon `poly` with the line ``\\ell``. Keeps the polygon on the _right_ side of the line.

Returns the input polygon if the line does not cut the polygon.
"""
function cut_poly_with_line(poly::ClockwiseOrientedPolygon{T}, ℓ) where {T}
    isect_pts = poly_line_intersections(poly, ℓ)
    if !_cut_yields_new_poly(poly, isect_pts)
        return ClosedPolygon(poly)
    end
    new_pts_right = Point{T}[]
    sizehint!(new_pts_right, num_vertices(poly) + 2)
    for point ∈ Iterators.flatten(zip(edge_starts(poly), isect_pts))
        all(isnan, point) && continue
        point_inside(poly, point) || continue
        if point_in_right_half_plane(ℓ, point)
            if (isempty(new_pts_right) || !is_in_neighborhood(point, last(new_pts_right)))
                # point is inside poly, on right side of the line, AND not at the end of the list
                push!(new_pts_right, point)
            end
        end
    end
    polyR = if isempty(new_pts_right)
        ClosedPolygon([Point(T(NaN), T(NaN))])
    elseif is_in_neighborhood(first(new_pts_right), last(new_pts_right))
        ClosedPolygon(new_pts_right[1:end-1])
    else
        ClosedPolygon(new_pts_right)
    end
    return polyR
end

function cut_poly_with_line(poly, ℓ)
    T = _numeric_dtype(poly)
    isect_pts = poly_line_intersections(poly, ℓ)
    if !_cut_yields_new_poly(poly, isect_pts)
        return collect(poly)
    end
    new_data = T[]
    #sizehint!(new_data, 2 * (num_vertices(poly)) + 2)
    for pt ∈ Iterators.flatten(zip(edge_starts(poly), isect_pts))
        all(isnan, pt) && continue
        point_inside(poly, pt) || continue
        if point_in_right_half_plane(ℓ, pt)
            if (
                isempty(new_data) ||
                !is_in_neighborhood(pt, Point(new_data[end-1], new_data[end]))
            )
                push!(new_data, pt...)
            end
        end
    end
    if isempty(new_data)
        push!(new_data, _POINT_DOES_NOT_EXIST(T)..., _POINT_DOES_NOT_EXIST(T)...)
    elseif is_in_neighborhood(
        Point(new_data[1], new_data[2]),
        Point(new_data[end-1], new_data[end]),
    )
        deleteat!(new_data, (lastindex(new_data) - 1, lastindex(new_data)))
    end
    return new_data
end

# 

"""
Check if `poly` is fully on the left side of `ax` and that `ax` does not coincide with any sides of `poly`

This is in order to solve the case of polygons that share an edge but have an intersection of size zero.

Returns 'true' if the above conditions are met.
"""
function _is_left_separating_axis(ax, poly)
    all_left = all(edge_starts(poly)) do pt
        return point_in_left_half_plane(ax, pt)
    end
    if !all_left
        return false
    end
    # make sure none of the edges of poly coincide with the proposed separating axis
    return all(edge_lines(poly)) do ℓ
        !lines_coincident(ℓ, ax)
    end
end

function _poly_image(ℓ, poly)
    return extrema(edge_starts(poly)) do pt
        projected_component(ℓ, pt)
    end
end

function _separating_axis(n, poly1, poly2)
    ℓ = SVector(zero(Point{eltype(n)})..., n...)
    p1_min, p1_max = _poly_image(ℓ, poly1)
    p2_min, p2_max = _poly_image(ℓ, poly2)

    # @show ℓ, (p1_L, p1_R), (p2_L, p2_R)
    if p1_max < p2_min || isapprox(p1_max, p2_min; atol = _HOW_CLOSE_IS_TOO_CLOSE)
        # largest projection of p1 is less than smallest projection of p2 
        # OR they are approx equal => axis separates
        return true
    elseif p2_max < p1_min || isapprox(p1_min, p2_max; atol = _HOW_CLOSE_IS_TOO_CLOSE)
        # largest projection of p2 is less than smallest projection of p1
        # AND they are not approx equal
        return true
    else
        # there is no gap between the projected intervals
        return false
    end
end

function are_polygons_intersecting(poly1, poly2;)
    has_separating_axis = any(
        n -> _separating_axis(n, poly1, poly2),
        Iterators.flatten((outward_edge_normals(poly1), outward_edge_normals(poly2))),
    )
    return !has_separating_axis
end

function _apply_cuts(poly1, poly2)
    res = poly1
    for ℓ ∈ edge_lines(poly2)
        res = cut_poly_with_line(res, ℓ)
    end
    return res
end

"""
    poly_intersection(poly1, poly2)

Returns the portion of `poly1` also contained by `poly2`.
"""
function poly_intersection(
    poly1::ClockwiseOrientedPolygon{T},
    poly2::ClockwiseOrientedPolygon{T},
) where {T}
    if !are_polygons_intersecting(poly1, poly2)
        return ClosedPolygon(_POINT_DOES_NOT_EXIST(T))
    end
    return _apply_cuts(poly1, poly2)
end

function poly_intersection(poly1, poly2)
    T = eltype(first(edge_starts(poly1)))
    if !are_polygons_intersecting(poly1, poly2)
        return [_POINT_DOES_NOT_EXIST(T)...]
    end
    return _apply_cuts(poly1, poly2)
end

function _renumber_edge_starts(poly, new_begin)
    e = edge_starts(poly)
    return Iterators.flatten((
        Iterators.drop(e, new_begin - 1),
        Iterators.take(e, new_begin),
    ))
end

function polygons_equal(p1, p2)
    if num_vertices(p1) != num_vertices(p2)
        return false
    end
    s1 = first(edge_starts(p1))
    idx = findfirst(pt -> is_in_neighborhood(s1, pt), edge_starts(p2))
    if isnothing(idx)
        return false
    end
    return all(zip(edge_starts(p1), _renumber_edge_starts(p2, idx))) do (p1, p2)
        return is_in_neighborhood(p1, p2)
    end
end