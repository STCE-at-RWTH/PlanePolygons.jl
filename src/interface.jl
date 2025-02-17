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
    vectors_parallel(u::Vec, v::Vec; atol=1.0e-12)

Test if vectors `u` and `v` are parallel up to `atol`.
"""
function vectors_parallel(u::Vec, v::Vec)
    scalarprod2 = (u ⋅ v)^2
    norm2prod = (u[1]^2 + u[2]^2) * (v[1]^2 + v[2]^2)
    return isapprox(scalarprod2, norm2prod; atol = _HOW_CLOSE_IS_TOO_CLOSE) ||
           isapprox(scalarprod2, -1 * norm2prod; atol = _HOW_CLOSE_IS_TOO_CLOSE)
end

function lines_coincident(ℓ1::Line, ℓ2::Line)
    return point_on_line(ℓ1, ℓ2.p) && vectors_parallel(ℓ1.dir, ℓ2.dir)
end

right_normal(ℓ::Line{T}) where {T} = Vec{T}(ℓ.dir[2], -ℓ.dir[1])
left_normal(ℓ::Line{T}) where {T} = Vec{T}(-ℓ.dir[2], ℓ.dir[1])

"""
    line_intersect(p, q, u, v; atol=1.0e-12); atol=atol

end
Computes the point of intersection between the lines ``\ell_1`` and ``\ell_2``.

Returns `Point(NaN, NaN)` if there is no intersection. May throw a SingularException.
"""
function line_intersect(ℓ1::Line{T}, ℓ2::Line{T};) where {T}
    if vectors_parallel(ℓ1.dir, ℓ2.dir)
        return _POINT_DOES_NOT_EXIST(T)
    end
    d = ℓ2.p - ℓ1.p
    A = hcat(ℓ1.dir, -1 * ℓ2.dir)
    (_, s) = A \ d
    return ℓ2.p + s * ℓ2.dir
end

"""
    line_intersect(p, q, u, v; atol=1.0e-12); atol=atol

end
Computes the point of intersection between the lines
``
    \vec p + t\vec q
``
and
``
    \vec u + s\vec v
``.
"""
function line_intersect(p::Point, q::Vec, u::Point, v::Vec;)
    return line_intersect(Line(p, q), Line(u, v))
end

function _right_half_plane(ℓ::Line, pt;)
    v1 = right_normal(ℓ) ⋅ pt
    v2 = right_normal(ℓ) ⋅ ℓ.p
    return (v1, v2)
end

"""
    point_in_right_half_plane(ℓ, pt; atol=1.0e-12)

Test if the point `point` is to the right of the hyperplane defined by the line `ℓ`.
"""
function point_in_right_half_plane(ℓ::Line, pt;)
    (v1, v2) = _right_half_plane(ℓ, pt)
    return v1 > v2 || isapprox(v1 - v2, 0; atol = _HOW_CLOSE_IS_TOO_CLOSE)
end

function point_in_right_half_plane_strict(ℓ::Line, pt;)
    (v1, v2) = _right_half_plane(ℓ, pt)
    return v1 > v2 && !isapprox(v1 - v2, 0; atol = _HOW_CLOSE_IS_TOO_CLOSE)
end

"""
    point_in_right_half_plane(ℓ, pt; atol=1.0e-12)

Test if the point `point` is to the left of the hyperplane defined by the line `ℓ`.
"""
function point_in_left_half_plane(ℓ::Line, pt;)
    return (point_on_line(ℓ, pt) || !point_in_right_half_plane(ℓ, pt))
end

function point_in_left_half_plane_strict(ℓ::Line, pt;)
    return !point_in_right_half_plane(ℓ, pt)
end

"""
    point_on_line(ℓ, point; atol=1.0e-12)

Test if the point `pt` is on the line `ℓ`.
"""
function point_on_line(ℓ::Line, pt;)
    return vectors_parallel(pt - ℓ.p, ℓ.dir)
end

function projected_component(ℓ::Line, pt)
    v = pt - ℓ.p
    return (v ⋅ ℓ.dir) / (ℓ.dir ⋅ ℓ.dir)
end

function project_point_onto(ℓ::Line, pt)
    return ℓ.p + ℓ.dir * projected_component(ℓ, pt)
end

"""
    lines_parallel(ℓ1, ℓ2; atol)

Test if two lines are parallel.
"""
function lines_parallel(ℓ1, ℓ2;)
    return vectors_parallel(ℓ1.dir, ℓ2.dir)
end

function edge_directions(poly::ClockwiseOrientedPolygon)
    return Iterators.map(zip(edge_starts(poly), edge_ends(poly))) do (p1, p2)
        return p2 - p1
    end
end

function edge_lines(poly::ClockwiseOrientedPolygon)
    return Iterators.map(zip(edge_starts(poly), edge_ends(poly))) do (a, b)
        return Line(a, b - a)
    end
end

edge_tangents(poly::ClockwiseOrientedPolygon) =
    Iterators.map(normalize, edge_directions(poly))

function inward_edge_normals(poly::ClockwiseOrientedPolygon)
    return Iterators.map(edge_tangents(poly)) do t̂
        return Point(t̂[2], -t̂[1])
    end
end

function outward_edge_normals(poly::ClockwiseOrientedPolygon)
    return Iterators.map(v -> -1 * v, inward_edge_normals(poly))
end

function point_inside(poly::ClockwiseOrientedPolygon, pt;)
    all(isnan, pt) && return false
    return all(edge_lines(poly)) do ℓ
        point_in_right_half_plane(ℓ, pt)
    end
end

"""
    poly_area(poly)

Computes the area of clockwise-oriented polygon `poly`.
"""
function poly_area(poly::ClockwiseOrientedPolygon{T}) where {T}
    twoA = zero(T)
    for (p1, p2) ∈ zip(edge_starts(poly), edge_ends(poly))
        twoA -= p1[1] * p2[2] - p1[2] * p2[1]
    end
    return twoA / 2
end

"""
    poly_line_intersections(poly, ℓ; atol)

Computes the intersections of each edge of polygon `poly` with the line ``\ell``.
"""
function poly_line_intersections(poly, ℓ;)
    isections = map(edge_lines(poly)) do ℓ1
        line_intersect(ℓ, ℓ1)
    end
    return isections
end

function poly_line_intersections(poly, p, q;)
    return poly_line_intersections(poly, Line(p, q))
end

"""
    cut_poly_with_line(poly, ℓ; atol)
    
Cuts polygon `poly` with the line ``\ell``. Keeps the polygon on the _right_ side of the line.

Returns the input polygon if the line does not cut the polygon.
"""
function cut_poly_with_line(poly::ClockwiseOrientedPolygon{T}, ℓ::Line{T}) where {T}
    isect_pts = poly_line_intersections(poly, ℓ)
    yields_no_new_poly = all(isect_pts) do pt
        all(isnan, pt) && return true
        point_inside(poly, pt) && return false
        return true
    end
    if yields_no_new_poly
        return ClosedPolygon(poly.pts)
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

# 

"""
Check if `poly` is fully on the left side of `ax` and that `ax` does not coincide with any sides of `poly`

This is in order to solve the case of polygons that share an edge but have an intersection of size zero.

Returns 'true' if the above conditions are met.
"""
function _is_left_separating_axis(ax, poly;)
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

function _poly_image(ℓ::Line, poly)
    return extrema(edge_starts(poly)) do pt
        projected_component(ℓ, pt)
    end
end

function _separating_axis(n::Vec, poly1, poly2;)
    ℓ = Line(Point(0.0, 0.0), Point(-n[2], n[1]))
    p1_L, p1_R = _poly_image(ℓ, poly1)
    p2_L, p2_R = _poly_image(ℓ, poly2)
    # @show ℓ, (p1_L, p1_R), (p2_L, p2_R)
    if p1_R < p2_L || isapprox(p1_R, p2_L; atol = _HOW_CLOSE_IS_TOO_CLOSE)
        # largest projection of p1 is less than smallest projection of p2 
        # AND they are not approx equal
        return true
    elseif p2_R < p1_L || isapprox(p1_L, p2_R; atol = _HOW_CLOSE_IS_TOO_CLOSE)
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

"""
    poly_intersection(poly1, poly2; atol=1.0e-12)

Returns the portion of `poly1` also contained by `poly2`.
"""
function poly_intersection(
    poly1::ClockwiseOrientedPolygon{T},
    poly2::ClockwiseOrientedPolygon{T},
) where {T}
    if !are_polygons_intersecting(poly1, poly2)
        return make_closed!([Point(T(NaN), T(NaN))])
    end
    res = poly1
    for ℓ ∈ edge_lines(poly2)
        res = cut_poly_with_line(res, ℓ)
    end
    return res
end

function _renumber_edge_starts(poly, new_begin)
    e = edge_starts(poly)
    return Iterators.flatten((
        Iterators.drop(e, new_begin - 1),
        Iterators.take(e, new_begin),
    ))
end

function polygons_equal(p1::ClockwiseOrientedPolygon, p2::ClockwiseOrientedPolygon;)
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