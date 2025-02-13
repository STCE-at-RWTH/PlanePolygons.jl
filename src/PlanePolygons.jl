module PlanePolygons
using StaticArrays
using LinearAlgebra

export Point, is_in_neighborhood
export Vec, vectors_parallel
export Line,
    right_normal,
    left_normal,
    point_in_right_half_plane,
    point_in_left_half_plane,
    point_on_line,
    lines_coincident,
    line_intersect,
    lines_parallel

export ClockwiseOrientedPolygon, SClosedPolygon, ClosedPolygon
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

"""

Inserted by extensions if AD is available.
"""
function unpack_polygon_tangent end

"""
    Point{T}

Alias for SVector{2, T}. Semantically represents a point in 2-D space.
"""
const Point{T} = SVector{2,T}

function is_in_neighborhood(p0::Point, p1::Point; atol = 1.0e-12)
    return isapprox(norm(p1 - p0), 0; atol = atol)
end

"""
    Vec{T}

Alias for SVector{2, T}. Semantically represents a vector in 2-D space.
"""
const Vec{T} = SVector{2,T}

"""
    vectors_parallel(u::Vec, v::Vec; atol=1.0e-12)

Test if vectors `u` and `v` are parallel up to `atol`.
"""
function vectors_parallel(u::Vec, v::Vec; atol = 1.0e-12)
    scalarprod2 = (u ⋅ v)^2
    norm2prod = (u[1]^2 + u[2]^2) * (v[1]^2 + v[2]^2)
    return isapprox(scalarprod2, norm2prod; atol = atol) ||
           isapprox(scalarprod2, -1 * norm2prod; atol = atol)
end

"""
    Line{T}

Represents a line in the plane that passes through `p` in direction `dir`.

Fields
---
- `p::Point{T}`
- `dir::Vec{T}`

Methods
---
- `right_normal(ℓ)`
- `left_normal(ℓ)`
- `line_intersect(ℓ1, ℓ2; atol)`
- `lines_parallel(ℓ1, ℓ2; atol)`
"""
struct Line{T}
    p::Point{T}
    dir::Vec{T}

    function Line(p::Point{T}, dir::Vec{T}) where {T}
        if isapprox(norm(dir), zero(T); atol = 1.0e-12)
            throw(ArgumentError("Cannot construct Line{T} with norm(dir) close to zero. "))
        end
        return new{T}(p, dir)
    end
end

function lines_coincident(ℓ1::Line, ℓ2::Line; atol = 1.0e-12)
    return point_on_line(ℓ1, ℓ2.p; atol = atol) &&
           vectors_parallel(ℓ1.dir, ℓ2.dir; atol = atol)
end

right_normal(ℓ::Line{T}) where {T} = Vec{T}(ℓ.dir[2], -ℓ.dir[1])
left_normal(ℓ::Line{T}) where {T} = Vec{T}(-ℓ.dir[2], ℓ.dir[1])

"""
    line_intersect(p, q, u, v; atol=1.0e-12); atol=atol

end
Computes the point of intersection between the lines ``\ell_1`` and ``\ell_2``.

Returns `Point(NaN, NaN)` if there is no intersection. May throw a SingularException.
"""
function line_intersect(ℓ1::Line{T}, ℓ2::Line{T}; atol = 1.0e-12) where {T}
    if vectors_parallel(ℓ1.dir, ℓ2.dir; atol = atol)
        return Point{T}(T(NaN), T(NaN))
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
function line_intersect(p::Point, q::Vec, u::Point, v::Vec; atol = 1.0e-12)
    return line_intersect(Line(p, q), Line(u, v); atol = atol)
end

function _right_half_plane(ℓ::Line, pt; atol = 1.0e-12)
    v1 = right_normal(ℓ) ⋅ pt
    v2 = right_normal(ℓ) ⋅ ℓ.p
    return (v1, v2)
end

"""
    point_in_right_half_plane(ℓ, pt; atol=1.0e-12)

Test if the point `point` is to the right of the hyperplane defined by the line `ℓ`.
"""
function point_in_right_half_plane(ℓ::Line, pt; atol = 1.0e-12)
    (v1, v2) = _right_half_plane(ℓ, pt; atol = atol)
    return v1 > v2 || isapprox(v1 - v2, 0; atol = atol)
end

function point_in_right_half_plane_strict(ℓ::Line, pt; atol = 1.0e-12)
    (v1, v2) = _right_half_plane(ℓ, pt; atol = atol)
    return v1 > v2 && !isapprox(v1 - v2, 0; atol = atol)
end

"""
    point_in_right_half_plane(ℓ, pt; atol=1.0e-12)

Test if the point `point` is to the left of the hyperplane defined by the line `ℓ`.
"""
function point_in_left_half_plane(ℓ::Line, pt; atol = 1.0e-12)
    return (
        point_on_line(ℓ, pt; atol = atol) || !point_in_right_half_plane(ℓ, pt; atol = atol)
    )
end

function point_in_left_half_plane_strict(ℓ::Line, pt; atol = 1.0e-12)
    return !point_in_right_half_plane(ℓ, pt; atol = atol)
end

"""
    point_on_line(ℓ, point; atol=1.0e-12)

Test if the point `pt` is on the line `ℓ`.
"""
function point_on_line(ℓ::Line, pt; atol = 1.0e-12)
    return vectors_parallel(pt - ℓ.p, ℓ.dir; atol = atol)
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
function lines_parallel(ℓ1, ℓ2; atol = 1.0e-12)
    return vectors_parallel(ℓ1.dir, ℓ2.dir; atol = atol)
end

"""
    ClockwiseOrientedPolygon{T}

A polygon, with its vertices listed in clockwise order (negative orientation).

Interface
---
- `num_vertices(poly::ClockwiseOrientedPolygon)`
- `edge_starts(poly::ClockwiseOrientedPolygon)`: Edge starting points
- `edge_ends(poly::ClockwiseOrientedPolygon)`: Edge endpoints
- `edge_directions(poly::ClockwiseOrientedPolygon)`: List of vectors parallel to the edges. Defined by default.
- `edge_lines(poly::ClockwiseOrientedPolygon)`: List of lines of which the edges of `poly` are segments.
- `edge_tangents(poly::ClockwiseOrientedPolygon)`: List of unit vectors parallel to the edges.
- `inward_edge_normals(poly::ClockwiseOrientedPolygon)`: List of inward normal unit vectors.
- `outward_edge_normals(poly::ClockwiseOrientedPolygon)`: List of outward normal unit vectors.
- `point_inside(poly::ClockwiseOrientedPolygon, point::Point; atol)`: Test if point `point` is inside `poly`. Defined by default.
- `poly_area(poly::ClockwiseOrientedPolygon)`: Area of `poly`.
- `cut_poly_with_line(poly::ClockwiseOrientedPolygon, line::Line; keep_right: true)`
- `are_polygons_intersecting(poly1::ClockwiseOrientedPolygon, poly2::ClockwiseOrientedPolygon)`
- `poly_intersection(poly1::ClockwiseOrientedPolygon, poly2::ClockwiseOrientedPolygon)`
"""
abstract type ClockwiseOrientedPolygon{T} end

edge_directions(poly::ClockwiseOrientedPolygon) = edge_ends(poly) - edge_starts(poly)

function edge_lines(poly::ClockwiseOrientedPolygon)
    return map(zip(edge_starts(poly), edge_ends(poly))) do (a, b)
        return Line(a, b - a)
    end
end

edge_tangents(poly::ClockwiseOrientedPolygon) = map(normalize, edge_directions(poly))

function inward_edge_normals(poly::ClockwiseOrientedPolygon)
    return map(edge_tangents(poly)) do t̂
        return Point(t̂[2], -t̂[1])
    end
end

function outward_edge_normals(poly::ClockwiseOrientedPolygon)
    return -1 .* inward_edge_normals(poly)
end

function point_inside(poly::ClockwiseOrientedPolygon, pt; atol = 1.0e-12)
    all(isnan, pt) && return false
    return all(edge_lines(poly)) do ℓ
        point_in_right_half_plane(ℓ, pt; atol = atol)
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
function poly_line_intersections(poly, ℓ; atol = 1.0e-12)
    isections = map(edge_lines(poly)) do ℓ1
        line_intersect(ℓ, ℓ1; atol = atol)
    end
    return isections
end

function poly_line_intersections(poly, p, q; atol = 1.0e-12)
    return poly_line_intersections(poly, Line(p, q); atol = atol)
end

"""
    cut_poly_with_line(poly, ℓ; atol)
    
Cuts polygon `poly` with the line ``\ell``. Keeps the polygon on the _right_ side of the line.

Returns the input polygon if the line does not cut the polygon.
"""
function cut_poly_with_line(poly::ClockwiseOrientedPolygon{T}, ℓ; atol = 1.0e-12) where {T}
    isect_pts = poly_line_intersections(poly, ℓ)
    yields_no_new_poly = all(isect_pts) do pt
        all(isnan, pt) && return true
        point_inside(poly, pt; atol = atol) && return false
        return true
    end
    if yields_no_new_poly
        return poly
    end
    new_pts_right = Point{T}[]
    sizehint!(new_pts_right, num_vertices(poly) + 2)
    for point ∈ Iterators.flatten(zip(edge_starts(poly), isect_pts))
        all(isnan, point) && continue
        point_inside(poly, point; atol = atol) || continue
        if point_in_right_half_plane(ℓ, point)
            if (
                isempty(new_pts_right) ||
                !is_in_neighborhood(point, last(new_pts_right); atol = atol)
            )
                # point is inside poly, on right side of the line, AND not at the end of the list
                push!(new_pts_right, point)
            end
        end
    end
    polyR = if isempty(new_pts_right)
        make_closed!([Point(T(NaN), T(NaN))])
    elseif is_in_neighborhood(first(new_pts_right), last(new_pts_right); atol = atol)
        ClosedPolygon(new_pts_right)
    else
        make_closed!(new_pts_right)
    end
    return polyR
end

# 

"""
Check if `poly` is fully on the left side of `ax` and that `ax` does not coincide with any sides of `poly`

This is in order to solve the case of polygons that share an edge but have an intersection of size zero.

Returns 'true' if the above conditions are met.
"""
function _is_left_separating_axis(ax, poly; atol = 1.0e-12)
    all_left = all(edge_starts(poly)) do pt
        return point_in_left_half_plane(ax, pt; atol = atol)
    end
    if !all_left
        return false
    end
    # make sure none of the edges of poly coincide with the proposed separating axis
    return all(edge_lines(poly)) do ℓ
        !lines_coincident(ℓ, ax; atol = atol)
    end
end

function _poly_image(ℓ::Line, poly)
    return extrema(edge_starts(poly)) do pt
        projected_component(ℓ, pt)
    end
end

function _separating_axis(n::Vec, poly1, poly2; atol = 1.0e-12)
    ℓ = Line(Point(0.0, 0.0), Point(-n[2], n[1]))
    p1_L, p1_R = _poly_image(ℓ, poly1)
    p2_L, p2_R = _poly_image(ℓ, poly2)
    # @show ℓ, (p1_L, p1_R), (p2_L, p2_R)
    if p1_R < p2_L || isapprox(p1_R, p2_L; atol = atol)
        # largest projection of p1 is less than smallest projection of p2 
        # AND they are not approx equal
        return true
    elseif p2_R < p1_L || isapprox(p1_L, p2_R; atol = atol)
        # largest projection of p2 is less than smallest projection of p1
        # AND they are not approx equal
        return true
    else
        # there is no gap between the projected intervals
        return false
    end
end

function are_polygons_intersecting(poly1, poly2; atol = 1.0e-12)
    has_separating_axis = any(
        n -> _separating_axis(n, poly1, poly2; atol = atol),
        Iterators.flatten((outward_edge_normals(poly1), outward_edge_normals(poly2))),
    )
    return !has_separating_axis
end

"""
    poly_intersection(poly1, poly2; atol=1.0e-12)

Returns the portion of `poly1` also contained by `poly2`.
"""
function poly_intersection(poly1::ClockwiseOrientedPolygon{T}, poly2::ClockwiseOrientedPolygon{T}; atol = 1.0e-12) where {T}
    if !are_polygons_intersecting(poly1, poly2; atol = atol)
        return make_closed!([Point(T(NaN), T(NaN))])
    end
    res = poly1
    for ℓ ∈ edge_lines(poly2)
        res = cut_poly_with_line(res, ℓ; atol = atol)
    end
    return res
end

function polygons_equal(
    p1::ClockwiseOrientedPolygon,
    p2::ClockwiseOrientedPolygon;
    atol = 1.0e-12,
)
    if num_vertices(p1) != num_vertices(p2)
        return false
    end
    s1 = first(edge_starts(p1))
    idx = findfirst(pt -> is_in_neighborhood(s1, pt; atol = atol), edge_starts(p2))
    if isnothing(idx)
        return false
    end
    i = 1
    for j ∈ idx:length(edge_starts(p2))
        if !is_in_neighborhood(edge_starts(p1)[i], edge_starts(p2)[j]; atol = atol)
            return false
        end
        i += 1
    end
    for j ∈ 1:(idx-1)
        if !is_in_neighborhood(edge_starts(p1)[i], edge_starts(p2)[j]; atol = atol)
            return false
        end
        i += 1
    end
    return true
end

"""
    ClosedPolygon{T}

Closed, clockwise oriented polygon.
"""
struct ClosedPolygon{T} <: ClockwiseOrientedPolygon{T}
    pts::Vector{Point{T}}
end

num_vertices(p::ClosedPolygon) = length(p.pts) - 1
edge_starts(p::ClosedPolygon) = @view p.pts[1:end-1]
edge_ends(p::ClosedPolygon) = @view p.pts[2:end]

function make_closed!(pts)
    push!(pts, pts[1])
    return ClosedPolygon(pts)
end

"""
    SClosedPolygon{NV, T}

Closed, clockwise oriented polygon with a fixed number of vertices.
"""
struct SClosedPolygon{NV,T} <: ClockwiseOrientedPolygon{T}
    pts::SVector{NV,Point{T}}
end

function SClosedPolygon(pts::Vararg{Point})
    return SClosedPolygon(SVector(pts...))
end

num_vertices(::SClosedPolygon{NV}) where {NV} = NV

function edge_starts(p::SClosedPolygon)
    return p.pts
end

function edge_ends(p::SClosedPolygon{NV}) where {NV}
    return p.pts[SVector(ntuple(i -> i + 1, NV - 1)..., 1)]
end

end
