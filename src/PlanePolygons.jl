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
    line_intersect

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
    cut_poly_with_line

"""
    Point{T}

Alias for SVector{2, T}. Semantically represents a point in 2-D space.
"""
Point{T} = SVector{2,T}

function is_in_neighborhood(p0::Point, p1::Point; atol = 1.0e-12)
    return isapprox(norm(p1 - p0), 0; atol = atol)
end

"""
    Vec{T}

Alias for SVector{2, T}. Semantically represents a vector in 2-D space.
"""
Vec{T} = SVector{2,T}

"""
    vectors_parallel(u::Vec, v::Vec; atol=1.0e-12)

Test if vectors `u` and `v` are parallel up to `atol`.
"""
function vectors_parallel(u::Vec, v::Vec; atol = 1.0e-12)
    scalarprod = u ⋅ v
    normprod = norm(u) * norm(v)
    return isapprox(scalarprod, normprod; atol = atol) ||
           isapprox(scalarprod, -1 * normprod; atol = atol)
end

"""
    Line{T}

Represents a line in the plane that passes through `p` in direction `dir`.

Fields
---
- `p::Point{T}`
- `dir::Vec{T}`  
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

right_normal(ℓ::Line{T}) where {T} = Vec{T}(ℓ.dir[2], -ℓ.dir[1])
left_normal(ℓ::Line{T}) where {T} = Vec{T}(-ℓ.dir[2], ℓ.dir[1])

"""
    line_intersect(p, q, u, v; atol=1.0e-12)

Computes the point of intersection between the lines ``\ell_1`` and ``\ell_2``.

Returns `nothing` if there is no intersection. May throw a SingularException.
"""
function line_intersect(ℓ1, ℓ2; atol = 1.0e-12)
    if vectors_parallel(ℓ1.dir, ℓ2.dir; atol = atol)
        return nothing
    end
    d = ℓ2.p - ℓ1.p
    A = hcat(ℓ1.dir, -1 * ℓ2.dir)
    (_, s) = A \ d
    return ℓ2.p + s * ℓ2.dir
end

"""
    line_intersect(p, q, u, v; atol=1.0e-12)

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

"""
    point_in_right_half_plane(ℓ, pt; atol=1.0e-12)

Test if the point `point` is to the right of the hyperplane defined by the line `ℓ`.
"""
function point_in_right_half_plane(ℓ::Line, pt; atol = 1.0e-10)
    v1 = right_normal(ℓ) ⋅ pt
    v2 = right_normal(ℓ) ⋅ ℓ.p
    return v1 > v2 || isapprox(v1, v2; atol = atol)
end

"""
    point_in_right_half_plane(ℓ, pt; atol=1.0e-12)

Test if the point `point` is to the left of the hyperplane defined by the line `ℓ`.
"""
function point_in_left_half_plane(ℓ::Line, pt; atol = 1.0e-10)
    v1 = left_normal(ℓ) ⋅ pt
    v2 = left_normal(ℓ) ⋅ ℓ.p
    return v1 > v2 || isapprox(v1, v2; atol = atol)
end

"""
    point_on_line(p0, dir, point; atol=1.0e-12)

Test if the point `pt` is on the line `ℓ`.
"""
function point_on_line(ℓ::Line, pt; atol = 1.0e-12)
    v1 = pt - ℓ.p
    v2 = v1 ./ ℓ.dir
    return isapprox(v2[1], v2[2]; atol = atol)
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
    
Cuts polygon `poly` with the line ``\ell``.

Returns a copy of `poly` if ``\ell`` does not intersect `poly`.
"""
function cut_poly_with_line(poly::ClockwiseOrientedPolygon{T}, ℓ; atol = 1.0e-12) where {T}
    isect_pts = poly_line_intersections(poly, ℓ)
    if all(isnothing, isect_pts)
        return copy(poly)
    end
    new_pts = SVector{2,T}[]
    sizehint!(new_pts, num_vertices(poly) + 2)
    for point ∈ Iterators.flatten(zip(edge_starts(poly), isect_pts))
        isnothing(point) && continue
        (!isempty(new_pts) && point == last(new_pts)) && continue
        if (
            point_in_right_half_plane(ℓ, point; atol = atol) &&
            point_inside(poly, point; atol = atol)
        )
            push!(new_pts, point)
        end
    end
    if first(new_pts) == last(new_pts)
        return ClosedPolygon(new_pts)
    end
    return make_closed!(new_pts)
end

function are_polygons_intersecting(poly1, poly2; atol = 1.0e-12)
    for ℓ ∈ edge_lines(poly1)
        # all points in poly1 on on the right side of u + sv
        # which means there is a separating axis iff all of the points
        # in poly2 are on the left side of u+sv
        if all(pt -> !point_in_right_half_plane(ℓ, pt; atol = atol), edge_starts(poly2))
            return false
        end
    end
    return true
end

function poly_intersection(poly1, poly2; atol = 1.0e-12)
    if !are_polygons_intersecting(poly1, poly2; atol = atol)
        return nothing
    end
    res = poly1
    for ℓ ∈ edge_lines(poly2)
        res = cut_poly_with_line(res, ℓ; atol = atol)
    end
    return res
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
    return p.pts[SVector(ntuple(i->i+1, NV-1)..., 1)]
end

end
