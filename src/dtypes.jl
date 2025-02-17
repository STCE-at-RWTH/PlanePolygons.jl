"""
    Point{T}

Alias for SVector{2, T}. Semantically represents a point in 2-D space.
"""
const Point{T} = SVector{2,T}

"""
    Vec{T}

Alias for SVector{2, T}. Semantically represents a vector in 2-D space.
"""
const Vec{T} = SVector{2,T}

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
end

"""
    ClockwiseOrientedPolygon{T}

A polygon, with its vertices listed in clockwise order (negative orientation).

Interface
---
- `num_vertices(poly::ClockwiseOrientedPolygon)`
- `edge_starts(poly::ClockwiseOrientedPolygon)`: Edge starting points
- `edge_ends(poly::ClockwiseOrientedPolygon)`: Edge endpoints
- `edge_directions(poly::ClockwiseOrientedPolygon)`: Iterable of vectors parallel to the edges. Defined by default.
- `edge_lines(poly::ClockwiseOrientedPolygon)`: Iterable of lines of which the edges of `poly` are segments.
- `edge_tangents(poly::ClockwiseOrientedPolygon)`: Iterable of unit vectors parallel to the edges.
- `inward_edge_normals(poly::ClockwiseOrientedPolygon)`: Iterable of inward normal unit vectors.
- `outward_edge_normals(poly::ClockwiseOrientedPolygon)`: Iterable of outward normal unit vectors.
- `point_inside(poly::ClockwiseOrientedPolygon, point::Point; atol)`: Test if point `point` is inside `poly`. Defined by default.
- `poly_area(poly::ClockwiseOrientedPolygon)`: Area of `poly`.
- `cut_poly_with_line(poly::ClockwiseOrientedPolygon, line::Line; keep_right: true)`
- `are_polygons_intersecting(poly1::ClockwiseOrientedPolygon, poly2::ClockwiseOrientedPolygon)`
- `poly_intersection(poly1::ClockwiseOrientedPolygon, poly2::ClockwiseOrientedPolygon)`
"""
abstract type ClockwiseOrientedPolygon{T} end

"""
    ClosedPolygon{T}

Closed, clockwise oriented polygon.
"""
struct ClosedPolygon{T} <: ClockwiseOrientedPolygon{T}
    pts::Vector{Point{T}}
end

num_vertices(p::ClosedPolygon) = length(p.pts)
edge_starts(p::ClosedPolygon) = (p.pts[i] for i = 1:num_vertices(p))
function edge_ends(p::ClosedPolygon)
    return (i == num_vertices(p) ? p.pts[1] : p.pts[i+1] for i = 1:num_vertices(p))
end

function ClosedPolygon(pts::Vararg{Point{T}}) where {T}
    return ClosedPolygon(collect(pts))
end

function ClosedPolygon(pts::AbstractArray{Point{T}}) where {T}
    return ClosedPolygon(collect(pts))
end

function ClosedPolygon(data::AbstractArray{T}) where {T}
    return ClosedPolygon(reinterpret(Point{T}, data))
end

abstract type SizedClockwiseOrientedPolygon{NV,T} <: ClockwiseOrientedPolygon{T} end

num_vertices(::SizedClockwiseOrientedPolygon{NV}) where {NV} = NV

"""
    SClosedPolygon{NV, T}

Closed, clockwise oriented polygon with a fixed number of vertices.
"""
struct SClosedPolygon{NV,T} <: SizedClockwiseOrientedPolygon{NV,T}
    pts::SVector{NV,Point{T}}
end

function SClosedPolygon(pts::Vararg{Point})
    return SClosedPolygon(SVector(pts...))
end

function SClosedPolygon(data::SVector{TWONV,T}) where {TWONV,T}
    NV = TWONV ÷ 2
    pts = SVector(ntuple(NV) do i
        Point(data[2 * i - 1], data[2*i])
    end)
    return SClosedPolygon(pts)
end

function edge_starts(p::SClosedPolygon)
    return p.pts
end

function edge_ends(p::SClosedPolygon{NV,T}) where {NV,T}
    return p.pts[SVector(ntuple(i -> i + 1, NV - 1)..., 1)]
end

struct MClosedPolygon{NV,T} <: SizedClockwiseOrientedPolygon{NV,T}
    pts::MVector{NV,Point{T}}
end

function MClosedPolygon(pts::Vararg{Point{T},N}) where {T,N}
    return MClosedPolygon{N,T}(MVector(pts...))
end

function MClosedPolygon(data::MVector{TWONV,T}) where {TWONV,T}
    NV = TWONV ÷ 2
    return MClosedPolygon(ntuple(i -> Point(data[2*i-1], data[2*i]), NV)...)
end

function edge_starts(p::MClosedPolygon{NV}) where {NV}
    return (p.pts[i] for i = 1:NV)
end

function edge_ends(p::MClosedPolygon{NV}) where {NV}
    return (i == num_vertices(p) ? p.pts[1] : p.pts[i+1] for i = 1:NV)
end

function _flatten(poly::ClockwiseOrientedPolygon{T}) where {T}
    return reduce(vcat, edge_starts(poly))
end

function _flatten(poly::MClosedPolygon{NV,T}) where {NV,T}
    return MVector{2 * NV,T}(reduce(vcat, edge_starts(poly)))
end

_reconstruct_polygon(data::SVector{TWONV,T}) where {TWONV,T} = SClosedPolygon(data)
_reconstruct_polygon(data::MVector{TWONV,T}) where {TWONV,T} = MClosedPolygon(data)

function _reconstruct_polygon(data::Vector{T}) where {T}
    ClosedPolygon(copy(reinterpret(Point{T}, data)))
end