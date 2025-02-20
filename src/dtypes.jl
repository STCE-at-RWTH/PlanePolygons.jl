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
- `point_on(ℓ)`: gets a point on `ℓ`
- `direction_of(ℓ)`: gets the direction of `ℓ`
- `right_normal(ℓ)`
- `left_normal(ℓ)`
- `line_intersect(ℓ1, ℓ2)`
- `lines_parallel(ℓ1, ℓ2)`
- `lines_coincident(ℓ1, ℓ2)`
"""
struct Line{T}
    p::Point{T}
    dir::Vec{T}
end

"""
    point_on(ℓ)

Get a point on the line `ℓ`.
"""
@inline function point_on(ℓ::Line{T}) where {T}
    return ℓ.p
end

"""
    direction_of(ℓ)

Get a vector parallel to `ℓ`.
"""
@inline function direction_of(ℓ::Line{T}) where {T}
    return ℓ.dir
end

"""
    _flatten(ℓ::Line{T}) where {T}

Flatten `ℓ` and convert it into a 4-element `SVector{T}`.
"""
function _flatten(ℓ::Line{T}) where {T}
    return SVector{4,T}(ℓ.p..., ℓ.dir...)
end

function point_on(ℓ::SVector{4,T}) where {T}
    return Point{T}(ℓ[1], ℓ[2])
end

function direction_of(ℓ::SVector{4,T}) where {T}
    return Vec{T}(ℓ[3], ℓ[4])
end

"""
    _unflatten_line(ℓ)

Unflatten `ℓ` and convert it back to a `Line{T}`.
"""
function _unflatten_line(ℓ::SVector{4,T}) where {T}
    return Line{T}(point_on(ℓ), direction_of(ℓ))
end

"""
    ClockwiseOrientedPolygon{T}

A polygon, with its vertices listed in clockwise order (negative orientation).

(Public) Interface
---
- `num_vertices(poly::ClockwiseOrientedPolygon)`: Number of vertices (and edges).
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

(Private) Interface
---
- `_numeric_dtype(poly)`, `_numeric_dtype(::Type{poly})`: Backing data type of a polygon.
"""
abstract type ClockwiseOrientedPolygon{T} end

_numeric_dtype(::ClockwiseOrientedPolygon{T}) where T = T
_numeric_dtype(::Type{ClockwiseOrientedPolygon{T}}) where T = T
_numeric_dtype(other_object) = eltype(other_object)

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

function ClosedPolygon(poly::ClockwiseOrientedPolygon{T}) where {T}
    return ClosedPolygon(collect(poly.pts))
end
 
function ClosedPolygon(poly::ClosedPolygon{T}) where {T}
    return ClosedPolygon(poly.pts)
end

function _flatten(poly::ClosedPolygon{T}) where {T}
    return copy(reinterpret(T, poly.pts))
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
        Point(data[2*i-1], data[2*i])
    end)
    return SClosedPolygon(pts)
end

function edge_starts(p::SClosedPolygon)
    return p.pts
end

function edge_ends(p::SClosedPolygon{NV,T}) where {NV,T}
    return p.pts[SVector(ntuple(i -> i + 1, NV - 1)..., 1)]
end

function _flatten(p::SClosedPolygon{NV,T}) where {NV,T}
    return reduce(vcat, p.pts)
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

##
# FLATTEN AND RECONSTRUCT
##

function _flatten(poly::MClosedPolygon{NV,T}) where {NV,T}
    return MVector{2 * NV,T}(reduce(vcat, edge_starts(poly)))
end

_unflatten_polygon(data::SVector{TWONV,T}) where {TWONV,T} = SClosedPolygon(data)
_unflatten_polygon(data::MVector{TWONV,T}) where {TWONV,T} = MClosedPolygon(data)

function _unflatten_polygon(data::Vector{T}) where {T}
    ClosedPolygon(copy(reinterpret(Point{T}, data)))
end

##
# WACKY NONSENSE
## 

function num_vertices(::SVector{TWONV,T}) where {TWONV,T}
    if isodd(TWONV)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $TWONV",
            ),
        )
    else
        return div(TWONV, 2)
    end
end

function edge_starts(poly_data::SVector{TWONV,T}) where {TWONV,T}
    if isodd(TWONV)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $TWONV",
            ),
        )
    else
        return SVector{div(TWONV, 2),Point{T}}((
            Point{T}(poly_data[i], poly_data[i+1]) for i ∈ 1:2:TWONV
        ))
    end
end

function edge_ends(poly_data::SVector{TWONV,T}) where {TWONV,T}
    if isodd(TWONV)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $TWONV",
            ),
        )
    else
        return SVector{div(TWONV, 2),Point{T}}((
            begin
                if i == TWONV - 1
                    Point{T}(poly_data[1], poly_data[2])
                else
                    Point{T}(poly_data[i+2], poly_data[i+3])
                end
            end for i ∈ 1:2:TWONV
        ))
    end
end

function num_vertices(::MVector{TWONV,T}) where {TWONV,T}
    if isodd(TWONV)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $TWONV",
            ),
        )
    else
        return div(TWONV, 2)
    end
end

function edge_starts(poly_data::MVector{TWONV,T}) where {TWONV,T}
    if isodd(TWONV)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $TWONV",
            ),
        )
    else
        return (Point{T}(poly_data[i], poly_data[i+1]) for i ∈ 1:2:TWONV)
    end
end

function edge_ends(poly_data::MVector{TWONV,T}) where {TWONV,T}
    if isodd(TWONV)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $TWONV",
            ),
        )
    else
        return (
            begin
                if i == TWONV - 1
                    Point{T}(poly_data[1], poly_data[2])
                else
                    Point{T}(poly_data[i+2], poly_data[i+3])
                end
            end for i ∈ 1:2:TWONV
        )
    end
end

function num_vertices(poly_data::Vector{T}) where {T}
    two_nv = length(poly_data)
    if isodd(two_nv)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $two_nv",
            ),
        )
    else
        return div(two_nv, 2)
    end
end

function edge_starts(poly_data::Vector{T}) where {T}
    two_nv = length(poly_data)
    if isodd(two_nv)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $two_nv",
            ),
        )
    else
        return (Point{T}(poly_data[i], poly_data[i+1]) for i ∈ 1:2:two_nv)
    end
end

function edge_ends(poly_data::Vector{T}) where {T}
    two_nv = length(poly_data)
    if isodd(two_nv)
        throw(
            ArgumentError(
                "Polygon as block array must have even number of entries... TWONV = $two_nv",
            ),
        )
    else
        return (
            begin
                if i == two_nv - 1
                    Point{T}(poly_data[1], poly_data[2])
                else
                    Point{T}(poly_data[i+2], poly_data[i+3])
                end
            end for i ∈ 1:2:two_nv
        )
    end
end