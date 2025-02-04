# PlanePolygons

[![Build Status](https://github.com/STCE-at-RWTH/PlanePolygons.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/STCE-at-RWTH/PlanePolygons.jl/actions/workflows/CI.yml?query=branch%3Amaster)

Offers simple primitives and algorithms for plane geometry.

## Usage

```julia
using PlanePolygons

pt1 = Point(0.0, 0.0)
dir1 = Vec(1.0, 2.0)

ell1 = Line(pt1, dir1)

poly1 = SClosedPolygon(Point(0.0, 1.0), Point(1.0, 2.0), Point(2.0, 1.0))
```

The polygon interface is listed under

```julia
help> ClockwiseOrientedPolygon
```
