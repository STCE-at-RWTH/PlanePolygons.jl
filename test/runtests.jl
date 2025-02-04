using PlanePolygons
using StaticArrays
using Test

@testset verbose = true "PlanePolygons.jl" begin
    p1 = Point(0.0, 0.0)
    p2 = Point(0.0, 1.0)
    p3 = Point(1.0, 1.0)
    p4 = Point(1.0, 0.0)
    p5 = Point(0.5, 0.5)

    l1 = Line(p1, p3)
    l2 = Line(p2, Vec(1.0, -1.0))

    @testset "Point" begin
        @test is_in_neighborhood(p1, p1)
        @test !is_in_neighborhood(p1, p2)
    end

    @testset "Vec" begin
        @test vectors_parallel(Vec(1.0, 1.0), Vec(2.0, 2.0))
        @test vectors_parallel(Vec(1.0, 2.0), Vec(-1.0, -2.0))
    end

    @testset "Line" begin
        @test_throws ArgumentError Line(p1, p1)
        @test point_on_line(l1, p1)
        A = line_intersect(l1, l2)
        @test !isnothing(A)
        @test is_in_neighborhood(p5, A)
        @test isnothing(line_intersect(Line(p1, p5), Line(p2, p3)))
    end

    @testset "SPoly" begin
        poly1 = SClosedPolygon(p1, p2, p4)
        poly2 = SClosedPolygon(p1, p3, p4)
        @test num_vertices(poly1) == 3
        @test poly_area(poly1) ≈ 0.5
        @test point_inside(poly1, Point(0.25, 0.25))
        @test are_polygons_intersecting(poly1, poly2)
    end

    @testset "Poly" begin
        using PlanePolygons: make_closed!
        poly1 = make_closed!([p1, p2, p4])
        poly2 = make_closed!([p1, p3, p4])
        @test num_vertices(poly1) == 3
        @test poly_area(poly1) ≈ 0.5
        @test point_inside(poly1, Point(0.25, 0.25))
        @test are_polygons_intersecting(poly1, poly2)
    end

    @testset "Cutting" begin
        poly = SClosedPolygon(p1, p2, p3, p4)
        poly1 = SClosedPolygon(p1, p3, p4)
        poly2 = cut_poly_with_line(poly, l1)
    end
end
