using LinearAlgebra
using StaticArrays
using Test


using PlanePolygons

@testset verbose = true "PlanePolygons.jl" begin
    p1 = Point(0.0, 0.0)
    p2 = Point(0.0, 1.0)
    p3 = Point(1.0, 1.0)
    p4 = Point(1.0, 0.0)
    p5 = Point(0.5, 0.5)

    l1 = Line(p1, p3)
    l2 = Line(p2, Vec(1.0, -1.0))
    l3 = Line(Point(-2., 0.), p4)

    @testset "Point" begin
        # are points points
        @test is_in_neighborhood(p1, p1)
        @test !is_in_neighborhood(p1, p2)
    end

    @testset "Vec" begin
        # are vectors vectors
        @test vectors_parallel(Vec(1.0, 1.0), Vec(2.0, 2.0))
        @test vectors_parallel(Vec(1.0, 2.0), Vec(-1.0, -2.0))
    end

    @testset "Line" begin
        @test_throws ArgumentError Line(p1, p1)
        @test point_on_line(l1, p1)
        # do lines intersect?
        A = line_intersect(l1, l2)
        @test !isnothing(A)
        @test is_in_neighborhood(p5, A)
        @test isnothing(line_intersect(Line(p1, p5), Line(p2, p3)))
        
        # are normals normal?
        @test right_normal(l1)⋅l1.dir ≈ 0
        @test left_normal(l1)⋅l1.dir ≈ 0

        @test point_in_right_half_plane(l2, Point(-1., -1.))
        @test point_in_left_half_plane(l2, Point(1., 1.))

        @test point_in_right_half_plane(l1, Point(0., 0.))
        @test point_in_left_half_plane(l1, Point(0., 0.))

        @test point_on_line(l1, Point(0.0, 0.0))
        @test point_on_line(l1, Point(10., 10.))
        @test point_on_line(l3, Point(-0.75, 0.))
    end

    @testset "SPoly" begin
        tri_low = SClosedPolygon(p1, p2, p4)
        tri_up = SClosedPolygon(p4, p1, p2)
        poly2 = SClosedPolygon(p1, p3, p4)
        poly3 = SClosedPolygon(p1, p2, p5)
        poly4 = SClosedPolygon(p1, p2, p3, p4)
        @test num_vertices(tri_low) == 3
        @test num_vertices(poly4) == 4
        @test !polygons_equal(tri_low, poly4)
        @test polygons_equal(tri_low, tri_up)
        @test !polygons_equal(tri_low, poly2)
        @test poly_area(tri_low) ≈ 0.5
        @test point_inside(tri_low, Point(0.25, 0.25))
        @test are_polygons_intersecting(tri_low, poly2)
    end

    @testset "Poly" begin
        using PlanePolygons: make_closed!
        
        tri_low = make_closed!([p1, p2, p4])
        tri_up = make_closed!([p2, p4, p1])
        poly1b = make_closed!([p4, p1, p2])
        @test polygons_equal(tri_low, tri_up)
        @test polygons_equal(tri_low, poly1b)
        @test polygons_equal(tri_up, poly1b)

        poly2 = make_closed!([p1, p3, p4])
        poly3 = make_closed!([p1, p2, p3])
        
        @test num_vertices(tri_low) == 3
        @test poly_area(tri_low) ≈ 0.5
        @test point_inside(tri_low, Point(0.25, 0.25))
        @test are_polygons_intersecting(tri_low, poly2)
        pts = (
            begin
                x = 0.02 + 0.95*rand(Float64)
                y = rand(Float64) * x
                Point(x, y)
            end for i = 1:1000
        )
        @test all(Base.Fix1(point_inside, poly2), pts)
        @test !any(Base.Fix1(point_inside, poly3), pts)
        @test all(Base.Fix1(point_inside, poly2), edge_starts(poly2))
    end


    @testset "Cutting" begin
        sqr = SClosedPolygon(p1, p2, p3, p4)
        tri_low = SClosedPolygon(p1, p3, p4)
        tri_up = SClosedPolygon(p1, p2, p3)
        cut_poly_1 = cut_poly_with_line(sqr, l1)
        reverse_l1 = Line(l1.p, -1 * l1.dir)
        @test polygons_equal(tri_low, cut_poly_1)
        poly2l = cut_poly_with_line(sqr, reverse_l1)
        @test polygons_equal(tri_up, poly2l)
    end

    @testset "Intersection" begin
        big = SClosedPolygon(Point(0.,0.), Point(0., 2.), Point(2., 2.), Point(2., 0.))
        low = SClosedPolygon(Point(0., 0.), Point(0., 2.), Point(0, -1.0))
        low_eps = SClosedPolygon(map(pt->pt+Vec(0.0, 1.0e-1), low.pts))
        @show low_eps
        tall = SClosedPolygon(Point(0.5, -0.5), Point(0.5, 0.5), Point(1.0, 1.0), Point(1.0, -0.5))
        res = SClosedPolygon(Point(0.5, 0.), Point(0.5, 0.5), Point(1.0, 1.0), Point(1.0, 0.))
        @test are_polygons_intersecting(big, tall)
        @test !are_polygons_intersecting(big, low)
        @test are_polygons_intersecting(big, low_eps)
        small = poly_intersection(big, tall)
        @test polygons_equal(res, small)
    end
end
