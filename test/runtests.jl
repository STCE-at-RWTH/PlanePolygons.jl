using LinearAlgebra
using StaticArrays
using Test

using PlanePolygons

@testset "PlanePolygons.jl" begin
    p1 = Point(0.0, 0.0)
    p2 = Point(0.0, 1.0)
    p3 = Point(1.0, 1.0)
    p4 = Point(1.0, 0.0)
    p5 = Point(0.5, 0.5)
    # Points are SVectors... let's assume that math works and we don't need to test it
    # same with Vecs
    @testset "Point/Neighborhood" begin
        # are points points
        @test is_in_neighborhood(p1, p1)
        @test !is_in_neighborhood(p1, p2)
    end

    @testset "Vec/Parallel" begin
        # are vectors vectors
        @test vectors_parallel(Vec(1.0, 1.0), Vec(2.0, 2.0))
        @test vectors_parallel(Vec(1.0, 2.0), Vec(-1.0, -2.0))
    end

    l1 = Line(p1, p3)
    l2 = Line(p2, Vec(1.0, -1.0))
    l3 = Line(Point(-2.0, 0.0), p4)

    @testset "Line" begin
        @test point_on_line(l1, p1)

        # are normals normal?
        @test right_normal(l1) ⋅ l1.dir ≈ 0
        @test left_normal(l1) ⋅ l1.dir ≈ 0

        @test point_in_right_half_plane(l2, Point(-1.0, -1.0))
        @test point_in_left_half_plane(l2, Point(1.0, 1.0))

        @test point_in_right_half_plane(l1, Point(0.0, 0.0))
        @test point_in_left_half_plane(l1, Point(0.0, 0.0))

        @test point_on_line(l1, Point(0.0, 0.0))
        @test point_on_line(l1, Point(10.0, 10.0))
        @test point_on_line(l3, Point(-0.75, 0.0))

        l4 = Line(Point(40.0, 40.0), Vec(-1.0, -1.0))
        @test lines_coincident(l1, l4)
        @test !lines_coincident(l1, l2)

        @testset "Intersection" begin
            # do lines intersect?
            A = line_intersect(l1, l2)
            @test all(!isnan, A)
            @test is_in_neighborhood(p5, A)
            @test all(isnan, line_intersect(Line(p1, p5), Line(p2, p3)))
        end
        
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
        tri_low = ClosedPolygon([p1, p2, p4])
        tri_up = ClosedPolygon([p2, p4, p1])
        poly1b = ClosedPolygon([p4, p1, p2])
        @test polygons_equal(tri_low, tri_up)
        @test polygons_equal(tri_low, poly1b)
        @test polygons_equal(tri_up, poly1b)

        poly2 = ClosedPolygon([p1, p3, p4])
        poly3 = ClosedPolygon([p1, p2, p3])

        @test num_vertices(tri_low) == 3
        @test poly_area(tri_low) ≈ 0.5
        @test point_inside(tri_low, Point(0.25, 0.25))
        @test are_polygons_intersecting(tri_low, poly2)
        pts = (
            begin
                x = 0.02 + 0.95 * rand(Float64)
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
        big_square = SClosedPolygon(
            Point(0.0, 0.0),
            Point(0.0, 2.0),
            Point(2.0, 2.0),
            Point(2.0, 0.0),
        )
        big_octagon = SClosedPolygon(
            Point(0.0, 0.0),
            Point(-1.0, 1.0),
            Point(-1.0, 2.0),
            Point(0.0, 3.0),
            Point(2.0, 3.0),
            Point(3.0, 2.0),
            Point(3.0, 1.0),
            Point(2.0, 0.0),
        )
        far_away = map(p -> p + Vec(100.0, 100.0), big_square.pts) |> SClosedPolygon
        tri_1 = SClosedPolygon(Point(0.0, 0.0), Point(2.0, 0.0), Point(0, -1.0))
        tri_1_eps = SClosedPolygon(map(pt -> pt + Vec(0.0, 1.0e-10), tri_1.pts))
        tri_2 = SClosedPolygon(Point(0.0, 1.0), Point(-1.0, 0.0), Point(-1.0, 2))
        tri_2_eps = map(pt -> pt + Vec(1.0e-6, 0.0), tri_2.pts) |> SClosedPolygon
        tri_3 = SClosedPolygon(Point(0.0, 0.0), Point(-1.0, 0.0), Point(-0.5, -1.0))
        tall_rectangle = SClosedPolygon(
            Point(0.5, -0.5),
            Point(0.5, 0.5),
            Point(1.0, 1.0),
            Point(1.0, -0.5),
        )
        big_tall_expected = SClosedPolygon(
            Point(0.5, 0.0),
            Point(0.5, 0.5),
            Point(1.0, 1.0),
            Point(1.0, 0.0),
        )
        @testset "Easy Cases" begin
            @test are_polygons_intersecting(big_square, tall_rectangle)
            @test are_polygons_intersecting(tall_rectangle, big_square)
            @test !are_polygons_intersecting(big_square, far_away)
            @test !are_polygons_intersecting(far_away, big_square)
            @test are_polygons_intersecting(big_square, big_octagon)
            @test are_polygons_intersecting(big_octagon, big_octagon)
            small = poly_intersection(big_square, tall_rectangle)
            @test polygons_equal(big_tall_expected, small)
        end
        @testset "Share Edge" begin
            @test !are_polygons_intersecting(big_square, tri_1)
            @test !are_polygons_intersecting(big_octagon, tri_1)
            @test are_polygons_intersecting(big_square, tri_1_eps)
        end

        @testset "Share Point" begin
            @test !are_polygons_intersecting(big_square, tri_2)
            @test !are_polygons_intersecting(tri_2, big_square)
            @test are_polygons_intersecting(big_square, tri_2_eps)
            @test !are_polygons_intersecting(tri_1, tri_3)
            @test !are_polygons_intersecting(tri_3, big_square)
            @test !are_polygons_intersecting(tri_3, big_octagon)
        end
    end
end

@testset "MooncakeExt" begin
    using Mooncake
    poly =
        SClosedPolygon(Point(0.0, 0.0), Point(0.0, 1.0), Point(1.0, 1.0), Point(1.0, 0.0))
    cache = Mooncake.prepare_gradient_cache(poly_area, poly)
    (res, tan) = Mooncake.value_and_gradient!!(cache, poly_area, poly)
    @test tan[1] == Mooncake.NoTangent()
    s = unpack_polygon_tangent(tan[2])
    @test res == 1.0
    @test length(s) == num_vertices(poly)
    @test s[1] == [-0.5, -0.5]
    @test s[2] == [-0.5, 0.5]
    @test s[3] == [0.5, 0.5]
    @test s[4] == [0.5, -0.5]
end
