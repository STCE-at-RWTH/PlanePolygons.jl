module MooncakeExt
    using PlanePolygons

    function PlanePolygons.unpack_polygon_tangent(tan)
        return map(tan.fields.pts.fields.data) do t
            t.fields.data |> Vec
        end
    end

    
end