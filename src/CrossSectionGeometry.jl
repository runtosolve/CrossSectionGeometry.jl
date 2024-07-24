module CrossSectionGeometry


using LinesCurvesNodes, LinearAlgebra, StaticArrays, LazySets, Statistics


@kwdef struct ThinWalledGeometry

    L::Vector{Float64}
    θ::Vector{Float64}
    n::Vector{Int}
    r::Vector{Float64}
    n_r::Vector{Int}
    t::Float64
    centerline_location::String
    offset::Vector{Float64}

    center::Vector{Vector{Float64}}
    left::Vector{Vector{Float64}}
    right::Vector{Vector{Float64}}

end


function generate_thin_walled(L, θ, n)

    #anchor points
    cross_section_nodes = lay_out_cross_section_nodes(L, θ)

    #no corners in this method
    corners = []
    flats = generate_straight_line_segments(cross_section_nodes, corners, n)
    
    cross_section = Array{Vector{Float64}}(undef, 0)
    
    #round 
    for i in eachindex(flats)
    
        flats[i] = [round.(flats[i][j], digits=5) for j in eachindex(flats[i])]
    
    end

    #combine flats
    for i in eachindex(flats)
    
        cross_section = vcat(cross_section, flats[i])
    
    end
    
    #remove negative zeros
    for i in eachindex(cross_section)
    
        if cross_section[i][1] === -0.0
    
            cross_section[i][1] = 0.0
    
        end
    
        if cross_section[i][2] === -0.0
    
            cross_section[i][2] = 0.0
    
        end
    
    end

    #remove repeats
    cross_section = unique(cross_section[2:end])
    cross_section = pushfirst!(cross_section, [0.0, 0.0]) #for closed sections 

    return cross_section

end


function generate_thin_walled(L, θ, n, r, n_r)

    cross_section_nodes = Geometry.lay_out_cross_section_nodes(L, θ)

    corners = Geometry.generate_cross_section_rounded_corners(cross_section_nodes, r, n_r)

    flats = Geometry.generate_straight_line_segments(cross_section_nodes, corners, n)

    num_flats = size(flats)[1]
    num_corners = size(corners)[1]

    num_segments = num_flats + num_corners

    cross_section = [flats; corners]

    cross_section_index = [1:2:num_segments; 2:2:num_segments]

    cross_section_index_sorted = sortperm(cross_section_index)

    cross_section = cross_section[cross_section_index_sorted]

    cross_section = vcat(cross_section...)

    for i in eachindex(cross_section)

        # if typeof(unit(L[1])) == Unitful.FreeUnits{(), NoDims, nothing}
            cross_section[i] = round.(cross_section[i], digits=5)  #check this 
        # else
            # cross_section[i] = round.(unit(cross_section[i][1]), cross_section[i], digits=5)
        # end

        cross_section[i] = Geometry.remove_negative_zeros(cross_section[i])

    end

    cross_section = unique(cross_section)

    return cross_section

end



function lay_out_cross_section_nodes(L, θ)

    num_segments = length(L)

    cross_section = Array{Vector{Float64}}(undef, num_segments)
    # cross_section = []

    for i in eachindex(L)

        if i == 1

            start_node = [0.0, 0.0]

        else

            start_node = cross_section[i-1]

        end

        # if typeof(unit(L[1])) == Unitful.FreeUnits{(), NoDims, nothing}
            cross_section[i] = round.(LinesCurvesNodes.transform_vector(L[i], start_node, θ[i]), digits=5)
        # else
            # cross_section = [cross_section; [round.(unit(L[i]), LinesCurvesNodes.transform_vector(L[i], start_node, θ[i]), digits=5)]]
        # end

        cross_section[i] = remove_negative_zeros(cross_section[i])

    end

    # cross_section = [[[0.0*unit(L[1]), 0.0*unit(L[1])]]; cross_section] #add start node at unity

    cross_section = [[[0.0, 0.0]]; cross_section]

    return cross_section

end



function generate_cross_section_rounded_corners(cross_section_nodes, r, n)

    corners = Array{Array{Vector{Float64}}}(undef, length(r))
    # corners = Array{Array{Vector{Any}}}(undef, length(r))

    for i in eachindex(r)

        if (i == length(r)) & (cross_section_nodes[1] == cross_section_nodes[end]) #closed section

            A = cross_section_nodes[end-1] 
            B = cross_section_nodes[1]
            C = cross_section_nodes[2]

        else #open section 

            A = cross_section_nodes[i] 
            B = cross_section_nodes[i+1]
            C = cross_section_nodes[i+2]

        end
        
        corners[i] = LinesCurvesNodes.generate_fillet(A, B, C, r[i], n[i])

    end

   return corners

end



function generate_straight_line_segments(cross_section_nodes, corners, n)

    segments = Array{Vector{Vector{Float64}}}(undef, length(n))

    if corners == []

        for i in eachindex(n)

            A = cross_section_nodes[i]
            B = cross_section_nodes[i+1]

            segments[i] = LinesCurvesNodes.discretize_vector(A, B, n[i])

        end

    else 

        if cross_section_nodes[1] != cross_section_nodes[end]  #open cross-sections 

            corner_index = 1

            for i in eachindex(n)

                if i == 1

                    A = cross_section_nodes[1]
                    B = corners[1][1]


                elseif i == length(n)

                    A = corners[corner_index][end]
                    B = cross_section_nodes[end]

                else
                    A = corners[corner_index][end]
                    corner_index += 1
                    B = corners[corner_index][1]

                end

                segments[i] = LinesCurvesNodes.discretize_vector(A, B, n[i])

            end

        else  #closed cross-sections

            corner_index = 1

            for i in eachindex(n)

                if i == 1

                    A = corners[end][end]
                    B = corners[1][1]


                elseif i == length(n)

                    A = corners[corner_index][end]
                    B = corners[end][1]

                else
                    A = corners[corner_index][end]
                    corner_index += 1
                    B = corners[corner_index][1]

                end

                segments[i] = LinesCurvesNodes.discretize_vector(A, B, n[i])

            end   
            
        end

    end

    return segments

end


function calculate_node_normal(A, B, C)

    BC = C-B
    BA = A-B

    BR = BA + BC/norm(BC) * norm(BA)

    return BR

end


function calculate_cross_section_unit_node_normals(cross_section)

    num_nodes = size(cross_section)[1]
    unit_node_normals = Array{Vector{Float64}}(undef, num_nodes)

    # cross_section = [ustrip.(cross_section[i]) for i in eachindex(cross_section)]


    for i=1:num_nodes

        if i == 1

            if cross_section[end] == cross_section[1]  #closed cross-sections

                A = cross_section[end-1]
                B = cross_section[1]
                C = cross_section[2]

                node_normal = calculate_node_normal(A, B, C)

                if isapprox(abs.(node_normal), [0.0, 0.0], atol=1e-8)
        
                    BA = A-B
                    unit_node_normals[i] = [-BA[2], BA[1]] / norm(BA)

                    unit_node_normals[i] = right_halfspace_normal_correction(A, B, unit_node_normals[i]) 

                else

                    unit_node_normals[i] = -node_normal / norm(node_normal)

                    unit_node_normals[i] = right_halfspace_normal_correction(A, B, unit_node_normals[i]) 

                end

            else  #open cross-sections 

                A = cross_section[1]
                B = cross_section[2]

                BA = A-B

                unit_node_normals[i] = [-BA[2], BA[1]] / norm(BA)

                unit_node_normals[i] = right_halfspace_normal_correction(A, B, unit_node_normals[i]) 

            end

        elseif i == num_nodes

            A = cross_section[end-1]
            B = cross_section[end]

            BA = A-B

            unit_node_normals[i] = [-BA[2], BA[1]] / norm(BA)

            unit_node_normals[i] = right_halfspace_normal_correction(A, B, unit_node_normals[i]) 

        else

            A = cross_section[i-1]
            B = cross_section[i]
            C = cross_section[i+1]

            node_normal = calculate_node_normal(A, B, C)

            if isapprox(abs.(node_normal), [0.0, 0.0], atol=1e-8)
    
                BA = A-B
                unit_node_normals[i] = [-BA[2], BA[1]] / norm(BA)

                unit_node_normals[i] = right_halfspace_normal_correction(A, B, unit_node_normals[i]) 

            else

                unit_node_normals[i] = -node_normal / norm(node_normal)

                unit_node_normals[i] = right_halfspace_normal_correction(A, B, unit_node_normals[i]) 

            end

        end


        unit_node_normals[i] = remove_negative_zeros(unit_node_normals[i])

    end


    return unit_node_normals

end


function get_coords_along_node_normals(cross_section, unit_node_normals, Δ)

    normal_cross_section = Array{Vector{Float64}}(undef, 0)

    for i in eachindex(cross_section)

        push!(normal_cross_section, cross_section[i] + unit_node_normals[i] * Δ)

    end

    return normal_cross_section

end


function right_halfspace_normal_correction(A, B, normal) 

    s = LineSegment(A, B)
    hs = halfspace_right(s) 
    dA = A + 1.0 * normal
    if !∈(dA, hs)  #is the node normal pointing left or right, if it is pointing left, reverse sign to make it point right
        normal = -normal
    end

    return normal

end


function remove_negative_zeros(coord)



    if coord[1] === -0.0

        coord[1] = 0.0

    end

    if coord[2] === -0.0

        coord[2] = 0.0

    end

    return coord

end



function calculate_axis_area(element_connectivity, element_thicknesses, node_geometry, axis_location, about_axis)

    num_elem = size(element_connectivity)[1]

    A_elements = zeros(Float64, num_elem)

    for i=1:num_elem

        node_i = trunc(Int, element_connectivity[i, 1])
        node_j = trunc(Int, element_connectivity[i, 2])

        x_i = node_geometry[node_i, 1]
        x_j = node_geometry[node_j, 1]
        y_i = node_geometry[node_i, 2]
        y_j = node_geometry[node_j, 2]

        length_i = norm([x_j, y_j] - [x_i, y_i])

        A_i = length_i * element_thicknesses[i]

        if about_axis == "y"

            z_i = mean([x_i, x_j])

        elseif about_axis == "x"

            z_i = mean([y_i, y_j])

        end

        if (z_i - axis_location) > 0  #Find which side of the axis A_i is on.

            A_elements[i] = A_i

        elseif (z_i - axis_location) < 0

            A_elements[i] = -A_i
            
        end

    end

    return A_elements

end

function create_thin_walled_cross_section_geometry(L, θ, n, r, n_r, t; centerline, offset)

    #bring in surface
    cross_section = Geometry.generate_thin_walled(L, θ, n, r, n_r)

    #calculate surface normals
    unit_node_normals = Geometry.calculate_cross_section_unit_node_normals(cross_section)

    #following along surface, is centerline to right or left?
    if centerline == "to left"
        increment = -t/2
    elseif centerline == "to right"
        increment = t/2
    elseif centerline == "at center"
        increment = 0.0
    end

    #calculate centerline geometry, surfaces
    center = Geometry.get_coords_along_node_normals(cross_section, unit_node_normals, increment)
    left = Geometry.get_coords_along_node_normals(center, unit_node_normals, -t/2)
    right = Geometry.get_coords_along_node_normals(center, unit_node_normals, t/2)

    #convert geometry to tuples, offset coordinates provided by users
    # convert_to_tuple(data) = [(data[i][1] + offset[1], data[i][2] + offset[2]) for i in eachindex(data)]
    # left = convert_to_tuple(left)
    # right = convert_to_tuple(right)
    # center = convert_to_tuple(center)

    convert_to_vectors(data) = [[data[i][1] + offset[1], data[i][2] + offset[2]] for i in eachindex(data)]
    left = convert_to_vectors(left)
    right = convert_to_vectors(right)
    center = convert_to_vectors(center)

    #collection up all the section info
    section_geometry = (center=center, left=left, right=right)

    return section_geometry

end


function create_thin_walled_cross_section_geometry(L, θ, n, t; centerline, offset)

    #bring in surface
    cross_section = Geometry.generate_thin_walled(L, θ, n)

    #calculate surface normals
    unit_node_normals = Geometry.calculate_cross_section_unit_node_normals(cross_section)

    #following along surface, is centerline to right or left?
    if centerline == "to left"
        increment = -t/2
    elseif centerline == "to right"
        increment = t/2
    elseif centerline == "at center"
        increment = 0.0
    end

    #calculate centerline geometry, surfaces
    center = Geometry.get_coords_along_node_normals(cross_section, unit_node_normals, increment)
    left = Geometry.get_coords_along_node_normals(center, unit_node_normals, -t/2)
    right = Geometry.get_coords_along_node_normals(center, unit_node_normals, t/2)

    #convert geometry to tuples, offset coordinates provided by users
    # convert_to_tuple(data) = [(data[i][1] + offset[1], data[i][2] + offset[2]) for i in eachindex(data)]
    # left = convert_to_tuple(left)
    # right = convert_to_tuple(right)
    # center = convert_to_tuple(center)

    convert_to_vectors(data) = [[data[i][1] + offset[1], data[i][2] + offset[2]] for i in eachindex(data)]
    left = convert_to_vectors(left)
    right = convert_to_vectors(right)
    center = convert_to_vectors(center)

    #collection up all the section info
    section_geometry = (center=center, left=left, right=right)

    return section_geometry

end





# """
#     wshape_nodes(shape_info, n)

# Accepts the Struct `shape_info` generated using CrossSection.AISC and the discretization Vector `n` and outputs the outline x-y coordinates of a W shape 'xcoords' and 'ycoords'.

# The Vector 'n' describes the number of segments in a quarter cross-section, i.e., `n = [half of outside flange face, flange thickness, half of inside flange face, flange-web radius, half of web]`.

# """


# function wshape_nodes(shape_info, n)

#     #from bottom of bottom flange, web centerline to left edge
#     xcoords = zeros(n[1]+1)
#     ycoords = zeros(n[1]+1)

#     flange_range = 0.0 : -shape_info.bf / 2 / n[1] : -shape_info.bf / 2
#     [xcoords[i] =  flange_range[i] for i in eachindex(flange_range)]
#     ycoords .= 0.0

#     #up along bottom flange thickness
#     flange_thickness_range = shape_info.tf/n[2]:shape_info.tf/n[2]:shape_info.tf
#     xcoords = [xcoords; ones(n[2])*xcoords[end]]
#     ycoords = [ycoords; flange_thickness_range]

#     #over to fillet radius at bottom flange - web intersection

#     # flange_flat = shape_info.bf/2 - shape_info.k1
#     flange_flat = shape_info.bf/2 - shape_info.tw/2 - (shape_info.kdes - shape_info.tf)

#     inside_flange_range = (xcoords[end] + flange_flat/n[3]) : flange_flat/n[3] : (xcoords[end] + flange_flat)

#     xcoords = [xcoords; inside_flange_range]
#     ycoords = [ycoords; ones(n[3])*ycoords[end]]

#     #go around the fillet
#     radius = -xcoords[end] - shape_info.tw/2
#     θ = (-π/2 + π/2/n[4]):π/2/n[4]: 0.0

#     xo = xcoords[end]
#     yo = ycoords[end] + radius

#     x_radius = xo .+ radius .* cos.(θ)
#     y_radius = yo .+ radius .* sin.(θ)

#     # plot(x_radius, y_radius, markershape = :o, linetype = :scatter)

#     xcoords = [xcoords; x_radius]
#     ycoords = [ycoords; y_radius]

#     #add web flat
#     web_flat = shape_info.d/2 - shape_info.tf - radius

#     web_flat_range = LinRange(ycoords[end] + web_flat/n[5], (ycoords[end] + web_flat), n[5])
#     # web_flat_range = (ycoords[end] + web_flat/n[5]): web_flat/n[5]: (ycoords[end] + web_flat)
#     xcoords = [xcoords; ones(n[5])*xcoords[end]]
#     ycoords = [ycoords; web_flat_range]

#     #mirror about horizontal axis
#     ycoords_horz_flip = ycoords .- ycoords[end]
#     ycoords_horz_flip = -ycoords_horz_flip
#     ycoords_horz_flip = ycoords_horz_flip .+ ycoords[end]

#     xcoords = [xcoords; reverse(xcoords)[2:end]]
#     ycoords = [ycoords; reverse(ycoords_horz_flip)[2:end]]

#     #mirror about vertical axis
#     xcoords_vert_flip = reverse(-xcoords)[2:end-1]

#     xcoords = [xcoords; xcoords_vert_flip]
#     ycoords = [ycoords; reverse(ycoords)[2:end-1]]


#     return xcoords, ycoords

# end




# function discretize_w_shape_centerline_model(shape, cross_section)

#     num_branches = size(shape)[1]

#     xcoords = []
#     ycoords = []

#     for i = 1:num_branches

#         ΔL = shape[i].magnitude
#         θ = shape[i].direction
#         n = shape[i].n
#         Δxy = Geometry.vector_components(ΔL, θ) 
#         anchor = shape[i].anchor

#         if i == 1
#             xcoords = range(0.0, Δxy[1], length = n + 1) .+ anchor[1]
#             ycoords = range(0.0, Δxy[2], length = n + 1) .+ anchor[2]

#         else
#             xcoords = [xcoords; range(0.0, Δxy[1], length = n + 1) .+ anchor[1]]
#             ycoords = [ycoords; range(0.0, Δxy[2], length = n + 1) .+ anchor[2]]

#         end

#     end

#     #Round here to help unique function.
#     xycoords = [(round(xcoords[i], digits = 3), round(ycoords[i], digits = 3)) for i = 1:length(xcoords)]

#     xycoords = unique(xycoords)

#     coord = [y[i] for y in xycoords, i in 1:2]

#     #Shift coordinates so that web is centered on x=0.
#     coord[:, 1] = coord[:, 1] .- cross_section.bf/2

#     #Shift coordinates so that bottom fiber is at y=0.
#     coord[:, 2] = coord[:, 2] .+ cross_section.tf/2


#     #Define element connectivity.

#     num_elem = sum([shape[i].n for i=1:num_branches])
    
#     node_start = 1
#     node_end = shape[1].n

#     node_i = node_start:node_end
#     node_j = node_i .+ 1

#     node_start = floor(Int, shape[1].n/2)+1
#     node_end = node_end + 2

#     node_i = [node_i; node_start]
#     node_j = [node_j; node_end]

#     node_start = node_end
#     node_end = node_end + shape[2].n - 2

#     node_i = [node_i; node_start:node_end]
#     node_j = [node_j; (node_start:node_end) .+ 1]

#     node_start = shape[1].n + shape[2].n + 2
#     node_end = node_start + floor(Int, shape[2].n/2) - 1
#     node_i_range = range(node_start, node_end-1)
#     node_j_range = node_i_range .+ 1

#     node_i = [node_i; node_i_range]
#     node_j = [node_j; node_j_range]

#     node_start = node_i[end] + 1
#     node_end = shape[1].n + shape[2].n + 1

#     node_i = [node_i; node_start]
#     node_j = [node_j; node_end]

#     node_start = node_j[end]
#     node_end = node_i[end] + 1

#     node_i = [node_i; node_start]
#     node_j = [node_j; node_end]

#     node_start = shape[1].n + shape[2].n + 2 + floor(Int, shape[3].n/2)
#     node_end = node_start + floor(Int, shape[3].n/2) - 1
#     node_i_range = range(node_start, node_end-1)
#     node_j_range = node_i_range .+ 1

#     node_i = [node_i; node_i_range]
#     node_j = [node_j; node_j_range]

#     t = [ones(Float64, shape[1].n)*cross_section.tf[1]; ones(Float64, shape[2].n)*cross_section.tw[1]; ones(Float64, shape[2].n)*cross_section.tf[1]]

#     ends = [node_i node_j t]

#     return coord, ends

# end

end # module CrossSectionGeometry
