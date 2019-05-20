include("util.jl")
using ProgressMeter

#A dictionary for different angles of view
camera_positions = Dict(
								"default" => (

        (0.2737621663770907,-0.24160876929304848,0.21589505869819867),
        (-0.10046632767779132,0.5972281587874799,0.7957543863242322),
		  (0.0,0.0,0.0))
		  )

"""
	write_xml(grid::Array{Float64,3},num_pixels::Int64,output_string::String)

	Write an yafaray .xml scene file to disk to be read by the yafaray-v3 raytracer. 
	Once written, one runs the command 
	yafaray-xml -f png -t {NUMBER_OF_CPU_CORES} {FILENAME}.xml
	to generate a .png of the scene
"""
function write_xml(grid::Array{Float64,3},num_pixels::Int64,output_string::String)
	#The header string
	#This was taken from another scene .xml
	#so a lot of stuff isn't used
	template_string_beginning = """<?xml version="1.0"?>
	<scene type="triangle">

	<material name="defaultMat">
		<type sval="shinydiffusemat"/>
	</material>

	<material name="Yellow.001">
		<IOR fval="1"/>
		<color r="0.55294118" g="0.81568627" b="0.9333333" a="1"/>
		<diffuse_reflect fval="1"/>
		<emit fval="0"/>
		<fresnel_effect bval="false"/>
		<mirror_color r="1" g="1" b="1" a="1"/>
		<specular_reflect fval="0"/>
		<translucency fval="0"/>
		<transmit_filter fval="1"/>
		<transparency fval="0"/>
		<type sval="shinydiffusemat"/>
	</material>

	<material name="Red.002">
		<IOR fval="1"/>
		<color r="0.83828" g="0.0722165" b="0.0187428" a="1"/>
		<diffuse_reflect fval="1"/>
		<emit fval="0"/>
		<fresnel_effect bval="false"/>
		<mirror_color r="1" g="1" b="1" a="1"/>
		<specular_reflect fval="0"/>
		<translucency fval="0"/>
		<transmit_filter fval="1"/>
		<transparency fval="0"/>
		<type sval="shinydiffusemat"/>
	</material>

	<material name="Green.002">
		<IOR fval="1"/>
		<color r="0.143054" g="0.711083" b="0.147195" a="1"/>
		<diffuse_reflect fval="1"/>
		<emit fval="0"/>
		<fresnel_effect bval="false"/>
		<mirror_color r="1" g="1" b="1" a="1"/>
		<specular_reflect fval="0"/>
		<translucency fval="0"/>
		<transmit_filter fval="1"/>
		<transparency fval="0"/>
		<type sval="shinydiffusemat"/>
	</material>

	<material name="Blue.002">
		<IOR fval="1"/>
		<color r="0.0678124" g="0.329362" b="0.762696" a="1"/>
		<diffuse_reflect fval="1"/>
		<emit fval="0"/>
		<fresnel_effect bval="false"/>
		<mirror_color r="1" g="1" b="1" a="1"/>
		<specular_reflect fval="0"/>
		<translucency fval="0"/>
		<transmit_filter fval="1"/>
		<transparency fval="0"/>
		<type sval="shinydiffusemat"/>
	</material>

	<material name="white">
		<IOR fval="1"/>
		<color r="1" g="1" b="1" a="1"/>
		<diffuse_reflect fval="1"/>
		<emit fval="0"/>
		<fresnel_effect bval="false"/>
		<mirror_color r="1" g="1" b="1" a="1"/>
		<specular_reflect fval="0"/>
		<translucency fval="0"/>
		<transmit_filter fval="1"/>
		<transparency fval="0"/>
		<type sval="shinydiffusemat"/>
	</material>

	<material name="Glass">
		<IOR fval="1.55"/>
		<absorption r="1" g="1" b="1" a="1"/>
		<absorption_dist fval="1"/>
		<dispersion_power fval="0"/>
		<fake_shadows bval="false"/>
		<filter_color r="1" g="1" b="1" a="1"/>
		<mirror_color r="1" g="1" b="1" a="1"/>
		<transmit_filter fval="1"/>
		<type sval="glossy"/>
	</material>

	<material name="Yellow">
		<IOR fval="1"/>
		<color r="1" g="1" b="1" a="1"/>
		<diffuse_reflect fval="1"/>
		<emit fval="0"/>
		<fresnel_effect bval="false"/>
		<mirror_color r="1" g="1" b="1" a="1"/>
		<specular_reflect fval="0"/>
		<translucency fval="0"/>
		<transmit_filter fval="1"/>
		<transparency fval="0"/>
		<type sval="shinydiffusemat"/>
	</material>

	<material name="Lamp.001">
		<IOR fval="1"/>
		<color r="1" g="1" b="1" a="1"/>
		<diffuse_reflect fval="1"/>
		<emit fval="0"/>
		<fresnel_effect bval="false"/>
		<mirror_color r="1" g="1" b="1" a="1"/>
		<power fval="5"/>
		<specular_reflect fval="0"/>
		<translucency fval="0"/>
		<transmit_filter fval="1"/>
		<transparency fval="0"/>
		<type sval="light_mat"/>
	</material>

	<light name="Lamp.001">
		<color r="1" g="1" b="1" a="1"/>
		<corner x="-0.25" y="-0.25" z="1.99646"/>
		<from x="0" y="0" z="1.99646"/>
		<point1 x="-0.25" y="0.25" z="1.99646"/>
		<point2 x="0.25" y="-0.25" z="1.99646"/>
		<power fval="5"/>
		<samples ival="16"/>
		<type sval="arealight"/>
	</light>

	<camera name="cam">
		<aperture fval="0"/>
		<bokeh_rotation fval="0"/>
		<bokeh_type sval="disk1"/>
		<dof_distance fval="0"/>
		<focal fval="1.37374"/>
		<from x="{FROM_X}" y="{FROM_Y}" z="{FROM_Z}"/>
		<resx ival="{NUM_PIXELS}"/>
		<resy ival="{NUM_PIXELS}"/>
		<to x="{TO_X}" y="{TO_Y}" z="{TO_Z}"/>
		<type sval="perspective"/>
		<up x="{UP_X}" y="{UP_Y}" z="{UP_Z}"/>
	</camera>
	"""
	template_string_end = """
	<background name="world_background">
		<color r="0.92156863" g="0.92156863" b="0.92156863" a="1"/>
		<power fval="1"/>
		<type sval="constant"/>
	</background>

	<integrator name="default">
		<bounces ival="3"/>
		<caustic_mix ival="5"/>
		<diffuseRadius fval="1"/>
		<fg_bounces ival="3"/>
		<fg_samples ival="32"/>
		<finalGather bval="true"/>
		<photons ival="800000"/>
		<raydepth ival="4"/>
		<search ival="150"/>
		<shadowDepth ival="2"/>
		<show_map bval="false"/>
		<transpShad bval="false"/>
		<type sval="photonmapping"/>
		<use_background bval="false"/>
	</integrator>

	<integrator name="volintegr">
		<type sval="none"/>
	</integrator>
	<render>
		<AA_inc_samples ival="2"/>
		<AA_minsamples ival="4"/>
		<AA_passes ival="2"/>
		<AA_pixelwidth fval="1.5"/>
		<AA_threshold fval="0.05"/>
		<background_name sval="world_background"/>
		<camera_name sval="cam"/>
		<clamp_rgb bval="true"/>
		<filter_type sval="mitchell"/>
		<gamma fval="2.2"/>
		<height ival="{NUM_PIXELS}"/>
		<integrator_name sval="default"/>
		<threads ival="1"/>
		<volintegrator_name sval="volintegr"/>
		<width ival="{NUM_PIXELS}"/>
		<xstart ival="0"/>
		<ystart ival="0"/>
		<z_channel bval="true"/>
	</render>
	</scene>
	"""

	#Add .xml if it is not already there
	if length(output_string) <=3 || output_string[end-3:end] != ".xml"
		output_string = string(output_string,".xml")
	end
	cam = "default"
	open(output_string, "w") do file
		write(file,insert_variables_in_template_string(template_string_beginning,num_pixels,cam))
		write_mesh_block(grid,file)
		write(file,insert_variables_in_template_string(template_string_end,num_pixels,cam))
	end
end

"""
	generate_xml_point_lines(grid,file)

Write the points to the .xml
"""
function generate_xml_point_lines(grid,file)
	GRID_SIZE::Int64 = size(grid)[1]
	@showprogress for i=1:GRID_SIZE, j=1:GRID_SIZE
		if j == 1
			sleep(0.03)
		end
		s = 
		X,Y,Z = grid[i,j,:]
		write(file,
				string("\t<p x=\"",X,"\" y=\"",Y,"\" z=\"",Z,"\"/>\n"))
	end


end

"""
	generate_xml_face_lines(grid,file)

Write which points form triangles to the .xml
"""
function generate_xml_face_lines(grid,file)
	GRID_SIZE::Int64 = size(grid)[1]
	@showprogress for i=1:GRID_SIZE, j=1:GRID_SIZE
		if j == 1
			sleep(0.03)
		end
		s = string("\t<f a=\"",
					  point_index(i,j,GRID_SIZE),
					  "\" b=\"",
					  point_index(i+1,j,GRID_SIZE),
					  "\" c=\"",
					  point_index(i,j+1,GRID_SIZE),
					  "\"/>\n")
		write(file,s)
		s = string("\t<f a=\"",
					  point_index(i+1,j+1,GRID_SIZE),
					  "\" b=\"",
					  point_index(i+1,j,GRID_SIZE),
					  "\" c=\"",
					  point_index(i,j+1,GRID_SIZE),
					  "\"/>\n")
		write(file,s)
	end
end

"""
	write_mesh_block(grid,file)

Write the mesh block to the scene .xml	
"""
function write_mesh_block(grid,file)
	num_vertices::Int64 = size(grid)[1]^2
	num_faces::Int64 = 2*num_vertices
	s = """<mesh vertices="{1}" faces="{2}" has_orco="false" has_uv="false" type="0">\n"""
	s = replace(s,"{1}"=>num_vertices)
	s = replace(s, "{2}"=>num_faces)
	write(file,s)
	generate_xml_point_lines(grid,file)
	write(file,"""\t<set_material sval="Yellow.001"/>\n""")
	generate_xml_face_lines(grid,file)
	write(file,"</mesh>\n")
end


function insert_variables_in_template_string(s::String,num_pixels,cam::String)
	FROM_X, FROM_Y, FROM_Z = camera_positions[cam][1]

	UP_X, UP_Y, UP_Z = camera_positions[cam][2]
	TO_X, TO_Y, TO_Z = camera_positions[cam][3]

	reduce(replace,[
			  "{FROM_X}"=>FROM_X,
				"{FROM_Y}"=>FROM_Y,
            "{FROM_Z}"=>FROM_Z,
            "{UP_X}"=>UP_X,
            "{UP_Y}"=>UP_Y,
            "{UP_Z}"=>UP_Z,
            "{TO_X}"=>TO_X,
            "{TO_Y}"=>TO_Y,
				"{TO_Z}"=>TO_Z,
				"{NUM_PIXELS}"=>num_pixels]
			 ,init = s)
end



