include("util.jl")
"""
	write_wrl(grid::Array{Float64,3},output_string::String)
Write a .wrl file with a regular square mesh of a map from the torus
to Euclidean 3-space.

PRECONDITION: grid is a square mesh
"""
function write_wrl(grid::Array{Float64,3},output_string::String)
	#A header
	template_string::String = """#VRML V2.0 utf8
Shape {
	appearance Appearance {
		material Material { diffuseColor 0.67 0.77 0.96}
	}
	geometry IndexedFaceSet {
		coord Coordinate {
			point [
""" 
	if length(output_string) <=3 || output_string[end-3:end] != ".wrl"
		output_string = string(output_string,".wrl")
	end
	open(output_string, "w") do file
		write(file,template_string)
		write_point_lines(grid,file)
		write(file,"\t\t\t]\n\t\t}\n\t\tcoordIndex [\n")
		write_face_lines(grid,file)
      write(file,"\t\t\t]\n\t\t}\n\t}\n}\n")
	end
end

"""
	write_face_lines(grid::Array{Float64,3},file)

write the lines in the "faces" block of the file,
which specifies which 4 points fit together to form 
a square polygon to draw.
"""
function write_face_lines(grid::Array{Float64,3},file)
	GRID_SIZE::Int64 = size(grid)[1]
	for i=1:GRID_SIZE, j=1:GRID_SIZE
		write(file,string(
					 "\t\t\t\t",
					 point_index(i,j,GRID_SIZE),
					 ", ",
					 point_index(i+1,j,GRID_SIZE),
					 ", ",
					 point_index(i+1,j+1,GRID_SIZE),
					 ", ",
					 point_index(i,j+1,GRID_SIZE),
					 ", -1,\n"
					 ))
	end
end

"""
	write_point_lines(grid::Array{Float64,3},file)

Write the lines in the "points" block of the file.
"""
function write_point_lines(grid::Array{Float64,3},file)
	GRID_SIZE::Int64 = size(grid)[1]
	for i=1:GRID_SIZE, j= 1:GRID_SIZE
		X::Float64, Y::Float64, Z::Float64 = grid[i,j,:]
		write(file,string("\t\t\t\t",X," ",Y," ",Z,",\n"))

	end
end
