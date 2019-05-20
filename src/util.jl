using SpecialFunctions: besselj0, besselj1
using StaticArrays
			

"""
	bessel0_inv(y::Float64)

Compute the inverse of the zeroth order Bessel
function of the first kind.
"""
function bessel0_inv(y::Float64)::Float64
	# We use Newton's method on f(x) = J0(x) -y
	a::Float64 = 2.4048
	b::Float64 = 0
	while true
		b = a + (besselj0(a) - y)/besselj1(a)
		if abs(b-a) <= 1e-10 #The tolerance
			break
		end
		a = b
	end
	return b
end


"""
	point_index(i,j,GRID_SIZE)

Compute the 0-indexed one-dimensional index based on the
two-dimensional 1-indexed indices

This function is used by the output file builders,
namely wrlBuilder and xmlBuilder
"""
function point_index(i,j,GRID_SIZE)
	#Because of 1-indexing we have to shift a bit
	(GRID_SIZE*(i-1) + (j-1)) % GRID_SIZE^2
end

"""
   standard_param(r_1,r_2)

Return a standard parametrization 
	f(p::SVector{2,Float64})::SVector{3,Float64}
	of the torus in euclidean 3-space. 
"""
function standard_param(r_1::Float64,r_2::Float64)
# We capture the variables for performance
let r_1::Float64 = r_1, r_2:Float64 = r_2 
	function f(p::SVector{2,Float64})::SVector{3,Float64}
		x::Float64,y::Float64 = p
		X::Float64 = 1/(2*pi)*(r_2 + r_1*cos(2*pi*x))*cos(2*pi*y)
		Y::Float64 = 1/(2*pi)*(r_2 + r_1*cos(2*pi*x))*sin(2*pi*y)
		Z::Float64 = 1/(2*pi)*(r_1*sin(2*pi*x))
		return SVector{3,Float64}(X,Y,Z)
	end
end
end

"""
   standard_param_jacobi(r_1,r_2)

Return the jacobian of the standard parametrization with radii
r_1 and r_2. The signature of the returned function is
jac(p::SVector{2,Float64})::MMAtrix{3,2,Float64}.
"""
function standard_param_jacobi(r_1::Float64,r_2::Float64)
# We capture the variables for performance
let r_1::Float64 = r_1, r_2:Float64 = r_2
	function jac(p::SVector{2,Float64})::MMatrix{3,2,Float64}
		x::Float64, y::Float64 = p
		MMatrix{3,2,Float64}(
			-r_1*sin(2*pi*x)*cos(2*pi*y), #f1_x
			-r_1*sin(2*pi*x)*sin(2*pi*y), #f2_x
			r_1*cos(2*pi*x), #f3_x
			-(r_2 + r_1*cos(2*pi*x))*sin(2*pi*y), #f1_y
			(r_2 + r_1*cos(2*pi*x))*cos(2*pi*y), #f2_y
			0.0 #f3_y
		)
	end
end
end


"""
	calculate_grid(f::Function,GRID_SIZE::Int64)::Array{Float64,3}

Compute a regular square grid given a map from the torus to 
Euclidean 3-space.
"""
	function calculate_grid(f::Function,GRID_SIZE::Int64)::Array{Float64,3}
	grid::Array{Float64,3} = Array{Float64,3}(undef,GRID_SIZE,GRID_SIZE,3)
	for i = 1:GRID_SIZE, j = 1:GRID_SIZE
		x = (i-1)/GRID_SIZE #We subtract 1 because of 1-indexing
		y = (j-1)/GRID_SIZE
		grid[i,j,:] = f(SVector{2,Float64}(x,y))
	end
	return grid
end


