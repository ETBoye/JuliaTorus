using StaticArrays



"""
A struct holding the data for the Interpolation of a map
from the flat square torus to euclidean n-space
"""
struct InterPolationData
	f_grid::Array{Float64,3}
	f_x::Array{Float64,3} #a grid with estimates of f derived w.r.t. x
	f_y::Array{Float64,3} #a grid with estimates of f derived w.r.t. x
	f_xy::Array{Float64,3} #a grid with estimates of the cross derivative
end


"""
	create_interp_data(f_grid::Array{Float64,3})::InterPolationData

Compute the interpoloation data from a grid of f
"""
function create_interp_data(f_grid::Array{Float64,3})::InterPolationData
	f_x::Array{Float64,3} = derive_first_direction(f_grid)
	f_y::Array{Float64,3} = derive_second_direction(f_grid)
	f_xy::Array{Float64,3} = derive_second_direction(f_x)
	InterPolationData(f_grid,f_x,f_y,f_xy)
end

"""
	derive_first_direction(f_grid::Array{Float64,3})::Array{Float64,3}

Compute a grid with estimates of the first partial derivative.
"""
function derive_first_direction(f_grid::Array{Float64,3})::Array{Float64,3}
	f_x::Array{Float64,3} = zeros(Float64,size(f_grid))
	size_i::Int64,size_j::Int64 = size(f_grid)
	for j=1:size_j, i=1:size_i
		#We use finite differences of order 4. See p. 48 in [B]
		f_x[i,j,:] = 1/12*(
			 f_grid[fix_idx(i-2,size_i),j,:] #We wrap the indices
			 -8*f_grid[fix_idx(i-1,size_i),j,:]
			 +8*f_grid[fix_idx(i+1,size_i),j,:]
			 -f_grid[fix_idx(i+2,size_i),j,:]
		)
	end
	f_x
end

"""
	derive_second_direction(f_grid::Array{Float64,3})::Array{Float64,3}

Compute a grid with estimates of the second partial derivative.
"""
function derive_second_direction(f_grid::Array{Float64,3})::Array{Float64,3}
	f_y::Array{Float64,3} = zeros(Float64,size(f_grid))
	size_i::Int64,size_j::Int64 = size(f_grid)
	for j=1:size_j, i=1:size_i
		#We use finite differences of order 4. See p. 48 in [B]
		f_y[i,j,:] = 1/12*(
			f_grid[i,fix_idx(j-2,size_j),:] #We wrap the indices
			-8*f_grid[i,fix_idx(j-1,size_j),:]
			+8*f_grid[i,fix_idx(j+1,size_j),:]
			-f_grid[i,fix_idx(j+2,size_j),:]
		)
	end
	f_y
end

"""
	create_interpolation_functions(f_grid::Array{Float64,3})

Return interpolation maps of f and the jacobian of f
given a regular square grid of f

PRECONDITION: f_grid is square
"""
function create_interpolation_functions(f_grid::Array{Float64,3})
#We capture f_grid for performance
let f_grid::Array{Float64,3}=f_grid
	GRID_SIZE::Int64 = size(f_grid)[1]
	interp_data::InterPolationData = create_interp_data(f_grid)
	
	#An index wrapper
	@inline function fix_idx(idx::Int64)::Int64
			mod(idx-1,GRID_SIZE)+1
	end

	#Below is an implementation of the algorithm 
	#described in pp.48-50 in [B].
	#This is very hand-coded for performance
	function f(p::SVector{2,Float64})::SVector{3,Float64}
		x::Float64,y::Float64 = p 
		i::Int64, j::Int64  = floor(Int64,GRID_SIZE*x), floor(Int64,GRID_SIZE*y) #(x,y) is inside the square p^(i,j) described p. 48 in [B]

		#The local coordinates
		x_i::Float64, y_j::Float64 = GRID_SIZE*x-i,GRID_SIZE*y-j 

		#Precompute the values of the hermite basis functions
		h_0_x::Float64 = hermite_0(x_i)
		h_1_x::Float64 = hermite_1(x_i)
		h_2_x::Float64 = hermite_2(x_i)
		h_3_x::Float64 = hermite_3(x_i)

		h_0_y::Float64 = hermite_0(y_j)
		h_1_y::Float64 = hermite_1(y_j)
		h_2_y::Float64 = hermite_2(y_j)
		h_3_y::Float64 = hermite_3(y_j)

		#adjust for 1-indexing
		i = fix_idx(i+1)
		j = fix_idx(j+1)

		i_p_1::Int64 = fix_idx(i+1)
		j_p_1::Int64 = fix_idx(j+1)


		#The computation of the estimate on f
		#	A comment on the macros:
		#		@inbounds means that Julia wont check if the indices
		#		are inside bounds
		#
		#		@views means that Julia does not create a copy of the 
		#		slice of the array we are accessing
		#
		#	For the rest, see eq. 4.5 at p. 49 in [B]
		q::SVector{3,Float64} = @inbounds @views  interp_data.f_grid[i,j,:] 
		result::SVector{3,Float64} =q * h_0_x * h_0_y


		q = @inbounds @views  interp_data.f_grid[i,j_p_1,:]
		result += q*h_0_x * h_1_y
		q = @inbounds @views  interp_data.f_grid[i_p_1,j,:]
		result += q*h_1_x * h_0_y
		q =@inbounds @views  interp_data.f_grid[i_p_1,j_p_1,:] 
		result += q*h_1_x * h_1_y

		q = @inbounds @views  interp_data.f_x[i,j,:]
		result += q*h_2_x*h_0_y
		q = @inbounds @views  interp_data.f_x[i,j_p_1,:]
		result += q*h_2_x * h_1_y
		q = @inbounds @views  interp_data.f_x[i_p_1,j,:]
		result += q*h_3_x * h_0_y
		q = @inbounds @views  interp_data.f_x[i_p_1,j_p_1,:]
		result += q*h_3_x * h_1_y

		q =  @inbounds @views  interp_data.f_y[i,j,:] 
		result +=q* h_0_x *h_2_y
		q = @inbounds @views  interp_data.f_y[i,j_p_1,:]
		result += q*h_0_x * h_3_y
		q = @inbounds @views  interp_data.f_y[i_p_1,j,:]
		result += q*h_1_x * h_2_y
		q = @inbounds @views  interp_data.f_y[i_p_1,j_p_1,:]
		result += q*h_1_x * h_3_y

		q = @inbounds @views  interp_data.f_xy[i,j,:] 
		result += q* h_2_x * h_2_y
		q = @inbounds @views  interp_data.f_xy[i,j_p_1,:]
		result += q*h_2_x * h_3_y
		q = @inbounds @views  interp_data.f_xy[i_p_1,j,:]
		result += q*h_3_x * h_2_y
		q = @inbounds @views  interp_data.f_xy[i_p_1,j_p_1,:]
		result += q*h_3_x * h_3_y

		return result
	end
	function f_jac(p::SVector{2,Float64})::MMatrix{3,2,Float64}
		x::Float64,y::Float64 = p
		i::Int64, j::Int64  = floor(Int64,GRID_SIZE*x), floor(Int64,GRID_SIZE*y) # (x,y) is inside the square p^(i,j) described p. 48 in [B]

		# The local coordinates
		x_i::Float64, y_j::Float64 = GRID_SIZE*x-i,GRID_SIZE*y-j

		#The computation of the estimate on f
		#	A comment on the macros:
		#		@inbounds means that Julia wont check if the indices
		#		are inside bounds
		#
		#		@views means that Julia does not create a copy of the 
		#		slice of the array we are accessing
		#
		#	For the rest, see eq. 4.5 at p. 49 in [B] 
		#	for the partial derivatives, we formally derive
		#	the interpolation map. This only means calculating
		#	the derivatives of the hermite basis maps
 
		h_0_x::Float64 = hermite_0(x_i)
		h_1_x::Float64 = hermite_1(x_i)
		h_2_x::Float64 = hermite_2(x_i)
		h_3_x::Float64 = hermite_3(x_i)

		h_0_x_p::Float64 = hermite_0_prime(x_i)
		h_1_x_p::Float64 = hermite_1_prime(x_i)
		h_2_x_p::Float64 = hermite_2_prime(x_i)
		h_3_x_p::Float64 = hermite_3_prime(x_i)


		h_0_y::Float64 = hermite_0(y_j)
		h_1_y::Float64 = hermite_1(y_j)
		h_2_y::Float64 = hermite_2(y_j)
		h_3_y::Float64 = hermite_3(y_j)

		h_0_y_p::Float64 = hermite_0_prime(y_j)
		h_1_y_p::Float64 = hermite_1_prime(y_j)
		h_2_y_p::Float64 = hermite_2_prime(y_j)
		h_3_y_p::Float64 = hermite_3_prime(y_j)

		#adjust for 1-indexing
		i = fix_idx(i+1)
		j = fix_idx(j+1)

		i_p_1::Int64 = fix_idx(i+1)
		j_p_1::Int64 = fix_idx(j+1)

		J::MMatrix{3,2,Float64} = MMatrix{3,2,Float64}(undef)

		q::SVector{3,Float64} = @inbounds @views  interp_data.f_grid[i,j,:]
		J[:,1] = q * h_0_x_p * h_0_y
		J[:,2] = q * h_0_x   * h_0_y_p
		q= @inbounds @views  interp_data.f_grid[i,j_p_1,:]
		J[:,1] += q * h_0_x_p * h_1_y
		J[:,2] += q * h_0_x   * h_1_y_p
		q= @inbounds @views  interp_data.f_grid[i_p_1,j,:]
		J[:,1] += q * h_1_x_p * h_0_y
		J[:,2] += q * h_1_x   * h_0_y_p
		q= @inbounds @views  interp_data.f_grid[i_p_1,j_p_1,:]
		J[:,1] += q * h_1_x_p * h_1_y
		J[:,2] += q * h_1_x   * h_1_y_p

		q= @inbounds @views  interp_data.f_x[i,j,:]
		J[:,1] += q * h_2_x_p * h_0_y
		J[:,2] += q * h_2_x   * h_0_y_p
		q= @inbounds @views  interp_data.f_x[i,j_p_1,:]
		J[:,1] += q * h_2_x_p * h_1_y
		J[:,2] += q * h_2_x   * h_1_y_p
		q= @inbounds @views  interp_data.f_x[i_p_1,j,:]
		J[:,1] += q * h_3_x_p * h_0_y
		J[:,2] += q * h_3_x   * h_0_y_p
		q= @inbounds @views  interp_data.f_x[i_p_1,j_p_1,:]
		J[:,1] += q * h_3_x_p * h_1_y
		J[:,2] += q * h_3_x   * h_1_y_p

		q= @inbounds @views  interp_data.f_y[i,j,:]
		J[:,1] += q * h_0_x_p * h_2_y
		J[:,2] += q * h_0_x   * h_2_y_p
		q= @inbounds @views  interp_data.f_y[i,j_p_1,:]
		J[:,1] += q * h_0_x_p * h_3_y
		J[:,2] += q * h_0_x   * h_3_y_p
		q= @inbounds @views  interp_data.f_y[i_p_1,j,:]
		J[:,1] += q * h_1_x_p * h_2_y
		J[:,2] += q * h_1_x   * h_2_y_p
		q= @inbounds @views  interp_data.f_y[i_p_1,j_p_1,:]
		J[:,1] += q * h_1_x_p * h_3_y
		J[:,2] += q * h_1_x   * h_3_y_p

		q= @inbounds @views  interp_data.f_xy[i,j,:]
		J[:,1] += q * h_2_x_p * h_2_y
		J[:,2] += q * h_2_x   * h_2_y_p
		q= @inbounds @views  interp_data.f_xy[i,j_p_1,:]
		J[:,1] += q * h_2_x_p * h_3_y
		J[:,2] += q * h_2_x   * h_3_y_p
		q= @inbounds @views  interp_data.f_xy[i_p_1,j,:]
		J[:,1] += q * h_3_x_p * h_2_y
		J[:,2] += q * h_3_x   * h_2_y_p
		q= @inbounds @views  interp_data.f_xy[i_p_1,j_p_1,:]
		J[:,1] += q * h_3_x_p * h_3_y
		J[:,2] += q * h_3_x   * h_3_y_p
		#Adjust for scaling
		@. J =  GRID_SIZE * J
		return J
	end
	return f,f_jac
end
end



@inline function fix_idx(idx::Int64, GRID_SIZE::Int64)::Int64
	mod(idx-1,GRID_SIZE)+1
end


# The Hermite basis functions 

@inline function hermite_0(t::Float64)::Float64
	(1+2*t)*(1-t)^2
end

@inline function hermite_1(t::Float64)::Float64
	t^2*(3-2*t)
end

@inline function hermite_2(t::Float64)::Float64
	t*(1-t)^2
end

@inline function hermite_3(t::Float64)::Float64
	t^2*(t-1)
end

#The derived Hermite basis functions
@inline function hermite_0_prime(t::Float64)::Float64
	6*t^2 - 6*t
end

@inline function hermite_1_prime(t::Float64)::Float64
	6*t-6*t^2
end

@inline function hermite_2_prime(t::Float64)::Float64
	3*t^2-4*t+1
end

@inline function hermite_3_prime(t::Float64)::Float64
	3*t^2-2*t
end
