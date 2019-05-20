using StaticArrays
using ProgressMeter

"""
	create_uniform_grid(glued_grid::Array{Float64,3},psi_grid::Array{Float64,2},V::SVector{2,Float64})

Create a grid with uniform samples form a grid of the ``f \\circ \\phi`` lines.

See pages 54-55 in [B] for the general idea.

PRECONDITION: glued_grid is square.
"""

function create_uniform_grid(glued_grid::Array{Float64,3},psi_grid::Array{Float64,2},V::SVector{2,Float64})

	GRID_SIZE::Int64 = size(glued_grid)[1]
	uniform_grid::Array{Float64,3} = Array{Float64,3}(undef,GRID_SIZE,GRID_SIZE,3)


	#Iterate over the line x -> u*U + x*V
	@showprogress "Creating a uniform grid... " for u = 1:GRID_SIZE
		indices, V_coeffs = get_V_coordinates(convert(SVector{2,Int64},V),u,GRID_SIZE)
		looking_for_idx::Int64 = 1
		for (t, v) in R_mod_Z_interpolator(psi_grid[:,u],V_coeffs)
			i::Int64,j::Int64 = indices[looking_for_idx,:]
			uniform_grid[i,j,:] = (1-t)*glued_grid[mod(v-1,GRID_SIZE)+1,u,:] + t*glued_grid[mod(v,GRID_SIZE)+1,u,:]
			looking_for_idx +=1
		end
	end
	uniform_grid
end
"""
	rep_just_below(low::Float64,mid::Float64,high::Float64)

Given three numbers low, mid, high with low < high and mid <high,
return representatives low',mid',high' of the same EC as low, mid, high in R/Z
with low' and mid' being representatives within distance 1 og high'
"""
@inline function rep_just_below(low::Float64,mid::Float64,high::Float64)::Tuple{Float64,Float64,Float64}
	low = mod(low,1.0)
	mid = mod(mid,1.0)
	high = mod(high,1.0)

	if low > high
		low -=1
	end
	if mid > high
		mid -=1
	end
	return (low,mid,high)
end	

"""
	t_calculator(looking_for,q_last,q_new)

given three representatives of numbers q_last, 
looking_for and q_new where looking_for = (1-t)*q_last +
t*q_last for t in [0,1), return t
"""
function t_calculator(looking_for::Float64,q_last::Float64,
							q_new::Float64)
	#We need the right representatives to find the correct t
	looking_for, q_last, q_new = rep_just_below(looking_for,q_last,q_new)
	(looking_for - q_last)/(q_new - q_last)
end

"""
	R_mod_Z_interpolator(fixed,search)

Given two sorted arrays fixed and search of reals, we look
for each n=1:length(search) for (t,i) such that

search[n] = (1-t)*fixed[i-1] + t*fixed[i] for t in [0,1)

We return the list of these (t,i) pairs
See page 54 in [B]

NB: This function could be refactored into a generator,
but the gains would be minimal - this is not the bottleneck,
the integration is.
"""
function R_mod_Z_interpolator(fixed::Array{Float64,1},search::Array{Float64,1})::Array{Tuple{Float64,Int64},1}

	#We are searching for an interval in which
	#search[looking_for_idx] is contained
	looking_for_idx::Int64 = 1
	looking_for::Float64 = search[looking_for_idx]

	i::Int64 = 1
	n::Int64 = length(fixed)

	#We make sure that we start looking for the interval in
	#the right place
	while !is_inside_half_open_interval(fixed[mod(i-1,n)+1],fixed[mod(i,n)+1],looking_for)
		i += 1
	end
	
	#Allocate the result
	result::Array{Tuple{Float64,Int64},1} = Array{Tuple{Float64,Int64},1}(undef,length(search))

	#We know that looking_for < fixed[i]
	interval_counter::Int64 = 0
	q_last::Float64 = fixed[mod(i-1,n)+1]
	num_found::Int64 = 0

	#We look until we have looked at all intervals
	while interval_counter <= n
		q_new::Float64 = fixed[mod(i,n)+1]
		#If we have found that looking_for is inside [q_last,q_new)
		#we put it in result
		while is_inside_half_open_interval(q_last,q_new,looking_for)
			
			t::Float64 = t_calculator(looking_for,q_last,q_new)
			result[num_found+1] = (t,i)

			num_found += 1
			
			#If we have found all in search, we break
			if num_found == length(search)
				return result
			end

			#We have found what we are looking for, so we update
			#what we are looking for
			looking_for_idx += 1
			looking_for = search[looking_for_idx]
		end
		interval_counter +=1
		q_last = q_new
		i +=1
		#we wrap around
		if i> n
			i -= n
		end
	end
	return result
end

"""
	is_inside_half_open_interval(low::Float64,high::Float64,x::Float64)

Given low, high, x which are representatives of numbers in R/Z,
return true if x is inside [low,high) in the sense that there exists
t in [0,1) such that x = (1-t)*low +t*high
"""
function is_inside_half_open_interval(low::Float64,high::Float64,x::Float64)

	low = mod(low,1.0)
	high = mod(high,1.0)
	x = mod(x,1.0)

	if	low < high
		return low <= x < high
	end

	0<=x<high || low <=x<=1
end

"""
	get_V_coordinates(V::SVector{2,Int64},u::Int64,n::Int64)

Given a vector V in Z^2 and ints u and n return a list 
of all indices (i,j) such that (i/n,j/n) = u/n*U + x*V
for some x in [0,1] and a list with the x values.

Since  julia is 1-indexed, we add 1 to both i and j.

See page 54 in [B].
"""
function get_V_coordinates(V::SVector{2,Int64},u::Int64,n::Int64)
	# The idea of this algorithm is that if
	#
	# (i/n,j/n) = u/n*U + x*V then 
	# (i,j) = u*U + n*x*V
	# and we know that then u = q*i - p*j (recall V=(p,q))
	# we use a Diophantine solver to find suitable i and j
	p::Int64, q::Int64 = V
	i::Int64, j::Int64, gcd::Int64 = diophantine_solver(q,p)
	i,j = u*i, -u*j
	#Now q*i-j*p == u 
	
	step_i::Int64 = div(p,gcd) #slightly unnecessary since gcd=1
	step_j::Int64 = div(q,gcd)

	indices::Array{Int64,2} = Array{Int64,2}(undef,n,2)
	for r = 1:n
		indices[r,1] = i + r*step_i
		indices[r,2] = j + r*step_j
	end	

	#Now indices contains all (i,j) that solves the diophantine
	#equation

	#We find out what the V component is
	V_coeffs::Array{Float64,1} = Array{Float64,1}(undef,n)
	for j = 1:n
		V_coeffs[j] = mod(1/(n*(p^2+ q^2))*
								(p*indices[j,1] + q*indices[j,2]),1.0)
	end

	#We wish to sort with respect to the coefficients
	perm = sortperm(V_coeffs)

	V_coeffs = V_coeffs[perm]
	#We sort so the indices match with the V_coeffs again
	indices[:,1] = indices[perm,1]
	indices[:,2] = indices[perm,2]
	#Fix indices for julia 
	for j = 1:n
		indices[j,1] = mod(indices[j,1],n) + 1
		indices[j,2] = mod(indices[j,2],n) + 1
	end

	indices,V_coeffs
end

"""
	diophantine_solver(a::Int64,b::Int64)::Tuple{Int64,Int64,Int64}

Find a solution to the Diophantine equation a*x+b*y = gcd(a,b).
Returns (x,y,gcd)
"""
function diophantine_solver(a::Int64,b::Int64)::Tuple{Int64,Int64,Int64}
	x_last::Int64 = 1
	y_last::Int64 = 0
	x_new::Int64 = 0
	y_new::Int64 = 1
	r_last::Int64 = a
	r_new::Int64 = b
	while r_new != 0
		q::Int64 = div(r_last,r_new)
		r_last, r_new =  r_new, r_last%r_new
		x_last, x_new = x_new, x_last - q*x_new
		y_last, y_new = y_new, y_last - q*y_new
	end
	#Signs the sign of the modulo operation is not like in 
	#math, we adjust a bit
	if r_last < 0
		x_last, y_last, r_last = -x_last, -y_last, -r_last
	end
	x_last,y_last,r_last
end
