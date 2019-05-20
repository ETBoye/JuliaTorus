using LinearAlgebra
using ProgressMeter
using StaticArrays
include("FlowCalculator.jl")
include("util.jl")


#The U vectors in direction 1, 2, and 3
const U_tuple = (SVector(1.0,0.0),1/5*SVector(1.0,2.0),1/5*SVector(1.0,-2.0))
#The V vectors in direction 1, 2, and 3
const V_tuple = (SVector(0.0,1.0),SVector(-2.0,1.0),SVector(2.0,1.0))

"""
	convex_integrate_flow_line!(...)::Nothing
	compute the convex integration of 
	``f \\circ \\phi(init_i/GRID_SIZE,\\cdot)`` with 
	respect to the direction j.
"""
function convex_integrate_flow_line!(
		psi_grid::Array{Float64,2}, integrated_grid::Array{Float64,3},
		init_i::Int64,
		GRID_SIZE::Int64,f,f_jac,F_0_jac,delta_k::Float64,
		N::Int64,j::Int64)::Nothing


	#Fetch the U,V vectors
	U::SVector{2,Float64} = U_tuple[j]	
	V::SVector{2,Float64} = V_tuple[j]	

	# The initial t value - we solve with initial value init_t*V
	init_t::Float64 = init_i/GRID_SIZE

	# Calculate psi_t(s) curve
	psi_t = psi_curve(init_t,U,V,f_jac,GRID_SIZE)

	function h(y::SVector{3,Float64},p,u::Float64)::SVector{3,Float64}
		# This function returns
		#h(init_t,u,{N*u}) as in p. 52 in [B]
		phi::SVector{2,Float64} = u*U + psi_t(u)*V
		JAC::MMatrix{3,2,Float64} = f_jac(phi)

		#Precompute the derivatives of f
		#in the directions of U and V
		f_dir_U::SVector{3,Float64} = JAC * U
		f_dir_V::SVector{3,Float64} = JAC * V

		#Compute W
		zeta::Float64 = -dot(f_dir_U,f_dir_V)/dot(f_dir_V,f_dir_V)
		W::SVector{2,Float64} = U + zeta*V

		#Compute the derivative of f 
		#in the direction of W
		f_dir_W::SVector{3,Float64} = JAC * W

		#Compute the gram matrix of D (See p. 53 in [B])
		D_gram_matrix::MMatrix{2,2,Float64} = 
			iso_default_gram_matrix(F_0_jac(phi),JAC,delta_k)

		#Compute the gram matrix of mu (See p. 53 in [B])
		mu_gram_matrix::MMatrix{2,2,Float64} = 
			pullback_metric_gram_matrix(JAC) + rho_j_decomp_gram_matrix(D_gram_matrix,j)

		#Compute \mu(W,W)
		mu_on_W_and_W::Float64 = transpose(W) * mu_gram_matrix * W

		r::Float64 = sqrt(mu_on_W_and_W)


		norm_f_dir_W::Float64 = norm(f_dir_W)
		#The tangent
		t::SVector{3,Float64} = f_dir_W/norm_f_dir_W

		f_dir_W_cross_f_dir_V::SVector{3,Float64} = cross(JAC[:,1],JAC[:,2])

		#A normal to the surface
		n::SVector{3,Float64} = f_dir_W_cross_f_dir_V/norm(f_dir_W_cross_f_dir_V)

		alpha::Float64 = bessel0_inv(norm_f_dir_W/r)
		trig_arg::Float64 = alpha*cos(2*pi*(N*u - floor(N*u)))

		#Finally we return h(t,u,{N*u})
		r * (cos(trig_arg)*t + sin(trig_arg)*n)
	end

	#We create the IVP
	initial_value_problem = ODEProblem(h,f(init_t*V),(0.0,1.0))

	#We solve it
	F_circ_phi_t = solve(
		initial_value_problem,
		Vern7(), #This is the recommended solver for these tolerances
		saveat=1/GRID_SIZE, 
		dense=false, #We do not need continuous output
		save_everystep = false, #We only want to save at the gridpoints
		abstol = 1e-10,
		reltol = 1e-10)


	#We glue
	f_phi_end::SVector{3,Float64} = f(U + psi_t(1.0)*V)

	for idx = 1:GRID_SIZE
		integrated_grid[init_i+1,idx,:] = F_circ_phi_t.u[idx] - w((idx-1)/GRID_SIZE)*(F_circ_phi_t.u[GRID_SIZE+1] - f_phi_end)
	end

	psi_grid[init_i+1,:] .= [psi_t(t) for t in 1/GRID_SIZE*(0:GRID_SIZE)]
	return
end


"""
	w(t::Float64)::Float64

This map is used to glue the integrated curves.
See page 53 in	[B]
"""
@inline function w(t::Float64)::Float64
	t^3.0*(6.0*t^2-15.0*t+10.0)
end


"""
	convex_integrate_grid(...)

Compute the	grid of the second component of the flow curves in the
basis of (U(j),V(j)) and compute the glued integrated grid.
"""
function convex_integrate_grid(
		GRID_SIZE::Int64,f,f_jac,F_0_jac,delta_k::Float64,
		N::Int64,dir::Int64)::Tuple{Array{Float64,2},Array{Float64,3}}

	# Allocate the grids
	integrated_grid::Array{Float64,3} = zeros(Float64,GRID_SIZE,GRID_SIZE,3)
	psi_grid::Array{Float64,2} = zeros(Float64,GRID_SIZE,GRID_SIZE+1)

	@showprogress "Integrating and gluing...  " for i=0:GRID_SIZE-1
		convex_integrate_flow_line!(psi_grid,integrated_grid,
				i,GRID_SIZE,f,f_jac,F_0_jac,delta_k,N,dir)
	end

	return psi_grid, integrated_grid
end

const EFG_to_rho_j_matrix =
	#Input is columns apparently, so this below represents
	#the transposed matrix of what it exactly is
	SMatrix{3,3,Float64}(1.0, 0.0, 0.0,
								1/5, 2/5,4/5,
								1/5	, -2/5, 4/5,
							  )
const EFG_to_rho_j_matrix_inv = inv(EFG_to_rho_j_matrix)

const ell_j_matrix = [
							 SMatrix{2,2,Float64}(1.0,0.0,0.0,0.0),
							 SMatrix{2,2,Float64}(1/5,2/5,2/5,4/5),
							 SMatrix{2,2,Float64}(1/5,-2/5,-2/5,4/5),
							]


"""
	rho_j_decomp_gram_matrix(D::MMatrix{2,2},j::Int64)

Compute the Gram matrix of ``rho_j(D)*\\ell(j) \\otimes  \\ell(j)``
See page 53 in [B]
"""
@inline function rho_j_decomp_gram_matrix(D::MMatrix{2,2},j::Int64)
	EFG_vec::SVector{3,Float64} = SVector{3,Float64}(D[1,1],D[1,2],D[2,2])
	rho_j = dot(EFG_to_rho_j_matrix_inv[j,:],EFG_vec)
	rho_j * ell_j_matrix[j]
end


"""
	g_k_gram_matrix(F_0_jac::MMatrix{3,2,Float64},delta_k::Float64)

Compute the Gram Matrix of g_k as in page 53 in [B]
"""
@inline function g_k_gram_matrix(F_0_jac::MMatrix{3,2,Float64},delta_k::Float64)
	(1.0-delta_k)*pullback_metric_gram_matrix(F_0_jac) + delta_k * SMatrix{2,2,Float64}(1.0,0.0,0.0,1.0)
end

"""
	iso_default_gram_matrix(F_0_JAC::MMatrix{3,2},f_jac::MMatrix{3,2},delta_k)

Compute the Gram matrix of the isometric default. This is
``D_{k,j}`` in page 53 in [B]
"""
@inline function iso_default_gram_matrix(F_0_JAC::MMatrix{3,2},f_jac::MMatrix{3,2},delta_k)
	g_k_gram_matrix(F_0_JAC,delta_k) - pullback_metric_gram_matrix(f_jac)
end


"""
	pullback_metric_gram_matrix(JAC::MMatrix{3,2,Float64})::MMatrix{2,2}

Compute the Gram matrix of the pullback metric given the jacobian
"""
@inline function pullback_metric_gram_matrix(JAC::MMatrix{3,2,Float64})::MMatrix{2,2}
	f_x::SVector{3,Float64} = JAC[:,1]
	f_y::SVector{3,Float64} = JAC[:,2]
	MMatrix{2,2,Float64}(dot(f_x,f_x), dot(f_x,f_y),
										 dot(f_x,f_y), dot(f_y,f_y)
   )
end
