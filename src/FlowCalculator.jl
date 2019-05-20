using StaticArrays
using DifferentialEquations
using Plots
using LinearAlgebra


"""
	psi_curve(init_psi::Float64,U::SVector{2,Float64},V::SVector{2,Float64},f_jac,GRID_SIZE::Int64)

Compute the second component of the integral curve of 
the flow described in [B] at page 23 with respect to 
the direct orthogonal basis (U,V) and with initial 
value init_psi. See equation 2.6.
"""
function psi_curve(init_psi::Float64,U::SVector{2,Float64},V::SVector{2,Float64},f_jac,GRID_SIZE::Int64)


	#psi derived with respect to s at the point
	#s*U  + psi*V
	#This is equal to zeta in [B] p. 23 
	function psi_deriv(psi::Float64,p,s::Float64)::Float64
			JAC::MMatrix{3,2,Float64} = f_jac(s*U+psi*V)
			f_dir_U::SVector{3,Float64} = JAC * U
			f_dir_V::SVector{3,Float64} = JAC * V
			-dot(f_dir_U,f_dir_V)/dot(f_dir_V,f_dir_V)
	end


	time_span::Tuple{Float64,Float64} = (0.0,1.0)

	#Create an ODEProblem  with initial value
	initial_value_problem = ODEProblem(psi_deriv,init_psi,time_span)
	return solve(initial_value_problem,
		OwrenZen3(), #We use this solver 
						 #since it uses 3rd degree interp
		reltol = 1e-10, 
		abstol = 1e-10,
		dense=true, #We want to be able to 
		            #call the solution as a function
	)

end

"""
	plot_flow(U,V,f_jac,GRID_SIZE)

Plot the flow curves of the cylinder described 
in [B] p. 23. Uses Plots
"""
function plot_flow(U,V,f_jac,GRID_SIZE)
	# Create a new plot
	Plots.plot()
	
	# Plot 100 curves
	@showprogress "Plotting the flow..." for t = range(0.0,1.0,length = 100)
		psi_t = psi_curve(t,U,V,f_jac,GRID_SIZE)
		# The points is put into an array
		XY = [s*U+psi_t(s)*V for s in 1/GRID_SIZE*range(0,GRID_SIZE,step=1)] 
		# We plot it
		Plots.plot!(XY,legend = false,color =:blue)
	end
	# Make the plot visible
	gui()
end
