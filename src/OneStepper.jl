include("ConvexIntegrator.jl")
include("UniformGrid.jl")
include("Interpolator.jl")
include("wrlBuilder.jl")
include("xmlBuilder.jl")


"""
	do_one_step(...)

This is the implementation equivalent to 
the One Step Theorem. We perform a convex integration
and return the psi_grid (Related to the flow,
see FlowCalculator.jl) and a uniform grid of the resulting map
"""
function do_one_step(
		GRID_SIZE::Int64,
		f::Function,f_jac::Function,F_0_jac::Function,
		delta_k::Float64,N::Int64, dir::Int64)::Tuple{Array{Float64,2},Array{Float64,3}}

	println("===  Initiating One Step Process ===")
	println("Direction = ", dir)
	println("delta     = ",delta_k)
	psi_grid, glued_grid = convex_integrate_grid(GRID_SIZE,f,f_jac,F_0_jac, delta_k, N, dir)
	uniform_grid = create_uniform_grid(glued_grid,psi_grid,V_tuple[dir])
	return psi_grid, uniform_grid
end


"""
	three_step_algorith(...)

Perform three corrugations and save to a specified file format.
Note that this is done in two "stages" as in [B], where we reset
the initial map F0 to be F1 in the first stage. This is	done
to lower the needed corrugation numbers.
"""
function three_step_algorithm(GRID_SIZE::Int64,F_0::Function,
										F_0_jac::Function, start_direction::Int64; save_every_step = true, save_wrl = true, save_xml = false, xml_resolution = 1000, should_plot_flow = false)
	if GRID_SIZE < 600
		println("WARNING: A grid size this low might make
				  the integration hang!")
	end

	println("====================================")
	println("===     Initiating Algorithm     ===")
	println("====================================")
	println("Grid size               = ", GRID_SIZE)
	println("Saving after every step = ", save_every_step)
	println("Saving wrl file format	= ", save_wrl)
	println("Saving yafaray-xml      = ", save_xml)
	println()

	N::Int64 = 12
	j::Int64 = start_direction
	delt::Float64 = 0.268371 #This was found running the code from [B]

	if  should_plot_flow
		plot_flow(U_tuple[j],V_tuple[j],F_0_jac,GRID_SIZE)
	end
	psi_grid, uni_1_grid = do_one_step(GRID_SIZE,F_0,F_0_jac,F_0_jac,delt,N,j)
	
	if save_every_step  && save_wrl 
		write_wrl(uni_1_grid,"torus_first_corrugation.wrl")
	end
	if save_every_step  && save_xml
		write_xml(uni_1_grid,xml_resolution,"torus_first_corrugation.xml")
	end


	F_1,F_1_jac = create_interpolation_functions(uni_1_grid)
	N = 80
	j = mod(j,3)+1
	delt = 0.8 #Also found running the code from [B]

	if should_plot_flow
		plot_flow(U_tuple[j],V_tuple[j],F_1_jac,GRID_SIZE)
	end
	psi_grid, uni_2_grid = do_one_step(GRID_SIZE,F_1,F_1_jac,F_1_jac,delt,N,j)

	if save_every_step && save_wrl
		write_wrl(uni_2_grid,"torus_second_corrugation.wrl")
	end
	if save_every_step && save_xml
		write_xml(uni_2_grid,xml_resolution,"torus_second_corrugation.xml	")
	end

	F_2,F_2_jac = create_interpolation_functions(uni_2_grid)
	N = 500
	j = mod(j,3)+1
  
	if should_plot_flow
		plot_flow(U_tuple[j],V_tuple[j],F_1_jac,GRID_SIZE)
	end
	psi_grid, uni_3_grid = do_one_step(GRID_SIZE,F_2,F_2_jac,F_1_jac,delt,N,j)

	
	if save_wrl
		write_wrl(uni_3_grid,"torus_third_corrugation.wrl")
	end
	if save_xml
		write_xml(uni_3_grid,xml_resolution,"torus_third_corrugation.wrl")
	end
end
