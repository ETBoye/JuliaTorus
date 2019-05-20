println("Loading the package, please wait a few minutes...")
include("util.jl")
include("OneStepper.jl")


three_step_algorithm(
	600, #grid size
	standard_param(0.2,0.5), #initial map
	standard_param_jacobi(0.2,0.5), #initial map jacobi
	1, #start direction
	save_every_step=true,
	save_xml = false,
	save_wrl = true,
	should_plot_flow = false) 
