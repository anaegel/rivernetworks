--------------------------------------------------------------------------------
-- Example of a problem on 1D line segments in 2D
-- (C) G-CSC, Uni Frankfurt, 2021
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")


local gridName = util.GetParam("-grid", "simple-river-y.ugx")
local charLength = 3.0


local ARGS = {}
ARGS.eps = 1.0 -- diffusion constant

ARGS.dim = util.GetParamNumber("-dim", 2, "dimension")
ARGS.numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
ARGS.numRefs    = util.GetParamNumber("-numRefs",    5, "number of refinements")
ARGS.startTime  = util.GetParamNumber("-start", 0.0, "start time")
ARGS.endTime    = util.GetParamNumber("-end", charLength*charLength/ARGS.eps, "end time")
ARGS.dt 		   = util.GetParamNumber("-dt", ARGS.endTime*0.05, "time step size")

util.CheckAndPrintHelp("Time-dependent problem setup example");


print(" Choosen Parameter:")
print("    numRefs      = " .. ARGS.numRefs)
print("    numPreRefs   = " .. ARGS.numPreRefs)
print("    startTime 	= " .. ARGS.startTime)
print("    endTime 		= " .. ARGS.endTime)
print("    dt 			= " .. ARGS.dt)
print("    grid         = " .. gridName)

-- choose algebra
InitUG(ARGS.dim, AlgebraType("CPU", 1));

-- Create, Load, Refine and Distribute Domain
local mandatorySubsets = {"River", "Sink"}
local dom = nil
if withRedist == true then
	dom = util.CreateDomain(gridName, 0, mandatorySubsets)
	balancer.RefineAndRebalanceDomain(dom, ARGS.numRefs)
else
	dom = util.CreateAndDistributeDomain(gridName, ARGS.numRefs, ARGS.numPreRefs, mandatorySubsets)
end

print("\ndomain-info:")
print(dom:domain_info():to_string())

-- create Approximation Space
print(">> Create ApproximationSpace")
local approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)

-- lets order indices using Cuthill-McKee
OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

print (">> Setting up Assembling")

local elemDisc = ConvectionDiffusionFV1("c", "River")
elemDisc:set_diffusion(ARGS.eps) -- Diffusion
-- elemDisc:set_upwind(FullUpwind())
-- elemDisc:set_velocity("Velocity")

local dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "c", "Source1")
dirichletBND:add(2.0, "c", "Source2")
dirichletBND:add(0.0, "c", "Sink")

local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print (">> Setting up Algebra Solver")

-- if the non-linear problem shall should be solved, one has to wrap the solver
-- inside a newton-solver. See 'solver_util.lua' for more details.
solverDesc = {
	type = "linear",
	precond = {
		type 		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
		approxSpace = approxSpace,
		smoother 	= {
			type		= "jac",-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
			damping 	= 0.8
		},
		cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
		preSmooth	= 2,		-- number presmoothing steps
		postSmooth 	= 2,		-- number postsmoothing steps
		baseSolver 	= "lu"
	},
	convCheck = {
		type		= "standard",
		iterations	= 100,
		absolute	= 1e-9,
		reduction	= 1e-12,
		verbose=true
	}
}

solver = util.solver.CreateSolver(solverDesc)

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

-- set initial value
print(">> Interpolation start values")
local u = GridFunction(approxSpace)
Interpolate(0.0, u, "c", ARGS.startTime)

-- perform time loop
--util.SolveNonlinearTimeProblem(u, domainDisc, solver, VTKOutput(), "Sol",
--							   "ImplEuler", 1, startTime, endTime, dt); 

util.SolveLinearTimeProblem(u, domainDisc, solver, VTKOutput(), "Sol",
							"ImplEuler", 1, ARGS.startTime, ARGS.endTime, ARGS.dt); 

print("Writing profile data")
WriteProfileData("profile_data.pdxml")
util.PrintProfile_TotalTime("main ")

-- end group app_convdiff
--[[!
\}
]]--
