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
ARGS.numRefs    = util.GetParamNumber("-numRefs",    2, "number of refinements")
ARGS.startTime  = util.GetParamNumber("-start", 0.0, "start time")
ARGS.endTime    = util.GetParamNumber("-end", 0.1*charLength*charLength/ARGS.eps, "end time")
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
InitUG(ARGS.dim, AlgebraType("CPU", 2));

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
approxSpace:add_fct("A", "Lagrange", 1)
approxSpace:add_fct("v", "Lagrange", 1)

-- lets order indices using Cuthill-McKee
OrderLex(approxSpace, "x");

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

print (">> Setting up Assembling")

function U(x,y,t)
  return 0.0
end

local upwind = FullUpwind()
local elemDisc = {}
elemDisc["A"] = ConvectionDiffusionFV1("A", "River")
elemDisc["v"] = ConvectionDiffusionFV1("v", "River")


local downStreamVec = EdgeOrientation(dom)

-- Gleichung für A
elemDisc["A"]:set_upwind(upwind)
elemDisc["A"]:set_mass_scale(1.0)                   -- \partial A / \partial t
elemDisc["A"]:set_velocity(elemDisc["v"]:value() * downStreamVec)   -- \partial_x (v*A)

--
local lambda = 0.03
local dhyd = 1.0

function sign(number)
    return number > 0 and 1 or (number == 0 and 0 or -1)
end



-- Hoehe (Rechteck)
local B  = 1.0 -- Breite
function Height(A) return A/B end
function Height_dA(A) return 1.0/B end

local height = LuaUserFunctionNumber("Height", 1)
height:set_input(0, elemDisc["A"]:value())
height:set_deriv(0, "Height_dA")


-- Rauheit
function Roughness(v) 
  return lambda/dhyd*0.5*math.abs(v) 
end

function Roughness_dv(v)  
  return lambda/dhyd*0.5*sign(v)
end

local roughness = LuaUserFunctionNumber("Roughness", 1)
roughness:set_input(0, elemDisc["v"]:value())
roughness:set_deriv(0, "Roughness_dv")




local nu_t = 1.0
local g = 9.81

function Sohle2d(x, y, t) return -x end
local sohle=LuaUserNumber2d("Sohle2d")

-- g*(h+zB)
local gefaelle = ScaleAddLinkerVector()
gefaelle:add(g*height,downStreamVec)
gefaelle:add(g*sohle,downStreamVec)

-- Fehlt: h, zB
elemDisc["v"]:set_upwind(upwind)
elemDisc["v"]:set_mass_scale(1.0) 
elemDisc["v"]:set_diffusion(nu_t)                          -- \partial_x (-nu_t \partial_x v)
elemDisc["v"]:set_velocity(0.5*elemDisc["v"]:value()*downStreamVec)      -- \partial_x (0.5*v*v)
elemDisc["v"]:set_flux(gefaelle)                           -- \partial_x (g *(h+zB)) 
elemDisc["v"]:set_reaction_rate(roughness)



local dirichletBND = DirichletBoundary()
dirichletBND:add(10.0, "A", "Source1")
dirichletBND:add(10.0, "A", "Source2")
dirichletBND:add(20.0, "A", "Sink")

dirichletBND:add(0.0, "v", "Source1")
dirichletBND:add(0.0, "v", "Source2")
dirichletBND:add(2.0, "v", "Sink")


local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc["A"])
domainDisc:add(elemDisc["v"])
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print (">> Setting up Algebra Solver")

-- if the non-linear problem shall should be solved, one has to wrap the solver
-- inside a newton-solver. See 'solver_util.lua' for more details.
solverDesc = 
{
  type = "newton",
  
  

  linSolver = {
	 type = "superlu",
	 
	 -- precond = "ilu",
	--[[ precond = {
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
	},--]]
	
	convCheck = {
		type		= "standard",
		iterations	= 100,
		absolute	= 1e-9,
		reduction	= 1e-12,
		verbose=true
	}
	
}
}

local dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

local solver = util.solver.CreateSolver(solverDesc)
solver:set_debug(dbgWriter)

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

-- set initial value
print(">> Interpolation start values")
local u = GridFunction(approxSpace)
Interpolate(1.0,  u, "v", ARGS.startTime)
Interpolate(10.0, u, "A", ARGS.startTime)

-- perform time loop
util.SolveNonlinearTimeProblem(u, domainDisc, solver, VTKOutput(), "Sol",
							   "ImplEuler", 1, ARGS.startTime, ARGS.endTime, ARGS.dt); 

--util.SolveLinearTimeProblem(u, domainDisc, solver, VTKOutput(), "Sol",
--							"ImplEuler", 1, ARGS.startTime, ARGS.endTime, ARGS.dt); 

print("Writing profile data")
WriteProfileData("profile_data.pdxml")
util.PrintProfile_TotalTime("main ")

-- end group app_convdiff
--[[!
\}
]]--
