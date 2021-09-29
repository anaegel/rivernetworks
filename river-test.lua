--------------------------------------------------------------------------------
-- Example of a problem on 1D line segments in 2D
-- (C) G-CSC, Uni Frankfurt, 2021
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")


--local gridName = util.GetParam("-grid", "grids/simple-river-y.ugx")
local gridName = util.GetParam("-grid", "grids/simple-river.ugx")
local charLength = 3.0


local ARGS = {}
ARGS.eps = 1.0 -- diffusion constant

ARGS.dim = util.GetParamNumber("-dim", 2, "dimension")
ARGS.numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
ARGS.numRefs    = util.GetParamNumber("-numRefs",    6, "number of refinements") 
ARGS.startTime  = util.GetParamNumber("-start", 0.0, "start time")
ARGS.endTime    = util.GetParamNumber("-end", 0.2*charLength*charLength/ARGS.eps, "end time")
ARGS.dt 		   = util.GetParamNumber("-dt", ARGS.endTime*5e-4, "time step size")

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

-- Lexicographic order of indices. 
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


function sign(number)
    return number > 0 and 1 or (number == 0 and 0 or -1)
end



-- Hoehe (Rechteck)
local B  = 1.0 -- Breite
local Banf=1.2
local Bend=1.0 --1.0
local xeng=0
--function Width(x,y,t) return (Banf+Bend) / 2 + (Bend-Banf)/ 2*math.tanh(50*(x - xeng)) end
function Width(x,y,t) if (x>-0.25 and x < 0.25) then return 1.2 else return 1.0 end end
local width=LuaUserNumber2d("Width")

function Height(A,width) return A/width end
function Height_dA(A,width) return 1.0/width end
function Height_dw(A,width) return -A/(width*width) end

local height = LuaUserFunctionNumber("Height", 2)
height:set_input(0, elemDisc["A"]:value())
height:set_input(1, width)
height:set_deriv(0, "Height_dA")
height:set_deriv(1, "Height_dw")


local lambda = 0.000001
local dhyd = 1.0

function Roughness(v,A,width) 
  return lambda/(4.0*width*(A/width)/(width+2*(A/width)))*0.5*math.abs(v) 
end

function Roughness_dv(v,A,width)  
  return lambda/(4.0*width*(A/width)/(width+2*(A/width)))*0.5*sign(v)
end

function Roughness_dA(v,A,width)  
  return lambda/(4.0*width*width*width/((2*A + width*width)*(2*A + width*width)))*0.5*math.abs(v)
end

local roughness = LuaUserFunctionNumber("Roughness", 3)
roughness:set_input(0, elemDisc["v"]:value())
roughness:set_input(1, elemDisc["A"]:value())
roughness:set_input(2, width)
roughness:set_deriv(0, "Roughness_dv")
roughness:set_deriv(1, "Roughness_dA")



local nu_t = 0.5
local g = 9.81

function Sohle2d(x, y, t) return -0.003*x end
local sohle=LuaUserNumber2d("Sohle2d")

-- g*(h+zB)
--FRAGE: height und sohle plotten, gibt es dafür funktion?
local gefaelle = g*(height+sohle)

-- Fehlt: h, zB
elemDisc["v"]:set_upwind(upwind)
elemDisc["v"]:set_mass_scale(1.0) 
elemDisc["v"]:set_diffusion(nu_t)                          -- \partial_x (-nu_t \partial_x v)
elemDisc["v"]:set_velocity(0.5*elemDisc["v"]:value()*downStreamVec)      -- \partial_x (0.5*v*v)
elemDisc["v"]:set_flux(gefaelle*downStreamVec)                           -- \partial_x (g *(h+zB)) 
elemDisc["v"]:set_reaction_rate(roughness)



-- Add inflow bnd cond.
-- local function AddInflowBC(domainDisc, Q, h, B, subsetID)
local function AddInflowBC(domainDisc, A0, v0, subsetID)
  local dirichletBND = DirichletBoundary()
  dirichletBND:add(A0, "A", subsetID)
  dirichletBND:add(v0, "v", subsetID) 
  domainDisc:add(dirichletBND)
end

-- Add outflow bnds.
local function AddOutflowBC(domainDisc, pointID, subsetID)
  local outflowBnd = {}

  outflowBnd["A"] = NeumannBoundaryFV1("A")      --  v*A
  outflowBnd["A"]:add(elemDisc["A"]:value()*elemDisc["v"]:value(), pointID, subsetID)

  outflowBnd["v"] = NeumannBoundaryFV1("v")      -- 0.5*v^2 + g*(h+z)
  outflowBnd["v"]:add(gefaelle, pointID, subsetID)
  outflowBnd["v"]:add(0.5*elemDisc["v"]:value()*elemDisc["v"]:value(), pointID, subsetID)

  domainDisc:add(outflowBnd["A"]) 
  domainDisc:add(outflowBnd["v"]) 
end


-- Create discretization.
local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc["A"])
domainDisc:add(elemDisc["v"])
AddInflowBC(domainDisc, 11.4, 1.12, "Source1")
AddOutflowBC(domainDisc, "Sink", "River")

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print (">> Setting up Algebra Solver")

-- if the non-linear problem shall should be solved, one has to wrap the solver
-- inside a newton-solver. See 'solver_util.lua' for more details.
solverDesc = 
{
  type = "newton",
  
 
  linSolver = 
  {
	 type = "bicgstab",
	 
	  precond = {
	   type = "ilu",
	   -- sort = true,
		},
	
		convCheck = {
			type		= "standard",
			iterations	= 100,
			absolute	= 1e-11,
			reduction	= 1e-13,
			verbose=true
		},
	
	},
	
  lineSearch =
  {
    type = "standard",
      maxSteps    = 10,
      lambdaStart   = 1.0,
      lambdaReduce  = 0.5,
      acceptBest    = false,
      checkAll    = false
  },	

}

local dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

local solver = util.solver.CreateSolver(solverDesc)
--solver:set_debug(dbgWriter)

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

-- Set initial value.
print(">> Interpolation start values")
local u = GridFunction(approxSpace)
Interpolate(1.0,  u, "v", ARGS.startTime)
Interpolate(10.0, u, "A", ARGS.startTime)

-- Configure VTK output.
local vtk = VTKOutput() 
vtk:select_element(height, "Height")
vtk:select_element(sohle, "Sohle")
vtk:select_element(gefaelle, "Gefaelle")
vtk:select("v", "velocity")
vtk:select("A", "area")

-- Perform time stepping loop.
util.SolveNonlinearTimeProblem(u, domainDisc, solver, vtk , "Sol_change_b",
							   "ImplEuler", 1, ARGS.startTime, ARGS.endTime, ARGS.dt, ARGS.dt * 1e-6); 

--util.SolveLinearTimeProblem(u, domainDisc, solver, VTKOutput(), "Sol",
--							"ImplEuler", 1, ARGS.startTime, ARGS.endTime, ARGS.dt); 

print("Writing profile data")
WriteProfileData("profile_data.pdxml")
util.PrintProfile_TotalTime("main ")

-- end group app_convdiff
--[[!
\}
]]--
