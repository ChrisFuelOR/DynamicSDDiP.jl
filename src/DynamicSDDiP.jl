# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

module DynamicSDDiP

#import JuMP
import SDDP
import Gurobi
import Revise
import TimerOutputs
import GAMS
#import SCIP
import MathOptInterface
using Printf
using Dates

import Reexport
Reexport.@reexport using JuMP

using Infiltrator

# Write your package code here.
include("state.jl")
include("JuMP.jl")

include("typedefs.jl")
include("logging.jl")

include("stopping.jl")
include("objective.jl")
include("bellman.jl")
include("cutSelection.jl")
include("solverHandling.jl")
include("binarization.jl")
include("binaryRefinement.jl")
include("regularizations.jl")
include("sigmaTest.jl")

include("algorithmMain.jl")
include("forwardPass.jl")
include("backwardPass.jl")
include("duals.jl")
include("lagrange.jl")

end
