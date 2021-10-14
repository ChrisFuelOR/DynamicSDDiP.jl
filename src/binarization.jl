# The function
# > "setup_state_backward",
# is derived from function "setup_state" in the 'SDDP.jl' package by
# Oscar Dowson and Lea Kapelevich released under the Mozilla Public License 2.0.
# The reproduced function and other functions in this file are also released
# under Mozilla Public License 2.0

# Copyright (c) 2021 Christian Fuellner <christian.fuellner@kit.edu>
# Copyright (c) 2021 Oscar Dowson, Lea Kapelevich.

# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

"""
Modifying the backward pass problem to include a binary expansion of the state
if such expansion is used.
"""
function changeStateSpace!(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    state::Dict{Symbol,Float64},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    ############################################################################
    # INITIALIZATION
    ############################################################################

    # The copy constraint is not modeled explicitly here. Instead,
    # the state variable is unfixed and takes the role of z in our paper.
    # It is then subtracted from the fixed value to obtain the so called slack.

    bw_data = node.ext[:backward_data]
    bw_data[:fixed_state_value] = Dict{Symbol,Float64}()
    bw_data[:bin_constraints] = JuMP.ConstraintRef[]
    bw_data[:bin_states] = Dict{Symbol,JuMP.VariableRef}()
    bw_data[:bin_x_names] = Dict{Symbol,Symbol}()
    bw_data[:bin_k] = Dict{Symbol,Int64}()

    number_of_states = 0
    binary_precision = state_approximation_regime.binary_precision

    ############################################################################
    # PREPARE NEW STATE VARIABLES
    ############################################################################
    for (state_name, value) in state
        # Get actual state from state_name
        state_comp = node.states[state_name]

        # Save fixed value of state
        fixed_value = JuMP.fix_value(state_comp.in)
        bw_data[:fixed_state_value][state_name] = fixed_value
        beta = binary_precision[state_name]

        # Set up state for backward pass using binary approximation
        setup_state_binarization!(subproblem, state_comp, state_name, beta, bw_data)
    end

    return
end


"""
Trivial modification of the backward pass problem if no binary approximation
is used.
"""
function changeStateSpace!(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    state::Dict{Symbol,Float64},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    return
end


"""
Modifying the backward pass problem to regain the original states
if a binary approximation was used.
"""
function rechangeStateSpace!(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    state::Dict{Symbol,Float64},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation
    )

    bw_data = node.ext[:backward_data]

    ############################################################################
    # FIX THE STATE VARIABLES AGAIN
    ############################################################################
    for (state_name, value) in state
        state_comp = node.states[state_name]

        JuMP.delete_lower_bound(state_comp.in)
        JuMP.delete_upper_bound(state_comp.in)

        # unset binary or integer type
        if JuMP.is_binary(state_comp.in)
            JuMP.unset_binary(state_comp.in)
        elseif JuMP.is_integer(state_comp.in)
            JuMP.unset_integer(state_comp.in)
        end

        JuMP.fix(state_comp.in, bw_data[:fixed_state_value][state_name])
    end

    ############################################################################
    # REPLACE THE NEW BY THE OLD OBJECTIVE
    ############################################################################
    """ Note that this is already done in Kelley's method """

    ############################################################################
    # DELETE ALL BINARY SPACE BASED VARIABLES AND CONSTRAINTS
    ############################################################################
    for (name, bin_state) in bw_data[:bin_states]
        JuMP.delete(subproblem, bin_state)
    end

    for constraint in bw_data[:bin_constraints]
        JuMP.delete(subproblem, constraint)
    end

    delete!(node.ext, :backward_data)

    return

end


"""
Trivial modification of the backward pass problem to regain the original states
if no binary approximation was used.
"""
function rechangeStateSpace!(
    node::SDDP.Node,
    subproblem::JuMP.Model,
    state::Dict{Symbol,Float64},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation
    )

    return
end


"""
Setting up the binary state variables.
"""
function setup_state_binarization!(
    subproblem::JuMP.Model,
    state_comp::State,
    state_name::Symbol,
    binary_precision::Float64,
    bw_data::Dict{Symbol,Any}
)

    # Get name of state variable in String representation
    name = JuMP.name(state_comp.in)

    ############################################################################
    # STATE IS ALREADY BINARY
    ############################################################################
    if state_comp.info.in.binary
        """ In this case, the variable must not be unfixed and, in principle,
        no new variables or constraints have to be introduced.

        Still, we introduce a new binary variable (binary_var), as this one
        is used when solving the Lagrangian dual.

        If we would not introduce a new variable, but just store the state
        itself per reference in binary_vars, then it would be deleted later.
        """

        ########################################################################
        # INTRODUCE ONE NEW BINARY VARIABLE TO THE PROBLEM
        ########################################################################
        binary_var = JuMP.@variable(
            subproblem,
            base_name = "bin_" * name,
        )
        subproblem[:binary_vars] = binary_var
        # store in list for later access and deletion
        sym_name = Symbol(JuMP.name(binary_var))
        bw_data[:bin_states][sym_name] = binary_var
        bw_data[:bin_x_names][sym_name] = state_name
        bw_data[:bin_k][sym_name] = 1

        ########################################################################
        # INTRODUCE BINARY EXPANSION CONSTRAINT TO THE PROBLEM
        ########################################################################
        binary_constraint = JuMP.@constraint(subproblem, state_comp.in == binary_var)
        # store in list for later access and deletion
        push!(bw_data[:bin_constraints], binary_constraint)

        ########################################################################
        # FIX NEW VARIABLE
        ########################################################################
        # Get fixed values from fixed value of original state
        fixed_binary_value = JuMP.fix_value(state_comp.in)
        # Fix binary variables
        #JuMP.unset_binary(binary_var.in)
        JuMP.fix(binary_var, fixed_binary_value)

        ########################################################################
        # UNFIX ORIGINAL STATE
        ########################################################################
        # Unfix the original state
        JuMP.unfix(state_comp.in)
        follow_state_unfixing!(state_comp)

    else
        if !isfinite(state_comp.info.in.upper_bound) || !state_comp.info.in.has_ub
            error("When using SDDiP, state variables require an upper bound.")
        end

        if state_comp.info.in.integer
            ####################################################################
            # STATE VARIABLE IS INTEGER
            ####################################################################
            """
            Note that we do not need to define the new "binary" variables as
            binary, as they are fixed (to 0 or 1) or relaxed to [0,1] anyway.
            Moreover, we do not need to define them as state variables.
            """

            ####################################################################
            # INTRODUCE BINARY VARIABLES TO THE PROBLEM
            ####################################################################
            num_vars = SDDP._bitsrequired(state_comp.info.in.upper_bound)

            binary_vars = JuMP.@variable(
                subproblem,
                [i in 1:num_vars],
                base_name = "bin_" * name,
            )
            subproblem[:binary_vars] = binary_vars

            # store in list for later access and deletion
            for i in 1:num_vars
               sym_name = Symbol(JuMP.name(binary_vars[i]))
               bw_data[:bin_states][sym_name] = binary_vars[i]
               bw_data[:bin_x_names][sym_name] = state_name
               bw_data[:bin_k][sym_name] = i
           end

            ####################################################################
            # INTRODUCE BINARY EXPANSION CONSTRAINT TO THE PROBLEM
            ####################################################################
            binary_constraint = JuMP.@constraint(
                subproblem,
                state_comp.in == SDDP.bincontract([binary_vars[i] for i in 1:num_vars])
            )
            # store in list for later access and deletion
            push!(bw_data[:bin_constraints], binary_constraint)

            ####################################################################
            # FIX NEW VARIABLES
            ####################################################################
            # Get fixed values from fixed value of original state
            fixed_binary_values = SDDP.binexpand(bw_data[:fixed_state_value][state_name], state_comp.info.in.upper_bound)
            # Fix binary variables
            for i in 1:num_vars
                #JuMP.unset_binary(binary_vars[i].in)
                JuMP.fix(binary_vars[i], fixed_binary_values[i])
            end

            ####################################################################
            # UNFIX ORIGINAL STATE
            ####################################################################
            # Unfix the original state
            JuMP.unfix(state_comp.in)
            follow_state_unfixing!(state_comp)

        else
            ####################################################################
            # STATE VARIABLE IS CONTINUOUS
            ####################################################################
            beta = binary_precision

            """ see comment above for integer case """

            ####################################################################
            # INTRODUCE BINARY VARIABLES TO THE PROBLEM
            ####################################################################
            num_vars = SDDP._bitsrequired(round(Int, state_comp.info.in.upper_bound / beta))

            binary_vars = JuMP.@variable(
                subproblem,
                [i in 1:num_vars],
                base_name = "bin_" * name,
            )

            subproblem[:binary_vars] = binary_vars
            # store in list for later access and deletion
            for i in 1:num_vars
                sym_name = Symbol(JuMP.name(binary_vars[i]))
                # Store binary state reference for later
                bw_data[:bin_states][sym_name] = binary_vars[i]
                bw_data[:bin_x_names][sym_name] = state_name
                bw_data[:bin_k][sym_name] = i
            end
            subproblem[:binary_vars] = binary_vars

            ####################################################################
            # INTRODUCE BINARY EXPANSION CONSTRAINT TO THE PROBLEM
            ####################################################################
            binary_constraint = JuMP.@constraint(
                subproblem,
                state_comp.in == SDDP.bincontract([binary_vars[i] for i in 1:num_vars], beta)
            )

            # store in list for later access and deletion
            push!(bw_data[:bin_constraints], binary_constraint)

            ####################################################################
            # FIX NEW VARIABLES
            ####################################################################
            # Get fixed values from fixed value of original state
            fixed_binary_values = SDDP.binexpand(bw_data[:fixed_state_value][state_name], state_comp.info.in.upper_bound, beta)
            # Fix binary variables
            for i in 1:num_vars
                #JuMP.unset_binary(binary_vars[i].in)
                JuMP.fix(binary_vars[i], fixed_binary_values[i])
            end

            ####################################################################
            # UNFIX ORIGINAL STATE
            ####################################################################
            # Unfix the original state
            JuMP.unfix(state_comp.in)
            follow_state_unfixing!(state_comp)
        end
    end

    return
end


"""
Determining the anchor points in the original space if BinaryApproximation
is used.
"""
function determine_anchor_states(
    node::SDDP.Node,
    outgoing_state::Dict{Symbol,Float64},
    state_approximation_regime::DynamicSDDiP.BinaryApproximation,
)

    anchor_states = Dict{Symbol,Float64}()
    for (name, value) in outgoing_state
        state_comp = node.states[name]
        beta = state_approximation_regime.binary_precision[name]
        (approx_state_value, )  = determine_anchor_state(state_comp, value, beta)
        anchor_states[name] = approx_state_value
    end

    return anchor_states
end


"""
Determining a single anchor state in the original space.
"""
function determine_anchor_state(
    state_comp::State,
    state_value::Float64,
    binaryPrecision::Float64,
)

    if state_comp.info.out.binary
        fixed_binary_values = state_value
        approx_state_value = state_value
    else
        if !isfinite(state_comp.info.out.upper_bound) || !state_comp.info.out.has_ub
            error("When using SDDiP, state variables require an upper bound.")
        end

        if state_comp.info.out.integer
            fixed_binary_values = SDDP.binexpand(state_value, state_comp.info.out.upper_bound)
            approx_state_value = SDDP.bincontract(fixed_binary_values)
        else
            fixed_binary_values = SDDP.binexpand(state_value, state_comp.info.out.upper_bound, binaryPrecision)
            approx_state_value = SDDP.bincontract(fixed_binary_values, binaryPrecision)
        end
    end
    return approx_state_value, fixed_binary_values
end


"""
Determining the anchor points in the original space if no state approximation
is used.
"""
function determine_anchor_states(
    node::SDDP.Node,
    outgoing_state::Dict{Symbol,Float64},
    state_approximation_regime::DynamicSDDiP.NoStateApproximation,
)

    anchor_states = Dict{Symbol,Float64}()
    return anchor_states
end
