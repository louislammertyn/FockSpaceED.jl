

function schrodinger_TD!(dψ, ψ, (fs, Ops), t)
    fill!(dψ, 0.0 + 0im)
    for (fk, Ok) in zip(fs, Ops)
        tmp_vec = Ok * ψ 
        dψ .+= -1im * fk(t) * tmp_vec
    end
    nothing
end

function schrodinger!(dψ, ψ, H, t)
    dψ[:] = -1im * H * ψ
    nothing
end

## This function takes as input the coefficients of the state we want to evolve
# as well as the Hamiltonian in the following time dependent form H=∑ᵢ fᵢ(t)⋅θᵢ where θᵢ is the appropriate operator in matrix form and fᵢ is a time dependent function
function Time_Evolution_TD(init::Vector{ComplexF64}, (fs, Ops), tspan; rtol = 1e-9, atol = 1e-9, solver=Vern7())
    prob = ODEProblem(schrodinger_TD!, init, tspan, (fs, Ops) )
    sol = solve(prob, solver; reltol=rtol, abstol=atol, save_everystep=false)
    return sol
end

function Time_Evolution(init::Vector{ComplexF64}, H, tspan; rtol = 1e-9, atol = 1e-9, solver=Vern7())
    prob = ODEProblem(schrodinger!, init, tspan,H)
    sol = solve(prob, solver; reltol=rtol, abstol=atol, save_everystep=false)
    return sol
end

