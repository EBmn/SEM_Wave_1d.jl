module MMS

export MMS_jet

mutable struct MMS_jet

    type::Array{Int64}
    dim::Int64
    coeff::Array{Float64}

    function MMS_jet(type, dim)

        coeff = ones(dim, 3*dim)

        new(type, dim, coeff)

    end


end


function MMSfun(x::Float64, idx::Int64, idim::Int, MMS)

    if MMS.type[idim] == 1
        # if trigonometric of type A*sin(B*x + C) in the direction given by idim

        i = mod(idx, 4)

        if (i == 0)

            u = MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*sin(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])

        elseif (i == 1)

            u = MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*cos(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])

        elseif (i == 2)

            u = -MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*sin(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])

        elseif (i == 3)

            u = -MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*cos(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])

        end


    elseif MMS.type[idim] == 2
        # if trigonometric of type A*cos(B*x + C) in the direction given by idim
    
        i = mod(idx, 4)

        if (i == 0)

            u = MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1] * cos(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])

        elseif (i == 1)

            u = -MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1] * sin(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])

        elseif (i == 2)
            
            u = -MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1] * cos(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])

        elseif (i == 3)

            u = MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1] * sin(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])

        end

    elseif MMS.type[idim] == 3
        # if of type A*(x-B)^C in the direction given by idim

        if (idx < MMS.coeff[idim, 3] + 1)

            #=
            println("dimension: " * string(idim))
            println(idx)
            
            println(MMS.coeff[idim, 1])
            println(MMS.coeff[idim, 2])
            println(MMS.coeff[idim, 3])
            =#
            
            u =  MMS.coeff[idim, 1] * (x - MMS.coeff[idim, 2])^(MMS.coeff[idim, 3] - idx) * (factorial(Int(round(MMS.coeff[idim, 3])))/factorial(Int(round(MMS.coeff[idim, 3])) - idx))

        else

            u = 0.0
        end
        
    elseif MMS.type[idim] == 4

        # if of type 1/2 + A * exp(-(x-B)^2/C)
        if (idx == 0)
            u = MMS.coeff[idim, 1]*exp(-(x - MMS.coeff[idim, 2])^2/MMS.coeff[idim, 3])
        elseif (idx == 1)
            u = (-2*(x - MMS.coeff[idim, 2])/MMS.coeff[idim, 3])*exp(-(x - MMS.coeff[idim, 2])^2/MMS.coeff[idim, 3])
        else 
            println("Warning! MMS derivative has not been calculated correctly")
            u = 0.0
        end

    end

    return u

end

function MMSfun(x::Float64, t::Float64, idx::Int64, idt::Int64, MMS)

    u = MMSfun(x, idx, 1, MMS)*MMSfun(t, idt, 2, MMS)        

    return u

end

function MMSfun(x::Float64, y::Float64, t::Float64, idx::Int64, idy::Int64, idt::Int64, MMS)

    u = MMSfun(x, idx, 1, MMS)*MMSfun(y, idy, 2, MMS)*MMSfun(t, idt, 3, MMS)        

    return u

end



end