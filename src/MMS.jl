module MMS

export MMS_jet

mutable struct MMS_jet

    type::Int64
    dim::Int64
    coeff::Array{Float64}

    function MMS_jet(type, dim)

        coeff = ones(dim, 3*dim)

        new(type, dim, coeff)

    end


end


function MMSfun(x::Float64, idx::Int64, idim::Int, MMS)

    if MMS.type == 1
        #if trigonometric

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

    end

    return u

end

function MMSfun(x::Float64, t::Float64, idx::Int64, idt::Int64, MMS)

    if MMS.type == 1

        u = MMSfun(x, idx, 1, MMS)*MMSfun(t, idt, 2, MMS)

    end

    return u

end


end