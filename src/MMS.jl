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
        #if trigonometric of type A*sin(B*x + C)

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


    elseif MMS.type == 2
            #if trigonometric of type A*cos(B*x + C)*D*cos(E*t + F)
    
            i = mod(idx, 4)
    
            if (i == 0)
    
                #println(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
                u = MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*cos(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
    
            elseif (i == 1)
    
                u = -MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*sin(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
    
            elseif (i == 2)
    
                u = -MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*cos(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
    
            elseif (i == 3)
    
                u = MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*sin(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
    
            end
    
        elseif MMS.type == 3
            #if of type u(x, t) = A*cos(B*x + C)*t^q
    
            i = mod(idx, 4)

            if idim == 1

    
                if (i == 0)
        
                    #println(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
                    u = MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*cos(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
        
                elseif (i == 1)
        
                    u = -MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*sin(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
        
                elseif (i == 2)
        
                    u = -MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*cos(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
        
                elseif (i == 3)
        
                    u = MMS.coeff[idim, 2]^idx * MMS.coeff[idim, 1]*sin(MMS.coeff[idim, 2]*x + MMS.coeff[idim, 3])
        
                end
            
            else

                if (idx < MMS.coeff[idim, 1] + 1)

                    u =  x^(MMS.coeff[idim, 1] - idx) * (factorial(Int(floor(MMS.coeff[idim, 1])))/factorial(Int(floor(MMS.coeff[idim, 1])) - idx))

                else

                    u = 0.0

                end
            
            end

        end


    return u

end

function MMSfun(x::Float64, t::Float64, idx::Int64, idt::Int64, MMS)

    if MMS.type == 1

        u = MMSfun(x, idx, 1, MMS)*MMSfun(t, idt, 2, MMS)

    elseif MMS.type == 2

        u = MMSfun(x, idx, 1, MMS)*MMSfun(t, idt, 2, MMS)

    elseif MMS.type == 3

        u = MMSfun(x, idx, 1, MMS)*MMSfun(t, idt, 2, MMS)        

    end

    return u

end


end