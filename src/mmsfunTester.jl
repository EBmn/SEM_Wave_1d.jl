
include("MMS.jl")

using .MMS

types = [3 3 2]

MMS_j = MMS.MMS_jet(types, 3) # 3-dimensional MMS of type 1

xVals = collect(LinRange(-pi, pi, 100))
yVals = collect(LinRange(-pi, pi, 100))
tVals = collect(LinRange(0, 2*pi, 200))

xVals = collect(LinRange(-1.0, 1.0, 100))
yVals = collect(LinRange(-1.0, 1.0, 100))
tVals = collect(LinRange(0, 1, 200))

# dimension two, evaluated at xVals:

MMS_j.coeff[1, 1] = 1
MMS_j.coeff[1, 2] = 0
MMS_j.coeff[1, 3] = 6

MMS_j.coeff[2, 1] = 1
MMS_j.coeff[2, 2] = 0
MMS_j.coeff[2, 3] = 6

MMS_j.coeff[3, 1] = 1
MMS_j.coeff[3, 2] = 2*pi
MMS_j.coeff[3, 3] = 0
#MMS_j.type = 2
println(MMS_j.coeff)

funVals = zeros(length(xVals), length(yVals))

anim = Animation()

for k = 1:length(tVals)

    for i = 1:length(xVals)
        for j = 1:length(yVals)
            funVals[i, j] = MMS.MMSfun(xVals[i], yVals[end+1-j], tVals[k], 0, 0, 0, MMS_j)
        end
    end

    surface(xVals, yVals[end:-1:1], funVals, zlims = (-2, 2), legend=:false)

    frame(anim)
    println(tVals[k])

end

println(maximum(funVals))
println(minimum(funVals))

gif(anim, "mmsfunTester.gif", fps=40)
