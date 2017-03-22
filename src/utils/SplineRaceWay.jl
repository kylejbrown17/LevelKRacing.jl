# module SplineRaceWay

# using AutomotiveDrivingModels
# using NearestNeighbors
# using SplineUtils

# export Raceway

type Raceway
    roadway::AutomotiveDrivingModels.Roadway
    tree::NearestNeighbors.KDTree
    obstacleMap::Dict
    models::Dict{Int,AutomotiveDrivingModels.DriverModel}
    Δs::Float64
    x::Array{Float64,1}
    y::Array{Float64,1}
    θ::Array{Float64,1}
    k::Array{Float64,1}
    s::Array{Float64,1}

    function Raceway(Pts::Array,degree::Int,num_points::Int,num_samples::Int,lane_width::Float64)
        retval = new()

        T, tt, rx, ry = ClosedB_Spline(Pts, degree, num_points)
        ṙx, ṙy = B_SplineDerivative(T,tt,Pts,degree)
        θ = atan2(ṙy,ṙx) # unit tangent vector
        s = zeros(size(rx))
        s[2:end] = cumsum(sqrt(diff(rx).^2 + diff(ry).^2))
        k = diff(θ)./diff(s) # curvature
        x,y,θ,s,k = ResampleSplineEven(rx,ry,θ,s,k,num_samples)
        retval.roadway, retval.tree = GenSplineRoadway(x,y,θ,s,k,lane_width)
        retval.x = x
        retval.y = y
        retval.θ = θ
        retval.s = s
        retval.k = k
        retval.Δs = s[2] - s[1]
        retval.obstacleMap = Dict()
        retval.models = Dict{Int, DriverModel}()

        retval
    end
end # type

# end # module
