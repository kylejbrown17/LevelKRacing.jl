# module SplineUtils
# using AutomotiveDrivingModels
# using NearestNeighbors
#
# import PyPlot

# export
#     ClosedB_Spline,
#     B_SplineDerivative,
#     ResampleSplineEven,
#     GenSplineRoadway,
#     PlotSplineRoadway


"""
Computes the coefficients of a B-spline of degree "degree" from the control points "Pts", and
samples "L_tt" points from that spline
"""
function ClosedB_Spline(Pts, degree, L_tt) # returns spline and derivative of spline
    k = degree
    p = k - 1 # degree of derivative
    n = size(Pts,2) # number of control points
    m = n + k + 1; # number of knots in T

    #knots
    T = 0.5*ones(m,1)
    T[1:k] = 0
    T[(k+1):end-(k)] = linspace(0,1,m-2*k)
    T[end-(k):end] = 1.0

    tt = linspace(0,1,L_tt)
    rx = zeros(size(tt))
    ry = zeros(size(tt))

    invd = Dict() # inverse denominator (to deal with zero denominators)
    for i = 1:k+1
        invd[i] = T[(i+1):end] - T[1:(end-(i))]
        # take care of empty knot spans
        for j = 1:size(invd[i])[1]
            if invd[i][j] != 0
                invd[i][j] = 1.0/invd[i][j]
            end
        end
    end

    N = Dict()
    for j = 1:size(tt)[1]
        t = tt[j]
        N[1] = 1.0*(t .>= T[1:end-1]) .* (t .< T[2:end])
        for i = 1:k+1
            if i > 1
                N[i] = (t - T[1:end-i]) .* (invd[i-1].*N[i-1])[1:end-1] + (T[i+1:end] - t) .*(invd[i-1].*N[i-1])[2:end]
            end
        end
        rx[j] = (Pts[1,:]'*N[k+1])[1]
        ry[j] = (Pts[2,:]'*N[k+1])[1]
    end

    rx[end] = Pts[1,end]
    ry[end] = Pts[2,end]

    return T, tt, rx, ry
end

#####################################################
"""
Computes the derivative of a spline
"""
function B_SplineDerivative(T, # knots
    tt, # sample locations t ∈ [0,1.0]
    Pts, # control points
    k) # degree of spline

    p = k - 1 # degree of derivative
    n = size(Pts,2) # number of control points
    m = n + k + 1; # number of knots in T

    invd = Dict() # inverse denominator (to deal with zero denominators)
    for i = 1:k+1
        invd[i] = T[(i+1):end] - T[1:(end-(i))]
        # take care of empty knot spans
        for j = 1:size(invd[i])[1]
            if invd[i][j] != 0
                invd[i][j] = 1.0/invd[i][j]
            end
        end
    end

    Q = zeros(size(Pts[:,1:end-1]))
    ṙx = zeros(size(tt)) # derivatives
    ṙy = zeros(size(tt)) # derivatives

    # Derivative
    Q[1,:] = p*invd[1][p+1:end-p].*(Pts[1,2:end] - Pts[1,1:end-1])
    Q[2,:] = p*invd[1][p+1:end-p].*(Pts[2,2:end] - Pts[2,1:end-1])

    N = Dict()
    for j = 1:size(tt)[1]
        t = tt[j]
        N[1] = 1.0*(t .>= T[1:end-1]) .* (t .< T[2:end])
        for i = 1:k+1
            if i > 1
                N[i] = (t - T[1:end-i]) .* (invd[i-1].*N[i-1])[1:end-1] + (T[i+1:end] - t) .*(invd[i-1].*N[i-1])[2:end]
            end
        end

        ṙx[j] = (Q[1,:]'*N[k][2:end-1])[1] # derivatives
        ṙy[j] = (Q[2,:]'*N[k][2:end-1])[1] # derivatives
    end

    ṙx[1] = ṙx[2]
    ṙx[end] = ṙx[end-1]
    ṙy[1] = ṙy[2]
    ṙy[end] = ṙy[end-1]

    return ṙx, ṙy
end

#################################################
"""
Evenly resamples spline at num_samples points along arc_length
"""
function ResampleSplineEven(rx,ry,θ,s,k,num_samples)
    s_span = linspace(0,s[end],num_samples+1)
    xP = zeros(num_samples)
    yP = zeros(num_samples)
    sP = zeros(num_samples)
    θP = zeros(num_samples)
    kP = zeros(num_samples)

    xP[1] = rx[1]
    # xP[end] = rx[1]
    yP[1] = ry[1]
    # yP[end] = ry[1]
    θP[1] = θ[1]
    # θP[end] = θ[1]
    sP[1] = s[1]
    # sP[end] = s[1]
    kP[1] = k[1]
    # kP[end] = k[1]

    j = 1
    for t in 2:length(s)
        if j > num_samples
            break
        end
        if (s[t] > s_span[j]) && (s[t-1] <= s_span[j])
            sP[j] = s_span[j]
            xP[j] = rx[t-1] + (rx[t] - rx[t-1])*(s_span[j] - s[t-1]) / (s[t] - s[t-1])
            yP[j] = ry[t-1] + (ry[t] - ry[t-1])*(s_span[j] - s[t-1]) / (s[t] - s[t-1])
            θP[j] = θ[t-1] + (θ[t] - θ[t-1])*(s_span[j] - s[t-1]) / (s[t] - s[t-1])
            kP[j] = k[t-1] + (k[t] - k[t-1])*(s_span[j] - s[t-1]) / (s[t] - s[t-1])
            j += 1
        end
    end

    return xP, yP, θP, sP, kP
end

##########################################
"""
Generates a roadway based on spline points (x,y,θ,s,k)
"""
function GenSplineRoadway(x,y,θ,s,k,lane_width)
    tree = KDTree([x';y'])
    nlanes = 1
    seg1 = RoadSegment(1, Array(Lane, nlanes))
    tag1=LaneTag(1,1)
    boundary_left = LaneBoundary(:solid, :white)
    boundary_right = LaneBoundary(:solid, :white)

    curvepts = Array(CurvePt, length(x))
    for i in 1:length(x)
        curvepts[i] = CurvePt(VecSE2(x[i],y[i],θ[i]), s[i], k[i], NaN)
    end

    curveind_lo = CurveIndex(1,0.0)
    curveind_hi = CurveIndex(length(curvepts)-1,1.0)

    seg1.lanes[1] = Lane(tag1, curvepts, width=lane_width,
                                  boundary_left=boundary_left, boundary_right=boundary_right,
                                  next = RoadIndex(curveind_lo, tag1),
                                  prev = RoadIndex(curveind_hi, tag1),
                                 )
    roadway = Roadway()
    push!(roadway.segments, seg1)

    return roadway, tree
end

#################################################
"""
Plots spline roadway points
"""
function PlotSplineRoadway(x,y,θ,lane_width)
    perp_lines1 = zeros(2,length(x))
    perp_lines2 = zeros(2,length(x))

    perp_lines1[1,:] = x + (lane_width/2.0)*sin(θ)
    perp_lines1[2,:] = y - (lane_width/2.0)*cos(θ)
    perp_lines2[1,:] = x - (lane_width/2.0)*sin(θ)
    perp_lines2[2,:] = y + (lane_width/2.0)*cos(θ)

    PyPlot.figure()
    PyPlot.scatter(x,y)
    PyPlot.plot(x,y)
    PyPlot.plot(perp_lines1[1,:],perp_lines1[2,:],color="green")
    PyPlot.plot(perp_lines2[1,:],perp_lines2[2,:],color="green")
    PyPlot.axis("equal")
    PyPlot.show()
end

# end # module
