# using AutomotiveDrivingModels
# using NearestNeighbors

# import AutomotiveDrivingModels: get_actions!, observe!, action_context, get_name
# import Base.rand
# import PyPlot

# export
#     HRHC,
#     curveDist,
#     wrap_to_π,
#     kdProject,
#     generateObstacleMap,
#     updateObstacleMap!,
#     generateMotionMap,
#     screenCollision,
#     tailgateAvoidance,
#     getSuccessorStates,
#     loopProjectionKD,
#     computeTrajectory,
#     screenTrajectory,
#     checkCollision,
#     calculateObjective,
#     plot_stϕ,
#     plotHRHCInfo,
#     plotObjectiveHorizon,
#     plotSplineRoadway

# utility functions
function curveDist(pt1::CurvePt, pt2::CurvePt)
    d = sqrt((pt1.pos.x - pt2.pos.x)^2 + (pt1.pos.y - pt2.pos.y)^2)
end
function wrap_to_π(θ)
    θ = θ - div(θ,2*Float64(π))*(2*Float64(π))
    θ = θ + (θ .< -π).*(2*Float64(π)) - (θ .> π).*(2*Float64(π))
end
function kdProject(x,y,θ,tree,roadway,hrhc)
    """
    project single (x,y,Θ) point to roadway spline using kdtree to find the nearest spline point
    """
    curve = roadway.segments[1].lanes[1].curve
    # Δs = roadway.segments[1].lanes[1].curve[2].s - roadway.segments[1].lanes[1].curve[1].s
    idx_list,dist = knn(tree,[x;y],1)
    idx = idx_list[1]
    idxA = idx-1
    idxB = idx+1

    if idx == length(curve)
        idxB = 1 # back to the beginning of the curve
    end
    if idx == 1
        idxA = length(curve)
    end
    dA = sqrt(sum(([curve[idxA].pos.x, curve[idxA].pos.y]-[x,y]).^2))
    dB = sqrt(sum(([curve[idxB].pos.x, curve[idxB].pos.y]-[x,y]).^2))
    if dA < dB
        idxB = idx
    else
        idxA = idx
    end

    # project
    vec1 = [curve[idxB].pos.x - curve[idxA].pos.x, curve[idxB].pos.y - curve[idxA].pos.y, 0]
    vec2 = [x - curve[idxA].pos.x, y - curve[idxA].pos.y, 0]
    idx_t = dot(vec2, vec1)/norm(vec1)^2

    pθ = curve[idxA].pos.θ + idx_t*(curve[idxB].pos.θ - curve[idxA].pos.θ)

    s = curve[idxA].s + idx_t*hrhc.Δs
    t = norm(vec2 - idx_t*vec1)*sign(sum(cross(vec1, vec2)))
    ϕ = wrap_to_π(θ - pθ)

    s,t,ϕ,idxA
end

# type Hierarchical Receding Horizion Controller
type HRHC <: DriverModel{AccelDesAng}
#     action_context::IntegratedContinuous
    car_ID::Int
    v_map
    δ_map
    motion_map # state changes associated with cmd = (v_command, δ_command) - this is the motion map
    v_cmd::Int # index of v in v_cmds
    δ_cmd::Int # index of δ in δ_cmds
    #car parameters with bicycle geometry model
    car_length::Float64 # wheel base
    wheel_base::Float64
    car_width::Float64
    ellipseA::Float64
    ellipseB::Float64

    a_step::Float64 # max acceleration (m/s)
    μ::Float64 # friction coefficient
    v_range # possible velocityies
    δ_range # possible steering angles
    successor_states # array of next states
    # current v, current δ
    v::Float64
    δ::Float64
    curve_ind::Int
    Δs::Float64
    # planning horizon
    h::Int
    Δt::Float64
    # logit level
    k::Int
    # maximum deviation from center of track (if |t| > T_MAX, car is out of bounds)
    T_MAX::Float64
    # Action = Next State
    action::VehicleState
    function HRHC(car_ID::Int,roadway;
        car_length::Float64=4.8,
        wheel_base::Float64=4.4,
        car_width::Float64=2.5,
        v::Float64=0.0,
        δ::Float64=0.0,
        h::Int=10,
        Δt::Float64=1.0/24,
        ΔV₊::Float64=1.55,
        ΔV₋::Float64=3.05,
        Δδ::Float64=Float64(π)/12,
        v_min::Float64=0.0,
        v_max::Float64=100.0,
        a_step::Float64=10.0,
        a_range=[-1,0,1],
        μ::Float64=20.0,
        g::Float64=9.81,
        δ_max::Float64=Float64(π)/6,
        δ_step::Float64=Float64(π)/128,
        k::Int=1
        )

        hrhc = new()

        hrhc.T_MAX=(roadway.segments[1].lanes[1].width - car_width)/2.0

        hrhc.car_ID = car_ID
        hrhc.car_length=car_length
        hrhc.wheel_base=wheel_base
        hrhc.car_width=car_width
        hrhc.h=h
        hrhc.Δt=Δt
        hrhc.a_step=a_step
        hrhc.μ=μ

        hrhc.k=k

        hrhc.motion_map, hrhc.v_map, hrhc.δ_map, hrhc.v_range, hrhc.δ_range = generateMotionMap(v_min,
            v_max,a_step,δ_max,δ_step,wheel_base,h,μ,Δt,a_range=a_range)

        hrhc.successor_states = zeros(size(hrhc.motion_map[1],1),size(hrhc.motion_map[1],2),3)

        hrhc.v=v
        hrhc.δ=δ
        hrhc.v_cmd = 1
        hrhc.δ_cmd = Int((length(hrhc.δ_range) - 1)/2)
        hrhc.curve_ind=1
        hrhc.Δs=roadway.segments[1].lanes[1].curve[2].s-roadway.segments[1].lanes[1].curve[1].s
#         hrhc.action_context=context
        hrhc.action = NextState(VehicleState(VecSE2(0,0,0),0.0))

        # calculate semimajor axes of bounding ellipse with minimal area (for collision checking)
        W = car_width/2.0
        L = car_length/2.0

        # use quadratic formula to find B^2
        a = 2
        b = -4*L - (W^2)*2*(L^2)
        c = 2*L^4

        B = sqrt((-b + sqrt(b^2 - 4*a*c))/(2*a))
        A = sqrt(((W)^2) / (1 - (L/(B))^2))

        hrhc.ellipseA = A
        hrhc.ellipseB = B

        hrhc
    end
end
function generateObstacleMap(scene, models)
    k = maximum([driver.k for (id, driver) in models])
    n = length(scene)
    h = maximum([driver.h for (id, driver) in models]) # h should be the same for all vehicles on the track
    obstacleDict = Dict()
    for level in 0:k
        obstacleDict[level] = Dict()
        for (id, driver) in models
            obstacleDict[level][id] = zeros(h,5) # x,y,θ,v,δ
        end
    end

    return obstacleDict
end
function updateObstacleMap!(obstacleMap, level, car_ID, trajectory)
    obstacleMap[level][car_ID][1:size(trajectory,1),1:size(trajectory,2)] = trajectory
end
function generateMotionMap(v_min,v_max,a_step,δ_max,δ_step,wheel_base,h,μ,Δt;g=9.81,a_range=[-1,0,1])

    v_step = a_step*Δt
    v_range = linspace(v_min,v_max,round(v_max/v_step))
    δ_range = linspace(-δ_max,δ_max,2*round((δ_max)/δ_step)+1)
    δ_range = (1/maximum(δ_range))*δ_range.*abs(δ_range) # concentrate toward the middle
    R_range = wheel_base./tan(δ_range)
    # a_range = [-1,0,1]
    motion_map = Dict()
    v_map = Dict()
    δ_map = Dict()

    for v_idx in 1:length(v_range)
        v = v_range[v_idx]
        v_map[v_idx] = zeros(Int,length(a_range),length(δ_range))
        δ_map[v_idx] = zeros(Int,length(a_range),length(δ_range))
        for a_idx in 1:length(a_range)
            Δv_idx = a_range[a_idx]
            tire_force = sqrt(((v^2)./abs(R_range)).^2 + (Δv_idx*a_step)^2)
            if v_idx + Δv_idx >= 1 && v_idx + Δv_idx <= length(v_range)
                v_map[v_idx][a_idx,:] = Int(v_idx + Δv_idx)# v_range[v_idx+Δv_idx]
            else
                v_map[v_idx][a_idx,:] = Int(v_idx) #v_range[v_idx]
            end
            if v_idx > 1
                for j in 1:length(δ_range)
                    if tire_force[j] >= μ*g
                        v_map[v_idx][a_idx,j] = Int(v_idx) - 1# v_range[v_idx-1] # decelerate
                    end
                end
            end
            #δ_map[v_idx][a_idx,:] = δ_range.*(tire_force .< μ*g) +
             #   maximum(abs(δ_range[tire_force .< μ*g])).*sign(δ_range).*(tire_force .>= μ*g)
            δ_map[v_idx][a_idx,:] = [i for i in 1:length(δ_range)].*(tire_force .< μ*g) +
            indmax(abs(δ_range).*(tire_force .< μ*g)).*(tire_force .>= μ*g).*-sign(δ_range) +
            (length(δ_range)+1).*(sign(δ_range) .> 0).*(tire_force .>= μ*g)
        end
    end


    # fill motion map
    for v_idx in 1:length(v_range)
        motion_map[v_idx] = zeros(length(a_range),length(δ_range),h,3)
        level_v_idx = copy(v_map[v_idx])
        level_δ_idx = copy(δ_map[v_idx])
        Δθ = zeros(length(a_range),length(δ_range))
        for i in 1:h
            radius = wheel_base./tan(δ_range[level_δ_idx])
            Δs = ((v_range[v_idx] + v_range[level_v_idx])./2.0)*Δt
            Δθ = Δs./radius
            ΔX = abs(radius) .* sin(abs(Δθ))
            ΔX[:,Int((length(δ_range)-1)/2)+1] = Δs[:,Int((length(δ_range)-1)/2)+1]
            ΔY = radius.*(1 - cos(Δθ))
            ΔY[:,Int((length(δ_range)-1)/2)+1] = 0
            if i == 1
                motion_map[v_idx][:,:,i,1] = ΔX
                motion_map[v_idx][:,:,i,2] = ΔY
                motion_map[v_idx][:,:,i,3] = Δθ
            else
                motion_map[v_idx][:,:,i,1] = motion_map[v_idx][:,:,i-1,1] +
                    ΔX.*cos(motion_map[v_idx][:,:,i-1,3]) - ΔY.*sin(motion_map[v_idx][:,:,i-1,3])
                motion_map[v_idx][:,:,i,2] = motion_map[v_idx][:,:,i-1,2] +
                    ΔX.*sin(motion_map[v_idx][:,:,i-1,3]) + ΔY.*cos(motion_map[v_idx][:,:,i-1,3])
                motion_map[v_idx][:,:,i,3] = Δθ + motion_map[v_idx][:,:,i-1,3]  # + motion_map[v_idx][:,:,i-1,3]
            end
            for j in 1:length(level_v_idx) # update to next velocity
                level_v_idx[j] = v_map[level_v_idx[j]][j]
                level_δ_idx[j] = δ_map[level_v_idx[j]][j]
            end
        end
    end
    return motion_map, v_map, δ_map, v_range, δ_range
end
function loopProjectionKD(hrhc,scene,roadway,tree)
    """
    projects all points in hrhc.successor_states to the kdtree representing
    the spline points along the centerline of roadway
    """
    curve = roadway.segments[1].lanes[1].curve

    s_grid = zeros(size(hrhc.successor_states,1),size(hrhc.successor_states,2))
    t_grid = zeros(size(s_grid))
    ϕ_grid = zeros(size(s_grid))
    idx_grid = zeros(Int,size(s_grid))

    pts = [reshape(hrhc.successor_states[:,:,1],size(hrhc.successor_states[:,:,1],1)*size(hrhc.successor_states[:,:,1],2),1)';
    reshape(hrhc.successor_states[:,:,2],size(hrhc.successor_states[:,:,2],1)*size(hrhc.successor_states[:,:,2],2),1)']
    idxs_list, _ = knn(tree,pts,1)
    idxs=reshape(idxs_list,size(hrhc.successor_states[:,:,2],1),size(hrhc.successor_states[:,:,2],2))


    for i in 1:size(s_grid,1)
        for j in 1:size(s_grid,2)
            idxA = idxs[i,j][1]-1
            idxB = idxs[i,j][1]+1
            if idxs[i,j][1] == length(curve)
                idxB = 1 # wrap to the beginning of the curve
            end
            if idxs[i,j][1] == 1
                idxA = length(curve) # wrap to the end of the curve
            end
            x = hrhc.successor_states[i,j,1]
            y = hrhc.successor_states[i,j,2]
            dA = sqrt(sum(([curve[idxA].pos.x, curve[idxA].pos.y]-[x,y]).^2))
            dB = sqrt(sum(([curve[idxB].pos.x, curve[idxB].pos.y]-[x,y]).^2))
            if dA < dB
                idxB = idxs[i,j][1]
            else
                idxA = idxs[i,j][1]
            end

            # project
            vec1 = [curve[idxB].pos.x - curve[idxA].pos.x, curve[idxB].pos.y - curve[idxA].pos.y, 0]
            vec2 = [x - curve[idxA].pos.x, y - curve[idxA].pos.y, 0]
            idx_t = dot(vec2, vec1)/norm(vec1)^2

            s_θ = curve[idxA].pos.θ + idx_t*(curve[idxB].pos.θ - curve[idxA].pos.θ)

            s_grid[i,j] = curve[idxA].s + idx_t*hrhc.Δs
            t_grid[i,j] = norm(vec2 - idx_t*vec1)*sign(sum(cross(vec1, vec2)))
            ϕ_grid[i,j] = wrap_to_π(hrhc.successor_states[i,j,3] - s_θ)
            idx_grid[i,j] = idxA
        end
    end
    # account for wrap-around
    s_grid[s_grid .< scene.vehicles[hrhc.car_ID].state.posF.s] += curve[end].s + hrhc.Δs

    return s_grid, t_grid, ϕ_grid
end
function screenTrajectory(trajectory, obstacleMap, scene, roadway, hrhc, tree, k_level)
    out_of_bounds = false
    collision_flag = false
    # check out of bounds
    for i in 1 : size(trajectory,1)
        x = trajectory[i,1]
        y = trajectory[i,2]
        θ = trajectory[i,3]
        s,t,ϕ = kdProject(x,y,θ,tree,roadway,hrhc)

        if abs(t) > hrhc.T_MAX
            out_of_bounds=true
            return out_of_bounds
        end
    end

    # collisionFlag = zeros(size(trajectory,1),1) # stores locations collisions
    threshold_dist = hrhc.car_length*4 # must be at least this close before we care to calculate collision cost

    # if k_level >= 1
    #     # for all other cars...
    #     for (id,traj) in obstacleMap[k_level - 1]
    #         if id != hrhc.car_ID
    #             state = scene.vehicles[hrhc.car_ID].state
    #             state2 = scene.vehicles[id].state
    #             diff = state.posG - state2.posG
    #             s1,_,_ = kdProject(state.posG.x,state.posG.y,state.posG.θ,tree,roadway,hrhc)
    #             s2,_,_ = kdProject(state2.posG.x,state2.posG.y,state2.posG.θ,tree,roadway,hrhc)
    #             if (norm([diff.x, diff.y]) < threshold_dist) && (s1 <= s2) # don't care if opponent is behind us
    #                 # R_idx = zeros(Int, size(trajectory,1),1) # to store the sorted indices of R
    #                 # ΔX = zeros(size(trajectory,1),1) # Δx, with opponent at origin
    #                 # ΔY = zeros(size(trajectory,1),1) # Δx, with opponent at origin
    #                 # Δθ = zeros(size(trajectory,1),1) # Δx, with opponent at origin
    #                 # for i in 1:hrhc.h
    #                     # pos = VecSE2(traj[i,1:3]) # x,y,θ of opponent at time step h
    #                 ΔX = trajectory[:,1] - traj[:,1] # Δx, with opponent at origin
    #                 ΔY = trajectory[:,2] - traj[:,2] # Δy with opponent at origin
    #                 Δθ = trajectory[:,3] - traj[:,3]  # Δθ with opponent at origin
    #                 # end
    #                 R = sqrt(ΔX.^2 + ΔY.^2)
    #                 R_idx = sortperm(R) # ordered by increasing distance
    #                 for idx in R_idx
    #                     # if collisionFlag[idx] == 1
    #                     #     continue
    #                     # end
    #                     if R[idx] > hrhc.car_length # no collision, and all other R values are greater
    #                         break
    #                     end
    #                     # R is less than hrhc.car_length
    #                     ψ = atan2(ΔY[idx],ΔX[idx]) - Δθ[idx]
    #                     r = W*L/sqrt((L*sin(ψ))^2 + (W*cos(ψ))^2) # radius of ellipse at given angle
    #                     if R[idx]-r < W*L/8 # collision
    #                         collision_flag = true
    #                         return collision_flag
    #                     end
    #                 end
    #             end
    #         end
    #     end
    # end

    return out_of_bounds || collision_flag
end
function getSuccessorStates(ΔXYθ, car_ID, h, scene::Scene)
    """ gets legal successor_states from motion primitives library """
    pos = scene.vehicles[car_ID].state.posG # global x,y,z of car

    ΔX = ΔXYθ[:,:,h,1] * cos(pos.θ) + ΔXYθ[:,:,h,2] * -sin(pos.θ)
    ΔY = ΔXYθ[:,:,h,1] * sin(pos.θ) + ΔXYθ[:,:,h,2] * cos(pos.θ)
    Δθ = ΔXYθ[:,:,h,3]

    successor_states = zeros(size(ΔXYθ[:,:,h,:]))
    successor_states[:,:,1] = ΔX + pos.x
    successor_states[:,:,2] = ΔY + pos.y
    successor_states[:,:,3] = Δθ + pos.θ

    return successor_states
end
function computeTrajectory(ΔXYθ, car_ID, scene, v_cmd, δ_cmd, h)
    pos = scene.vehicles[car_ID].state.posG

    trajectory = zeros(h,3)
    trajectory[:,1] = pos.x + ΔXYθ[v_cmd,δ_cmd,1:h,1]*cos(pos.θ) + ΔXYθ[v_cmd,δ_cmd,1:h,2]*-sin(pos.θ)
    trajectory[:,2] = pos.y + ΔXYθ[v_cmd,δ_cmd,1:h,1]*sin(pos.θ) + ΔXYθ[v_cmd,δ_cmd,1:h,2]*cos(pos.θ)
    trajectory[:,3] = pos.θ + ΔXYθ[v_cmd,δ_cmd,1:h,3]

    return trajectory
end
function checkCollisionElliptic(R,Δx,Δy,Δθ,L,W)
    b = W # buffer
    ζ = atan2(Δy,Δx)
    ψ = ζ - Δθ
    r1 = W*L/sqrt((L*sin(ζ))^2 + (W*cos(ζ))^2) # radius of ellipse at given angle
    r2 = W*L/sqrt((L*sin(ψ))^2 + (W*cos(ψ))^2)
    if R < r1 + r2 + b
        return true
    else
        return false
    end
end
function screenCollision(hrhc, obstacleMap, tree, roadway, scene, k_level)
    L = hrhc.ellipseB+2
    W = hrhc.ellipseA+2

    collisionFlag = zeros(size(hrhc.successor_states[:,:,1])) # stores locations collisions
    threshold_dist = hrhc.car_length*4 # must be at least this close before we care to calculate collision cost

    if k_level >= 1
        # for all other cars...
        for (id,trajectory) in obstacleMap[k_level - 1]
            if id != hrhc.car_ID
                state = scene.vehicles[hrhc.car_ID].state
                state2 = scene.vehicles[id].state
                diff = state.posG - state2.posG
                s1,_,_ = kdProject(state.posG.x,state.posG.y,state.posG.θ,tree,roadway,hrhc)
                s2,_,_ = kdProject(state2.posG.x,state2.posG.y,state2.posG.θ,tree,roadway,hrhc)
                if (norm([diff.x, diff.y]) < threshold_dist) && (s1 <= s2) # don't care if opponent is behind us
                    R_idx = zeros(Int, size(hrhc.successor_states,1)*size(hrhc.successor_states,2)) # to store the sorted indices of R
                    for i in 1:hrhc.h
                        pos = VecSE2(trajectory[i,1:3]) # x,y,θ of opponent at time step h
                        successor_states = getSuccessorStates(hrhc.motion_map[hrhc.v_cmd],hrhc.car_ID,i,scene)
                        ΔX = successor_states[:,:,1] - pos.x # Δx, with opponent at origin
                        ΔY = successor_states[:,:,2] - pos.y # Δy with opponent at origin
                        Δθ = successor_states[:,:,3] - pos.θ # Δθ with opponent at origin
                        R = sqrt(ΔX.^2 + ΔY.^2)
                        R_idx = sortperm(reshape(R,size(R,1)*size(R,2),1)[:])
                        for idx in R_idx
                            if collisionFlag[idx] == 1
                                continue
                            end
                            if R[idx] > hrhc.car_length # no collision, and all other R values are greater
                                break
                            end
                            # R is less than hrhc.car_length
                            if checkCollisionElliptic(R[idx],ΔX[idx],ΔY[idx],Δθ[idx],L,W)
                                collisionFlag[idx] = 1
                            end
                            # ψ = atan2(ΔY[idx],ΔX[idx]) - Δθ[idx]
                            # r = W*L/sqrt((L*sin(ψ))^2 + (W*cos(ψ))^2) # radius of ellipse at given angle
                            # if R[idx]-r < W*L/8 # collision
                            #     collisionFlag[idx] = 1
                            # end
                        end
                    end
                end
            end
        end
    end
    return collisionFlag
end
function tailgateAvoidance(hrhc, obstacleMap, tree, roadway, scene, k_level)
    L = 2*hrhc.car_length
    W = hrhc.car_width

    collisionFlag = zeros(size(hrhc.successor_states[:,:,1])) # stores locations collisions
    threshold_dist = hrhc.car_length*4 # must be at least this close before we care to calculate collision cost

    R = zeros(size(hrhc.successor_states[:,:,1]))
    ψ = zeros(size(hrhc.successor_states[:,:,1]))
    cost = zeros(size(hrhc.successor_states[:,:,1]))

    if k_level >= 1
        # for all other cars...
        for (id,trajectory) in obstacleMap[k_level - 1]
            if id != hrhc.car_ID && id == 3
                state = scene.vehicles[hrhc.car_ID].state
                state2 = scene.vehicles[id].state
                diff = state.posG - state2.posG
                s1,_,_ = kdProject(state.posG.x,state.posG.y,state.posG.θ,tree,roadway,hrhc)
                s2,_,_ = kdProject(state2.posG.x,state2.posG.y,state2.posG.θ,tree,roadway,hrhc)
                if (norm([diff.x, diff.y]) < threshold_dist) && (s1 <= s2) # don't care if opponent is behind us
                    pos = VecSE2(trajectory[hrhc.h,1:3]) # x,y,θ of opponent at time step h
                    successor_states = getSuccessorStates(hrhc.motion_map[hrhc.v_cmd],hrhc.car_ID,hrhc.h,scene)
                    ΔX = successor_states[:,:,1] - pos.x # Δx, with opponent at origin
                    ΔY = successor_states[:,:,2] - pos.y # Δy with opponent at origin
                    Δθ = successor_states[:,:,3] - pos.θ # Δθ with opponent at origin
                    R = sqrt(ΔX.^2 + ΔY.^2)
                    ψ = atan2(ΔY,ΔX) - Δθ

                    cost = (1./(W*R.*cos(ψ).^2 + 2*L*R.*sin(ψ).^2 + 1)).*(-cos(ψ).^3 + 1.1) + 1

                end
            end
        end
    end
    return cost
end
function calculateObjective(car_ID,s₀,s,t,ϕ,T_MAX;ϕ_MAX=Float64(π),s_factor=1.0)
    s_norm = (s-s₀)/maximum(s)
    t_norm = t/T_MAX
    ϕ_norm = ϕ/ϕ_MAX

    #costs
    t_cost = abs(t_norm).^6
    ϕ_cost = abs(ϕ_norm).^6
    s_factor = 1
    s_cost = s_factor*(1-s_norm)
    A = [1 .5; #  [ϕ t] [a1 a2] [ϕ]
        .5 0] #         [a2 a3] [t]
    tϕ_cost = A[1,1]*(ϕ_norm).^2 + (A[1,2]+A[2,1])*(ϕ_norm).*(t_norm) + A[2,2]*(t_norm).^2

    objective = 1+s_cost+t_cost+ϕ_cost+tϕ_cost
    return objective
end
function AutomotiveDrivingModels.observe!(hrhc::HRHC, scene::Scene, roadway::Roadway, egoid::Int, tree::KDTree, obstacleMap, k_level)
    """
    Observe the current environment and select optimal action to apply at next
    time step
    """

    state = scene.vehicles[hrhc.car_ID].state
    hrhc.curve_ind = state.posF.roadind.ind.i
    v = state.v # current v
    hrhc.v = v
    trajectory = zeros(hrhc.h,3)
    action_selected = false
    δ_cmd = 1
    a_cmd = 1

    i = 0
    for i in 0:(hrhc.h-1)
        if action_selected
            break # out of for loop
        end
        # calculate successor states
        hrhc.successor_states = getSuccessorStates(hrhc.motion_map[hrhc.v_cmd], hrhc.car_ID, hrhc.h-i, scene)
        # project successor states onto track
        s,t,ϕ = loopProjectionKD(hrhc, scene, roadway, tree)
        # optimization objective
        objective = calculateObjective(hrhc.car_ID, scene.vehicles[hrhc.car_ID].state.posF.s,
        s, t, ϕ, hrhc.T_MAX)
        tailgateCost = tailgateAvoidance(hrhc, obstacleMap, tree, roadway, scene, k_level)
        objective = objective + tailgateCost
        collisionFlag = screenCollision(hrhc, obstacleMap, tree, roadway, scene, k_level)
        objective[collisionFlag .> 0] = Inf

        while (action_selected==false) && (minimum(objective) != Inf)
            index = indmin(objective) # find get a better method of optimizing this
            a_cmd, δ_cmd = ind2sub(s, index)
            hrhc.δ_cmd = δ_cmd
            # compute full trajectory up to horizon
            trajectory = computeTrajectory(hrhc.motion_map[hrhc.v_cmd], hrhc.car_ID, scene, a_cmd, δ_cmd, hrhc.h-i)
            # screen trajectory for collisions / validity
            out_of_bounds = screenTrajectory(trajectory, obstacleMap, scene, roadway, hrhc, tree, k_level)
            if out_of_bounds
                objective[index] = Inf
            else
                action_selected=true
                updateObstacleMap!(obstacleMap, k_level, hrhc.car_ID, trajectory)
                if k_level != hrhc.k # only assign commands for the given k_level
                    return # do not update the action
                end
            end
        end
    end

    hrhc.δ_cmd = δ_cmd
    # hrhc.v_cmd = max(hrhc.v_cmd + a_cmd - 2,1) # assumes 3 options for acceleration
    hrhc.v_cmd = hrhc.v_map[hrhc.v_cmd][a_cmd,δ_cmd]

    hrhc.δ = hrhc.δ_range[hrhc.δ_cmd] # next δ
    hrhc.v = hrhc.v_range[hrhc.v_cmd] # next v

    next_state = VehicleState(VecSE2(trajectory[1,1:3]),roadway,hrhc.v)
    hrhc.action = next_state # action
end
AutomotiveDrivingModels.get_name(::HRHC) = "HRHC"
# AutomotiveDrivingModels.action_context(driver::HRHC) = driver.action_context # AutomotiveDrivingModels.action_context
Base.rand(hrhc::HRHC) = hrhc.action

# Plotting functions
function plotSplineRoadway(x,y,θ,lane_width)
    perp_lines1 = zeros(2,length(x))
    perp_lines2 = zeros(2,length(x))

    perp_lines1[1,:] = x + (lane_width/2.0)*sin(θ)
    perp_lines1[2,:] = y - (lane_width/2.0)*cos(θ)
    perp_lines2[1,:] = x - (lane_width/2.0)*sin(θ)
    perp_lines2[2,:] = y + (lane_width/2.0)*cos(θ)

    # PyPlot.figure()
    # PyPlot.scatter(x,y)
    PyPlot.plot(x,y)
    PyPlot.plot(perp_lines1[1,:],perp_lines1[2,:],color="green")
    PyPlot.plot(perp_lines2[1,:],perp_lines2[2,:],color="green")
    PyPlot.axis("equal")
    # PyPlot.show()
end
function plotObjectiveHorizon(hrhc,scene,roadway,tree,trajectory,obstacleMap,xR,yR,θR)
    lo=hrhc.curve_ind
    if :V_MAX in fieldnames(hrhc)
        hi = hrhc.curve_ind + Int(1+div(hrhc.V_MAX*hrhc.Δt*hrhc.h,hrhc.Δs))
    else
        hi = hrhc.curve_ind + Int(1+div(hrhc.v_range[end]*hrhc.Δt*hrhc.h,hrhc.Δs))
    end
    lane_width = roadway.segments[1].lanes[1].width

    x = zeros(hrhc.h,size(hrhc.successor_states,1),size(hrhc.successor_states,2))
    y = zeros(size(x))
    Θ = zeros(size(x))
    s = zeros(size(x))
    t = zeros(size(x))
    ϕ = zeros(size(x))
    objective = zeros(size(x))

    for i in 1:hrhc.h
        getLegalMoves!(hrhc, scene, h=i)
        getSuccessorStates!(hrhc, scene)
        x[i,:,:] = copy(hrhc.successor_states[:,:,1])
        y[i,:,:] = copy(hrhc.successor_states[:,:,2])
        Θ[i,:,:] = copy(hrhc.successor_states[:,:,3])
        s[i,:,:], t[i,:,:], ϕ[i,:,:] = loopProjectionKD(hrhc,scene,roadway,tree)
        objective[i,:,:] = calculateObjective(hrhc,scene, roadway, tree,s[i,:,:],t[i,:,:],ϕ[i,:,:],obstacleMap,hrhc.k,hrhc.h)
    end

    PyPlot.figure(figsize=[12,4])

    PyPlot.subplot(141) # ϕ
    plotSplineRoadway(xR[lo:hi],yR[lo:hi],θR[lo:hi],lane_width)
    PyPlot.scatter(x,y,c=ϕ,edgecolor="none")
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.axis("off")
    PyPlot.title("|phi|")

    PyPlot.subplot(142) # s
    plotSplineRoadway(xR[lo:hi],yR[lo:hi],θR[lo:hi],lane_width)
    PyPlot.scatter(x,y,c=s,edgecolor="none")
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.axis("off")
    PyPlot.title("s")

    PyPlot.subplot(143) # t
    plotSplineRoadway(xR[lo:hi],yR[lo:hi],θR[lo:hi],lane_width)
    PyPlot.scatter(x,y,c=t,edgecolor="none")
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.axis("off")
    PyPlot.title("t")

    PyPlot.subplot(144) # objective
    plotSplineRoadway(xR[lo:hi],yR[lo:hi],θR[lo:hi],lane_width)
    PyPlot.scatter(x,y,c=log(objective),edgecolor="none")
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.axis("off")
    PyPlot.title("log objective")
end
function plot_stϕ(hrhc,roadway,scene,x,y,θ,trajectory,s,t,ϕ,objective)
    lo=hrhc.curve_ind
    hi=hrhc.curve_ind + Int(1+2*div(hrhc.V_MAX*hrhc.Δt*hrhc.h,hrhc.Δs))
    if hi > length(roadway.segments[1].lanes[1].curve)
        lo = length(roadway.segments[1].lanes[1].curve)
        hi=hrhc.curve_ind + Int(1+2*div(hrhc.V_MAX*hrhc.Δt*hrhc.h,hrhc.Δs))
    end
    lane_width = roadway.segments[1].lanes[1].width

    PyPlot.figure(figsize=[12,4])
    PyPlot.subplot(141)
    plotSplineRoadway(x[lo:hi],y[lo:hi],θ[lo:hi],lane_width)
    # PyPlot.scatter(Pts[1,:],Pts[2,:],color="red")
    PyPlot.scatter(hrhc.successor_states[:,:,1],hrhc.successor_states[:,:,2],c=abs(ϕ),edgecolor="none")
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.scatter(scene.vehicles[hrhc.car_ID].state.posG.x, scene.vehicles[hrhc.car_ID].state.posG.y, c="k", edgecolors="none",s=40)
    PyPlot.axis("off")
    PyPlot.title("|phi|")

    PyPlot.subplot(142)
    plotSplineRoadway(x[lo:hi],y[lo:hi],θ[lo:hi],lane_width)
    # PyPlot.scatter(Pts[1,:],Pts[2,:],color="red")
    PyPlot.scatter(hrhc.successor_states[:,:,1],hrhc.successor_states[:,:,2],c=abs(t),edgecolor="none")
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.scatter(scene.vehicles[hrhc.car_ID].state.posG.x, scene.vehicles[hrhc.car_ID].state.posG.y, c="k", edgecolors="none",s=40)
    PyPlot.axis("off")
    PyPlot.title("|t|")

    PyPlot.subplot(143)
    plotSplineRoadway(x[lo:hi],y[lo:hi],θ[lo:hi],lane_width)
    # PyPlot.scatter(Pts[1,:],Pts[2,:],color="red")
    PyPlot.scatter(hrhc.successor_states[:,:,1],hrhc.successor_states[:,:,2],c=s,edgecolor="none")
    PyPlot.scatter(scene.vehicles[hrhc.car_ID].state.posG.x, scene.vehicles[hrhc.car_ID].state.posG.y, c="k", edgecolors="none",s=40)
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.axis("off")
    PyPlot.title("s")


    PyPlot.subplot(144)
    plotSplineRoadway(x[lo:hi],y[lo:hi],θ[lo:hi],lane_width)
    # PyPlot.scatter(Pts[1,:],Pts[2,:],color="red")
    PyPlot.scatter(hrhc.successor_states[:,:,1],hrhc.successor_states[:,:,2],c=log(objective),edgecolor="none")
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.scatter(scene.vehicles[hrhc.car_ID].state.posG.x, scene.vehicles[hrhc.car_ID].state.posG.y, c="k", edgecolors="none",s=40)
    PyPlot.axis("off")
    PyPlot.title("objective")
end
function plotHRHCInfo(hrhc,models,scene,roadway,trajectory,cmd,x,y,Θ,s,t,ϕ,objective)
    lo = hrhc.curve_ind
    hi = hrhc.curve_ind + Int(1+2*div(hrhc.V_MAX*hrhc.Δt*hrhc.h,hrhc.Δs))
    lane_width = roadway.segments[1].lanes[1].width
    if hi > length(roadway.segments[1].lanes[1].curve)
        lo = length(roadway.segments[1].lanes[1].curve)
        hi=hrhc.curve_ind + Int(1+2*div(hrhc.V_MAX*hrhc.Δt*hrhc.h,hrhc.Δs))
    end
    PyPlot.figure(figsize=[12,10])
    # Plot Raceway
    PyPlot.subplot(221)
    # plotSplineRoadway(x[lo:hi],y[lo:hi],θ[lo:hi],lane_width)
    plotSplineRoadway(x,y,Θ,lane_width)
    PyPlot.scatter(hrhc.successor_states[:,:,1],hrhc.successor_states[:,:,2],color="red")
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.scatter(roadway.segments[1].lanes[1].curve[hrhc.curve_ind].pos.x, roadway.segments[1].lanes[1].curve[hrhc.curve_ind].pos.y, c="white", s=40)
    for (id,car) in models
        if id == hrhc.car_ID
            PyPlot.scatter(scene.vehicles[id].state.posG.x,scene.vehicles[id].state.posG.y,c="red",s=20)
        else
            PyPlot.scatter(scene.vehicles[id].state.posG.x,scene.vehicles[id].state.posG.y,c="blue",s=20)
        end
    end
    PyPlot.axis("off")
    PyPlot.title("Raceway with Motion Primitives")

    PyPlot.subplot(222)
    plotSplineRoadway(x[lo:hi],y[lo:hi],Θ[lo:hi],lane_width)
    PyPlot.scatter(scene.vehicles[hrhc.car_ID].state.posG.x, scene.vehicles[hrhc.car_ID].state.posG.y, c="red", edgecolors="none",s=40)
    PyPlot.scatter(hrhc.successor_states[:,:,1], hrhc.successor_states[:,:,2],c=log(objective),edgecolors="none")
    PyPlot.scatter(hrhc.successor_states[cmd[1],cmd[2],1], hrhc.successor_states[cmd[1],cmd[2],2],c="white",s=40)
    PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
    PyPlot.axis("off")
    PyPlot.title("Log Objective Function")
end
