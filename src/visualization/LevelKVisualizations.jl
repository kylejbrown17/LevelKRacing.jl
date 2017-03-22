using AutomotiveDrivingModels
using AutoViz
using NearestNeighbors

import PyPlot

# export
#     DemoState,
#     DemoMotionPrimitives,
#     DemoIncreasedHorizonBehavior,
#     DemoBuildingBlocks,
#     DemoObjective,
#     DemoTailgateAvoidance,
#     DemoObserveActObjective,
#     DemoObstacleMap


function DemoState(track,scene,carcolors;targetid=1,zoom=5.0)
    x = scene.vehicles[1].state.posG.x
    y = scene.vehicles[1].state.posG.y
    θ = scene.vehicles[1].state.posG.θ
    v = scene.vehicles[1].state.v

    s,t,ϕ = kdProject(x,y,θ,track.tree,track.roadway,track.models[1])
    println("Vehicle 1: ")
    @show x
    @show y
    @show θ
    @show v
    @show s
    @show t
    @show ϕ

    render(scene, track.roadway, cam=CarFollowCamera(targetid, zoom), car_colors=carcolors, canvas_height=300)
end
function DemoMotionPrimitives(track,context)
    demoCar = HRHC(1,track.roadway,context,h=50,v_max=120.0,μ=30.0,a_step=12.0,a_range=[-1,0,1],k=2)

    v_idx = [50,150,220]
    δ_idx = [18,14,40]
    a_idx = [3,1,2]
    titles = ["V = 50", "V = 150", "V = 220"]
    h=12
    motion_map = demoCar.motion_map

    PyPlot.figure(figsize=[12,4])
    for i in 1:3
        PyPlot.subplot(1,3,i)
        PyPlot.axis("tight")
        PyPlot.axis("off")
        ΔXYθ = motion_map[v_idx[i]][:,:,h,:]
        trajectory = motion_map[v_idx[i]][a_idx[i],:,1:h,:]
        for j in 1:size(trajectory,1)
            PyPlot.plot(trajectory[j,:,2],trajectory[j,:,1],c="red")
        end
        PyPlot.plot(trajectory[δ_idx[i],:,2],trajectory[δ_idx[i],:,1],c="lime",linewidth=2)
        PyPlot.scatter(ΔXYθ[:,:,2],ΔXYθ[:,:,1],edgecolor="none")

        PyPlot.axis("equal")
        PyPlot.ylim([0,80])
        PyPlot.xlim([-40,40])

        PyPlot.title(titles[i])
    end
end
function DemoIncreasedHorizonBehavior(track,context)
    demoCar = HRHC(1,track.roadway,context,h=50,v_max=120.0,μ=30.0,a_step=12.0,a_range=[-1,0,1],k=2)
    v_idx = 220
    δ_idx = 40
    a_idx = 3
    h=demoCar.h
    motion_map = demoCar.motion_map
    PyPlot.figure()
    ΔXYθ = motion_map[v_idx][:,:,h,:]
    trajectory = motion_map[v_idx][a_idx,:,1:h,:]
    for j in 1:size(trajectory,1)
        PyPlot.plot(trajectory[j,:,2],trajectory[j,:,1],c="red")
    end
    PyPlot.plot(trajectory[δ_idx,:,2],trajectory[δ_idx,:,1],c="lime",linewidth=2)
    PyPlot.scatter(ΔXYθ[:,:,2],ΔXYθ[:,:,1],edgecolor="none")

    PyPlot.axis("tight")
    PyPlot.axis("off")
    PyPlot.axis("equal")
    PyPlot.title("Increased Horizon Behavior (V = 220)")
end
function DemoBuildingBlocks()
    t1 = linspace(-13,13,200)
    T_MAX = 10
    t_cost = exp(abs(t1/T_MAX).^6)

    ϕ1 = linspace(-Float64(π),Float64(π),200)
    ϕ_MAX = Float64(π)/2
    ϕ_cost = exp(abs(ϕ1/ϕ_MAX).^2)

    s1 = linspace(0,1,20)
    s_factor = 5
    s_cost = s_factor*(1-s1/maximum(s1))

    PyPlot.figure(figsize=[12,4])
    PyPlot.subplot(1,4,1)
    PyPlot.plot(t1/T_MAX,t_cost)
    PyPlot.ylim([0,50])
    PyPlot.title("cost(t)")
    PyPlot.subplot(1,4,2)
    PyPlot.plot(ϕ1/ϕ_MAX,ϕ_cost)
    PyPlot.ylim([0,50])
    PyPlot.title("cost(phi)")
    PyPlot.subplot(1,4,3)
    PyPlot.plot(s1,s_cost)
    PyPlot.ylim([0,50])
    PyPlot.title("cost(s)")
end
function DemoObjective(t_shift)
    θ = linspace(-Float64(π),Float64(π),360)
    r = linspace(1.0,20.0,20)
    x = [cos(i) for i in θ]*[i for i in r]'
    y = [sin(i) for i in θ]*[i for i in r]'
    R = [1 for i in θ]*[i for i in r]'
    Θ = [i for i in θ]*[1 for i in r]'
    s = x
    t = y - t_shift
    ϕ = Θ

    car_ID = 1
    s₀ = -21 # scene.vehicles[car_ID].state.posF.s

    T_MAX= 10
    ϕ_MAX = Float64(π)

    s_norm = s/maximum(s)
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

    objective = calculateObjective(car_ID,s,s₀,t,ϕ,T_MAX)

    PyPlot.figure(figsize=[12,3])

    PyPlot.subplot(1,4,1)
    PyPlot.plot([-20,20],[T_MAX+t_shift, T_MAX+t_shift],c="black",linewidth=2)
    PyPlot.plot([-20,20],[-T_MAX+t_shift, -T_MAX+t_shift],c="black",linewidth=2)
    PyPlot.scatter(x,y,c=s_cost,edgecolor="none",s=4)
    PyPlot.axis("off")
    PyPlot.axis("equal")
    PyPlot.title("cost(s)")
    PyPlot.subplot(1,4,2)
    PyPlot.plot([-20,20],[T_MAX+t_shift, T_MAX+t_shift],c="black",linewidth=2)
    PyPlot.plot([-20,20],[-T_MAX+t_shift,-T_MAX+t_shift],c="black",linewidth=2)
    PyPlot.scatter(x,y,c=t_cost,edgecolor="none",s=4)
    PyPlot.axis("off")
    PyPlot.axis("equal")
    PyPlot.title("cost(t)")
    PyPlot.subplot(1,4,3)
    PyPlot.plot([-20,20],[T_MAX+t_shift, T_MAX+t_shift],c="black",linewidth=2)
    PyPlot.plot([-20,20],[-T_MAX+t_shift,-T_MAX+t_shift],c="black",linewidth=2)
    PyPlot.scatter(x,y,c=tϕ_cost,edgecolor="none",s=4)
    PyPlot.axis("off")
    PyPlot.axis("equal")
    PyPlot.title("cost(t,phi)")
    PyPlot.subplot(1,4,4)
    PyPlot.plot([-20,20],[T_MAX+t_shift, T_MAX+t_shift],c="black",linewidth=2)
    PyPlot.plot([-20,20],[-T_MAX+t_shift,-T_MAX+t_shift],c="black",linewidth=2)
    PyPlot.scatter(x,y,c=log(objective),edgecolor="none",s=4)
    PyPlot.scatter(x[indmin(objective)],y[indmin(objective)],c="white",s=30)
    PyPlot.axis("off")
    PyPlot.axis("equal")
    PyPlot.title("objective = cost(s,t,phi)")
end
function DemoTailgateAvoidance()
    θ = linspace(-Float64(π),Float64(π),360)
    r = linspace(1.0,20.0,20)
    x = [cos(i) for i in θ]*[i for i in r]'
    y = [sin(i) for i in θ]*[i for i in r]'
    R = [1 for i in θ]*[i for i in r]'
    ψ = [i for i in θ]*[1 for i in r]'

    #costs
    R_cost = 1./(2.5*R.*cos(ψ).^2 + 2*4.8*R.*sin(ψ).^2 + 1)
    # ψ_cost = 1./(R.*(sin(ψ).^2 + 1))
    ψ_cost = -cos(ψ).^3 + 1.1
    Rψ_cost = R_cost.*ψ_cost + 1

    PyPlot.figure(figsize=[12,4])

    PyPlot.subplot(1,3,1)
    PyPlot.scatter(x,y,c=R_cost,edgecolor="none",s=8)
    PyPlot.axis("equal")
    PyPlot.axis("off")
    PyPlot.title("cost(r)")
    PyPlot.subplot(1,3,2)
    PyPlot.scatter(x,y,c=ψ_cost,edgecolor="none",s=8)
    PyPlot.axis("equal")
    PyPlot.axis("off")
    PyPlot.title("cost(psi)")
    PyPlot.subplot(1,3,3)
    PyPlot.scatter(x,y,c=Rψ_cost,edgecolor="none",s=8)
    PyPlot.axis("equal")
    PyPlot.axis("off")
    PyPlot.title("cost(r,psi)")

    return x,y,Rψ_cost
end
function DemoObserveActObjective(track,scene;steps=50)
    lane_width = track.roadway.segments[1].lanes[1].width
    actions = Array(DriveAction, length(scene))
    for i in 1:steps
        k_level = 0 # needs to be updated into a loop
        for k_level in 0:maximum([model.k for (id,model) in track.models])
            for (i, veh) in enumerate(scene)
                model = track.models[veh.def.id]
                observe!(model, scene, track.roadway, veh.def.id, track.tree, track.obstacleMap, k_level)
                actions[i] = rand(model)
            end
        end
        for (veh, action) in zip(scene, actions)
            model = track.models[veh.def.id]
            context = action_context(model)
            veh.state = propagate(veh, action, context, track.roadway)
        end
    end

    PyPlot.figure(figsize=[12,8])

    lo = Inf
    hi = 0

    for (id,hrhc) in track.models
        lo = Int(min(lo,hrhc.curve_ind))
        hi = Int(max(hi,hrhc.curve_ind + Int(1+div(hrhc.v_range[end]*hrhc.Δt*hrhc.h,hrhc.Δs))))
        lane_width = track.roadway.segments[1].lanes[1].width

        x = zeros(hrhc.h,size(hrhc.successor_states,1),size(hrhc.successor_states,2))
        y = zeros(size(x))
        Θ = zeros(size(x))
        s = zeros(size(x))
        t = zeros(size(x))
        ϕ = zeros(size(x))
        objective = zeros(size(x))
        avoidanceCost = zeros(size(x))
        collisionFlag = zeros(size(x))

        s₀ = scene.vehicles[hrhc.car_ID].state.posF.s
        for i in 1:hrhc.h
            successor_states = getSuccessorStates(hrhc.motion_map[hrhc.v_cmd],hrhc.car_ID, i, scene)
            x[i,:,:] = copy(successor_states[:,:,1])
            y[i,:,:] = copy(successor_states[:,:,2])
            Θ[i,:,:] = copy(successor_states[:,:,3])
            s[i,:,:], t[i,:,:], ϕ[i,:,:] = loopProjectionKD(hrhc,scene,track.roadway,track.tree)
            objective[i,:,:] = calculateObjective(hrhc.car_ID,s[i,:,:],s₀,t[i,:,:],ϕ[i,:,:],hrhc.T_MAX)
            avoidanceCost[i,:,:] = tailgateAvoidance(hrhc, track.obstacleMap, track.tree, track.roadway, scene, hrhc.k)
            objective[i,:,:] = objective[i,:,:] + avoidanceCost[i,:,:]
            collisionFlag[i,:,:] = screenCollision(hrhc, track.obstacleMap, track.tree, track.roadway, scene, hrhc.k)
        end
        objective[collisionFlag .> 0] = Inf
        trajectory = track.obstacleMap[hrhc.k][hrhc.car_ID]

        PyPlot.scatter(x,y,c=log(objective),edgecolor="none")
        PyPlot.plot(trajectory[:,1],trajectory[:,2],color="red")
        PyPlot.axis("off")
        PyPlot.title("log objective functions of vehicles on track")

    end
    plotSplineRoadway(track.x[lo:hi],track.y[lo:hi],track.θ[lo:hi],lane_width)
end
function DemoObstacleMap(track,k_level)
    lo = Inf
    hi = 0
    lane_width = track.roadway.segments[1].lanes[1].width
    for (id,hrhc) in track.models
        lo = Int(min(lo,hrhc.curve_ind))
        hi = Int(max(hi,hrhc.curve_ind + Int(1+div(hrhc.v_range[end]*hrhc.Δt*hrhc.h,hrhc.Δs))))
        # lane_width = track.roadway.segments[1].lanes[1].width
    end
    plotSplineRoadway(track.x[lo:hi],track.y[lo:hi],track.θ[lo:hi],lane_width)

    for (car_id,trajectory) in track.obstacleMap[k_level]
        PyPlot.plot(trajectory[:,1],trajectory[:,2],c="r",linewidth=1)
        PyPlot.scatter(trajectory[end,1],trajectory[end,2],c="r",edgecolor="none",s=10)
    end
    PyPlot.axis("equal")
end
