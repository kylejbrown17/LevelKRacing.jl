module LevelKRacing

using AutomotiveDrivingModels
using AutoViz
using NearestNeighbors

import AutomotiveDrivingModels: get_actions!, observe!, get_name
import Base.rand
import PyPlot

export
    # HierarchicalRecedingHorizonController
    HRHC,
    curveDist,
    wrap_to_π,
    kdProject,
    generateObstacleMap,
    updateObstacleMap!,
    generateMotionMap,
    screenCollision,
    tailgateAvoidance,
    getSuccessorStates,
    loopProjectionKD,
    computeTrajectory,
    screenTrajectory,
    checkCollision,
    calculateObjective,
    plot_stϕ,
    plotHRHCInfo,
    plotObjectiveHorizon,
    plotSplineRoadway,

    # LevelKVisualizations
    DemoState,
    DemoMotionPrimitives,
    DemoIncreasedHorizonBehavior,
    DemoBuildingBlocks,
    DemoObjective,
    DemoTailgateAvoidance,
    DemoObserveActObjective,
    DemoObstacleMap,

    # SplineRaceWay
    Raceway,

    # SplineUtils
    ClosedB_Spline,
    B_SplineDerivative,
    ResampleSplineEven,
    GenSplineRoadway,
    PlotSplineRoadway

include(Pkg.dir("LevelKRacing","src","controllers","HierarchicalRecedingHorizonController.jl"))
include(Pkg.dir("LevelKRacing","src","visualization","LevelKVisualizations.jl"))
include(Pkg.dir("LevelKRacing","src","utils","SplineRaceWay.jl"))
include(Pkg.dir("LevelKRacing","src","utils","SplineUtils.jl"))

end
