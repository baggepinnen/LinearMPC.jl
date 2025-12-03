# lb <= Au uk + Ax xk <= ub for k ∈ ks
# (additional terms Ar rₖ, Aw wₖ, Ad dₖ, Aup u⁻ₖ)

struct Constraint{KST <: AbstractVector{Int64}}
    Au::FixedSizeMatrixDefault{Float64}
    Ax::FixedSizeMatrixDefault{Float64}
    Ar::FixedSizeMatrixDefault{Float64}
    Aw::FixedSizeMatrixDefault{Float64}
    Ad::FixedSizeMatrixDefault{Float64}
    Aup::FixedSizeMatrixDefault{Float64}
    ub::Vector{Float64}
    lb::Vector{Float64}
    ks::KST
    soft::Bool
    binary::Bool
    prio::Int
end

# Weights used to define the objective function of the OCP
struct MPCWeights
    Q::FixedSizeMatrixDefault{Float64}
    R::FixedSizeMatrixDefault{Float64}
    Rr::FixedSizeMatrixDefault{Float64}
    S::FixedSizeMatrixDefault{Float64}
    Qf::FixedSizeMatrixDefault{Float64}
    Qfx::FixedSizeMatrixDefault{Float64}
end

function MPCWeights(nu,nx,nr)
    return MPCWeights(fs(1.0*I(nr)),fs(1.0*I(nu)),fszeros(nu,nu),
                      fszeros(nx,nu),fszeros(nr,nr),fszeros(nx,nx))
end

"""
MPC controller settings.

# Fields
- `reference_condensation::Bool = false`: Collapse reference trajectory to setpoint 
- `reference_tracking::Bool = true`: Enable reference tracking
- `reference_preview::Bool = false`: Enable time-varying reference preview
- `soft_weight::Float64 = 1e6`: Penalty weight for soft constraint violations
- `solver_opts::Dict{Symbol,Any}`: Additional solver options
"""
Base.@kwdef mutable struct MPCSettings
    move_block_foh::Bool= false
    reference_condensation::Bool= false
    reference_tracking::Bool= true
    reference_preview::Bool = false
    soft_weight::Float64= 1e6
    solver_opts::Dict{Symbol,Any} = Dict()
    traj2setpoint::FixedSizeMatrixDefault{Float64} = fszeros(0,0)
end

# MPC controller
mutable struct MPC 

    model::Model

    # parameters 
    nr::Int
    nuprev::Int

    # Horizons 
    Np::Int # Prediction
    Nc::Int # Control

    ## 
    weights::MPCWeights

    # lb <= u <=ub
    umin::FixedSizeVectorDefault{Float64}
    umax::FixedSizeVectorDefault{Float64}
    binary_controls::FixedSizeVectorDefault{Int64}

    # General constraints 
    constraints::Vector{<:Constraint}

    # Settings
    settings::MPCSettings

    #Optimization problem
    mpQP

    # DAQP optimization model
    opt_model::DAQPBase.Model

    # Prestabilizing feedback
    K::FixedSizeMatrixDefault{Float64}

    # Move blocks
    move_blocks::FixedSizeVectorDefault{Int}

    mpqp_issetup::Bool

    uprev::FixedSizeVectorDefault{Float64}

    traj2setpoint::FixedSizeMatrixDefault{Float64}

    state_observer

    Δx0::FixedSizeVectorDefault{Float64}
end

function MPC(model::Model;Np=10,Nc=Np)
    MPC(model,0,0,Np,Nc,
        MPCWeights(model.nu,model.nx,model.ny),
        fszeros(0),fszeros(0),FixedSizeVector(Int[]),
        Constraint[],MPCSettings(),nothing,
        DAQP.Model(),fszeros(model.nu,model.nx),FixedSizeVector(Int[]),false, fszeros(model.nu),fszeros(0,0),
       nothing,fszeros(model.nx))
end

function MPC(F,G;Gd=fszeros(0,0), C=fszeros(0,0), Dd= fszeros(0,0), offset=fszeros(0), Ts= -1.0, Np=10, Nc = Np)
    MPC(Model(fs(F),fs(G);Gd=fs(Gd),offset=fs(offset),C=fs(C),Dd=fs(Dd),Ts=Ts);Np,Nc);
end

function MPC(A,B,Ts::Float64; Bd = fszeros(0,0), offset=fszeros(0), C = fszeros(0,0), Dd = fszeros(0,0), Np=10, Nc=Np)
    MPC(Model(fs(A),fs(B),Ts; Bd=fs(Bd), offset=fs(offset),C=fs(C),Dd=fs(Dd));Np,Nc)
end

function MPC(sys; Ts=1.0, Np=10, Nc=Np)
    MPC(Model(sys;Ts);Np,Nc)
end

struct ParameterRange
    xmin::Vector{Float64}
    xmax::Vector{Float64}

    rmin::Vector{Float64}
    rmax::Vector{Float64}

    dmin::Vector{Float64}
    dmax::Vector{Float64}

    umin::Vector{Float64}
    umax::Vector{Float64}
end


function ParameterRange(mpc::MPC)

    nx,nr,nd,nuprev = get_parameter_dims(mpc);

    xmin,xmax = -100*ones(nx),100*ones(nx)
    rmin,rmax = -100*ones(nr),100*ones(nr)
    dmin,dmax = -100*ones(nd),100*ones(nd)
    if(nuprev > 0)
        nmin,nmax = length(mpc.umin),length(mpc.umax)
        nb = max(nmin,nmax)
        umin = [mpc.umin;-100*ones(nb-nmin)]
        umax = [mpc.umax;+100*ones(nb-nmax)]
    else
        umin,umax = zeros(0),zeros(0)
    end

    return ParameterRange(xmin,xmax,
                          rmin,rmax,
                          dmin,dmax,
                          umin,umax)
end
