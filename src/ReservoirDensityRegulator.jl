"Defaults for the optional charge-neutral upstream density regulator."
default_reservoir_density_regulator_config() = (
    enabled=false,
    target_density_e=1.0,
    target_density_i=1.0,
    relaxation_time_omega_pe=200.0,
    maximum_fractional_correction_per_step=1.0e-3,
    injection_side=:upstream,
    inject_charge_neutral_pairs=true,
)

function reservoir_density_regulator_config(raw=nothing)
    defaults=default_reservoir_density_regulator_config()
    isnothing(raw) && return defaults
    return merge(defaults,raw)
end

"Checkpointable remainder prevents bias when the requested correction is sub-macroparticle."
mutable struct ReservoirDensityRegulatorState
    pair_remainder::Float64
    cumulative_pairs::Int
end
ReservoirDensityRegulatorState()=ReservoirDensityRegulatorState(0.0,0)

function regulator_injection_side(config,flow_speed)
    config.injection_side === :upstream && return flow_speed >= 0 ? :left : :right
    config.injection_side in (:left,:right) ||
        error("reservoir_density_regulator.injection_side must be :upstream, :left, or :right")
    return config.injection_side
end

function inject_regulated_particle!(species,side,cfg,dt,rng)
    thermal_speed=sqrt(species.name === :electrons ?
        cfg.electron_temperature : cfg.ion_temperature/cfg.mass_ratio)
    drift=cfg.flow_sound_speed/sqrt(cfg.mass_ratio)
    velocity=sample_flux_velocity(rng,drift,thermal_speed,side)
    position=side === :left ? rand(rng)*abs(velocity)*dt :
        cfg.domain_length-rand(rng)*abs(velocity)*dt
    push!(species.x,clamp(position,0.0,cfg.domain_length))
    push!(species.v,velocity)
    return species.charge*species.weight
end

"Inject the common electron/ion density deficit as neutral upstream pairs."
function apply_reservoir_density_regulator!(electrons,ions,cfg,config,state,dt,rng_e,rng_i)
    config.enabled || return (0.0,0,0.0)
    isapprox(electrons.charge*electrons.weight+
             ions.charge*ions.weight,0.0;rtol=1e-12,atol=1e-14) ||
        error("charge-neutral density regulation requires equal electron and ion macroparticle charge magnitudes")
    target_e=config.target_density_e*cfg.domain_length/electrons.weight
    target_i=config.target_density_i*cfg.domain_length/ions.weight
    pair_deficit=max(0.0,min(target_e-length(electrons.x),
                             target_i-length(ions.x)))
    requested=min(dt/config.relaxation_time_omega_pe*pair_deficit,
        config.maximum_fractional_correction_per_step*min(target_e,target_i))
    state.pair_remainder += requested
    number=max(0,floor(Int,state.pair_remainder))
    state.pair_remainder -= number
    side=regulator_injection_side(config,cfg.flow_sound_speed)
    injected_charge=0.0
    for _ in 1:number
        injected_charge += inject_regulated_particle!(electrons,side,cfg,dt,rng_e)
        injected_charge += inject_regulated_particle!(ions,side,cfg,dt,rng_i)
    end
    state.cumulative_pairs += number
    return injected_charge,number,pair_deficit
end

reservoir_density_regulator_history_path(history_file)=
    splitext(history_file)[1]*"_density_regulator.tsv"

function append_reservoir_density_regulator_history!(path,step,time,electrons,ions,
        cfg,state,injected_pairs,pair_deficit;new_file=false)
    mkpath(dirname(path))
    mode=new_file ? "w" : "a"
    open(path,mode) do io
        new_file && println(io,"step\ttime\tmean_ne\tmean_ni\tpair_deficit\tinjected_pairs\tcumulative_pairs\tpair_remainder")
        mean_ne=length(electrons.x)*electrons.weight/cfg.domain_length
        mean_ni=length(ions.x)*ions.weight/cfg.domain_length
        println(io,"$step\t$time\t$mean_ne\t$mean_ni\t$pair_deficit\t$injected_pairs\t$(state.cumulative_pairs)\t$(state.pair_remainder)")
    end
    return nothing
end
