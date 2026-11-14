"""Optional external debris-charge forcing.

`:prescribed_charge_offset` retains the historical Poisson-only offset.
`:sinusoidal_current` injects charge into the dynamic debris state so that the
resulting potential feeds back on particle collection. Disabled is an exact
zero and consumes no random numbers.
"""
default_external_debris_charge_driver_config() = (
    enabled=false,
    model=:prescribed_charge_offset,
    amplitude_e=0.0,
    current_amplitude_e_per_omega_pe=0.0,
    maximum_absolute_driven_charge_e=Inf,
    frequency_omega_pe=0.019,
    phase=0.0,
    ramp_start_time=100.0,
    ramp_duration=100.0,
)

mutable struct ExternalDebrisCurrentState
    cumulative_driven_charge_e::Float64
end
ExternalDebrisCurrentState()=ExternalDebrisCurrentState(0.0)

function external_debris_charge_driver_config(cfg)
    defaults=default_external_debris_charge_driver_config()
    isnothing(cfg) ? defaults : merge(defaults,cfg)
end

function external_charge_ramp(time,cfg)
    time<=cfg.ramp_start_time && return (value=0.0,derivative=0.0)
    cfg.ramp_duration==0 && return (value=1.0,derivative=0.0)
    u=(time-cfg.ramp_start_time)/cfg.ramp_duration
    u>=1 && return (value=1.0,derivative=0.0)
    u=clamp(u,0.0,1.0)
    return (value=u^2*(3-2u),
            derivative=6u*(1-u)/cfg.ramp_duration)
end

function external_debris_charge_sample(cfg,time,normalized_charge_per_e)
    cfg.enabled || return (charge_e=0.0,charge_normalized=0.0,
        current_e_per_omega_pe=0.0,envelope=0.0,phase=0.0,
        cumulative_driven_charge_e=0.0)
    ramp=external_charge_ramp(time,cfg)
    omega=2pi*cfg.frequency_omega_pe
    phase=omega*time+cfg.phase
    if cfg.model===:prescribed_charge_offset
        charge_e=cfg.amplitude_e*ramp.value*cos(phase)
        current_e=cfg.amplitude_e*(ramp.derivative*cos(phase)-
                                   ramp.value*omega*sin(phase))
    elseif cfg.model===:sinusoidal_current
        charge_e=0.0
        current_e=cfg.current_amplitude_e_per_omega_pe*ramp.value*sin(phase)
    else
        error("Unsupported external debris charge model=$(cfg.model)")
    end
    return (;charge_e,charge_normalized=charge_e*normalized_charge_per_e,
        current_e_per_omega_pe=current_e,envelope=ramp.value,phase,
        cumulative_driven_charge_e=0.0)
end

"Advance the external-current contribution with trapezoidal time integration."
function advance_external_debris_current!(debris,state,cfg,time,dt,charge_per_e)
    sample=external_debris_charge_sample(cfg,time,charge_per_e)
    (!cfg.enabled || cfg.model!==:sinusoidal_current) && return sample,0.0
    previous=external_debris_charge_sample(cfg,time-dt,charge_per_e)
    proposed=0.5dt*(previous.current_e_per_omega_pe+
                    sample.current_e_per_omega_pe)
    lower=-cfg.maximum_absolute_driven_charge_e-state.cumulative_driven_charge_e
    upper= cfg.maximum_absolute_driven_charge_e-state.cumulative_driven_charge_e
    delta_e=clamp(proposed,lower,upper)
    state.cumulative_driven_charge_e += delta_e
    debris.total_charge_e += delta_e
    delta_normalized=delta_e*charge_per_e
    debris.total_charge_normalized += delta_normalized
    enriched=merge(sample,
        (cumulative_driven_charge_e=state.cumulative_driven_charge_e,))
    return enriched,delta_normalized
end

external_debris_charge_history_path(history_file)=
    joinpath(dirname(history_file),"external_debris_charge_driver.tsv")

function append_external_debris_charge_history!(path,step,time,sample,
        self_charge_e,debris_potential;new_file=false)
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        new_file && println(io,
            "step\ttime_omega_pe\tself_charge_e\texternal_charge_e\tdeposited_charge_e\texternal_current_e_per_omega_pe\tramp_envelope\tphase_radians\tdebris_potential\texternal_power_normalized\tcumulative_current_driven_charge_e")
        deposited=self_charge_e+sample.charge_e
        power=debris_potential*sample.current_e_per_omega_pe
        println(io,join((step,time,self_charge_e,sample.charge_e,deposited,
            sample.current_e_per_omega_pe,sample.envelope,sample.phase,
            debris_potential,power,sample.cumulative_driven_charge_e),'\t'))
    end
    nothing
end
