"""Independent deterministic illumination multiplier for debris photoemission.

This represents a modulated UV source or periodic illumination of rotating
debris.  The disabled path is exactly unity and consumes no random numbers.
"""
default_periodic_photoemission_driver_config() = (
    enabled=false,
    waveform=:sinusoidal,
    modulation_depth=0.0,
    frequency_omega_pe=0.024,
    phase=0.0,
    ramp_start_time=100.0,
    ramp_duration=100.0,
)

function periodic_photoemission_driver_config(cfg)
    defaults=default_periodic_photoemission_driver_config()
    isnothing(cfg) ? defaults : merge(defaults,cfg)
end

"C1 smoothstep ramp, zero before start and unity after the ramp interval."
function periodic_photoemission_ramp(time,cfg)
    time<=cfg.ramp_start_time && return 0.0
    cfg.ramp_duration==0 && return 1.0
    u=clamp((time-cfg.ramp_start_time)/cfg.ramp_duration,0.0,1.0)
    u^2*(3-2u)
end

function periodic_photoemission_sample(cfg,time)
    cfg.enabled || return (multiplier=1.0,envelope=0.0,phase=0.0)
    cfg.waveform===:sinusoidal ||
        error("Unsupported periodic photoemission waveform=$(cfg.waveform)")
    envelope=periodic_photoemission_ramp(time,cfg)
    phase=2pi*cfg.frequency_omega_pe*time+cfg.phase
    multiplier=1+envelope*cfg.modulation_depth*cos(phase)
    return (;multiplier,envelope,phase)
end

periodic_photoemission_history_path(history_file)=
    joinpath(dirname(history_file),"periodic_photoemission_driver.tsv")

function append_periodic_photoemission_history!(path,step,time,sample;
        new_file=false)
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        new_file && println(io,
            "step\ttime_omega_pe\tmultiplier\tramp_envelope\tphase_radians")
        println(io,join((step,time,sample.multiplier,sample.envelope,
                         sample.phase),'\t'))
    end
    nothing
end
