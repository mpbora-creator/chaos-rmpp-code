"Deterministic, opt-in small-signal probe for SEE-to-plasma response."
default_see_response_probe_config() = (
    enabled=false,
    waveform=:log_chirp,
    amplitude=0.02,
    start_frequency=0.002,
    end_frequency=0.20,
    start_time=100.0,
    end_time=1000.0,
)

function see_response_probe_config(cfg)
    defaults=default_see_response_probe_config()
    isnothing(cfg) ? defaults : merge(defaults,cfg)
end

function see_response_probe_sample(cfg,time)
    (!cfg.enabled || time<cfg.start_time || time>cfg.end_time) &&
        return (multiplier=1.0,frequency=0.0,active=false)
    duration=cfg.end_time-cfg.start_time
    u=clamp((time-cfg.start_time)/duration,0.0,1.0)
    ratio=cfg.end_frequency/cfg.start_frequency
    frequency=cfg.start_frequency*ratio^u
    phase=2pi*cfg.start_frequency*duration/log(ratio)*(ratio^u-1)
    return (multiplier=1+cfg.amplitude*sin(phase),frequency=frequency,
            active=true)
end

see_response_probe_history_path(history_file)=
    joinpath(dirname(history_file),"see_response_probe.tsv")

function append_see_response_probe_history!(path,step,time,sample;new_file=false)
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        new_file && println(io,
            "step\ttime_omega_pe\tyield_multiplier\t" *
            "instantaneous_frequency\tactive")
        println(io,join((step,time,sample.multiplier,sample.frequency,
                         sample.active ? 1 : 0),'\t'))
    end
    return nothing
end
