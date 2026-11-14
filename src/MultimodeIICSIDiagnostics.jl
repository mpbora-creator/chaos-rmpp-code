"""Read-only, pluggable diagnostics for simultaneous resolved ion-density modes.

This module never changes fields, particles, charging, emission, or RNG state.
With `enabled=false` it performs no sampling and creates no file.
"""

default_multimode_iicsi_diagnostics_config() = (
    enabled=false,
    modes=(1,2),
    window_fraction=(0.0,1.0),
    every=5,
)

multimode_iicsi_diagnostics_config(raw=nothing) = isnothing(raw) ?
    default_multimode_iicsi_diagnostics_config() :
    merge(default_multimode_iicsi_diagnostics_config(),raw)

function multimode_iicsi_sample(grid,cfg)
    left,right=cfg.window_fraction
    selected=findall(x->left*grid.length<=x<=right*grid.length,grid.x)
    isempty(selected)&&error("multimode IICSI diagnostic window is empty")
    density=grid.ni[selected]
    density=density.-sum(density)/length(density)
    values=Float64[]
    for mode in cfg.modes
        phase=2pi*mode.*grid.x[selected]./grid.length
        normalization=2/length(selected)
        cosine=normalization*sum(density.*cos.(phase))
        sine=-normalization*sum(density.*sin.(phase))
        append!(values,(cosine,sine,hypot(cosine,sine)))
    end
    values
end

multimode_iicsi_history_path(history_file)=
    joinpath(dirname(history_file),"multimode_iicsi.tsv")

function append_multimode_iicsi_history!(path,step,time,grid,cfg;new_file=false)
    cfg.enabled||return nothing
    (new_file || step%cfg.every==0)||return nothing
    values=multimode_iicsi_sample(grid,cfg)
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        if new_file
            columns=String["step","time_omega_pe"]
            for mode in cfg.modes
                append!(columns,("mode$(mode)_cos","mode$(mode)_sin","mode$(mode)_amplitude"))
            end
            println(io,join(columns,'\t'))
        end
        println(io,join((step,time,values...), '\t'))
    end
    nothing
end
