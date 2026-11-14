module SLChargingReduction

using FFTW, LinearAlgebra, Statistics, Printf, Serialization

export fit_sl_charging_reduction, fit_sl_charging_global

function load_history(path)
    lines=readlines(path)
    header=split(lines[1],'\t')
    rows=[split(line,'\t') for line in lines[2:end]]
    Dict(name => [parse(Float64,row[j]) for row in rows]
         for (j,name) in enumerate(header))
end

function fft_filter(x,dt,low,high)
    n=length(x); spectrum=rfft(x)
    frequency=(0:length(spectrum)-1)./(n*dt)
    spectrum[(frequency.<low).|(frequency.>high)].=0
    irfft(spectrum,n)
end

# The oscillator quadrature y=dx/dtau is constructed spectrally after the
# narrow-band projection. Regression itself uses integral increments below.
function spectral_derivative(x,dt)
    n=length(x); spectrum=rfft(x)
    frequency=(0:length(spectrum)-1)./(n*dt)
    spectrum .*= im.*(2pi.*frequency)
    iseven(n)&&(spectrum[end]=0)
    irfft(spectrum,n)
end

function detrend(x)
    t=collect(1.0:length(x)); X=hcat(ones(length(x)),t)
    x-X*(X\x)
end

function dominant_frequency(x,dt,band)
    power=abs2.(rfft(detrend(x)))
    frequency=(0:length(power)-1)./(length(x)*dt)
    ids=findall((frequency.>=band[1]).&(frequency.<=band[2]))
    isempty(ids)&&error("empty IICSI frequency band")
    frequency[ids[argmax(power[ids])]]
end

function ridge_fit(X,target,ridge)
    scale=vec(sqrt.(sum(abs2,X;dims=1)./size(X,1)))
    scale[scale.<eps()].=1
    Xs=X./reshape(scale,1,:)
    ((Xs'Xs+ridge*I)\(Xs'target))./scale
end

"Ridge regression with one coefficient constrained to be nonnegative."
function ridge_fit_nonnegative(X,target,ridge,positive_column)
    coefficient=ridge_fit(X,target,ridge)
    coefficient[positive_column]>=0 && return coefficient
    keep=setdiff(axes(X,2),positive_column)
    constrained=zeros(size(X,2))
    constrained[keep]=ridge_fit(X[:,keep],target,ridge)
    constrained
end

"Cumulative trapezoidal integral on a uniformly sampled grid."
function cumulative_trapezoid(x,step)
    result=zeros(length(x))
    @inbounds for i in 2:length(x)
        result[i]=result[i-1]+step*(x[i-1]+x[i])/2
    end
    result
end

function interpolate_linear(source_time,source_values,target_time)
    result=similar(target_time,Float64)
    @inbounds for (j,t) in pairs(target_time)
        i=searchsortedlast(source_time,t)
        if i<1
            result[j]=source_values[1]
        elseif i>=length(source_time)
            result[j]=source_values[end]
        else
            weight=(t-source_time[i])/(source_time[i+1]-source_time[i])
            result[j]=(1-weight)*source_values[i]+weight*source_values[i+1]
        end
    end
    result
end

"Load the leading coherent ion-density mode from serialized Vlasov frames."
function add_ion_density_mode!(history,path,start_time,end_time)
    isfile(path)||error("Vlasov snapshots required for ion_density_mode: $path")
    times=Float64[]; profiles=Vector{Float64}[]
    open(path,"r") do io
        while !eof(io)
            frame=deserialize(io)
            start_time<=frame.time<=end_time || continue
            push!(times,Float64(frame.time))
            push!(profiles,vec(sum(frame.ion_f;dims=2)))
        end
    end
    length(times)>=20||error("too few Vlasov frames in the SL fitting window")
    order=sortperm(times); times=times[order]; profiles=profiles[order]
    density=reduce(hcat,profiles)'
    density .-= mean(density;dims=1)
    factor=svd(density;full=false)
    score=factor.U[:,1].*factor.S[1]
    # Fix the otherwise arbitrary SVD sign for reproducible coefficients.
    sum(score.*collect(1.0:length(score)))<0 && (score .*= -1)
    history["ion_density_mode"]=interpolate_linear(times,score,
        history["time_omega_pe"])
    return factor.S[1]^2/sum(abs2,factor.S)
end

r2(y,prediction)=1-sum(abs2,y.-prediction)/max(sum(abs2,y.-mean(y)),eps())
nrmse(y,prediction)=sqrt(mean(abs2,y.-prediction))/max(std(y),eps())

function fit_observable(history,name,cfg,ids,dt)
    raw=history[String(name)][ids]
    f0=dominant_frequency(raw,dt,cfg.frequency_band)
    df=1/(length(raw)*dt)
    halfwidth=max(cfg.carrier_halfwidth_fraction*f0,2df)
    x=fft_filter(detrend(raw),dt,max(cfg.frequency_band[1],f0-halfwidth),
                 min(cfg.frequency_band[2],f0+halfwidth))
    xscale=sqrt(2)*std(x)
    xscale>eps()||error("$name is constant in the selected carrier band")
    x./=xscale
    omega0=2pi*f0
    y=spectral_derivative(x,dt)./omega0

    q=detrend(history["Qd_over_e"][ids])
    q=fft_filter(q,dt,0.0,cfg.charge_lowpass_frequency)
    z=q./cfg.charge_scale_e

    edge=max(100,round(Int,cfg.edge_fraction*length(x)))
    keep=(edge+1):(length(x)-edge)
    x,y,z=x[keep],y[keep],z[keep]
    time=history["time_omega_pe"][ids][keep]

    # Integrate the exact centered forms of Eqs. (15)-(17) over short windows:
    # Delta z = integral(C*x+D*y-Lambda*z) d tau
    # Delta y + integral(x)d tau = integral[eps(1-x^2)y-H*x*z+K*z]d tau.
    # This changes only the estimator, not the theoretical model library.
    dτ=omega0*dt
    window=max(2,round(Int,cfg.integral_window_omega_pe/dt))
    window<length(x)-20||error("integral fitting window is too long")
    charge_terms=(x,y,-z)
    oscillator_terms=((1 .-x.^2).*y,-x.*z,z)
    charge_integrals=[cumulative_trapezoid(v,dτ) for v in charge_terms]
    oscillator_integrals=[cumulative_trapezoid(v,dτ) for v in oscillator_terms]
    x_integral=cumulative_trapezoid(x,dτ)
    starts=1:(length(x)-window); stops=(1+window):length(x)
    Xcharge=hcat([(v[stops].-v[starts]) for v in charge_integrals]...)
    Xosc=hcat([(v[stops].-v[starts]) for v in oscillator_integrals]...)
    target_charge=z[stops].-z[starts]
    target_osc=y[stops].-y[starts].+x_integral[stops].-x_integral[starts]
    interval_time=time[stops]
    split=clamp(floor(Int,cfg.training_fraction*length(starts)),10,length(starts)-10)
    train=1:split
    test=split+1:length(starts)
    charge=ridge_fit_nonnegative(Xcharge[train,:],target_charge[train],cfg.ridge,3)
    oscillator=ridge_fit_nonnegative(Xosc[train,:],target_osc[train],cfg.ridge,1)
    charge_prediction=Xcharge*charge
    oscillator_prediction=Xosc*oscillator
    (;name,f0,omega0,xscale,time=interval_time,x=x[stops],y=y[stops],z=z[stops],
      target_charge,target_osc,charge,oscillator,window,
      Xcharge,Xosc,
      charge_train_r2=r2(target_charge[train],charge_prediction[train]),
      charge_test_r2=r2(target_charge[test],charge_prediction[test]),
      oscillator_train_r2=r2(target_osc[train],oscillator_prediction[train]),
      oscillator_test_r2=r2(target_osc[test],oscillator_prediction[test]),
      charge_test_nrmse=nrmse(target_charge[test],charge_prediction[test]),
      oscillator_test_nrmse=nrmse(target_osc[test],oscillator_prediction[test]),
      charge_prediction,oscillator_prediction,split)
end

"Global in-sample SINDy-like regression with no train/test partition."
function fit_sl_charging_global(parameters)
    hasproperty(parameters,:sl_charging_reduction)||return nothing
    cfg=parameters.sl_charging_reduction
    history_file=parameters.checkpoint.history_file
    isfile(history_file)||error("SL reduction history not found: $history_file")
    history=load_history(history_file)
    t=history["time_omega_pe"]
    finish=isfinite(cfg.end_time_omega_pe) ? cfg.end_time_omega_pe : maximum(t)
    modal_variance=nothing
    if :ion_density_mode in cfg.observables
        modal_variance=add_ion_density_mode!(history,
            parameters.sindy_export.vlasov.file,cfg.start_time_omega_pe,finish)
    end
    ids=findall((t.>=cfg.start_time_omega_pe).&(t.<=finish))
    length(ids)>=1000||error("SL global-fit interval contains too few samples")
    dt=median(diff(t[ids]))
    prepared=[fit_observable(history,name,cfg,ids,dt) for name in cfg.observables]
    results=map(prepared) do r
        charge=ridge_fit_nonnegative(r.Xcharge,r.target_charge,cfg.ridge,3)
        oscillator=ridge_fit_nonnegative(r.Xosc,r.target_osc,cfg.ridge,1)
        charge_fit=r.Xcharge*charge
        oscillator_fit=r.Xosc*oscillator
        merge(r,(charge=charge,oscillator=oscillator,
            charge_global_r2=r2(r.target_charge,charge_fit),
            oscillator_global_r2=r2(r.target_osc,oscillator_fit),
            charge_global_nrmse=nrmse(r.target_charge,charge_fit),
            oscillator_global_nrmse=nrmse(r.target_osc,oscillator_fit)))
    end
    score(r)=min(r.charge_global_r2,r.oscillator_global_r2)
    best=results[argmax(score.(results))]
    output_directory=dirname(history_file)
    table=joinpath(output_directory,"sl_reduction_global_fit.tsv")
    open(table,"w") do io
        println(io,"observable\tf0\tomega0\tx_scale\tcharge_scale_e\tC\tD\tLambda\tepsilon\tH\tK\tcharge_global_r2\toscillator_global_r2\tcharge_global_nrmse\toscillator_global_nrmse")
        for r in results
            C,D,Lambda=r.charge; epsilon,H,K=r.oscillator
            println(io,join((r.name,r.f0,r.omega0,r.xscale,cfg.charge_scale_e,
                C,D,Lambda,epsilon,H,K,r.charge_global_r2,
                r.oscillator_global_r2,r.charge_global_nrmse,
                r.oscillator_global_nrmse),'\t'))
        end
    end
    summary=joinpath(output_directory,"sl_reduction_global_summary.txt")
    open(summary,"w") do io
        C,D,Lambda=best.charge; epsilon,H,K=best.oscillator
        println(io,"Global in-sample SINDy-like SL charging fit")
        println(io,"window = $(cfg.start_time_omega_pe) .. $finish")
        println(io,"partition = none; all usable samples enter one regression")
        println(io,"best_observable = $(best.name)")
        println(io,"f0 = $(best.f0)")
        println(io,"omega0 = $(best.omega0)")
        println(io,"epsilon = $epsilon")
        println(io,"H = $H")
        println(io,"K = $K")
        println(io,"C = $C")
        println(io,"D = $D")
        println(io,"Lambda = $Lambda")
        println(io,"charge_global_r2 = $(best.charge_global_r2)")
        println(io,"oscillator_global_r2 = $(best.oscillator_global_r2)")
        println(io,"charge_global_nrmse = $(best.charge_global_nrmse)")
        println(io,"oscillator_global_nrmse = $(best.oscillator_global_nrmse)")
        modal_variance===nothing||println(io,
            "ion_density_mode_variance_fraction = $modal_variance")
        println(io,"These are descriptive in-sample residual scores, not predictive validation.")
    end
    @printf("Global SL charging fit: best=%s charge R2=% .3f oscillator R2=% .3f\n",
        String(best.name),best.charge_global_r2,best.oscillator_global_r2)
    println("  wrote $table")
    (;best,results,table,summary)
end

function fit_sl_charging_reduction(parameters)
    hasproperty(parameters,:sl_charging_reduction)||return nothing
    cfg=parameters.sl_charging_reduction
    cfg.enabled||return nothing
    history_file=parameters.checkpoint.history_file
    isfile(history_file)||error("SL reduction history not found: $history_file")
    history=load_history(history_file)
    t=history["time_omega_pe"]
    finish=isfinite(cfg.end_time_omega_pe) ? cfg.end_time_omega_pe : maximum(t)
    modal_variance=nothing
    if :ion_density_mode in cfg.observables
        hasproperty(parameters,:sindy_export)&&parameters.sindy_export.vlasov.enabled||
            error("ion_density_mode requires sindy_export.vlasov.enabled=true")
        modal_variance=add_ion_density_mode!(history,
            parameters.sindy_export.vlasov.file,cfg.start_time_omega_pe,finish)
    end
    ids=findall((t.>=cfg.start_time_omega_pe).&(t.<=finish))
    length(ids)>=1000||error("SL reduction interval contains too few samples")
    dt=median(diff(t[ids]))

    results=[fit_observable(history,name,cfg,ids,dt) for name in cfg.observables]
    score(r)=min(r.charge_test_r2,r.oscillator_test_r2)
    best=results[argmax(score.(results))]
    output_directory=dirname(history_file)
    table=joinpath(output_directory,"sl_reduction_fit.tsv")
    open(table,"w") do io
        println(io,"observable\tf0\tomega0\tx_scale\tcharge_scale_e\tC\tD\tLambda\tepsilon\tH\tK\tcharge_train_r2\tcharge_test_r2\toscillator_train_r2\toscillator_test_r2\tcharge_test_nrmse\toscillator_test_nrmse")
        for r in results
            C,D,Lambda=r.charge; epsilon,H,K=r.oscillator
            println(io,join((r.name,r.f0,r.omega0,r.xscale,cfg.charge_scale_e,
                C,D,Lambda,epsilon,H,K,r.charge_train_r2,r.charge_test_r2,
                r.oscillator_train_r2,r.oscillator_test_r2,
                r.charge_test_nrmse,r.oscillator_test_nrmse),'\t'))
        end
    end

    trajectory=joinpath(output_directory,"sl_reduction_best_trajectory.tsv")
    open(trajectory,"w") do io
        println(io,"interval_end_time_omega_pe\ttau\txi\ty\tz\tdelta_z\tdelta_z_fit\tdelta_y_plus_integral_x\tdelta_y_plus_integral_x_fit\tpartition")
        for i in eachindex(best.time)
            partition=i<=best.split ? "train" : "test"
            println(io,join((best.time[i],best.omega0*best.time[i],best.x[i],
                best.y[i],best.z[i],best.target_charge[i],best.charge_prediction[i],
                best.target_osc[i],best.oscillator_prediction[i],partition),'\t'))
        end
    end

    summary=joinpath(output_directory,"sl_reduction_summary.txt")
    open(summary,"w") do io
        C,D,Lambda=best.charge; epsilon,H,K=best.oscillator
        println(io,"Stuart-Landau charging reduction")
        println(io,"window = $(cfg.start_time_omega_pe) .. $finish")
        println(io,"best_observable = $(best.name)")
        println(io,"f0 = $(best.f0)")
        println(io,"omega0 = $(best.omega0)")
        println(io,"charge_scale_e = $(cfg.charge_scale_e)")
        println(io,"integral_window_omega_pe = $(cfg.integral_window_omega_pe)")
        modal_variance===nothing||println(io,"ion_density_mode_variance_fraction = $modal_variance")
        println(io,"epsilon = $epsilon")
        println(io,"H = $H")
        println(io,"K = $K")
        println(io,"C = $C")
        println(io,"D = $D")
        println(io,"Lambda = $Lambda")
        println(io,"charge_test_r2 = $(best.charge_test_r2)")
        println(io,"oscillator_test_r2 = $(best.oscillator_test_r2)")
        println(io,"The fit is diagnostic only; poor held-out R2 rejects the reduction.")
    end
    @printf("SL charging reduction: best=%s f0=%.5f charge R2test=% .3f oscillator R2test=% .3f\n",
            String(best.name),best.f0,best.charge_test_r2,best.oscillator_test_r2)
    println("  wrote $table")
    (;best,results,table,trajectory,summary)
end

end
