# =============================================================================
# OPTIONAL FLOW-CORRECTED CHARGE-DEPENDENT COLLECTION
# =============================================================================
# The 1D PIC field already accelerates particles longitudinally. This plug-in
# supplies only the unresolved transverse collection-area correction: an
# attractive debris potential focuses trajectories, whereas a repulsive
# surface barrier reduces their probability of reaching the material.
# enabled=false returns the original MCC probability bit-for-bit and consumes
# no random numbers.

default_charge_dependent_collection_config() = (
    enabled=false,
    model=:flow_corrected_oml,
    sampling_offset_debye=1.0,
    barrier_smoothing_normalized=0.20,
    electron_coupling_strength=0.10,
    ion_coupling_strength=0.02,
    maximum_cross_section_multiplier=1.20,
    minimum_cross_section_multiplier=0.50,
)

charge_dependent_collection_config(raw=nothing) = isnothing(raw) ?
    default_charge_dependent_collection_config() :
    merge(default_charge_dependent_collection_config(),raw)

struct ChargeDependentCollectionSample
    phi_debris::Float64
    phi_left_reference::Float64
    phi_right_reference::Float64
end

mutable struct ChargeDependentCollectionState
    electron_multiplier_sum::Float64
    electron_candidates::Int
    ion_multiplier_sum::Float64
    ion_candidates::Int
    cumulative_electron_candidates::Int
    cumulative_ion_candidates::Int
end
ChargeDependentCollectionState()=
    ChargeDependentCollectionState(0.0,0,0.0,0,0,0)

function begin_charge_dependent_collection_step!(state)
    state.electron_multiplier_sum=0.0; state.electron_candidates=0
    state.ion_multiplier_sum=0.0; state.ion_candidates=0
    return nothing
end

function charge_dependent_collection_sample(grid,pic_cfg,cfg)
    cfg.enabled || return ChargeDependentCollectionSample(0.0,0.0,0.0)
    center=pic_cfg.debris.center_fraction*pic_cfg.domain_length
    radius=pic_cfg.debris.size_debye
    offset=cfg.sampling_offset_debye
    phi_d=scalar_at_position(grid.phi,center,grid,pic_cfg)
    phi_l=scalar_at_position(grid.phi,max(0.0,center-radius-offset),grid,pic_cfg)
    phi_r=scalar_at_position(grid.phi,
        min(pic_cfg.domain_length,center+radius+offset),grid,pic_cfg)
    ChargeDependentCollectionSample(phi_d,phi_l,phi_r)
end

"Apply the unresolved orbital/barrier correction to an MCC probability."
function charge_dependent_capture_probability(base_probability,kinetic_energy,
        species_charge,face,sample,cfg,state)
    cfg.enabled || return base_probability
    phi_reference=face<0 ? sample.phi_left_reference : sample.phi_right_reference
    potential_energy=species_charge*(sample.phi_debris-phi_reference)
    coupling=species_charge<0 ? cfg.electron_coupling_strength :
        cfg.ion_coupling_strength
    denominator=max(kinetic_energy+cfg.barrier_smoothing_normalized,eps())
    ratio=potential_energy/denominator
    multiplier=ratio>=0 ? exp(-coupling*ratio) : 1+coupling*(-ratio)
    multiplier=clamp(multiplier,cfg.minimum_cross_section_multiplier,
        cfg.maximum_cross_section_multiplier)
    if species_charge<0
        state.electron_multiplier_sum+=multiplier
        state.electron_candidates+=1
        state.cumulative_electron_candidates+=1
    else
        state.ion_multiplier_sum+=multiplier
        state.ion_candidates+=1
        state.cumulative_ion_candidates+=1
    end
    clamp(base_probability*multiplier,0.0,1.0)
end

charge_dependent_collection_history_path(history_file)=
    joinpath(dirname(history_file),"charge_dependent_collection.tsv")

function append_charge_dependent_collection_history!(path,step,time,state,sample;
        new_file=false)
    me=state.electron_candidates==0 ? 1.0 :
        state.electron_multiplier_sum/state.electron_candidates
    mi=state.ion_candidates==0 ? 1.0 :
        state.ion_multiplier_sum/state.ion_candidates
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        new_file&&println(io,"step\ttime_omega_pe\tphi_debris\t" *
            "phi_left_reference\tphi_right_reference\t" *
            "electron_candidates\tion_candidates\t" *
            "mean_electron_multiplier\tmean_ion_multiplier\t" *
            "cumulative_electron_candidates\tcumulative_ion_candidates")
        println(io,join((step,time,sample.phi_debris,sample.phi_left_reference,
            sample.phi_right_reference,state.electron_candidates,
            state.ion_candidates,me,mi,state.cumulative_electron_candidates,
            state.cumulative_ion_candidates),'\t'))
    end
    nothing
end
