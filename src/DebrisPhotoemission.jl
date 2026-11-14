# =============================================================================
# OPTIONAL DEBRIS PHOTOEMISSION
# =============================================================================
# The disabled path is an immediate no-op and consumes no random numbers, so
# photoemission=false preserves the original IICSI trajectory exactly.

mutable struct PhotoemissionState
    macroparticle_remainder::Float64
    emitted_macroparticles::Int
    emitted_physical_electrons::Float64
end

PhotoemissionState() = PhotoemissionState(0.0, 0, 0.0)

default_photoemission_config() = (
    enabled=false, photon_flux_m2_s=0.0, photon_energy_ev=20.0,
    quantum_yield=0.0, illumination_fraction=0.0, sunward_side=:left,
    reference_density_m3=1.0e11, reference_temperature_ev=1.0,
    effective_emitting_area_m2=nothing, cold_temperature_ev=2.7,
    hot_temperature_ev=10.0, hot_fraction=0.05,
)

const PHOTO_ELECTRON_CHARGE_C = 1.602176634e-19
const PHOTO_ELECTRON_MASS_KG = 9.1093837139e-31
const PHOTO_VACUUM_PERMITTIVITY = 8.8541878128e-12

photoelectron_omega_pe(density_m3) = sqrt(density_m3 * PHOTO_ELECTRON_CHARGE_C^2 /
    (PHOTO_ELECTRON_MASS_KG * PHOTO_VACUUM_PERMITTIVITY))

photoelectron_debye_length(temperature_ev, density_m3) =
    sqrt(PHOTO_VACUUM_PERMITTIVITY * temperature_ev /
         (density_m3 * PHOTO_ELECTRON_CHARGE_C))

function effective_photoemitting_area(photo, debris_radius_debye)
    !isnothing(photo.effective_emitting_area_m2) &&
        return photo.effective_emitting_area_m2
    radius_m = debris_radius_debye *
        photoelectron_debye_length(photo.reference_temperature_ev,
                                   photo.reference_density_m3)
    return pi * radius_m^2 # projected sunlit area of a spherical grain
end

function sample_photoelectron_energy_ev(rng, photo)
    temperature = rand(rng) < photo.hot_fraction ?
                  photo.hot_temperature_ev : photo.cold_temperature_ev
    # Flux-emitted Maxwellian normal energy.
    return -temperature * log(max(rand(rng), eps(Float64)))
end

"Emit charge-conserving electron macroparticles from the illuminated debris face."
function emit_debris_photoelectrons!(electrons, debris, debris_cfg, photo,
                                     state::PhotoemissionState, dt, center,
                                     radius, domain_length, rng;
                                     rate_multiplier=1.0)
    photo.enabled || return 0
    rate_multiplier >= 0 || error("photoemission rate_multiplier must be nonnegative")
    omega_pe = photoelectron_omega_pe(photo.reference_density_m3)
    area = effective_photoemitting_area(photo, radius)
    physical_electrons = rate_multiplier * photo.photon_flux_m2_s * photo.quantum_yield *
        photo.illumination_fraction * area * dt / omega_pe
    state.macroparticle_remainder +=
        physical_electrons / debris_cfg.physical_particles_per_macroparticle
    number = floor(Int, state.macroparticle_remainder)
    state.macroparticle_remainder -= number
    number == 0 && return 0

    launch_offset = max(1.0e-9, 1.0e-6 * radius)
    for _ in 1:number
        # :both represents uniform emission from the two faces in this 1D
        # geometry. Each macroparticle independently chooses either face with
        # equal probability and is launched outward from that face.
        side_sign = photo.sunward_side === :left ? -1.0 :
                    photo.sunward_side === :right ? 1.0 :
                    (rand(rng) < 0.5 ? -1.0 : 1.0)
        launch_position = clamp(center + side_sign * (radius + launch_offset),
                                0.0, domain_length)
        energy_ev = sample_photoelectron_energy_ev(rng, photo)
        speed = sqrt(2 * energy_ev / photo.reference_temperature_ev)
        push!(electrons.x, launch_position)
        push!(electrons.v, side_sign * speed)
    end

    # Internal transfer: plasma gains -q and debris gains +q.
    debris.total_charge_normalized += number * electrons.weight
    physical_charge_e = number * debris_cfg.physical_particles_per_macroparticle
    debris.total_charge_e += physical_charge_e
    state.emitted_macroparticles += number
    state.emitted_physical_electrons += physical_charge_e
    return number
end
