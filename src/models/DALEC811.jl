struct DALEC811{FT} <: AbstractDALECModel{FT}
    # Meterological Forcings
    TIME::Vector{FT}
    T_MIN::Vector{FT}
    T_MAX::Vector{FT}
    RADIATION::Vector{FT}
    ATMOSPHERIC_CO2::Vector{FT}
    DOY::Vector{FT}
    BURNED_AREA::Vector{FT}
    VPD::Vector{FT}
    PRECIPITATION::Vector{FT}
    
    # Constant Variables
    LATITUDE::FT
    DELTA_T::Int
    MEAN_TEMP::FT
    MEAN_PRECIP::FT
    NODAYS::Int
    
    # Trainable DALEC Parameters
    decomposition_rate::FT
    f_gpp::FT
    f_fol::FT
    f_root::FT
    leaf_lifespan::FT
    tor_wood::FT
    tor_root::FT
    tor_litter::FT
    tor_som::FT
    Q10::FT
    canopy_efficiency::FT
    Bday::FT
    f_lab::FT
    clab_release_period::FT
    Fday::FT
    leaf_fall_period::FT
    LMCA::FT
    Clab::FT
    Cfol::FT
    Croot::FT
    Cwood::FT
    Clitter::FT
    Csom::FT
    IWUE::FT
    runoff_focal_point::FT
    wilting_point::FT
    initial_water::FT
    foliar_cf::FT
    ligneous_cf::FT
    dom_cf::FT
    resilience::FT
    lab_lifespan::FT
    moisture_factor::FT 
end

ClimaLSM.name(model::DALEC811) = :dalec811;

ClimaLSM.prognostic_vars(::DALEC811) = (:LAI, :GPP, :ET, :temperate, :respiration_auto, :leaf_production, :labile_production, 
    :root_production, :wood_production, :lff, :lrf, :labile_release, :leaf_litter, :wood_litter,
      :root_litter, :respiration_hetero_litter, :respiration_hetero_som, :litter_to_som, :runoff,
       :labile_fire_combust, :foliar_fire_combust, :root_fire_combust, :wood_fire_combust,
        :litter_fire_combust, :som_fire_combust, :labile_fire_transfer, :foliar_fire_transfer,
         :root_fire_transfer, :wood_fire_transfer, :litter_fire_transfer, :total_fire_combust,
          :nee, :next_labile_pool, :next_foliar_pool, :next_root_pool, :next_wood_pool, :next_litter_pool,
           :next_som_pool, :next_water_pool);

# currently the model runs at single point site.
ClimaLSM.Domains.coordinates(model::DALEC811{FT}) where {FT} = FT.([0.0]);

ClimaLSM.prognostic_types(::DALEC811{FT}) where {FT} = (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT,
    FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT);

    ClimaLSM.auxiliary_vars(::DALEC811) = (:LAI, :GPP, :ET, :temperate, :respiration_auto, :leaf_production, :labile_production, 
    :root_production, :wood_production, :lff, :lrf, :labile_release, :leaf_litter, :wood_litter,
      :root_litter, :respiration_hetero_litter, :respiration_hetero_som, :litter_to_som, :runoff,
       :labile_fire_combust, :foliar_fire_combust, :root_fire_combust, :wood_fire_combust,
        :litter_fire_combust, :som_fire_combust, :labile_fire_transfer, :foliar_fire_transfer,
         :root_fire_transfer, :wood_fire_transfer, :litter_fire_transfer, :total_fire_combust,
          :nee, :next_labile_pool, :next_foliar_pool, :next_root_pool, :next_wood_pool, :next_litter_pool,
           :next_som_pool, :next_water_pool);

ClimaLSM.auxiliary_types(::DALEC811{FT}) where {FT} = (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT,
    FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT);
    
"""
    ClimaLSM.make_rhs(model::DALEC811{FT}) where {FT}

The method extends ClimaLSM.make_rhs. The ode function is discrete and is expected to run on weekly to monthly temporal resolution.
To intrage the model, the function should be cast in a DiscreteProblem and intetrated using FunctionMap{true} of the OrdinaryDiffEq module.
"""
function ClimaLSM.make_rhs(model::DALEC811{FT}) where {FT}
    function rhs!(dY, Y, p, t)
        
        # normalize foliar pool by LCMA (leaf carbon mass per area) to obtain LAI. 
        p.dalec811.LAI[1] = Y.dalec811.next_foliar_pool[1] / model.LMCA
        
        # compute gpp FLUXES[0]
        p.dalec811.GPP[1] = ACM(lat = model.LATITUDE, doy = model.DOY[t], t_max = model.T_MAX[t],
            t_min = model.T_MIN[t], lai = p.dalec811.LAI[1], rad = model.RADIATION[t],
            ca = model.ATMOSPHERIC_CO2[t], ce = model.canopy_efficiency, FT = FT) * min(Y.dalec811.next_water_pool[1] / model.wilting_point, 1)
        
        # compute ET_flux FLUXES[28]
        p.dalec811.ET[1] = p.dalec811.GPP[1] * model.VPD[t] / model.IWUE 
        
        # compute temperate FLUXES[1]
        p.dalec811.temperate[1] = exp(model.Q10 * (FT(0.5) * (model.T_MAX[t] + model.T_MIN[t]) - model.MEAN_TEMP)) * ((model.PRECIPITATION[t] / model.MEAN_PRECIP - FT(1)) * model.moisture_factor + FT(1))
        
        # compute autotrophic respiration FLUXES[2]
        p.dalec811.respiration_auto[1] = model.f_gpp * p.dalec811.GPP[1]
        
        # compute leaf production FLUXES[3]
        p.dalec811.leaf_production[1] = (p.dalec811.GPP[1] - p.dalec811.respiration_auto[1]) * model.f_fol
        
        # compute labile production FLUXES[4]
        p.dalec811.labile_production[1] = (p.dalec811.GPP[1] - p.dalec811.respiration_auto[1] - p.dalec811.leaf_production[1]) * model.f_lab
        
        # compute root production FLUXES[5]
        p.dalec811.root_production[1] = (p.dalec811.GPP[1] - p.dalec811.respiration_auto[1] - p.dalec811.leaf_production[1] - p.dalec811.labile_production[1]) * model.f_root
    
        # compute wood production FLUXES[6]
        p.dalec811.wood_production[1] = p.dalec811.GPP[1] -  p.dalec811.respiration_auto[1] - p.dalec811.leaf_production[1] - p.dalec811.labile_production[1] - p.dalec811.root_production[1]
        
        # compute leaf fall factor FLUXES[8]
        p.dalec811.lff[1] = leaf_fall_factor(model.TIME[t], model.leaf_lifespan, model.leaf_fall_period, model.Fday, FT)
        
        # compute labile release factor FLUXES[15]
        p.dalec811.lrf[1] = lab_release_factor(model.TIME[t], model.lab_lifespan, model.clab_release_period, model.Bday, FT)
        
        # compute labile release FLUXES[7]
        p.dalec811.labile_release[1] = Y.dalec811.next_labile_pool[1] * (FT(1) - (FT(1) - p.dalec811.lrf[1]) ^ model.DELTA_T) / model.DELTA_T
        
        # compute leaf litter production FLUXES[9]
        p.dalec811.leaf_litter[1] = Y.dalec811.next_foliar_pool[1] * (FT(1) - (FT(1) - p.dalec811.lff[1]) ^ model.DELTA_T) / model.DELTA_T
        
        # compute wood litter production FLUXES[10]
        p.dalec811.wood_litter = Y.dalec811.next_wood_pool[1] * (FT(1) - (FT(1) - model.tor_wood)^ model.DELTA_T) / model.DELTA_T

        # compute root litter production FLUXES[11]
        p.dalec811.root_litter[1] = Y.dalec811.next_root_pool[1] * (FT(1) - (FT(1) - model.tor_root) ^ model.DELTA_T) / model.DELTA_T
        
        # compute respiration heterotrophic litter FLUXES[12]
        p.dalec811.respiration_hetero_litter[1] = Y.dalec811.next_litter_pool[1] * (FT(1) - (FT(1) - p.dalec811.temperate[1] * model.tor_litter) ^ model.DELTA_T) / model.DELTA_T
        
        # compute respiration heterotrophic SOM FLUXES[13]
        p.dalec811.respiration_hetero_som[1] = Y.dalec811.next_som_pool[1] * (FT(1) - (FT(1) - p.dalec811.temperate[1] * model.tor_som) ^ model.DELTA_T) / model.DELTA_T

        # compute litter to som flux
        p.dalec811.litter_to_som[1] = Y.dalec811.next_litter_pool[1] * (FT(1) - (FT(1) - p.dalec811.temperate[1] * model.decomposition_rate) ^ model.DELTA_T) / model.DELTA_T
        
        # compute runoff flux
        p.dalec811.runoff[1] = Y.dalec811.next_water_pool[1] ^ FT(2) / model.runoff_focal_point / model.DELTA_T
        
        if Y.dalec811.next_water_pool[1] > model.runoff_focal_point / FT(2)
            p.dalec811.runoff[1] = (Y.dalec811.next_water_pool[1] - model.runoff_focal_point / FT(4)) / model.DELTA_T
        end
        
        # all pools before including fire.
        p.dalec811.next_labile_pool[1] = Y.dalec811.next_labile_pool[1] + (p.dalec811.labile_production[1] - p.dalec811.labile_release[1]) * model.DELTA_T
        p.dalec811.next_foliar_pool[1] = Y.dalec811.next_foliar_pool[1] + (p.dalec811.leaf_production[1] - p.dalec811.leaf_litter[1] + p.dalec811.labile_release[1]) * model.DELTA_T
        p.dalec811.next_root_pool[1] = Y.dalec811.next_root_pool[1] + (p.dalec811.root_production[1] - p.dalec811.root_litter[1]) * model.DELTA_T
        p.dalec811.next_wood_pool[1] = Y.dalec811.next_wood_pool[1] + (p.dalec811.wood_production[1] - p.dalec811.wood_litter[1]) * model.DELTA_T
        p.dalec811.next_litter_pool[1] = Y.dalec811.next_litter_pool[1] + (p.dalec811.leaf_litter[1] + p.dalec811.root_litter[1] - p.dalec811.respiration_hetero_litter[1] - p.dalec811.litter_to_som[1]) * model.DELTA_T
        p.dalec811.next_som_pool[1] = Y.dalec811.next_som_pool[1] + (p.dalec811.litter_to_som[1] - p.dalec811.respiration_hetero_som[1] + p.dalec811.wood_litter[1]) * model.DELTA_T
        p.dalec811.next_water_pool[1] = Y.dalec811.next_water_pool[1] - p.dalec811.runoff[1] * model.DELTA_T + model.PRECIPITATION[t] * model.DELTA_T - p.dalec811.ET[1] * model.DELTA_T
    
        
        p.dalec811.labile_fire_combust[1] = p.dalec811.next_labile_pool[1] * model.BURNED_AREA[t] * model.ligneous_cf / model.DELTA_T
        p.dalec811.foliar_fire_combust[1] = p.dalec811.next_foliar_pool[1] * model.BURNED_AREA[t] * model.foliar_cf / model.DELTA_T
        p.dalec811.root_fire_combust[1] = p.dalec811.next_root_pool[1] * model.BURNED_AREA[t] * model.ligneous_cf / model.DELTA_T
        p.dalec811.wood_fire_combust[1] = p.dalec811.next_wood_pool[1] * model.BURNED_AREA[t] * model.ligneous_cf / model.DELTA_T
        p.dalec811.litter_fire_combust[1] = p.dalec811.next_litter_pool[1] * model.BURNED_AREA[t] * (model.ligneous_cf + model.foliar_cf) * FT(0.5) / model.DELTA_T
        p.dalec811.som_fire_combust[1] = p.dalec811.next_som_pool[1] * model.BURNED_AREA[t] * model.dom_cf / model.DELTA_T
        
        p.dalec811.labile_fire_transfer[1] = p.dalec811.next_labile_pool[1] * model.BURNED_AREA[t] *(FT(1) - model.ligneous_cf) * (FT(1) - model.resilience) / model.DELTA_T
        p.dalec811.foliar_fire_transfer[1] = p.dalec811.next_foliar_pool[1] * model.BURNED_AREA[t] *(FT(1)- model.foliar_cf) * (FT(1)- model.resilience) / model.DELTA_T
        p.dalec811.root_fire_transfer[1] = p.dalec811.next_root_pool[1] * model.BURNED_AREA[t] *(FT(1) - model.ligneous_cf) * (FT(1) - model.resilience) / model.DELTA_T
        p.dalec811.wood_fire_transfer[1] = p.dalec811.next_wood_pool[1] * model.BURNED_AREA[t] * (FT(1) - model.ligneous_cf)*(FT(1) - model.resilience) / model.DELTA_T
        p.dalec811.litter_fire_transfer[1] = p.dalec811.next_litter_pool[1] * model.BURNED_AREA[t] * (FT(1) -(model.ligneous_cf + model.foliar_cf) * FT(0.5)) * (FT(1) - model.resilience) / model.DELTA_T
        
        # include fire combust and fire transfer fluxes to the carbon pools/
        p.dalec811.next_labile_pool[1] = p.dalec811.next_labile_pool[1] - (p.dalec811.labile_fire_combust[1] + p.dalec811.labile_fire_transfer[1]) * model.DELTA_T
        p.dalec811.next_foliar_pool[1] = p.dalec811.next_foliar_pool[1] - (p.dalec811.foliar_fire_combust[1] + p.dalec811.foliar_fire_transfer[1]) * model.DELTA_T
        p.dalec811.next_root_pool[1] = p.dalec811.next_root_pool[1] - (p.dalec811.root_fire_combust[1] + p.dalec811.root_fire_transfer[1]) * model.DELTA_T
        p.dalec811.next_wood_pool[1] = p.dalec811.next_wood_pool[1] - (p.dalec811.wood_fire_combust[1] + p.dalec811.wood_fire_transfer[1]) * model.DELTA_T
        p.dalec811.next_litter_pool[1] = p.dalec811.next_litter_pool[1] + (p.dalec811.labile_fire_transfer[1] + p.dalec811.foliar_fire_transfer[1] + p.dalec811.root_fire_transfer[1] - p.dalec811.litter_fire_combust[1] - p.dalec811.litter_fire_transfer[1]) * model.DELTA_T
        p.dalec811.next_som_pool[1] = p.dalec811.next_som_pool[1] + (p.dalec811.wood_fire_transfer[1] + p.dalec811.litter_fire_transfer[1] - p.dalec811.som_fire_combust[1]) * model.DELTA_T
    
        p.dalec811.total_fire_combust[1] = p.dalec811.labile_fire_combust[1] + p.dalec811.foliar_fire_combust[1] + p.dalec811.root_fire_combust[1] + p.dalec811.wood_fire_combust[1] + p.dalec811.litter_fire_combust[1] + p.dalec811.som_fire_combust[1]
    
        p.dalec811.nee[1] = -p.dalec811.GPP[1] + p.dalec811.respiration_auto[1] + p.dalec811.respiration_hetero_litter[1] + p.dalec811.respiration_hetero_som[1] + p.dalec811.total_fire_combust[1]
        
        dY.dalec811.LAI[1] = p.dalec811.LAI[1] - Y.dalec811.LAI[1]
        dY.dalec811.GPP[1] = p.dalec811.GPP[1] - Y.dalec811.GPP[1]
        dY.dalec811.ET[1] = p.dalec811.ET[1] - Y.dalec811.ET[1]
        dY.dalec811.temperate[1] = p.dalec811.temperate[1] - Y.dalec811.temperate[1]
        dY.dalec811.respiration_auto[1] = p.dalec811.respiration_auto[1] - Y.dalec811.respiration_auto[1]
        dY.dalec811.leaf_production[1] = p.dalec811.leaf_production[1] - Y.dalec811.leaf_production[1]
        dY.dalec811.labile_production[1] = p.dalec811.labile_production[1] - Y.dalec811.labile_production[1]
        dY.dalec811.root_production[1] = p.dalec811.root_production[1] - Y.dalec811.root_production[1]
        dY.dalec811.wood_production[1] = p.dalec811.wood_production[1] - Y.dalec811.wood_production[1]
        dY.dalec811.lff[1] = p.dalec811.lff[1] - Y.dalec811.lff[1]
        dY.dalec811.lrf[1] = p.dalec811.lrf[1] - Y.dalec811.lrf[1]
        dY.dalec811.labile_release[1] = p.dalec811.labile_release[1] - Y.dalec811.labile_release[1]
        dY.dalec811.leaf_litter[1] = p.dalec811.leaf_litter[1] - Y.dalec811.leaf_litter[1]
        dY.dalec811.wood_litter[1] = p.dalec811.wood_litter[1] - Y.dalec811.wood_litter[1]
        dY.dalec811.root_litter[1] = p.dalec811.root_litter[1] - Y.dalec811.root_litter[1]
        dY.dalec811.respiration_hetero_litter[1] = p.dalec811.respiration_hetero_litter[1] - Y.dalec811.respiration_hetero_litter[1]
        dY.dalec811.respiration_hetero_som[1] = p.dalec811.respiration_hetero_som[1] - Y.dalec811.respiration_hetero_som[1]
        dY.dalec811.litter_to_som[1] = p.dalec811.litter_to_som[1] - Y.dalec811.litter_to_som[1]
        dY.dalec811.runoff[1] = p.dalec811.runoff[1] - Y.dalec811.runoff[1]
        dY.dalec811.labile_fire_combust[1] = p.dalec811.labile_fire_combust[1] - Y.dalec811.labile_fire_combust[1]
        dY.dalec811.foliar_fire_combust[1] = p.dalec811.foliar_fire_combust[1] - Y.dalec811.foliar_fire_combust[1]
        dY.dalec811.root_fire_combust[1] = p.dalec811.root_fire_combust[1] - Y.dalec811.root_fire_combust[1]
        dY.dalec811.wood_fire_combust[1] = p.dalec811.wood_fire_combust[1] - Y.dalec811.wood_fire_combust[1]
        dY.dalec811.litter_fire_combust[1] = p.dalec811.litter_fire_combust[1] - Y.dalec811.litter_fire_combust[1]
        dY.dalec811.som_fire_combust[1] = p.dalec811.som_fire_combust[1] - Y.dalec811.som_fire_combust[1]
        dY.dalec811.labile_fire_transfer[1] = p.dalec811.labile_fire_transfer[1] - Y.dalec811.labile_fire_transfer[1]
        dY.dalec811.foliar_fire_transfer[1] = p.dalec811.foliar_fire_transfer[1] - Y.dalec811.foliar_fire_transfer[1]
        dY.dalec811.root_fire_transfer[1] = p.dalec811.root_fire_transfer[1] - Y.dalec811.root_fire_transfer[1]
        dY.dalec811.wood_fire_transfer[1] = p.dalec811.wood_fire_transfer[1] - Y.dalec811.wood_fire_transfer[1]
        dY.dalec811.litter_fire_transfer[1] = p.dalec811.litter_fire_transfer[1] - Y.dalec811.litter_fire_transfer[1]
        dY.dalec811.total_fire_combust[1] = p.dalec811.total_fire_combust[1] - Y.dalec811.total_fire_combust[1]
        dY.dalec811.nee[1] = p.dalec811.nee[1] - Y.dalec811.nee[1]
        dY.dalec811.next_labile_pool[1] = p.dalec811.next_labile_pool[1] - Y.dalec811.next_labile_pool[1]
        dY.dalec811.next_foliar_pool[1] = p.dalec811.next_foliar_pool[1] - Y.dalec811.next_foliar_pool[1]
        dY.dalec811.next_root_pool[1] = p.dalec811.next_root_pool[1] - Y.dalec811.next_root_pool[1]
        dY.dalec811.next_wood_pool[1] = p.dalec811.next_wood_pool[1] - Y.dalec811.next_wood_pool[1]
        dY.dalec811.next_litter_pool[1] = p.dalec811.next_litter_pool[1] - Y.dalec811.next_litter_pool[1]
        dY.dalec811.next_som_pool[1] = p.dalec811.next_som_pool[1] - Y.dalec811.next_som_pool[1]
        dY.dalec811.next_water_pool[1] = p.dalec811.next_water_pool[1] - Y.dalec811.next_water_pool[1]   
    end
end


"""
    dalec_811_parmin()

Return minimum parameter bounds for the dalec 811 model.
"""
function dalec_811_parmin()
    parmin=[0.00001, 0.2, 0.01, 0.01, 1.001, 0.000025, 0.0001, 0.0001,
    0.0000001, 0.018, 5, 365.25, 0.01, 365.25/12, 365.25, 365.25/12,
     5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10, 1, 1, 1, 0.01, 0.01, 0.01, 0.01, 1.001, 0.01]
     return parmin
end

"""
    dalec_811_parmax()

Return maximum parameter bounds for the dalec 811 model.
"""
function dalec_811_parmax()
   parmax = [0.01, 0.8, 0.5, 1, 8, 0.001, 0.01, 0.01, 0.001, 0.08, 50, 365.25*4,
    0.5, 100, 365.25 * 4, 150, 200, 2000.0, 2000.0, 2000.0, 100000.0, 2000.0,
     200000.0, 50, 100000, 10000, 10000, 1, 1, 1, 1, 8, 1]
   return parmax
end

"""
    dalec_811_parnames()

Return trainable parameter names for the dalec 811 model.
"""
function dalec_811_parnames()
   parnames = (:decomposition_rate, :f_gpp, :f_fol, :f_root, :leaf_lifespan, :tor_wood, :tor_root,
   :tor_litter, :tor_som, :Q10, :canopy_efficiency, :Bday, :f_lab, :clab_release_period,
    :Fday, :leaf_fall_period, :LMCA, :Clab, :Cfol, :Croot, :Cwood, :Clitter, :Csom,
     :IWUE, :runoff_focal_point, :wilting_point, :initial_water, :foliar_cf, :ligneous_cf,
      :dom_cf, :resilience, :lab_lifespan, :moisture_factor)
   return parnames
end