"""
    load_dalec_811(CBF_PATH::String, PARAM_PATH::String, ::Type{FT}) where {FT <: AbstractFloat}

Load meterological forcing inputs (cbf file in NetCDf format) and dalec parameters (in csv format) to
a DALEC811 struct prior to running the simulation.
"""
function load_dalec_811(CBF_PATH::String, PARAM_PATH::String, ::Type{FT}) where {FT <: AbstractFloat}
    ID=convert(Int, ncread(CBF_PATH, "ID")[1])
    if ID != 811
        error("CBF driver file does not mach ID==811, check if the correct load function is used.")
    end
    
    TIME = convert.(FT, ncread(CBF_PATH, "time"))
    T_MIN = convert.(FT, ncread(CBF_PATH, "T2M_MIN"))
    T_MAX = convert.(FT, ncread(CBF_PATH, "T2M_MAX"))
    RADIATION = convert.(FT, ncread(CBF_PATH, "SSRD"))
    ATMOSPHERIC_CO2 = convert.(FT, ncread(CBF_PATH, "CO2"))
    DOY = convert.(FT, ncread(CBF_PATH, "DOY"))
    BURNED_AREA = convert.(FT, ncread(CBF_PATH, "BURNED_AREA"))
    VPD = convert.(FT, ncread(CBF_PATH, "VPD"))
    PRECIPITATION = convert.(FT, ncread(CBF_PATH, "TOTAL_PREC"))

    LATITUDE=convert(FT, ncread(CBF_PATH, "LAT")[1])
    DELTA_T = convert(Int, round(TIME[2] - TIME[1]))
    MEAN_TEMP = sum((T_MIN + T_MAX)/2) / length(T_MIN)
    MEAN_PRECIP = sum(PRECIPITATION) / length(PRECIPITATION)
    
    NODAYS = convert(Int, length(TIME))

    dalec_params = FT.(CSV.read(PARAM_PATH, DataFrame).dalec811_params);
    check_dalec_811_parameter_bounds(dalec_params)

    dalec = DALEC811(TIME, T_MIN, T_MAX, RADIATION, ATMOSPHERIC_CO2, DOY, BURNED_AREA, VPD,
        PRECIPITATION, LATITUDE, DELTA_T, MEAN_TEMP, MEAN_PRECIP, NODAYS, dalec_params...)

    return dalec
end

"""
    check_dalec_811_parameter_bounds(params::Vector{FT}) where {FT <: AbstractFloat}

Validate if a parameter vector is within the parameter bounds of the dalec 811 model.
"""
function check_dalec_811_parameter_bounds(params::Vector{FT}) where {FT <: AbstractFloat}
    parnames = dalec_811_parnames()
    parmax = dalec_811_parmax()
    parmin = dalec_811_parmin()
    
    npar = length(params)
    
    if npar != length(parnames)
        error("The length of the param vector should be 33 for dalec 811 model!")
    end
    
    for i in 1:npar
        if (params[i] < parmin[i] || params[i] > parmax[i])
            error("Parameter no. " * string(i) * ": " * string(parnames[i]) * " is out of bound with min: " * string(parmin[i]) * " and max: " * string(parmax[i]) * ".") 

        end
    end
end

"""
    load_initial_condition!(model::DALEC811{FT}, Y::AbstractVector{FT}) where {FT <: AbstractFloat}

Load initia parameters into the Y vector.
"""
function load_initial_condition!(model::DALEC811{FT}, Y::AbstractVector{FT}) where {FT <: AbstractFloat}
    Y.dalec811.next_labile_pool[1] = model.Clab
    Y.dalec811.next_foliar_pool[1] = model.Cfol
    Y.dalec811.next_root_pool[1] = model.Croot
    Y.dalec811.next_wood_pool[1] = model.Cwood
    Y.dalec811.next_litter_pool[1] = model.Clitter
    Y.dalec811.next_som_pool[1] = model.Csom
    Y.dalec811.next_water_pool[1] = model.initial_water
    return nothing
end