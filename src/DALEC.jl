module DALEC

using ClimaCore
using ClimaLSM
using ClimaLSM.Domains
using NetCDF
using CSV
using DataFrames

import ClimaLSM: name, make_rhs, prognostic_vars, prognostic_types
import ClimaLSM.Domains: coordinates


abstract type AbstractDALECModel{FT} <: AbstractModel{FT} end

include("auxi/ACM.jl")
include("auxi/seasonality.jl")
include("models/DALEC811.jl")
include("util/dalec_io.jl")

export DALEC811, load_dalec_811, load_initial_condition!

end # module
