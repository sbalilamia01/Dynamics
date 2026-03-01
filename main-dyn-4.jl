using Pkg
Pkg.activate(".")
Pkg.instantiate()

using MPSGE
using DataFrames
using PlotlyJS
import JuMP.Containers: DenseAxisArray, @container
import JuMP
using Ipopt

include("Ramsey-dyn-4.jl")

#
#  DYN4 uses balanced growth

d4_data  = IndexedData();
d4_model = dynamic_dyn4_model(d4_data);

solve!(d4_model, cumulative_iteration_limit = 0)
