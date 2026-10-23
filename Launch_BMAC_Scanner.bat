@echo off
title BMAC Lab Scanner
cd C:\bmac-julia
start julia +1.12 --project=C:\BindingAndCatalysis.jl -e "include(\"C:\\\\bmac-julia\\\\bmac_final_tool.jl\")"
timeout /t 8
start http://localhost:8270