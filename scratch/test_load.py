import os
import sys

project_path = os.path.abspath(".")
print(f"Project path: {project_path}")

from juliacall import Main as jl
jl.eval(f'push!(LOAD_PATH, "{project_path.replace("\\\\", "/")}")')
print("LOAD_PATH updated.")

jl.eval('using SchroedingerEquation')
print("Successfully loaded SchroedingerEquation!")
