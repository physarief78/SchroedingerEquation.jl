using SchroedingerEquation
using Test
using SparseArrays
using LinearAlgebra
using Unitful
using Base.Threads

@testset "SchroedingerEquation.jl Tests" begin
    @testset "Basis and Hamiltonian Construction" begin
        basis = RealSpaceGrid1D(-5.0, 5.0, 100)
        potential = CustomPotential(x -> 0.5 * x^2)
        ham = build_hamiltonian(basis, potential)
        
        @test ham.basis.N == 100
        @test size(ham.H_matrix) == (100, 100)
    end

    @testset "TISE Solver - Harmonic Oscillator" begin
        basis = RealSpaceGrid1D(-10.0, 10.0, 500)
        potential = HarmonicPotential(1.0) # omega=1.0
        ham = build_hamiltonian(basis, potential)
        
        energies, wavefunctions = solve_tise(ham, 3)
        
        # Ground state energy should be approx 0.5 (hbar=1, m=1, omega=1)
        @test energies[1] ≈ 0.5 atol=1e-2
        @test energies[2] ≈ 1.5 atol=1e-2
        @test energies[3] ≈ 2.5 atol=1e-2

        # Test expectation value of x for ground state (should be 0 due to symmetry)
        @test expectation_value(wavefunctions[1], x -> x) ≈ 0.0 atol=1e-5
    end

    @testset "TDSE Solver" begin
        basis = RealSpaceGrid1D(-5.0, 5.0, 100)
        potential = CustomPotential(x -> 0.0) # Free particle
        ham = build_hamiltonian(basis, potential)
        
        # Gaussian wavepacket
        psi0_vals = exp.(-basis.x.^2) .+ 0im
        psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
        normalize!(psi0)
        
        history = propagate_tdse(psi0, ham, 0.1, 0.01; save_stride=2)
        
        @test length(history) == 6 # 0, 2, 4, 6, 8, 10 step indices (actually i starts at 1, so steps 2, 4, 6, 8, 10 are pushed)
        # i=1: no, i=2: yes, ..., i=10: yes. Plus initial state. Total 6.
        @test abs(sum(probability_density(history[end])) * basis.dx - 1.0) < 1e-10 # Norm preservation
    end

    @testset "Momentum Operator" begin
        basis = RealSpaceGrid1D(-10.0, 10.0, 500)
        # Ground state of HO is exp(-x^2/2) / pi^(1/4)
        # Momentum expectation value should be 0
        potential = HarmonicPotential(1.0)
        ham = build_hamiltonian(basis, potential)
        _, wfs = solve_tise(ham, 1)
        
        P = build_momentum_operator(basis)
        @test real(expectation_value(wfs[1], P)) ≈ 0.0 atol=1e-5
    end

    @testset "Time-Dependent TDSE" begin
        basis = RealSpaceGrid1D(-10.0, 10.0, 200)
        # Static part: Harmonic Oscillator
        static_pot = HarmonicPotential(1.0)
        # Dynamic part: Oscillating field
        dynamic_pot = OscillatingField(0.1, 1.0)
        combined = CombinedPotential(static_pot, dynamic_pot)
        
        # Kinetic energy matrix
        # Need to build it manually or use build_hamiltonian and extract H_matrix - V_matrix
        ham_static = build_hamiltonian(basis, static_pot)
        V_static = spdiagm(0 => [static_pot(xi) for xi in basis.x])
        T_matrix = ham_static.H_matrix - V_static
        
        psi0_vals = exp.(-basis.x.^2) .+ 0im
        psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
        normalize!(psi0)
        
        history = propagate_tdse(psi0, basis, T_matrix, combined, 0.5, 0.05)
        
        @test length(history) == 11
        @test abs(sum(probability_density(history[end])) * basis.dx - 1.0) < 1e-10
    end

    @testset "Observables Tracking" begin
        basis = RealSpaceGrid1D(-5.0, 5.0, 100)
        potential = HarmonicPotential(1.0)
        ham = build_hamiltonian(basis, potential)
        
        psi0_vals = exp.(-(basis.x .- 1.0).^2) .+ 0im # Displaced Gaussian
        psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
        normalize!(psi0)
        
        # Track position <x> using NamedTuple
        obs = (x = wf -> expectation_value(wf, x -> x),)
        
        history, obs_history = propagate_tdse(psi0, ham, 0.5, 0.05; observables=obs)
        
        @test hasproperty(obs_history, :x)
        @test length(obs_history.x) == 11 # 0.5/0.05 + 1
        # Initial position should be around 1.0
        @test obs_history.x[1] ≈ 1.0 atol=1e-2
    end

    @testset "Memory Efficient TDSE (No History)" begin
        basis = RealSpaceGrid1D(-5.0, 5.0, 100)
        potential = HarmonicPotential(1.0)
        ham = build_hamiltonian(basis, potential)
        
        psi0_vals = exp.(-basis.x.^2) .+ 0im
        psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
        normalize!(psi0)
        
        obs = (norm = wf -> sum(abs2, wf.psi) * wf.basis.dx,)
        
        # keep_history=false should return only obs_history
        result = propagate_tdse(psi0, ham, 0.2, 0.05; observables=obs, keep_history=false)
        
        @test result isa NamedTuple
        @test hasproperty(result, :norm)
        @test length(result.norm) == 5 # 0.2/0.05 + 1
        @test result.norm[end] ≈ 1.0 atol=1e-10
    end

    @testset "Complex Absorbing Potential (CAP)" begin
        basis = RealSpaceGrid1D(-10.0, 10.0, 200)
        # CAP at the edges
        cap = AbsorbingBoundary(2.0)
        potential = HarmonicPotential(1.0)
        
        # Free particle wavepacket moving towards the right boundary
        psi0_vals = exp.(-(basis.x .+ 5.0).^2 ./ 2.0) .* exp.(im * 5.0 .* basis.x)
        psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
        normalize!(psi0)
        
        # Static Hamiltonian with CAP
        ham = build_hamiltonian(basis, potential, cap)
        
        # Propagate
        history = propagate_tdse(psi0, ham, 2.0, 0.1)
        
        # Norm should decrease as the wavepacket hits the CAP
        final_norm = sum(abs2, history[end].psi) * basis.dx
        @test final_norm < 1.0 # Significant absorption
    end

    @testset "Unitful Integration" begin
        # Define in nanometers
        basis = RealSpaceGrid1D(-5.0u"nm", 5.0u"nm", 500)
        # Harmonic potential in eV/nm^2
        potential = HarmonicPotential(1.0u"eV/nm^2")
        
        ham = build_hamiltonian(basis, potential; m=1.0u"me", hbar=1.0u"hbar")
        
        energies, wavefunctions = solve_tise(ham, 1)
        
        # The solver works in atomic units internally
        # Result should be consistent with physical units
        @test energies[1] > 0.0
    end

    @testset "SSFM Solver" begin
        basis = RealSpaceGrid1D(-10.0, 10.0, 512) # Power of 2 for FFT speed
        potential = HarmonicPotential(1.0)
        
        # Ground state
        ham = build_hamiltonian(basis, potential)
        _, wfs = solve_tise(ham, 1)
        psi0 = wfs[1]
        
        # Propagate with SSFM
        # Ground state should be stationary (only phase change)
        history = propagate_ssfm(psi0, potential, 1.0, 0.01) # Smaller dt
        
        final_psi2 = probability_density(history[end])
        initial_psi2 = probability_density(psi0)
        
        diff = sum(abs.(final_psi2 .- initial_psi2)) * basis.dx
        println("SSFM Density diff: ", diff)
        @test diff < 1e-3 # Looser tolerance for finite-diff vs spectral mismatch
    end

    @testset "Parallel Thread Safety" begin
        basis = RealSpaceGrid1D(-10.0, 10.0, 256)
        potential = HarmonicPotential(1.0)
        ham = build_hamiltonian(basis, potential)
        _, wfs = solve_tise(ham, 1)
        psi0 = wfs[1]
        
        n_runs = 8
        results = Vector{Vector{ComplexF64}}(undef, n_runs)
        
        # Run in parallel
        @threads for i in 1:n_runs
            ws = TDSEWorkspace(ham, 0.1)
            history = propagate_tdse(psi0, ham, 0.5, 0.1; workspace=ws, keep_history=true)
            results[i] = history[end].psi
        end
        
        # Verify all results are identical
        for i in 2:n_runs
            @test results[1] ≈ results[i]
        end
    end
end
