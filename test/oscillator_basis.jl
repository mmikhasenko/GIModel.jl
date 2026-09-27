@testset "The oscillator basis has a fixed phase" begin
    # LAPACK's QR picks its own column signs, and they varied with nbasis, beta
    # and the mesh (flips at n = 10 and n = 18 for beta = 0.65, nbasis = 24).
    # Invisible numerically -- a sign flip is unitary -- but fatal once a closed
    # form enters, since analytic elements are written in ho_reduced_radial's
    # convention. Wiring exact p^2 into an unfixed basis put charmonium 1S
    # 9.4 MeV BELOW the finite-difference answer, which a variational
    # calculation in a finite basis cannot do.
    r, h = GIModel.radial_grid(450, 24.0)
    for β in (0.25, 0.65, 1.35, 2.35), nb in (6, 12, 24)
        B = GIModel.ho_basis_matrix(0, β, r, nb)
        U = GIModel.orthonormalize_physical_basis(B, h)
        # every column keeps the sign of the basis function it came from
        @test all(sum(view(U, :, j) .* view(B, :, j)) > 0 for j in 1:nb)
        # and it is still orthonormal under the physical inner product
        @test maximum(abs, h * (transpose(U) * U) - I) < 1e-10
    end

    # With the phase fixed, the projected p^2 agrees with the closed form up to
    # mesh error, and that error is a property of the mesh vs beta: fine at the
    # variational optimum, poor where the basis is too compact or too diffuse
    # for the grid (the resolution wall and the beta railing, as numbers).
    ana(β, nb) = Matrix(ho_p2_matrix(0, β, nb))
    proj(β, nb) = Matrix(GIModel.projected_matrix(
        GIModel.orthonormalize_physical_basis(GIModel.ho_basis_matrix(0, β, r, nb), h),
        h, GIModel.p2_operator(1.0, 0, r, h)))
    rel(β, nb) = maximum(abs, proj(β, nb) - ana(β, nb)) / maximum(abs, ana(β, nb))
    @test rel(0.65, 12) < 0.01          # near the variational optimum
    @test rel(0.65, 24) < 0.01
    @test rel(2.35, 24) > 0.05          # too compact for h: the resolution wall
    @test rel(0.25, 24) > 1.0           # too diffuse for rmax: the beta railing
end

@testset "Oscillator matrix elements validate themselves" begin
    # These are a self-contained mathematics problem: matrix elements of
    # operators in the 3D oscillator basis. They can be — and here are —
    # validated with no mesh, no reference data and no reference to the GI
    # model at all. Doing that BEFORE any table comparison means a later
    # disagreement with the paper is about the physics, not about whether the
    # algebra is right.

    # 1. The strongest check available: p^2 and r^2 must RECONSTRUCT the
    #    oscillator Hamiltonian. H = p^2/2mu + (1/2) mu w^2 r^2 is diagonal with
    #    eigenvalue (2n+L+3/2)w, and beta^2 = mu*w, so
    #        p^2/2mu + (beta^4/2mu) r^2
    #    must be EXACTLY diagonal with exactly those eigenvalues. An error in
    #    either operator's magnitude, sign, or power of beta destroys the
    #    cancellation, so this tests both operators and their relative
    #    normalization simultaneously.
    for L in (0, 1, 2), β in (0.35, 0.65, 1.35, 2.35)
        nb, μ = 24, 0.7                      # any mu; beta^2 = mu*w fixes w
        ω = β^2 / μ
        H = Matrix(ho_p2_matrix(L, β, nb)) ./ (2μ) .+
            (β^4 / (2μ)) .* Matrix(ho_r2_matrix(L, β, nb))
        @test maximum(abs, H - Diagonal(diag(H))) < 1e-12          # exactly diagonal
        @test maximum(abs, diag(H) .- [(2n + L + 1.5) * ω for n = 0:(nb-1)]) < 1e-12
    end

    # 2. Virial theorem for the oscillator: <T> = <V> in every eigenstate.
    for L in (0, 2), β in (0.45, 1.85)
        nb, μ = 16, 1.3
        T = Matrix(ho_p2_matrix(L, β, nb)) ./ (2μ)
        V = (β^4 / (2μ)) .* Matrix(ho_r2_matrix(L, β, nb))
        for n in 1:nb
            @test isapprox(T[n, n], V[n, n]; rtol = 1e-12)
        end
    end

    # 3. Dimensional scaling: p^2 ~ beta^2, r^2 ~ 1/beta^2, so their product is
    #    beta-independent — a check no single operator can provide alone.
    for L in (0, 1)
        A = Matrix(ho_p2_matrix(L, 0.5, 10)) * Matrix(ho_r2_matrix(L, 0.5, 10))
        B = Matrix(ho_p2_matrix(L, 1.9, 10)) * Matrix(ho_r2_matrix(L, 1.9, 10))
        @test isapprox(A, B; rtol = 1e-10)
    end

    # 4. Structure: symmetric, tridiagonal, and r^2 positive definite.
    for op in (ho_p2_matrix(1, 0.8, 9), ho_r2_matrix(1, 0.8, 9))
        M = Matrix(op)
        @test M ≈ transpose(M)
        @test all(iszero, [M[i, j] for i in 1:9, j in 1:9 if abs(i - j) > 1])
        @test all(>(0), eigvals(Symmetric(M)))
    end
    @test_throws ArgumentError ho_r2_matrix(0, -1.0, 4)
    @test_throws ArgumentError ho_r2_matrix(0, 0.5, 0)
end

@testset "A17 position side: Gauss-Laguerre in Golub-Welsch form" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))

    # The Jacobi matrix for weight x^(L+1/2) e^{-x}, x = (beta r)^2, IS
    # beta^2 * ho_r2_matrix. Diagonalizing it gives nodes as eigenvalues and
    # sqrt(w_i) p_n(x_i) as eigenvector entries -- so the 1e19 polynomial value
    # and the 1e-95 weight never exist separately. That is the whole trick:
    # forming them apart loses the digits before any summation happens, which no
    # compensated summation can undo.
    for L in (0, 1, 2), β in (0.25, 0.65, 2.35)
        # g = 1 must be the identity, because Z is orthogonal
        @test maximum(abs, Matrix(ho_operator_matrix(L, β, 24, r -> 1.0)) - I) < 1e-13
        # g = r^2 must reconstruct the very Jacobi matrix that generated the rule
        @test maximum(abs, Matrix(ho_operator_matrix(L, β, 24, r -> r^2)) -
                           Matrix(ho_r2_matrix(L, β, 24))) < 1e-10
    end

    # Higher pure moments too, against the Gaussian-moment closed form
    # <n|r^0|n> = 1 and the r^2 case above; r^4 is checked for consistency
    # between two independent quadrature sizes rather than a closed form.
    for L in (0, 1), β in (0.45, 1.35)
        a = Matrix(ho_operator_matrix(L, β, 12, r -> r^4; nq = 64))
        b = Matrix(ho_operator_matrix(L, β, 12, r -> r^4; nq = 256))
        @test maximum(abs, a - b) < 1e-9
    end

    # The real operator: the Appendix-A smeared Coulomb + confinement, checked
    # against independent adaptive quadrature (QuadGK), not against itself.
    mc = mq["c"]
    g(r) = GIModel.smeared_coulomb_G_closed(params, mc, mc, r) +
           GIModel.smeared_confinement_S_closed(params, mc, mc, r)
    L, β, nb = 0, 0.65, 24
    M = Matrix(ho_operator_matrix(L, β, nb, g))
    for (a, b) in ((1, 1), (1, 2), (3, 3), (5, 8), (12, 12))
        ref, _ = quadgk(r -> GIModel.ho_reduced_radial(a - 1, L, β, r) * g(r) *
                             GIModel.ho_reduced_radial(b - 1, L, β, r),
                        0, Inf; rtol = 1e-13, order = 21)
        @test isapprox(M[a, b], ref; atol = 1e-11)
    end

    # Symmetry, and independence of the quadrature size once converged.
    @test M ≈ transpose(M)
    @test maximum(abs, M - Matrix(ho_operator_matrix(L, β, nb, g; nq = 512))) < 1e-10
    @test_throws ArgumentError ho_operator_matrix(0, -1.0, 4, r -> 1.0)
    @test_throws ArgumentError ho_operator_matrix(0, 0.5, 0, r -> 1.0)

    # The rule itself carries no beta. In x = (beta r)^2 the Jacobi matrix is
    # beta^2 * ho_r2_matrix(L, beta, n) = ho_r2_matrix(L, 1, n), every beta
    # cancelling, so nodes and DVR weights are functions of (L, nq) alone and
    # beta enters only as r_i = sqrt(x_i)/beta. `gauss_laguerre_dvr` memoizes on
    # that fact; if it ever stopped holding, the cache would silently hand one
    # beta's rule to another.
    for L in (0, 1, 2), β in (0.25, 0.65, 2.35)
        @test Matrix(β^2 * ho_r2_matrix(L, β, 40)) ≈ Matrix(ho_r2_matrix(L, 1, 40))
    end
    # The memo is capped by BYTES, not entry count: entries span 24x64 to
    # 24x8192, so 256 small ones cost less than one large one. Flushing can only
    # cost time -- a rebuilt entry is bit-identical, which is what makes the
    # whole memo safe to discard at any moment.
    let c = GIModel._GAUSS_LAGUERRE_DVR
        a = GIModel.gauss_laguerre_dvr(0, 24, 1024)
        empty!(c)
        b = GIModel.gauss_laguerre_dvr(0, 24, 1024)
        @test a[1] == b[1] && a[2] == b[2]
        @test GIModel._dvr_cache_bytes() == sizeof(b[1]) + sizeof(b[2])
    end

    sqrt_x, Z = GIModel.gauss_laguerre_dvr(1, 12, 64)
    @test Z * transpose(Z) ≈ I            # orthogonality of the leading rows
    @test issorted(sqrt_x)                # nodes come out ordered
    @test length(sqrt_x) == 64 && size(Z) == (12, 64)
end

@testset "Exact oscillator p^2 (A17 momentum side)" begin
    # p^2 = 2*mu*H_osc - beta^4 r^2 is tridiagonal in the oscillator basis:
    #   <n|p^2|n>   = beta^2 (2n + L + 3/2)
    #   <n|p^2|n+1> = beta^2 sqrt((n+1)(n + L + 3/2))
    # Checked against the mesh projection it is meant to replace -- and the
    # check is CONVERGENCE, not a tolerance: the difference must shrink as the
    # mesh refines, which is what shows the formula is exact and the mesh is
    # the approximation.
    for (L, β, nb) in ((0, 0.55, 6), (1, 0.75, 6), (2, 0.45, 5))
        ana = Matrix(ho_p2_matrix(L, β, nb))
        errs = Float64[]
        for (ng, rm) in ((450, 24.0), (2000, 40.0), (8000, 56.0))
            r, h = GIModel.radial_grid(ng, rm)
            U = GIModel.orthonormalize_physical_basis(GIModel.ho_basis_matrix(L, β, r, nb), h)
            num = Matrix(GIModel.projected_matrix(U, h, GIModel.p2_operator(1.0, L, r, h)))
            push!(errs, maximum(abs, num .- ana))
        end
        @test errs[1] > errs[2] > errs[3]        # monotone convergence to the formula
        @test errs[3] < 1e-3
    end

    # Structure: symmetric, tridiagonal, positive definite (p^2 is).
    M = Matrix(ho_p2_matrix(0, 0.6, 8))
    @test M ≈ transpose(M)
    @test all(iszero, [M[i, j] for i in 1:8, j in 1:8 if abs(i - j) > 1])
    @test all(>(0), eigvals(Symmetric(M)))
    # beta scaling: p^2 has dimensions of beta^2.
    @test Matrix(ho_p2_matrix(0, 1.2, 5)) ≈ 4 .* Matrix(ho_p2_matrix(0, 0.6, 5))
    @test_throws ArgumentError ho_p2_matrix(0, -1.0, 4)
    @test_throws ArgumentError ho_p2_matrix(0, 0.5, 0)
end

@testset "Continuum momentum projection is not a finite p² matrix function" begin
    L, beta, nbasis = 1, 0.9, 4
    m = 0.22
    f(p) = sqrt(m^2 + p^2)
    op = GIModel._ho_momentum_operator_matrix(L, beta, nbasis, f)
    for (i, j) in ((1, 1), (1, 3), (3, 4))
        exact, _ = quadgk(0.0, 20.0; rtol = 1e-11) do p
            (-1.0)^(i+j-2) * GIModel.ho_reduced_radial(i-1, L, inv(beta), p) *
            GIModel.ho_reduced_radial(j-1, L, inv(beta), p) * f(p)
        end
        @test op[i, j] ≈ exact atol = 1e-9
    end
    @test GIModel._ho_momentum_operator_matrix(L, beta, nbasis, p -> p^2) ≈
          Matrix(GIModel.ho_p2_matrix(L, beta, nbasis)) atol = 1e-11
end
