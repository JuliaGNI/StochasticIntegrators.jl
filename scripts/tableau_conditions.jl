#
# Verification of the tableau coefficients of the strong (mean-square) methods against the
# conditions of
#
#   M. Kraus and T. M. Tyranowski,
#   "Variational integrators for stochastic dissipative Hamiltonian systems",
#   IMA Journal of Numerical Analysis (2021).
#
# Two independent sets of conditions are checked, both purely algebraic in the coefficients:
#
#   * the Lagrange-d'Alembert (symplecticity) conditions, Eq. (Symplectic conditions for SPRK),
#     which are what make a scheme a variational integrator;
#   * the conditions for mean-square convergence of order 1.0, Eq. (SPRK order conditions),
#     Theorem (Mean-square convergence theorem).
#
# This is deliberately separate from the test suite's energy-error checks. Those exercise the
# integrator; this exercises the coefficients, so a transcription error in a tableau is isolated
# from an error in the code that consumes it. Run it before the integrator tests, not after.
#
# Only four of the tableaus in this package come from that paper (Störmer-Verlet, the 2-stage
# DIRK, SRKw1 and SRKw2). The rest are checked *against its conditions* rather than against a
# result stated there, which is the stronger check: the expected outcomes below, including the
# two deliberate failures, are asserted rather than merely printed.
#
# Usage:  julia --project=. scripts/tableau_conditions.jl

using Printf
using RungeKutta: Tableau
using StochasticIntegrators
using Test

#
# The coefficients of the paper, extracted from a tableau.
#
# In the notation of Eq. (SPRK tableau) the layout is (a | ā | â | b | b̄ | b̂ ; α α α̂ β β β̂):
# `a`, `α` act on ∂H/∂p, `ā`, `b̄` on the Hamiltonian part of the p equation, and `â`, `b̂`,
# `α̂`, `β̂` on the forcing terms F and f.
#

struct Coefficients{T}
    name::Symbol
    a::Matrix{T}
    ā::Matrix{T}
    â::Matrix{T}
    b::Matrix{T}
    b̄::Matrix{T}
    b̂::Matrix{T}
    α::Vector{T}
    α̂::Vector{T}
    β::Vector{T}
    β̂::Vector{T}
end

mat(t::Tableau) = Matrix(t.a)
vec_(t::Tableau) = Vector(t.b)

"A split partitioned tableau carries all six matrices and all four weight vectors explicitly."
function Coefficients(tab::TableauSISPRK)
    Coefficients(tab.name,
        mat(tab.qdrift), mat(tab.pdrift1), mat(tab.pdrift2),
        mat(tab.qdiff), mat(tab.pdiff1), mat(tab.pdiff2),
        vec_(tab.qdrift), vec_(tab.pdrift2), vec_(tab.qdiff), vec_(tab.pdiff2))
end

"""
An unsplit partitioned tableau applies the *same* coefficients to the Hamiltonian and the forcing
parts of the p equation, so `â = ā`, `b̂ = b̄`, `α̂ = α`, `β̂ = β`.
"""
function Coefficients(tab::TableauSIPRK)
    Coefficients(tab.name,
        mat(tab.qdrift), mat(tab.pdrift), mat(tab.pdrift),
        mat(tab.qdiff), mat(tab.pdiff), mat(tab.pdiff),
        vec_(tab.qdrift), vec_(tab.pdrift), vec_(tab.qdiff), vec_(tab.pdiff))
end

"A non-partitioned tableau uses one drift and one diffusion tableau throughout."
function Coefficients(tab::TableauSIRK)
    Coefficients(tab.name,
        mat(tab.qdrift), mat(tab.qdrift), mat(tab.qdrift),
        mat(tab.qdiff), mat(tab.qdiff), mat(tab.qdiff),
        vec_(tab.qdrift), vec_(tab.qdrift), vec_(tab.qdiff), vec_(tab.qdiff))
end

#
# Eq. (Symplectic conditions for SPRK), for all i, j = 1 … s.
#
# Note the asymmetry: (1)-(4) pair ā, b̄ at index ij against a, b at index ji, and never mention
# â, b̂; (5)-(8) pair â, b̂ at ij against a, b at ji, and never mention ā, b̄.
#

const LDA_CONDITIONS = (
    ("(1)  αᵢ āᵢⱼ + αⱼ aⱼᵢ = αᵢ αⱼ",
        (c, i, j) -> c.α[i] * c.ā[i, j] + c.α[j] * c.a[j, i] -
                     c.α[i] * c.α[j]),
    ("(2)  βᵢ b̄ᵢⱼ + βⱼ bⱼᵢ = βᵢ βⱼ",
        (c, i, j) -> c.β[i] * c.b̄[i, j] + c.β[j] * c.b[j, i] -
                     c.β[i] * c.β[j]),
    ("(3)  βᵢ āᵢⱼ + αⱼ bⱼᵢ = βᵢ αⱼ",
        (c, i, j) -> c.β[i] * c.ā[i, j] + c.α[j] * c.b[j, i] -
                     c.β[i] * c.α[j]),
    ("(4)  αᵢ b̄ᵢⱼ + βⱼ aⱼᵢ = αᵢ βⱼ",
        (c, i, j) -> c.α[i] * c.b̄[i, j] + c.β[j] * c.a[j, i] -
                     c.α[i] * c.β[j]),
    ("(5)  αᵢ âᵢⱼ + α̂ⱼ aⱼᵢ = αᵢ α̂ⱼ",
        (c, i, j) -> c.α[i] * c.â[i, j] + c.α̂[j] * c.a[j, i] -
                     c.α[i] * c.α̂[j]),
    ("(6)  αᵢ b̂ᵢⱼ + β̂ⱼ aⱼᵢ = αᵢ β̂ⱼ",
        (c, i, j) -> c.α[i] * c.b̂[i, j] + c.β̂[j] * c.a[j, i] -
                     c.α[i] * c.β̂[j]),
    ("(7)  βᵢ âᵢⱼ + α̂ⱼ bⱼᵢ = βᵢ α̂ⱼ",
        (c, i, j) -> c.β[i] * c.â[i, j] + c.α̂[j] * c.b[j, i] -
                     c.β[i] * c.α̂[j]),
    ("(8)  βᵢ b̂ᵢⱼ + β̂ⱼ bⱼᵢ = βᵢ β̂ⱼ",
        (c, i, j) -> c.β[i] * c.b̂[i, j] + c.β̂[j] * c.b[j, i] -
                     c.β[i] * c.β̂[j])
)

"Residuals of the eight Lagrange-d'Alembert conditions; one entry per condition."
function lda_residuals(c::Coefficients)
    s = length(c.α)
    [maximum(abs(f(c, i, j)) for i in 1:s, j in 1:s) for (_, f) in LDA_CONDITIONS]
end

#
# Eq. (SPRK order conditions). Note that the drift matrices a, ā, â do not appear at all — only
# the four weight sums and the six β/β̂ × b/b̄/b̂ products.
#

const ORDER_CONDITIONS = (
    ("Σ αᵢ = 1", c -> sum(c.α) - 1),
    ("Σ α̂ᵢ = 1", c -> sum(c.α̂) - 1),
    ("Σ βᵢ = 1", c -> sum(c.β) - 1),
    ("Σ β̂ᵢ = 1", c -> sum(c.β̂) - 1),
    ("βᵀ b e  = 1/2",
        c -> sum(c.β[i] * c.b[i, j] for i in axes(c.b, 1), j in axes(c.b, 2)) -
             1 // 2),
    ("βᵀ b̄ e  = 1/2",
        c -> sum(c.β[i] * c.b̄[i, j] for i in axes(c.b̄, 1), j in axes(c.b̄, 2)) -
             1 // 2),
    ("βᵀ b̂ e  = 1/2",
        c -> sum(c.β[i] * c.b̂[i, j] for i in axes(c.b̂, 1), j in axes(c.b̂, 2)) -
             1 // 2),
    ("β̂ᵀ b e  = 1/2",
        c -> sum(c.β̂[i] * c.b[i, j] for i in axes(c.b, 1), j in axes(c.b, 2)) -
             1 // 2),
    ("β̂ᵀ b̄ e  = 1/2",
        c -> sum(c.β̂[i] * c.b̄[i, j] for i in axes(c.b̄, 1), j in axes(c.b̄, 2)) -
             1 // 2),
    ("β̂ᵀ b̂ e  = 1/2",
        c -> sum(c.β̂[i] * c.b̂[i, j] for i in axes(c.b̂, 1), j in axes(c.b̂, 2)) -
             1 // 2)
)

"Residuals of the ten mean-square order conditions."
order_residuals(c::Coefficients) = [abs(f(c)) for (_, f) in ORDER_CONDITIONS]

const ATOL = 1e-14

satisfies_lda(c) = all(<(ATOL), lda_residuals(c))
satisfies_order(c) = all(<(ATOL), order_residuals(c))

function report(c::Coefficients)
    lda = lda_residuals(c)
    ord = order_residuals(c)

    @printf("\n%s  (s = %d)\n", c.name, length(c.α))
    println("  Lagrange-d'Alembert:")
    for (k, (label, _)) in enumerate(LDA_CONDITIONS)
        @printf("    %-32s %-4s  max residual %.3e\n", label,
            lda[k] < ATOL ? "ok" : "FAIL", lda[k])
    end
    println("  Mean-square order 1.0:")
    for (k, (label, _)) in enumerate(ORDER_CONDITIONS)
        @printf("    %-32s %-4s  residual %.3e\n", label,
            ord[k] < ATOL ? "ok" : "FAIL", ord[k])
    end
end

#
# The schemes, with the outcome each is expected to produce.
#

const SCHEMES = [
    # (tableau, satisfies Lagrange-d'Alembert, satisfies order 1.0)
    (TableauStochasticGLRK(1), true, true),
    (TableauStochasticGLRK(2), true, true),
    (TableauStochasticGLRK(3), true, true),
    (TableauStochasticDIRK(0.0), true, true),
    (TableauStochasticDIRK(0.5), true, true),
    (TableauStochasticDIRK(1.0), true, true),
    (TableauStochasticDIRK(0.3), true, true),
    (TableauStochasticStoermerVerlet(), true, true),
    (TableauStochasticSymplecticEuler(), true, false),
    (TableauModifiedStochasticStoermerVerlet(0.0), true, true),
    (TableauModifiedStochasticStoermerVerlet(0.5), true, true),
    (TableauModifiedStochasticStoermerVerlet(1.0), true, true),
    (TableauStochasticLobattoIIIABD2(), false, true)
]

function main()
    for (tab, _, _) in SCHEMES
        report(Coefficients(tab))
    end

    println("\n", "="^78, "\n")

    @testset "Tableau conditions of Kraus & Tyranowski" begin
        for (tab, lda, ord) in SCHEMES
            c = Coefficients(tab)
            @testset "$(c.name)" begin
                @test satisfies_lda(c) == lda
                @test satisfies_order(c) == ord
            end
        end

        # The two deliberate failures, pinned down precisely rather than left as "not all pass".
        # A script that reported everything satisfied would be reporting a bug in itself.

        @testset "LobattoIIIABD2 fails exactly the forcing conditions" begin
            r = lda_residuals(Coefficients(TableauStochasticLobattoIIIABD2()))
            # (1)-(4) involve only the IIIA/IIIB pair and hold exactly
            @test all(iszero, r[1:4])
            # (5)-(8) are the conditions coupling the forcing tableaus â, b̂ to a, b. All four
            # fail, each by 1/8: with â = b̂ = LobattoIIID(2) = [1/4 -1/4; 3/4 1/4] and weights
            # (1/2, 1/2), at i = j = 1 the left side is 1/2·1/4 + 1/2·0 = 1/8 against 1/4.
            @test all(≈(1 / 8), r[5:8])
        end

        @testset "SymplecticEuler fails exactly the diffusion order conditions" begin
            c = Coefficients(TableauStochasticSymplecticEuler())
            r = order_residuals(c)
            # all four weight sums are 1
            @test all(<(ATOL), r[1:4])
            # but βᵀ b e = 1 ≠ 1/2 and βᵀ b̄ e = 0 ≠ 1/2
            @test sum(c.β[i] * c.b[i, j] for i in axes(c.b, 1), j in axes(c.b, 2)) ≈ 1
            @test sum(c.β[i] * c.b̄[i, j] for i in axes(c.b̄, 1), j in axes(c.b̄, 2)) ≈ 0
        end

        # The modified Störmer-Verlet family satisfies both sets across its whole parameter
        # range: c enters the order conditions solely through Σ α̂ᵢ = c + (1-c) = 1, and the
        # Lagrange-d'Alembert residuals are identically zero in c.
        @testset "ModifiedStoermerVerlet holds across its parameter range" begin
            for c in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
                coeff = Coefficients(TableauModifiedStochasticStoermerVerlet(c))
                @test satisfies_lda(coeff)
                @test satisfies_order(coeff)
            end
        end

        # Nothing in either set of conditions restricts c to [0,1] — the residuals above are
        # independent of it. `TableauModifiedStochasticStoermerVerlet` nevertheless asserts
        # 0 ≤ c ≤ 1, which keeps the collocation nodes c_drift2 = (c, c) inside the time step.
        # That is a restriction on the *method*, not one the conditions impose; it is recorded
        # here so that the discrepancy is deliberate rather than forgotten.
        @testset "the constructor restricts c beyond what the conditions require" begin
            @test_throws AssertionError TableauModifiedStochasticStoermerVerlet(-0.5)
            @test_throws AssertionError TableauModifiedStochasticStoermerVerlet(1.5)
        end

        # At c = 1/2 it reduces to the plain stochastic Störmer-Verlet of the paper.
        @testset "ModifiedStoermerVerlet(1/2) reduces to StoermerVerlet" begin
            m = Coefficients(TableauModifiedStochasticStoermerVerlet(0.5))
            v = Coefficients(TableauStochasticStoermerVerlet())
            @test m.a ≈ v.a
            @test m.ā ≈ v.ā
            @test m.â ≈ v.â
            @test m.b ≈ v.b
            @test m.b̄ ≈ v.b̄
            @test m.b̂ ≈ v.b̂
            @test m.α ≈ v.α
            @test m.α̂ ≈ v.α̂
            @test m.β ≈ v.β
            @test m.β̂ ≈ v.β̂
        end
    end
end

main()
