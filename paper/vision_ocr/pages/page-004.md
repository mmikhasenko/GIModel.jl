<!-- PDF page 4; model gpt-4.1 -->

192

FIG. 2. The saturating $\alpha_s(Q^2)$ [Eq. (12), solid curve] compared to lowest-order QCD with $\Lambda=200$ MeV [Eq. (11), dashed curve]; to allow for thresholds we let $N_f$ in (11) be the number of flavors with $4m_f^2 < Q^2$ but demanded that $\alpha_s$ be continuous ($\Lambda$ refers to the $N_f=2$ regime); the fitted function is $\alpha_s(Q^2) = 0.25 \exp(-Q^2) + 0.15 \exp(-Q^2 / 10) + 0.20 \exp(-Q^2 / 1000)$, with $Q$ in GeV.

We have solved for mesons with the Hamiltonian (1) in three stages. In the first two of these stages we treat the Hamiltonian

\[
\tilde{H}_1 = (p^2 + m_1^2)^{1/2} + (p^2 + m_2^2)^{1/2} + \tilde{H}_{12}^{\text{conf}}
+ \tilde{H}_{12}^{\text{hyp}} + \tilde{H}_{12}^{\text{so}}
\tag{14}
\]

($\tilde{H}$ denotes an operator that has been modified by the relativistic effects described above and detailed in Appendix A) by directly diagonalizing in a large harmonic-oscillator sectors. This diagonalization is first performed in $| jm; ls \rangle$ sectors where ${\bf L} = {\bf r} \times {\bf p}$, ${\bf S} = {\bf S}_1 + {\bf S}_2$, and ${\bf J} = {\bf L} + {\bf S}$. The off-diagonal effects of $\tilde{H}_{12}^{\text{tensor}}$ [the tensor part of (4) which can cause $^3L_J \leftrightarrow {}^3L'_J$ mixing] and of $\tilde{H}_{[12]}^{\text{so}}$ (the antisymmetric piece of the spin-orbit interaction which arises only if the quark masses are unequal, in which circumstance it can cause $^3L_J \leftrightarrow {}^1L_J$ mixing) are then treated perturbatively by diagonalizing the mass matrix in the basis of eigenvectors of the $| jm; ls \rangle$ sectors. At both stages the basis used is expanded until we find convergence.

For most states the solution of our Hamiltonian problem is complete at this point, but for self-conjugate isoscalar mesons we must also consider the effects$^5$ of $H_A$.

\[
A(^{2S+1}L_J)_{ji} = 4\pi (2L+1) \left\{
A(^{2S+1}L_J) \left[ \frac{\alpha_s(M_j^2) \alpha_s(M_i^2)}{\pi^2} \right]^{n/2} S_L(\Psi_j) S_L(\Psi_i)
\over m_i m_j
\right\},
\]

where $A(^{2S+1}L_J)_{ji}$ depends on the unperturbed annihilation channel masses $M_j$ and $M_i$, $n$ is as above, and where $S_L(\Psi_i)$ is a smearing of the $q_i \bar{q}_i$ wave function at the origin:


---

In mesons, single-gluon annihilation is forbidden by color conservation, but annihilation via multiple gluons is expected. For heavy quarks where the annihilation is controlled by a small $\alpha_s$, this process will (at least in the absence of anomalies) be dominated by the minimum number of gluons allowed: two for even and three for odd charge-conjugation states. On general grounds we expect this effect to lead to a contribution to the mass matrix with diagonal entries of the form

\[
A_{Q\bar{Q} \to Q\bar{Q}} \approx \alpha_s^n \frac{|\Psi_{Q\bar{Q}}(0)|^2}{M_Q^2},
\tag{15}
\]

where $n=2$ or 3 as $C=+$ or $-$. Since $\Psi_{Q\bar{Q}}(0) \simeq 0$ if $L>0$, we may further expect this effect to be very small in heavy-quark systems unless $L=0$ [it will not be exactly zero both because the annihilation actually occurs over a region of size $m_Q^{-1}$ and because there will be relativistic smearing of the quarks over a region of size $m_Q^{-1}$. Even in $S$ waves, however, this effect should be quite small in the triplet states as can be seen by comparing to $H^{\text{hyp}}$ and noting that $\alpha_s^3 \ll \alpha_s$ even for $c\bar{c}$. Thus in heavy quark systems we can anticipate that the only place where annihilation might be noticeable is in the states $n\,{}^1S_0$.

In light-quark systems we must, on the other hand, expect $H_A$ to play a more important role. In the absence of a calculation of the annihilation amplitudes we must then treat $H_A$ phenomenologically and consequently the predictive power of our model is reduced for such light-isoscalar mesons. This weakness is somewhat alleviated by two factors: (1) Even in the light mesons $H_A$ is usually small, and there is considerable phenomenological evidence to reinforce one's expectation that it becomes weaker as a meson system becomes more excited. Thus in practice $H_A$ can often simply be ignored. (2) All of the self-conjugate isoscalar mesons in a given $^{2S+1}L_J$ sector can, if it is necessary to consider annihilation at all, be described by the introduction of a single new annihilation parameter $A(^{2S+1}L_J)$, and, with the exception of the pseudoscalar mesons which we will discuss extensively below, this description is insensitive to uncertainties in how the effects of $H_A$ should be implemented.

With the exception of the pseudoscalar mesons, our prescription for gluon annihilation mixing is adapted directly from Eq. (15) above with relativistic modifications motivated by the observations of Appendix A: for the annihilation amplitude from $q_i \bar{q}_i \to q_j \bar{q}_j$ in the channel $^{2S+1}L_J$ we take

\[
A(^{2S+1}L_J)_{ji} = 4\pi (2L + 1) \left\{
A (^{2S+1}L_J) \left[ \frac{\alpha_s(M_j^2)\alpha_s(M_i^2)}{\pi^2} \right]^{n/2} \frac{S_L(\Psi_j) S_L(\Psi_i)}{m_i m_j}
\right\},
\tag{16}
\]

\[
S_L(\Psi_i) \equiv \frac{1}{(2\pi)^{3/2}}
\int d^3p \frac{1}{\sqrt{4\pi}} \Phi_i(p)
\left[ \frac{p}{E_i} \right]^L \frac{m_i}{E_i}.
\tag{17}
\]

Here $\phi_i(p) = \Phi_i(p) Y_{LM}(\theta_p, \phi_p)$ is the full normalized
