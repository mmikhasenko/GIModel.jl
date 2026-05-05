<!-- PDF page 36; model gpt-4.1 -->

224                                                                                                        32

STEPHEN GODFREY AND NATHAN ISGUR

\[
G(Q^2) = \frac{-4\alpha_s(Q^2)}{3} \frac{4\pi}{Q^2}
\]

and a long-range $1\otimes 1$ linear confining interaction $S(Q^2)$ suggested by lattice QCD calculations.$^{41}$ We include the effects of asymptotic freedom by using the Ansatz (12) for $\alpha_s(Q^2)$.

By this method we find for $|\mathbf{p}| = |\mathbf{p}'| \equiv p$ that

\[
\chi_s^{\dagger} \chi_{\bar{s}}^{\dagger} V_{\text{eff}}(\mathbf{P}, \mathbf{r}) \chi_s \chi_{\bar{s}}
\]
\[
= \frac{1}{(2\pi)^3} \int d^3Q \, e^{i\mathbf{Q} \cdot \mathbf{r}} \overline{U} (\mathbf{p}', s') \overline{V}(-\mathbf{p}, \bar{s}) I(Q^2)
\times U(\mathbf{p}, s) V(-\mathbf{p}', \bar{s}'), 
\tag{A1}
\]
where

\[
I(Q^2) = G(Q^2)(\gamma^{\mu})_q (\gamma_{\mu})_{\bar{q}} - S(Q^2)(1)_q (1)_{\bar{q}}
\tag{A2}
\]

with $\mathbf{Q} = \mathbf{p}' - \mathbf{p}$ and $\mathbf{P} = (\mathbf{p} + \mathbf{p}')/2$. With

\[
f(\mathbf{P}, \mathbf{r}) \equiv \frac{1}{(2\pi)^3} \int d^3 Q \, e^{i \mathbf{Q}\cdot\mathbf{r}} f(\mathbf{P}, \mathbf{Q})
\tag{A3}
\]

we find

\[
V_{\text{eff}}(\mathbf{P}, \mathbf{r}) = G_{\text{eff}}(\mathbf{P}, \mathbf{r}) + S_{\text{eff}}(\mathbf{P}, \mathbf{r}),
\tag{A4}
\]

where

\[
G_{\text{eff}}(\mathbf{P}, \mathbf{r}) = \frac{1}{(2\pi)^3} \int d^3 Q \, e^{i \mathbf{Q} \cdot \mathbf{r}} G(Q^2)
\left\{
\left[
1 - \frac{Q^2}{4E(E+m)} + \frac{i\mathbf{Q} \times \mathbf{P} \cdot \mathbf{S}_q}{E(E+m)}
\right]
\right.
\]
\[
\left.
+ \left[
\frac{\mathbf{P} - 2i\mathbf{Q} \times \mathbf{S}_q}{2E}
\right]
\left[
1 - \frac{Q^2}{4\overline{E}(\overline{E} + \overline{m})} + \frac{i\mathbf{Q} \times \mathbf{P} \cdot \mathbf{S}_q}{\overline{E}(\overline{E} + \overline{m})}
\right]
\right\}
\tag{A5}
\]

and

\[
S_{\text{eff}}(\mathbf{P},\mathbf{r}) = \frac{1}{(2\pi)^3} \int d^3 Q \, e^{i \mathbf{Q} \cdot \mathbf{r}} S(Q^2)
\left[
\frac{m}{E} + \frac{Q^2}{4E(E+m)} - \frac{i\mathbf{Q} \times \mathbf{P} \cdot \mathbf{S}_q}{E(E+m)}
\right.
\]
\[
\left.
\times
\left(
\frac{\overline{m}}{\overline{E}} + \frac{Q^2}{4\overline{E}(\overline{E}+\overline{m})} - \frac{i\mathbf{Q} \times \mathbf{P} \cdot \mathbf{S}_q}{\overline{E}(\overline{E}+\overline{m})}
\right)
\right]
\tag{A6}
\]

where $m$ and $\overline{m}$ are the quark and antiquark masses, $\mathbf{S}_q$ and $\mathbf{S}_{\bar{q}}$ are their spins, and $E = (p^2 + m^2)^{1/2}$, $\overline{E} = (p^2 + \overline{m}^2)^{1/2}$.

For on-shell $q\bar{q}$ scattering at c.m. momentum $p$ these results are exact, i.e., these potentials will reproduce exactly the scattering amplitudes we have assumed. There are several reasons, however, for not using these potentials directly in (1). The first is that as they stand they do not adequately reflect the full expected momentum dependence of the potential $V(\mathbf{P}, \mathbf{r})$ in (1) since it will in general have off-energy-shell behavior not considered in (A1). In field theory the Schrödinger equation (1) arises in the $q\bar{q}$ sector of Fock space by integrating over more complex components of Fock space such as $| q\bar{q}g \rangle$, and this integration will introduce additional $\mathbf{P}$ dependence in $V$ not seen in (A1). Related to this deficiency is the fact that potentials like (A5) and (A6) have the usual ambiguity in the ordering of the classical quantities $E$ and $\overline{E}$ into a quantum operator: a matrix element involving $E$ will in general depend on $\mathbf{p}$ and $\mathbf{p}'$. Another reason for not using the above potentials without modification has specifically to do with (A6). It seems to us very unlikely that the confinement potential is a simple $1 \times 1$ interaction: the picture that emerges from studies of lattice QCD indicates that it is spin-independent, but that it arises from a distortion of a Coulomb interaction. If in fact we use one-dimensional QED as a guide, then we would expect the confinement potential to not only be spin-independent, but also $p$-independent.

We respond to this situation by treating (A5) and (A6) as a framework on which to build a semiquantitative model of relativistic effects. We roughly classify these effects into three categories: (a) the strengths of the various interactions will depend on the c.m. momentum of the interacting quarks, (b) the interactions will, since they depend on both $\mathbf{P}$ and $\mathbf{Q}$, be nonlocal, and (c) the interactions will, through $\mathbf{Q}$ dependence, take on new $\mathbf{r}$ dependences. Based on (b) and (c) we introduce a smearing function for a meson $q_i \bar{q}_j$

\[
\rho_{ij}(\mathbf{r}-\mathbf{r}') = \frac{\sigma_{ij}^3}{\pi^{3/2}} e^{-\sigma_{ij}^2(\mathbf{r}-\mathbf{r}')^2}
\tag{A7}
\]

which we apply to our basic potentials $G(r)$ and $S(r)$ to obtain smeared potentials $\widetilde{G}(r)$ and $\widetilde{S}(r)$ via

\[
\widetilde{f}_{ij}(r) \equiv \int d^3 r' \, \rho_{ij}(\mathbf{r}-\mathbf{r}') f(r')
\tag{A8}
\]
with the prescription

\[
\sigma_{ij}^2 = \sigma_0^2
\left[
\frac{1}{2} + \frac{1}{2}
\left(
\frac{4m_i m_j}{(m_i + m_j)^2}
\right)^4
\right]
+s^2
\left[
\frac{2m_i m_j}{m_i + m_j}
\right]^2
\tag{A9}
\]
