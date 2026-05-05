<!-- PDF page 2; model gpt-4.1 -->

190                                                                                   STEPHEN GODFREY AND NATHAN ISGUR                                                        32

masses and wave functions, we then took another step which, in previous treatments, has at best been only partly done: as a test of our results we performed a very extensive analysis of the couplings (strong, electromagnetic, and weak) of our model states. In fact, as far as we are aware, our calculations (which embrace not only all known mesons but also many predicted ones) represent the most wide-ranging and complete set of such calculations ever done. Given the success of our model in understanding the properties of most known states, it is our hope and belief that such a unified treatment of all mesons and their couplings can provide a useful guide to experimenters in their searches for new states. Finally, in addition to the new elements of our model and our extensive analysis of meson couplings, we also make a number of phenomenological observations (on, e.g., the $1\,^3D_1$–$2\,^3S_1$ ambiguity and the scalar-meson problem in light-meson systems) which are new. Many of these observations, while made within the context of our model, have a more general validity.

While the primary impetus for this work was, as indicated, to understand mesons, we had some secondary motivations. One of these was to provide a reasonably reliable model of the meson “background” against which one hopes to see some of the more exotic hadrons (pure glue states and hybrids) expected in QCD. Another was to use mesons as a testing ground for ideas on the relativization of the quark model before applying those ideas to the richer and experimentally better known baryons.

---

**II. SOFT QCD AND THE MESON HAMILTONIAN**

Soft QCD, as we define it here, is based on the hypothesis that hadrons may be approximately described in terms of rest-frame valence-quark configurations, the dynamics of which are governed by a Hamiltonian with one-gluon exchange dominant at short distances and with confinement implemented by a flavor-independent Lorentz-scalar interaction.\textsuperscript{4}

We take as our basic equation the (not manifestly covariant but relativistic) rest-frame Schrödinger-type equation

\[
H \,|\, \Psi \rangle = (H_0 + V) \,|\, \Psi \rangle = E \,|\, \Psi \rangle
\tag{1a}
\]

where

\[
H_0 = (p^2 + m_1^2)^{1/2} + (p^2 + m_2^2)^{1/2},
\tag{1b}
\]

$V = V(\mathbf{p}, \mathbf{r})$ is a momentum-dependent potential, $\mathbf{p} = \mathbf{p}_1 = -\mathbf{p}_2$ is a center-of-mass momentum, and where $\mathbf{r}$ becomes the usual spatial coordinate in the nonrelativistic limit. The derivation of this equation and of the potential $V(\mathbf{p}, \mathbf{r})$ is given below and in Appendix A, but before proceeding we comment briefly on its status. In a Fock-space representation, appropriate to a field-theoretic description of bound states, it is always possible to use the Schrödinger equation $H\Psi = E\Psi$, where $H$ is the Hamiltonian of the field theory and $\Psi$ a superposition of the states of the theory. (Strictly speaking, this equation is only well defined in the infinite-momentum frame, but this technicality is easily circumvented.) In this form the effects of, for example, transverse-gluon exchange on the $q\bar q$ component of the total wave function are felt in terms of mixing-matrix elements to $q\bar qg$ states. By integrating out the effects of all higher Fock components in the wave function, one can from this starting point always arrive at an equation of the form of our equation (1). Our key assumptions are that (i) with QCD cut off at some small scale $\mu$ of the order of the appropriate constituent quark mass, the $q\bar q$ wave functions described by (1) will dominate the total Fock-space wave functions so that their normalizations can be taken to be approximately unity and (ii) $V(\mathbf{p}, \mathbf{r})$ is a variant of the usual one-gluon-exchange-plus-linear-confinement potential with modifications reflecting various expected relativistic effects to be discussed below. While we will return to the general case momentarily, for orientation we first note that in the nonrelativistic limit this equation becomes the familiar nonrelativistic Schrödinger equation with

\[
H_0 \rightarrow \sum_{i=1}^2 \left( m_i + \frac{p^2}{2m_i} \right)
\tag{2a}
\]
and
\[
V_{ij}(\mathbf{p}, \mathbf{r}) \rightarrow H_{ij}^{\text{conf}} + H_{ij}^{\text{hyp}} + H_{ij}^{\text{so}} + H_A
\tag{2b}
\]
where
\[
H_{ij}^{\text{conf}} = -\left[ \frac{3}{4}c + \frac{3}{4}br - \frac{\alpha_s(r)}{r} \right] \mathbf{F}_i \cdot \mathbf{F}_j
\tag{3}
\]
includes the spin-independent linear confinement and Coulomb-type interactions,

\[
H_{ij}^{\text{hyp}} = -\frac{\alpha_s(r)}{m_i m_j} \left[ \frac{8\pi}{3} \mathbf{S}_i \cdot \mathbf{S}_j \delta^3(\mathbf{r}) + \frac{1}{r^3} \left( \frac{3 \mathbf{S}_i \cdot \mathbf{r} \, \mathbf{S}_j \cdot \mathbf{r}}{r^2} - \mathbf{S}_i \cdot \mathbf{S}_j \right) \right] \mathbf{F}_i \cdot \mathbf{F}_j
\tag{4}
\]
is the color hyperfine interaction, and

\[
H_{ij}^{\text{so}} = H_{ij}^{\text{so(cm)}} + H_{ij}^{\text{so(tp)}}
\tag{5}
\]
is the spin-orbit interaction with

\[
H_{ij}^{\text{so(cm)}} = -\frac{\alpha_s(r)}{r^3} \left[ \left( \frac{1}{m_i} + \frac{1}{m_j} \right) \left( \frac{\mathbf{S}_i}{m_i} + \frac{\mathbf{S}_j}{m_j} \right) \right] \cdot \mathbf{L} (\mathbf{F}_i \cdot \mathbf{F}_j),
\tag{6}
\]
its color-magnetic piece and with

\[
H_{ij}^{\text{so(tp)}} = -\frac{1}{2r} \frac{\partial H_{ij}^{\text{conf}}}{\partial r} \left[ \frac{\mathbf{S}_i}{m_i^2} + \frac{\mathbf{S}_j}{m_j^2} \right] \cdot \mathbf{L}
\tag{7}
\]
being the Thomas-precession term. In these formulas $\mathbf{L} = \mathbf{r} \times \mathbf{p}$,

\[
\mathbf{F}_i =
\begin{cases}
\frac{\lambda_i}{2} & \text{for quarks}, \\
\frac{\lambda_i^c}{2} = -\frac{\lambda_i^*}{2} & \text{for antiquarks},
\end{cases}
\tag{8}
\]
