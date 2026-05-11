---
title: "Godfrey-Isgur 1985 Vision OCR"
format:
  html: default
  pdf:
    filters:
      - tag-display-math.lua
---

# Godfrey-Isgur 1985 Vision OCR

Generated from rendered page/column images with a vision model. Use the PDF and saved crops as the authority for final equation audits.

## PDF Page 1

PHYSICAL REVIEW D                     VOLUME 32, NUMBER 1                          1 JULY 1985

Mesons in a relativized quark model with chromodynamics

Stephen Godfrey and Nathan Isgur  
Department of Physics, University of Toronto, Toronto, MSS 1A7 Canada  
(Received 12 December 1983; revised manuscript received 10 May 1985)

We show that mesons—from the $\pi$ to the $\Upsilon$—can be described in a unified quark model with chromodynamics. The key ingredient of the model is a universal one-gluon-exchange-plus-linear-confinement potential motivated by QCD, but it is crucial to the success of the description to take into account relativistic effects. The spectroscopic results of the model are supported by an extensive analysis of strong, electromagnetic, and weak meson couplings.

***

# I. INTRODUCTION

The discovery and exploration of the charmonium system and the parallel development of quantum chromodynamics (QCD) have revolutionized hadron physics. It is becoming clear that heavy-quark systems can be well described by nonrelativistic potential models and that many of their properties reflect the dynamics expected from QCD.$^1$

However, until confinement is better understood it is unlikely that we will be able to rigorously compare quarkonia with the predictions of QCD. This is illustrated in Table I which shows the fraction of the $2S$-$1S$ splitting which arises from confinement in a typical "Coulomb-plus-linear" fit to heavy quarkonia: clearly the properties of the confinement potential will continue to play an important role in the foreseeable future. Therefore, to study these systems we must for the present rely on models which are a mixture of "true" QCD and phenomenological treatments of confinement which are motivated by QCD: we call such models "soft QCD" to remind ourselves both of their rigor and their region of applicability.

Despite this appellation, and notwithstanding the many possible criticisms, such models have been successfully applied to the charmonium $(c\bar{c})$ and more recently $b$-quarkonium $(b\bar{b})$ families,$^2$ so that the value of such a picture is now widely accepted. Moreover, the success of soft QCD in these sectors at least raises the question of where, as the quark mass is decreased, such models become useless. Even though it seems certain that they will become inaccurate if small masses are involved, one might hope that the dynamics of light quarks can be at least qualitatively understood on the basis of these same models.$^3$

In this paper we present the results of a study of light and heavy mesons in soft QCD. *We have found that all mesons—from the pion to the upsilon—can be described in a unified framework.* We substantiate this conclusion by first calculating meson spectra and then performing an extensive analysis of meson couplings. Section II describes the model, while Secs. III and IV deal with spectroscopy and decays, respectively. In Sec. V we discuss our results. Some conclusions and comments are given in the final section.

Since most of the elements of our model have appeared in one form or another elsewhere, some general comments on its relationship to earlier work in this area seem to us to be mandatory; we will make more specific comments in the appropriate sections below. Almost all quark potential models are based on some variant of the Coulomb-plus-linear potential expected from QCD and ours turned out to be no exception (we tried and rejected several alternatives). Many models have also included some form of the running coupling constant of QCD, but we know of no other work in which the effects of $\alpha_s(Q^2)$ have been treated in such a consistent and complete manner as is the case here. Relativistic effects have also often been discussed, but their treatment has normally been a somewhat patchwork affair. Here we have attempted to identify all possible types of relativistic effects, including smearing, nonlocality, and momentum-dependent effective potentials, and to then treat them for all mesons in a unified and physically motivated way. We are not satisfied with our relativization of the quark model, but we believe it to be a step forward. Aside from such fundamental differences in the framework of our model (we have mentioned the most important new features, but there are others), we believe the work presented here is also distinguishable from earlier work by its breadth of application. We have not only compared the results of our unified model to all known mesons simultaneously, but we have also made predictions for hundreds of the as-yet-unseen low-lying excitations of the various meson flavor sectors. As explained in the text, our calculations were not only extensive, but they were also accurately done: we did not rely on dubious perturbative treatments of various terms in the Hamiltonian. Once in possession of predictions for meson

***

TABLE I. The importance of confinement in $Q\bar{Q}$.

$$
\begin{array}{cccc}
m_Q\ (\text{GeV}) & \text{Typical}\ \alpha_s & a_0 = \left(\frac{2}{3} \alpha_s m_Q \right)^{-1}\ (\text{fm}) & \text{Approximate \% of 2S-1S from confinement} \\
\hline
1.5 & 0.34 & 0.58 & 50 \\
5 & 0.21 & 0.28 & 35 \\
25 & 0.16 & 0.07 & 20 \\
50 & 0.14 & 0.04 & 15 \\
100 & 0.12 & 0.02 & 10 \\
\end{array}
$$

## PDF Page 2

190                                                                                   STEPHEN GODFREY AND NATHAN ISGUR                                                        32

masses and wave functions, we then took another step which, in previous treatments, has at best been only partly done: as a test of our results we performed a very extensive analysis of the couplings (strong, electromagnetic, and weak) of our model states. In fact, as far as we are aware, our calculations (which embrace not only all known mesons but also many predicted ones) represent the most wide-ranging and complete set of such calculations ever done. Given the success of our model in understanding the properties of most known states, it is our hope and belief that such a unified treatment of all mesons and their couplings can provide a useful guide to experimenters in their searches for new states. Finally, in addition to the new elements of our model and our extensive analysis of meson couplings, we also make a number of phenomenological observations (on, e.g., the $1\,^3D_1$–$2\,^3S_1$ ambiguity and the scalar-meson problem in light-meson systems) which are new. Many of these observations, while made within the context of our model, have a more general validity.

While the primary impetus for this work was, as indicated, to understand mesons, we had some secondary motivations. One of these was to provide a reasonably reliable model of the meson “background” against which one hopes to see some of the more exotic hadrons (pure glue states and hybrids) expected in QCD. Another was to use mesons as a testing ground for ideas on the relativization of the quark model before applying those ideas to the richer and experimentally better known baryons.

***

**II. SOFT QCD AND THE MESON HAMILTONIAN**

Soft QCD, as we define it here, is based on the hypothesis that hadrons may be approximately described in terms of rest-frame valence-quark configurations, the dynamics of which are governed by a Hamiltonian with one-gluon exchange dominant at short distances and with confinement implemented by a flavor-independent Lorentz-scalar interaction.\textsuperscript{4}

We take as our basic equation the (not manifestly covariant but relativistic) rest-frame Schrödinger-type equation

$$
H \,|\, \Psi \rangle = (H_0 + V) \,|\, \Psi \rangle = E \,|\, \Psi \rangle
\tag{1a}
$$

where

$$
H_0 = (p^2 + m_1^2)^{1/2} + (p^2 + m_2^2)^{1/2},
\tag{1b}
$$

$V = V(\mathbf{p}, \mathbf{r})$ is a momentum-dependent potential, $\mathbf{p} = \mathbf{p}_1 = -\mathbf{p}_2$ is a center-of-mass momentum, and where $\mathbf{r}$ becomes the usual spatial coordinate in the nonrelativistic limit. The derivation of this equation and of the potential $V(\mathbf{p}, \mathbf{r})$ is given below and in Appendix A, but before proceeding we comment briefly on its status. In a Fock-space representation, appropriate to a field-theoretic description of bound states, it is always possible to use the Schrödinger equation $H\Psi = E\Psi$, where $H$ is the Hamiltonian of the field theory and $\Psi$ a superposition of the states of the theory. (Strictly speaking, this equation is only well defined in the infinite-momentum frame, but this technicality is easily circumvented.) In this form the effects of, for example, transverse-gluon exchange on the $q\bar q$ component of the total wave function are felt in terms of mixing-matrix elements to $q\bar qg$ states. By integrating out the effects of all higher Fock components in the wave function, one can from this starting point always arrive at an equation of the form of our equation (1). Our key assumptions are that (i) with QCD cut off at some small scale $\mu$ of the order of the appropriate constituent quark mass, the $q\bar q$ wave functions described by (1) will dominate the total Fock-space wave functions so that their normalizations can be taken to be approximately unity and (ii) $V(\mathbf{p}, \mathbf{r})$ is a variant of the usual one-gluon-exchange-plus-linear-confinement potential with modifications reflecting various expected relativistic effects to be discussed below. While we will return to the general case momentarily, for orientation we first note that in the nonrelativistic limit this equation becomes the familiar nonrelativistic Schrödinger equation with

$$
H_0 \rightarrow \sum_{i=1}^2 \left( m_i + \frac{p^2}{2m_i} \right)
\tag{2a}
$$
and
$$
V_{ij}(\mathbf{p}, \mathbf{r}) \rightarrow H_{ij}^{\text{conf}} + H_{ij}^{\text{hyp}} + H_{ij}^{\text{so}} + H_A
\tag{2b}
$$
where
$$
H_{ij}^{\text{conf}} = -\left[ \frac{3}{4}c + \frac{3}{4}br - \frac{\alpha_s(r)}{r} \right] \mathbf{F}_i \cdot \mathbf{F}_j
\tag{3}
$$
includes the spin-independent linear confinement and Coulomb-type interactions,

$$
H_{ij}^{\text{hyp}} = -\frac{\alpha_s(r)}{m_i m_j} \left[ \frac{8\pi}{3} \mathbf{S}_i \cdot \mathbf{S}_j \delta^3(\mathbf{r}) + \frac{1}{r^3} \left( \frac{3 \mathbf{S}_i \cdot \mathbf{r} \, \mathbf{S}_j \cdot \mathbf{r}}{r^2} - \mathbf{S}_i \cdot \mathbf{S}_j \right) \right] \mathbf{F}_i \cdot \mathbf{F}_j
\tag{4}
$$
is the color hyperfine interaction, and

$$
H_{ij}^{\text{so}} = H_{ij}^{\text{so(cm)}} + H_{ij}^{\text{so(tp)}}
\tag{5}
$$
is the spin-orbit interaction with

$$
H_{ij}^{\text{so(cm)}} = -\frac{\alpha_s(r)}{r^3} \left[ \left( \frac{1}{m_i} + \frac{1}{m_j} \right) \left( \frac{\mathbf{S}_i}{m_i} + \frac{\mathbf{S}_j}{m_j} \right) \right] \cdot \mathbf{L} (\mathbf{F}_i \cdot \mathbf{F}_j),
\tag{6}
$$
its color-magnetic piece and with

$$
H_{ij}^{\text{so(tp)}} = -\frac{1}{2r} \frac{\partial H_{ij}^{\text{conf}}}{\partial r} \left[ \frac{\mathbf{S}_i}{m_i^2} + \frac{\mathbf{S}_j}{m_j^2} \right] \cdot \mathbf{L}
\tag{7}
$$
being the Thomas-precession term. In these formulas $\mathbf{L} = \mathbf{r} \times \mathbf{p}$,

$$
\mathbf{F}_i =
\begin{cases}
\frac{\lambda_i}{2} & \text{for quarks}, \\
\frac{\lambda_i^c}{2} = -\frac{\lambda_i^*}{2} & \text{for antiquarks},
\end{cases}
\tag{8}
$$

## PDF Page 3

and $\alpha_s(r)$ is the running coupling constant of QCD which we will discuss below in more detail. Since

$$
\langle \mathbf{F}_i \cdot \mathbf{F}_j \rangle =
\begin{cases}
  -\frac{4}{3} & \text{in a meson} \,, \\
  -\frac{2}{3} & \text{in a baryon} \,,
\end{cases}
\tag{9}
$$

the dynamics of these two confined systems are very closely related by soft QCD. Multiquark systems are also related by the same basic dynamics, although with many possible internal color states their characteristics are not expected to be very closely related to those of mesons and baryons. The Hamiltonian (3) also has the property that it only allows the creation of color singlets as isolated hadrons.

Finally, the term $H_A$ is the annihilation interaction of Fig. 1, which must be taken into account in any sector where $q\bar{q}$ annihilation via gluons can occur.$^5$ While $H_A$ is in principle calculable, at the present time it must normally be parametrized. Fortunately, it can only contribute in isoscalar channels so in many cases its effects can be avoided. We discuss this situation more extensively below.

The Hamiltonian (2), though derived from the well behaved equation (1), is actually inconsistent as it stands: the spin-dependent terms $H^{\mathrm{hyp}}$ and $H^{so}$ are more singular than $r^{-2}$ and are therefore illegal operators in the Schrödinger equation. The resolution of this paradox requires that we return to the more general case of Eq. (1) as discussed in Appendix A. It is shown there that the relativistic potential $V(\mathbf{p}, \mathbf{r})$ differs from its nonrelativistic limit in two qualitatively important ways: (i) the coordinate $\mathbf{r}$ (which in the nonrelativistic limit is the relative coordinate $\mathbf{r}_{12} = \mathbf{r}_1 - \mathbf{r}_2$) becomes smeared out over distances of the order of the inverse quark masses and (ii) the coefficients of the various potentials [which in the nonrelativistic limit have the strengths shown in Eqs. (3) to (7)] become dependent on the momentum of the interacting quarks. The smearing of the potentials has the consequence of taming all of their singularities, making them legal operators in (2) and, more directly relevant for our purposes, in (1) (which demands potentials less singular than $r^{-1}$). The details of our implementation of this smearing, which we accomplish via a smearing function

$$
\rho_{ij} (\mathbf{r}' - \mathbf{r}) = \frac{\sigma_{ij}^3}{\pi^{3/2}} e^{-\sigma_{ij}^2 (\mathbf{r}' - \mathbf{r})^2} ,
\tag{10}
$$

are relegated to Appendix A. The momentum dependence of the potentials, as well as some technical issues related to performing calculations with such potentials, is also discussed in detail in Appendix A, but as an illustration of such a dependence we note that in a relativistic treatment, in general factors of $m^{-1}$ can become factors like $(p^2 + m^2)^{-1/2}$. Thus, for example, the hyperfine interaction of a light quark should not blow up like $m^{-1}$ but rather should have (as in the bag model) a finite limit as $m \to 0$ determined by $\langle p^{-1} \rangle$, which is in turn controlled by the radius of confinement. Of course such modifications play a significant role only in light-quark systems.

While both of these requirements are semiquantitatively defined by the considerations of Appendix A, the method we have chosen for implementing them is very coarse. Each type of interaction would in principle, for example, have a distinct smearing function as well as more complicated energy-dependent factors than those we assume. We should also stress that these two types of effects are intimately connected: they are together describing a momentum dependence of our potentials which we cannot readily impose on the usual spatial Schrödinger equation. Despite these shortcomings, we believe that our method correctly portrays the main characteristics of these relativistic effects.

We now turn to a discussion of the running coupling constant $\alpha_s$. With $N_f$ quark flavors with masses much less than $Q^2$, in lowest-order QCD

$$
\alpha_s(Q^2) = \frac{12\pi}{(33 - 2N_f) \ln(Q^2 / \Lambda^2)} .
\tag{11}
$$

For $100 \leq \Lambda \leq 300$ MeV and for $3 \leq N_f \leq 5$, $\alpha_s$ is always around $0.2$ when $Q = 10$ GeV and it varies very slowly from $Q = 5$ to 20 GeV. On the other hand, as $Q \to \Lambda$ the perturbative formula (11) diverges, a behavior commonly taken to be a signal of confinement. Since we are interested in this soft regime, we cannot avoid this divergence; rather we assume that $\alpha_s$ saturates at some value $\alpha_s^{\mathrm{critical}}$ for low $Q^2$ as confinement emerges. We parametrize this qualitative behavior in the convenient form

$$
\alpha_s(Q^2) = \sum_k a_k e^{-Q^2 / 4 \gamma_k^2} ,
\tag{12}
$$

where

$$
\alpha_s^{\mathrm{critical}} \equiv \sum_k a_k
$$

is a free parameter, but where the remaining parameters are used to fit (12) to a QCD curve for $\alpha_s(Q^2)$ as shown in Fig. 2. The form (12) is convenient both because it is easily transformed into

$$
\alpha_s(r) = \sum_k a_k \frac{2}{\sqrt{\pi}} \int_0^{\gamma_k r} e^{-x^2} dx
\tag{13}
$$

and because the resulting color-charge distribution is easily convoluted with the relativistic smearing (10). These details are also discussed in Appendix A.

***

**FIG. 1.** The origin of the annihilation term $H_A$: a typical graph.

## PDF Page 4

192

FIG. 2. The saturating $\alpha_s(Q^2)$ [Eq. (12), solid curve] compared to lowest-order QCD with $\Lambda=200$ MeV [Eq. (11), dashed curve]; to allow for thresholds we let $N_f$ in (11) be the number of flavors with $4m_f^2 < Q^2$ but demanded that $\alpha_s$ be continuous ($\Lambda$ refers to the $N_f=2$ regime); the fitted function is $\alpha_s(Q^2) = 0.25 \exp(-Q^2) + 0.15 \exp(-Q^2 / 10) + 0.20 \exp(-Q^2 / 1000)$, with $Q$ in GeV.

We have solved for mesons with the Hamiltonian (1) in three stages. In the first two of these stages we treat the Hamiltonian

$$
\tilde{H}_1 = (p^2 + m_1^2)^{1/2} + (p^2 + m_2^2)^{1/2} + \tilde{H}_{12}^{\text{conf}}
+ \tilde{H}_{12}^{\text{hyp}} + \tilde{H}_{12}^{\text{so}}
\tag{14}
$$

($\tilde{H}$ denotes an operator that has been modified by the relativistic effects described above and detailed in Appendix A) by directly diagonalizing in a large harmonic-oscillator sectors. This diagonalization is first performed in $| jm; ls \rangle$ sectors where $\mathbf{L} = \mathbf{r} \times \mathbf{p}$, $\mathbf{S} = \mathbf{S}_1 + \mathbf{S}_2$, and $\mathbf{J} = \mathbf{L} + \mathbf{S}$. The off-diagonal effects of $\tilde{H}_{12}^{\text{tensor}}$ [the tensor part of (4) which can cause $^3L_J \leftrightarrow {}^3L'_J$ mixing] and of $\tilde{H}_{[12]}^{\text{so}}$ (the antisymmetric piece of the spin-orbit interaction which arises only if the quark masses are unequal, in which circumstance it can cause $^3L_J \leftrightarrow {}^1L_J$ mixing) are then treated perturbatively by diagonalizing the mass matrix in the basis of eigenvectors of the $| jm; ls \rangle$ sectors. At both stages the basis used is expanded until we find convergence.

For most states the solution of our Hamiltonian problem is complete at this point, but for self-conjugate isoscalar mesons we must also consider the effects$^5$ of $H_A$.

In mesons, single-gluon annihilation is forbidden by color conservation, but annihilation via multiple gluons is expected. For heavy quarks where the annihilation is controlled by a small $\alpha_s$, this process will (at least in the absence of anomalies) be dominated by the minimum number of gluons allowed: two for even and three for odd charge-conjugation states. On general grounds we expect this effect to lead to a contribution to the mass matrix with diagonal entries of the form

$$
A_{Q\bar{Q} \to Q\bar{Q}} \approx \alpha_s^n \frac{|\Psi_{Q\bar{Q}}(0)|^2}{M_Q^2},
\tag{15}
$$

where $n=2$ or 3 as $C=+$ or $-$. Since $\Psi_{Q\bar{Q}}(0) \simeq 0$ if $L>0$, we may further expect this effect to be very small in heavy-quark systems unless $L=0$ [it will not be exactly zero both because the annihilation actually occurs over a region of size $m_Q^{-1}$ and because there will be relativistic smearing of the quarks over a region of size $m_Q^{-1}$. Even in $S$ waves, however, this effect should be quite small in the triplet states as can be seen by comparing to $H^{\text{hyp}}$ and noting that $\alpha_s^3 \ll \alpha_s$ even for $c\bar{c}$. Thus in heavy quark systems we can anticipate that the only place where annihilation might be noticeable is in the states $n\,{}^1S_0$.

In light-quark systems we must, on the other hand, expect $H_A$ to play a more important role. In the absence of a calculation of the annihilation amplitudes we must then treat $H_A$ phenomenologically and consequently the predictive power of our model is reduced for such light-isoscalar mesons. This weakness is somewhat alleviated by two factors: (1) Even in the light mesons $H_A$ is usually small, and there is considerable phenomenological evidence to reinforce one's expectation that it becomes weaker as a meson system becomes more excited. Thus in practice $H_A$ can often simply be ignored. (2) All of the self-conjugate isoscalar mesons in a given $^{2S+1}L_J$ sector can, if it is necessary to consider annihilation at all, be described by the introduction of a single new annihilation parameter $A(^{2S+1}L_J)$, and, with the exception of the pseudoscalar mesons which we will discuss extensively below, this description is insensitive to uncertainties in how the effects of $H_A$ should be implemented.

With the exception of the pseudoscalar mesons, our prescription for gluon annihilation mixing is adapted directly from Eq. (15) above with relativistic modifications motivated by the observations of Appendix A: for the annihilation amplitude from $q_i \bar{q}_i \to q_j \bar{q}_j$ in the channel $^{2S+1}L_J$ we take

$$
A(^{2S+1}L_J)_{ji} = 4\pi (2L + 1) \left\{
A (^{2S+1}L_J) \left[ \frac{\alpha_s(M_j^2)\alpha_s(M_i^2)}{\pi^2} \right]^{n/2} \frac{S_L(\Psi_j) S_L(\Psi_i)}{m_i m_j}
\right\},
\tag{16}
$$

where $A(^{2S+1}L_J)_{ji}$ depends on the unperturbed annihilation channel masses $M_j$ and $M_i$, $n$ is as above, and where $S_L(\Psi_i)$ is a smearing of the $q_i \bar{q}_i$ wave function at the origin:


$$
S_L(\Psi_i) \equiv \frac{1}{(2\pi)^{3/2}}
\int d^3p \frac{1}{\sqrt{4\pi}} \Phi_i(p)
\left[ \frac{p}{E_i} \right]^L \frac{m_i}{E_i}.
\tag{17}
$$

Here $\phi_i(p) = \Phi_i(p) Y_{LM}(\theta_p, \phi_p)$ is the full normalized

## PDF Page 5

32 MESONS IN A RELATIVISTIC QUARK MODEL . . . 193

Fourier transform of the wave function $\Psi_i(r)$ and
$$
E_i = (m_i^2 + p^2)^{1/2}.
$$

In the pseudoscalar mesons, (16) and (17) fail. We believe that the behavior of this channel is related to the U(1) problem of QCD in which the mass of the ninth pseudoscalar meson is lifted from zero in the chiral limit by *nonperturbative* annihilation amplitudes.$^{6}$ These nonperturbative effects actually have the opposite sign to the perturbative two-gluon annihilation amplitude, implying that the pseudoscalar amplitude $A(^1S_0)$ must have a complicated dependence on the annihilation channel mass $M$. We can consequently offer no compelling description of annihilation in this channel, but we have found two examples of possible behavior for $A(^1S_0)$ with interesting phenomenological consequences. The simplest possibility (P1) is that there is a large positive nonperturbative annihilation amplitude in this channel which dies away exponentially with a scale of $m_\eta^2$ and which is to be added directly to the perturbative piece (which is, as already implied above, known in this case). Thus in P1, in place of the bracketed factor in (16), we take

$$
\begin{equation}
\left[ A_{np} e^{-(m_i^2 + m_j^2)/m_\eta^2} \right.
\left. + \frac{2\pi}{3} (\ln 2 - 1) \left( \frac{\alpha_s(M_j^2) \alpha_s(M_i^2)}{\pi^2} \right) \right]. \tag{18a}
\end{equation}
$$

This simple possibility has many attractive features, but, as we shall see below, it would be ruled out if the experimental indication for an isoscalar pseudoscalar meson at around 1275 MeV is confirmed. We consequently consider a second more exotic possibility (P2). At large annihilation-channel invariant mass $M^2$, the perturbative calculation must be correct so we know that in this region
$$
A(^1S_0)_{ii} \propto \left[ \frac{\alpha_s(M^2)}{\pi} \right]^2
$$
is small and negative. As $M^2$ is decreased it is possible that it becomes (nonperturbatively) large and negative before becoming large and positive at $M^2 \rightarrow 0$ as it is constrained to do.$^6$ In this picture the simple exponential in (18a) is replaced by a nonperturbative contribution which rapidly changes sign at some $M_0^2$: the bracketed factor in (16) becomes

$$
\begin{equation}
\bigg[ A_{np} \Big\{ 1 - \left( \frac{M}{M_0} \right)^4 \Big\} e^{-(m_i^2 + m_j^2)/M_0^2 - M^4/4M_0^4}
+ \frac{2\pi}{3} (\ln 2 - 1) \left( \frac{\alpha_s(M^2)}{\pi} \right)^2 \bigg]. \tag{18b}
\end{equation}
$$

The reader will note that, apart from satisfying the general features required for P2, we allowed ourselves great freedom in parametrizing the dependence of $A(^1S_0)$ on the annihilation-channel invariant mass $M$. With this elaborately defined (but perhaps not implausible) model we are able to accommodate an isoscalar in this sector near or below the $\pi'$.

The results of applying (16) and (18) are given in Sec. III. As already mentioned, our results are very stable under any reasonable modification of (16). On the other hand, we find ourselves unable to draw definite conclusions about pseudoscalar mixing. The consequences of the models (18) will, however, be discussed in Sec. V A as two possible scenarios for these states.

This completes the specification of our model and leaves us ready to determine its parameters. A search in this parameter space converged more or less uniquely to the physically reasonable values quoted in Table II which form the basis for the results of this paper. An exception to the uniqueness of these parameters is that we found that our results were relatively insensitive (after readjustment of our other parameters) to the addition to our quark masses of an overall constant that kept $\frac{1}{2}(m_u + m_d)$ in the range from zero to about 250 MeV. We initially hoped to take advantage of this freedom to assign our quarks current quark masses. However, this proved to be inconsistent: once isospin splittings were taken into account; $^7$ indeed, observed isospin splittings along with the constraint $m_d - m_u \sim 5$ MeV led us to the typical relativized constituent quark mass$^8$ $\frac{1}{2}(m_u + m_d) = 220$ MeV. In retrospect this seems inevitable: to maintain the dominance of the $q\bar{q}$ sector of Fock space for light quarks, we clearly must choose a cutoff $\mu$ which produces an extended (constituent) quark with "gluonic mass" of order $\Lambda$.

We close this section by commenting on the level of accuracy we expect from the model we have just described. A source of error common to all mesons is our restriction to the simplest sectors of Fock space. This means, in particular, that we are not considering the mass shifts and mixings that will arise from the interaction between these mesons and their decay channels (both open and closed).

***
TABLE II. The parameters of soft-QCD spectroscopy.

|                |                        |
|----------------|------------------------|
| **Masses$^a$** |                        |
| $\frac{1}{2}(m_u + m_d)$ | 220 MeV                |
| $m_s$          | 419 MeV                |
| $m_c$          | 1628 MeV               |
| $m_b$          | 4977 MeV               |
| **Potentials** |                        |
| $b$            | 0.18 GeV$^2$           |
| $\alpha_s^{\text{critical}}$ | 0.60                 |
| $\Lambda$      | 200 MeV                |
| $c$            | $-253$ MeV             |

Relativistic effects  
(see Appendix A)  
Smearing: $\sigma_0 = 1.80$ GeV  
$s = 1.55$  
$m \leftrightarrow E$ ambiguity:  
$\epsilon_c = -0.168$  
$\epsilon_f = +0.025$  
$\epsilon_{\text{sol}(V)} = -0.035$  
$\epsilon_{\text{sol}(S)} = +0.055$  

$^a$Note that we ignore isospin violation here and throughout this paper unless otherwise indicated.

## PDF Page 6

Since such mass shifts are of the order of magnitude of 10 MeV, this provides a limit to our expected accuracy. It should be noted that such shifts must be expected even in $\Upsilon$ spectroscopy: in that case 

$$
\Delta M \sim \Gamma_s^2 \langle E_{\text{threshold}} - M \rangle^{-1},
$$

where $\Gamma_s$ is a typical strong-interaction width and $E_{\text{threshold}} - M$ is the distance of the state in question from $B\bar{B}$ threshold. This uncertainty, which we should also apply to level splittings, will in fact (as it must according to decoupling theorems) disappear as $M_Q \to \infty$, but as can be seen from the above formula this decoupling will not occur with high accuracy until $M_Q$ is of the order of 100 GeV. The other principal source of error in these calculations is associated with the relativistic corrections. As we have stressed, our implementation of these effects is rather schematic. Taken together we therefore believe that we can expect only 25-MeV average accuracy for mesons containing light quarks and 10-MeV average accuracy for heavy-quark systems. We measure the quality of the results which follow against these expectations.

***

### III. THE MESON SPECTRUM OF SOFT QCD

The meson spectra predicted by the dynamics of the previous section are shown in Figs. 3 to 11. The compositions of most states are given in the figures or figure captions; for isoscalars see Table III. We show separately the isovector mesons ($-ud$, $2^{-1/2}(u\bar{u} - d\bar{d}), d\bar{u}$) (Fig. 3), the strange mesons ($-u\bar{s}$, $-d\bar{s}$, $-s\bar{d}$, $+s\bar{u}$) (Fig. 4), the light-isoscalar mesons (dominantly mixtures of $u\bar{u}$, $d\bar{d}$, and $s\bar{s}$) (Fig. 5), the $J/\psi$ family ($c\bar{c}$) (Fig. 6), and the $\Upsilon$ family ($b\bar{b}$) (Fig. 8). In Fig. 7 we show the charmed ($-c\bar{d}, c\bar{u}$) and charm-strange ($c\bar{s}$) mesons, and in Fig. 9 the various states containing a single $b$ quark.

***

![FIG. 3]

FIG. 3. The isovector mesons $[-ud, \sqrt{1/2}(u\bar{u}-d\bar{d}), d\bar{u}]$. The dominant spectral composition and predicted masses of states in GeV are shown near solid bars representing their masses. Shaded areas correspond to the experimental masses and their uncertainties, normally taken from the Particle Data Group (1984). The comparison of the $1^{--}$ and $0^{++}$ sectors with experiment requires special consideration: see Secs. V A and V D, respectively. Significant spectroscopic mixing in this sector: $1^{--}(1.45) \simeq 1.00(2\,^3S_1) + 0.04(1\,^3D_1)$.

## PDF Page 7

FIG. 4. The strange mesons ( $-u\bar{s}, -d\bar{s}$ ). The legend is as for Fig. 3. Significant spectroscopic mixing in this sector:

(a) With 
$$
\left[ \begin{array}{c} Q_{\text{low}} \\ Q_{\text{high}} \end{array} \right]
\simeq
\left[ \begin{array}{cc}
\cos \theta_{nL} & \sin \theta_{nL} \\
-\sin \theta_{nL} & \cos \theta_{nL}
\end{array} \right]
\left[ \begin{array}{c}
n\,^1 LL \\
n\,^3 LL
\end{array} \right]
$$
we find $\theta_{1P} \simeq 34^\circ$, $\theta_{1D} \simeq 33^\circ$, $\theta_{2P} \simeq 15^\circ$, $\theta_{1F} \simeq 32^\circ$, $\theta_{2D} \simeq 25^\circ$, $\theta_{1G} \simeq 33^\circ$;

(b)
$$
1^{--}(1.58) \simeq 1.00 (2\,^3S_1) + 0.04 (1\,^3D_1).
$$

($-b\bar{d}, b\bar{u}, b\bar{s}, b\bar{c}$). Finally, in Figs. 10 and 11 we show the spectrum of a charge $+\frac{2}{3}$ $t$ quark of mass 35 GeV in $t\bar{t}$ and in states containing a single $t$ quark.

It would be premature to make a detailed comparison between the model and experiment until the decay analysis of the next section has been presented, but some general observations are already possible if we anticipate the support of this analysis. The most important of these are the following.

(1) The gross spectrum of mesons is determined by a single universal potential whose main features are quark confinement at large distances and a Coulomb-type attraction at short distances. A good fit to all spectra requires, in the context of this model, a strong coupling constant which evolves along the lines expected from QCD. The universal potential is shown in Fig. 12 with the rms values of the interquark separations of various representative mesons shown to illustrate the region over which the potential can be considered to be tested.$^{10}$

(2) The existence of such a universal potential is revealed only if the relativistic effects of Sec. II are properly taken into account.

(3) The expected spin-dependent interactions are present with strengths which are correctly correlated to those of the spin-independent potentials, and the expected mass and spatial dependences of these spin-dependent interactions is borne out.

(4) Finally, the central message of these results is that there are no qualitative changes in the behavior of meson systems as the quark masses decrease. For example, features present in the $c\bar{c}$ system persist in the isovectors, the principle difference between the two spectra being that relativistic corrections (including hyperfine and fine splittings) have become more prominent in the latter.

## PDF Page 8

![Figure 5 - isoscalar meson spectrum]

FIG. 5. The isoscalar mesons (mainly $u\bar{u}$, $d\bar{d}$, $s\bar{s}$). The legend is as for Fig. 3. The assignment is given in Table III. The comparison of the $0^{-+}$ isoscalars with experiment requires special consideration: see Sec. V A. For the $\eta'$ meson see Ref. 9.

![Figure 6 - charmonium spectrum]

FIG. 6. The charmonia ($c\bar{c}$). The legend is as for Fig. 3. Significant spectroscopic mixing in this sector: 
$$
1^{--}(3.82) \simeq 1.00(1\,^3D_1) + 0.01(1\,^3S_1) - 0.03(2\,^3S_1) - 0.01(3\,^3S_1).
$$

## PDF Page 9

**FIG. 7.** The charmed mesons ( $c\bar{d}, c\bar{u}, c\bar{s}$ ). The legend is as for Fig. 3. Significant spectroscopic mixing in these sectors:

With

$$
\begin{bmatrix}
Q^{cq}_{\text{low}} \\[1ex]
Q^{cq}_{\text{high}}
\end{bmatrix}
\simeq
\begin{bmatrix}
\cos{\theta^{cq}_{NL}} & \sin{\theta^{cq}_{NL}} \\
-\sin{\theta^{cq}_{NL}} & \cos{\theta^{cq}_{NL}}
\end{bmatrix}
\begin{bmatrix}
n\, ^1 L_L \\[1ex]
n\, ^3 L_L
\end{bmatrix}
$$

we find $\theta^{c\bar{u}}_{1P} \simeq -41^\circ$, $\theta^{c\bar{u}}_{1D} \simeq -39^\circ$, $\theta^{c\bar{s}}_{1P} \simeq -44^\circ$, $\theta^{c\bar{s}}_{1D} \simeq -39^\circ$.

The $2^-$ states are not shown in the figure, but they are all within 20 MeV of their respective $3^-$ partners.

**FIG. 8.** The $b$-quarkonia ($b\bar{b}$). The legend is as for Fig. 3.

## PDF Page 10

![Left: Mass spectrum for b-flavored mesons, right: t-quarkonium hypothetical spectrum.]

**FIG. 9.** The $b$-flavored mesons. The legend is as for Fig. 3. Significant spectroscopic mixing in this sector:

With

$$
\begin{pmatrix}
Q_{\text{low}}^{b\bar{q}} \\
Q_{\text{high}}^{b\bar{q}}
\end{pmatrix}
\simeq
\begin{pmatrix}
\cos \theta_{nL}^{b\bar{q}} & \sin \theta_{nL}^{b\bar{q}} \\
-\sin \theta_{nL}^{b\bar{q}} & \cos \theta_{nL}^{b\bar{q}}
\end{pmatrix}
\begin{pmatrix}
n\,{}^1L_L \\
n\,{}^3L_L
\end{pmatrix}
$$

we find $\theta_{1P}^{b\bar{u}} \simeq -43^\circ$, $\theta_{1P}^{b\bar{s}} \simeq -45^\circ$, $\theta_{1P}^{b\bar{c}} \simeq -53^\circ$.

These $1^+$ states are not shown in the figure, but they are all within 40 MeV of their respective $2^+$ partners.

**FIG. 10.** Some hypothetical $t$-quarkonium ($t\bar{t}$) with $m_t = 35$ GeV.

(b) $b\bar{s}$  (c) $b\bar{c}$

Significant spectroscopic mixing in this sector:

$\theta_{1P}^{b\bar{u}} \simeq -43^\circ$, $\theta_{1P}^{b\bar{s}} \simeq -45^\circ$, $\theta_{1P}^{b\bar{c}} \simeq -53^\circ$.

within 40 MeV of their respective $2^+$ partners.


IV. AN ANALYSIS OF MESON COUPLINGS

A. Introduction

A successful model of hadrons must address not only hadronic spectra, but also the internal structures of hadronic systems. Of course a spectrum and its associated quantum numbers (spin, parity, isospin, etc.) depend on this structure, but the most sensitive measures of the internal compositions of hadrons are their couplings to other hadrons and to electromagnetic and weak currents. In this section we accordingly discuss the results of an extensive analysis of meson couplings based on the states found in the previous section.

B. Decays via pseudoscalar-meson emission

Our approach to describing meson decay via the emission of a pseudoscalar meson is almost identical to a recent analysis of baryon decays in the soft QCD model$^{11}$ in terms of the usual quark-pion-emission model.$^{12}$ In our variation on this classic approach, the decay of a meson is assumed to proceed through a single-quark transition as depicted in Fig. 13. In mesons there is an element to these calculations which does not appear in baryons: there is a danger of double counting as a consequence of approximating the emission of the pseudoscalar meson $P$ by a pointlike emission by a single quark. This is illustrated in Fig. 14 which demonstrates that one should not sum the amplitude for emission from both the quark and antiquark as one would for photon emission. It is therefore

## PDF Page 11

Labels under diagrams:  

**FIG. 11.** Some hypothetical t-flavored mesons with $m_t = 35$ GeV.


TABLE III. The approximate composition of some mixed isoscalars. These are amplitude decompositions in terms of the eigenstates in the absence of annihilation. The annihilation parameters are $A(^{3}S_1) = +2.5$, $A(^{3}P_2) = -0.8$, and in the case of the pseudoscalars in model P1, $A_{np} = +0.5$, while in model P2 we used $A_{np} = +0.55$ and $M_0 = 1.17$ GeV. Other states, in the absence of compelling evidence to the contrary, have been assumed for now to be ideally mixed. Also shown are the predicted and observed splittings of these states from their isovector companion states. We have denoted $(1/\sqrt{2})(u\bar{u}+d\bar{d})$ by ns; $n$ in the column labeling is the radial quantum number. Note that since (18b) is mass dependent, poles of the inverse propagator are not orthogonal.

| State          | (Name)      | Model | ns    | ss    | cc     | bb       | ns    | ss    | cc     | $\Delta m^{\mathrm{theory}}$ (MeV) | $\Delta m^{\mathrm{expt}}$ (MeV) |
|----------------|-------------|-------|-------|-------|--------|----------|-------|-------|--------|----------|----------------------|--------------------|
| $1\,^1S_0$     | $\eta(548)$         | P1    | +0.67 | -0.73 | +0.001  | +2$\times 10^{-4}$ | +0.11 | +0.042 | $-5\times10^{-4}$ | +370                | +410               |
|                |                        | P2    | +0.68 | -0.73 | -0.005  | +3$\times 10^{-4}$ | +0.09 | +0.051 | +0.002          | +340                |                    |
| $1\,^1S_0$     | $\eta'(958)$         | P1    | +0.58 | +0.62 | +0.004  | +5$\times 10^{-4}$ | +0.47 | +0.13 | -0.002          | +810                | +820               |
|                |                        | P2    | +0.48 | +0.78 | -0.008  | +9$\times 10^{-4}$ | +0.34 | +0.14 | +0.004          | +780                |                    |
| $1\,^1S_0$     | $\eta_c(2980)$     | P1    | -0.008 | -0.005 | +1.000 | +2$\times 10^{-4}$ | +0.004 | +0.003 | -0.002          |                     |                    |
|                |                        | P2    | -0.002 | -0.001 | +1.000 | +2$\times 10^{-4}$ | +0.001 | +0.001 | -0.002          |                     |                    |
| $2\,^1S_0$     | $\eta_{?}(?)$             | P1    | -0.26 | -0.17 | -0.003  | -3$\times 10^{-4}$ | +0.79 | -0.44 | +0.001          | see Sec. VA         |                    |
|                |                        | P2    | +0.07 | +0.08 | -0.005  | -2$\times 10^{-4}$ | +0.99 | +0.07 | +0.002          | see Sec. VA         |                    |
| $2\,^1S_0$     | $\eta'_{?}(?)$            | P1    | -0.17 | -0.10 | -0.003  | -3$\times 10^{-4}$ | +0.26 | +0.86 | +0.001          | see Sec. VA         |                    |
|                |                        | P2    | +0.09 | +0.08 | -0.006  | -1$\times 10^{-4}$ | -0.16 | +0.97 | +0.003          | see Sec. VA         |                    |
| $1\,^3S_1$     | $\omega(783)$             |         | +0.999| -0.02 | -0.001  | -1$\times 10^{-4}$|        |        |                 | +10                 | +15 $\pm$ 10       |
|                | $\phi(1019)$             |         | +0.02 | +0.999| -0.0006 | -7$\times 10^{-5}$|        |        |                 | +250                | +250 $\pm$ 10      |
|                | $\psi(3097)$             |         | +0.0009| +0.0006| +1.000| -3$\times 10^{-5}$|        |        |                 | -40                 | -45 $\pm$ 15       |
| $1\,^3P_2$     | $f_2(1270)$             |         | +0.997| +0.06 | +0.007  | +9$\times 10^{-4}$|        |        |                 | +215                | +200 $\pm$ 20      |
|                | $f_2'(1515)$             |         | -0.07 | +0.997| +0.004  | +4$\times 10^{-4}$|        |        |                 |                     |                    |

## PDF Page 12

![FIG. 12 illustration]

**FIG. 12.** The universal $q\bar{q}$ potential $H^{\text{conf}}(r)$ in a color-singlet meson; also shown is the rms $q\bar{q}$ separation in some representative mesons calculated by analogy to $r_{E}^{2}$ defined in Table VII.


![FIG. 13 illustration]

**FIG. 13.** A single quark transition $q \rightarrow qP$.

always necessary to bear in mind that amplitudes like that of Fig. 13 are being used as approximations to a pair-creation amplitude and that the correct amplitude is obtained by summing over distinct pair creation plus rearrangement processes. This awkward feature is a symptom of the fact that the spectator model we are using cannot be strictly valid; we comment further on this point below.

To proceed with the model in the simplest way we add two other assumptions (both of which are supported by our results): (1) pair creation of $u, d,$ and $s$ quarks is approximately SU(3) symmetric and (2) violations of Zweig's rule proceed mainly via meson wave functions. By this latter assumption we mean that in a decay like $K^* \rightarrow K\eta$, only diagrams like Fig. 14 are taken into account, it being assumed that the effects of disconnected diagrams like Fig. 15(a) are mostly taken into account by using the Zweig-rule-violating $\eta$ wave function which arises from the annihilation mixing of Fig. 15(b) as discussed in Sec. II. With these approximations, the amplitudes for pseudoscalar emission from a quark (antiquark) take the form


$$
A_{q(\bar{q})}\langle M^*(\mathbf{k},s)\rightarrow M(\mathbf{k}',s')P^i(q)\rangle = \pm i \frac{(4\omega\omega')^{1/2}}{(2\pi)^{9/2}}\langle M(s') \mid (\mathbf{g}\sigma_{q(\bar{q})}\cdot \mathbf{q} \pm h\sigma_{q(\bar{q})}\cdot \mathbf{p}')e^{\mp i\mathbf{q}'\cdot \mathbf{r}} X^{i}_{q(\bar{q})}|M^*(s)\rangle,
\tag{19}
$$


where $\sigma_{q(\bar{q})}/2$ and $\mathbf{r}_{q(\bar{q})}$ are the spin and position of the quark (antiquark), $\mathbf{p}' = -i\mathbf{\nabla}$ is a gradient acting on the final-state wave function, and the upper (lower) sign refers to the $q\,(\bar{q})$ case. The $X^{i}_{q(\bar{q})}$ are flavor operators defined and detailed in Appendix B. The calculations are most readily performed by taking $\mathbf{q} = q\hat{\mathbf{z}}$ thereby calculating helicity amplitudes $H_{m}$ where $m = s' = s$, and then transforming to the usual partial-wave basis; details of this process are given in Appendix C.

We apply this crude decay model by mimicking completely the previous baryon analysis,$^{{11}}$ forsaking our full wave functions for the harmonic-oscillator wave functions of the SU(6) limit. This allows us to calculate the amplitudes analytically to reveal their basic simplicity and the intrinsic relations between them. We discuss the effect of using more realistic wave functions and the possibility of using a more realistic decay model below.

In our approach we find that all of the states in a given SU(6)$\times$O(3) multiplet often share (apart from individual strength factors) common partial-wave amplitudes in-





![FIG. 14 illustration]

**FIG. 14.** On the danger of double counting in meson decays to two pseudoscalar mesons.


![FIG. 15 illustration]

**FIG. 15.** (a) A disconnected diagram. (b) A diagram contributing to $\eta$-$\eta'$ mixing.

## PDF Page 13

dependent of their flavor, internal spins, or total angular momentum. These "universal" amplitudes are displayed in Table IV from which one can see that, as in baryons, the amplitudes fall into two classes. The first class, called "structure independent," consists of $A, A', A'', A_c$, and $A_0$, which have only the momentum dependence dictated by angular momentum considerations along with the gentle "elastic form factor" $e^{-q^2/16\beta^2}$. The second class of amplitudes, called "structure dependent," consists of $S$, $D$, $P$, and $S_c$ since they have additional polynomial momentum dependences which are highly sensitive to the structure of the states. We follow Ref. 11 at this point and forego attempting to calculate these various reduced amplitudes in terms of $g$ and $h$. Instead we allow a new constant for each such amplitude; in principle this means that our harmonic-oscillator decay amplitudes would be described by an expanded number of parameters instead of just two. However, we note from Table IV that we might expect $A \simeq A' \simeq A'' \simeq A_c \simeq A_0$ and $S \simeq D \simeq P \simeq S_c$ so in practice we tentatively employ such a two-parameter fit to the decay amplitudes. This fit can easily be relaxed as more accurate data on highly excited meson decays become available (since we calculate with $A$ and $S$ but explicitly display factors of $A'/A$, $D/S$, etc.); for now it provides an adequate guide. Recalling that our main objective here is to test the model for meson structure, this relaxation of what is obviously a very rudimentary decay model seems to us both sensible and prudent. With this generalization, the model has much in common with more algebraic approaches.$^{13}$ The values of the reduced partial-wave amplitudes (i.e., the amplitudes in square brackets in Table IV) which we used in our calculations


***

TABLE IV. The reduced partial-wave amplitudes. The full amplitudes have in addition a factor of $q^L e^{-q^2/16\beta^2}$; $\beta_c$ is defined in Table V.

| Amplitude | Representative decays | $L$ |
|:----------|:----------------------|:---:|
| $A=\left[g+\frac{1}{4} h\right]\beta$ | $1\,{}^3S_1 \to 1\,{}^1S_0 + P$ | 1 |
| | $1\,{}^3P_2 \to 1\,{}^1S_0 + P \quad 1\,{}^3S_1 + P$ | 2 |
| | $1\,{}^3P_1 \to 1\,{}^3S_1 + P$ | 2 |
| | $1\,{}^1P_1 \to 1\,{}^3S_1 + P$ | 2 |
| | $1\,{}^3D_3 \to 1\,{}^1S_0 + P \quad 1\,{}^3S_1 + P$ | 3 |
| $A'=\left[g-\frac{1}{4} h\right]\beta$ | $1\,{}^1P_1 \to 1\,{}^3P_0 + P$ | 1 |
| $A''=\left[g+\frac{1}{8} h\right]\beta$ | $1\,{}^3D_3 \to 1\,{}^1P_1 + P$ | 2 |
| $A_c=\left[g+\frac{1}{2}\left(\frac{m_d}{m_d+m_c}\right)h\right]\beta$ | as $A$, but for charmed-meson decays | |
| $A_0=\{g\}\beta$ | $1\,{}^3P_1 \to 1\,{}^3P_0 + P$ | 1 |
| $S=\left[3h-\frac{1}{2}\left(g+\frac{1}{4} h\right)\frac{q^2}{\beta^2}\right]\beta$ | $1\,{}^3P_1 \to 1\,{}^3S_1 + P$ | 0 |
| | $1\,{}^1P_1 \to 1\,{}^3S_1 + P$ | 0 |
| | $1\,{}^3P_0 \to 1\,{}^1S_0 + P$ | 0 |
| $D=\left[3h-\frac{3}{10}\left(g+\frac{1}{4} h\right)\frac{q^2}{\beta^2}\right]\beta$ | $1\,{}^3D_1 \to 1\,{}^1S_0 + P \quad 1\,{}^3S_1 + P$ | 1 |
| | $1\,{}^3D_3 \to 1\,{}^1S_0 + P \quad 1\,{}^3S_1 + P$ | 1 |
| $P=\left[3h-\frac{3}{4}\left(g+\frac{1}{4} h\right)\frac{q^2}{\beta^2}\right]\beta$ | $2\,{}^1S_0 \to 1\,{}^3S_1 + P$ | 1 |
| | $2\,{}^3S_1 \to 1\,{}^1S_0 + P \quad 1\,{}^3S_1 + P$ | 1 |
| $S_c=\left[3h-\frac{m_c A_c}{(m_d+m_c)\beta}\frac{q^2}{\beta_c^2}\right]\beta_c$ | as $S$, but for charmed-meson decays | |

## PDF Page 14

TABLE V. The strong decay amplitudes for $M^* \to M + P$. The amplitudes tabulated here are for the process $M^* \to M + P(q)$ where $P(q)$ is the pseudoscalar meson shown recoiling with momentum $q$ from the decay of $M^*$ at rest; the decay rate is the square of the listed numerical amplitude. The amplitude formula shown has been abbreviated by suppressing a factor of $\langle q/2\pi \rangle^{1/2} \exp(-q^2/16\beta^2)$ as well as a factor of $+i$ in all $P$ and $F$ wave amplitudes. Isospin coefficients are to be calculated in the MP order. Note also that a factor of $1/\sqrt{2}$ has already been included for amplitudes involving identical particles. The two decays $\rho \to \pi\pi$ and $B \to [\omega\pi]_S$ are used as input to the fit. We have allowed for isoscalar mixing in our numerical results, but results in the formula column are for ideal mixing [as in (B10) and (B11)] in every nonet except $1^1S_0$ where we show "perfect-mixing" formulas from (B14) and (B15). The isoscalar mixings are taken from Table III (using P1 for pseudoscalars) in every nonet except $1^1S_0$ where we simply use (B14) and (B15). A mixing angle not given explicitly in Table III is assumed, for now, to be zero. Note that the table gives amplitudes to $K^*\bar{K}$; the full rate to $K^*K$ is (approximately) twice that to $K^*\bar{K}$.

***

#### Harmonic-oscillator amplitude

*Columns (broken OCR layout): Decay; $(\bar{q} \equiv q / \beta)$; (MeV$^{1/2}$); $u$–$d$–$s$ mesons; “Realistic” factor; Experiment; References.*

$1^3S_1$  

$\rho \to \pi\pi$  
$\phi \to KK$  
$K^* \to K\pi$  

- $+ \left( \frac{4}{3} \right)^{1/2} A \bar{q}$
- $- \left( \frac{4}{3} \right)^{1/2} A \bar{q}$
- $+ A\bar{q}$

- $+12.4$
- $-2.4$
- $+7.9$

- 12.4 ; fit
- 1.9$\pm$0.1 ; 
- 7.1$\pm$0.1 ; 

$1^3P_2$

$A_2 \to (\pi\pi)\rho\pi$  
$A_2 \to \eta\pi$  
$A_2 \to KK$  
$A_2 \to \eta'\pi$  
$f \to \pi\pi$  
$f \to KK$  
$f' \to \eta\eta$  
$f' \to (K\pi)K^*\bar{K}$  
$K^* \to K^*\eta$  
$K^* \to K^*\eta'$  

- $- \left( \frac{1}{3} \right)^{1/2} A \bar{q}^2$
- $+ \left( \frac{1}{30} \right)^{1/2} A \bar{q}^2$
- $- \left( \frac{1}{30} \right)^{1/2} A \bar{q}^2$
- $+ \left( \frac{1}{30} \right)^{1/2} A \bar{q}^2$
- $- \left( \frac{1}{10} \right)^{1/2} A \bar{q}^2$
- $- \left( \frac{1}{30} \right)^{1/2} A \bar{q}^2$
- $+ \left( \frac{1}{120} \right)^{1/2} A \bar{q}^2$
- $0$
- $- \left( \frac{1}{15} \right)^{1/2} A \bar{q}^2$
- $+ \left( \frac{1}{60} \right)^{1/2} A \bar{q}^2$
- $- \left( \frac{1}{30} \right)^{1/2} A \bar{q}^2$
- $- \left( \frac{1}{10} \right)^{1/2} A \bar{q}^2$
- $+ \left( \frac{1}{10} \right)^{1/2} A \bar{q}^2$
- $+ \left( \frac{1}{20} \right)^{1/2} A \bar{q}^2$
- $- \left( \frac{3}{40} \right)^{1/2} A \bar{q}^2$
- $- \left( \frac{3}{40} \right)^{1/2} A \bar{q}^2$
- $+ \left( \frac{3}{40} \right)^{1/2} A \bar{q}^2$
- $- \dfrac{\sqrt{2}-1}{\sqrt{120}} A \bar{q}^2$
- $+ \dfrac{\sqrt{2}+1}{\sqrt{120}} A \bar{q}^2$

- $-5.6$
- $+4.5$
- $-2.8$
- $+1.0$
- $-11$
- $-2.5$
- $+0.8$
- $+1.1$
- $-7.1$
- $+2.9$
- $-0.04$
- $+1.7$
- $+7.7$
- $-3.7$
- $-2.1$
- $+1.2$
- $-0.8$
- below threshold

- (1.5)
- 8.8$\pm$0.3 ; f
- 4.0$\pm$0.1 ; 
- 2.3$\pm$0.1 ; 
- $<$1.5 ; 
- 12$\pm$1 ; 
- 2.3$\pm$0.2 ; 
- 1.0$\pm$0.2 ; g
- 0.8$\pm$0.4 ; 
- $\approx$8$\pm$1 ; 
- $<$6 ; 
- (1.4) $<$5 ; f
- 6.7$\pm$0.5 ; 
- (1.4) 5.0$\pm$0.5 ; f
- (1.4) 3.0$\pm$0.4 ; f
- (1.4) 2.0$\pm$0.4 ; 
- 2.2$^{+1.1}_{-2.2}$ ; 

$1^3P_1$ (nonstrange)

$A_1 \to [(\pi\pi)\rho\pi]_S$  
$A_1 \to [(\pi\pi)\rho\pi]_D$  
$A_1 \to (\pi\pi)\epsilon\pi$  

- $+\left( \frac{8}{9} \right)^{1/2} S$
- $- \left( \frac{1}{9} \right)^{1/2} A \bar{q}^2$
- $- \left( \frac{8}{9} \right)^{1/2} A_0 \bar{q}$

- $+20$
- $-3.0$
- $-4.3 \left( \frac{A_0}{A} \right)$

- (1.5)
- 18$\pm$2 ; f
- f
- 3.6$\pm$0.4 ; f

## PDF Page 15

TABLE V. (Continued).

| Decay | Harmonic-oscillator amplitude $(\vec{q} \equiv q/\beta)$ | (MeV$^{1/2}$) | "Realistic" factor | Experiment$^a$ (MeV$^{1/2}$) | References, footnotes |
|:------|:--------------------------------------------------------|:-------------:|:------------------:|:-----------------------------:|:----------------------:|
| $D \rightarrow [(K\pi)_{K^*} \overline{K}]_S$ | $+ ( \frac{2}{9})^{1/2} S$ | $+0.69$ |  |  | f |
| $D \rightarrow [(K\pi)_{K^*} \overline{K}]_D$ | $- (\frac{1}{36})^{1/2} A \vec{q}^{\,2}$ | $-0.0$ |  |  | f |
| $D \rightarrow (\eta\pi)_{\delta_2} \pi$ | $+ ( \frac{8}{3})^{1/2} A_0 \vec{q}$ | $+4.6 \left( \frac{A_0}{A} \right)$ | $4.0 \pm 0.8$ | e,f |
| $D \rightarrow (\pi\pi)_{\epsilon} \eta $ | $- ( \frac{4}{9})^{1/2} A_0 \vec{q}$ | $-3.4 \left( \frac{A_0}{A} \right)$ | $1.8 \pm 0.8$ | f |
| $E \rightarrow [(K\pi)_{K^*} \overline{K}]_S$ | $- ( \frac{4}{9})^{1/2} S$ | $-10$ | $7^{+4}_{-1}$ | f, h |
| $E \rightarrow [(K\pi)_{K^*} \overline{K}]_D$ | $+ ( \frac{1}{18})^{1/2} A \vec{q}^2$ | $+0.7$ | | f |
| $E \rightarrow \delta_2 \pi $ | $0$ | $\sim 0$ | | e,f |
| $E \rightarrow \epsilon \eta$ | $0$ | $\sim 0$ | | f |

$ ^{1}P_1 $ and $ ^{1}P_1 $ (strange)
| $Q_1 \rightarrow [(K\pi)_{K^*} \pi]_S$ | $- ( \frac{1}{6})^{1/2} S $ | $-0.3$ |  | $2.7 \pm 1.5$ | b,f,i,j |
| $Q_1 \rightarrow [(K\pi)_{K^*} \pi]_D$ | $- ( \frac{1}{12})^{1/2}A \vec{q}^2$ | $-3.0$ | (1.4) | $2.6 \pm 0.5$ | b,f,i |
| $Q_1 \rightarrow [(\pi\pi)_\rho K]_S $ | $+ (\frac{1}{6})^{1/2} S$ | $+9.9$ | (1.4) | $6.2 \pm 1.0$ | b,f,i |
| $Q_1 \rightarrow [(\pi\pi)_\rho K]_D $ | $+ ( \frac{1}{12})^{1/2} A \vec{q}^2$ | $+0.5$ | | not seen | b,f,i |
| $Q_1 \rightarrow [\omega K]_S$ | $- ( \frac{1}{18})^{1/2} S$ | $-7.0$ | | $3.1 \pm 0.5$ | b,i |
| $Q_1 \rightarrow [\omega K]_D $ | $- ( \frac{1}{36})^{1/2} A \vec{q}^2$ | $-0.2$ | | not seen | b,i |
| $Q_1 \rightarrow (K\pi)_\kappa \pi$ | $- ( \frac{1}{3})^{1/2} A' \vec{q}$ | $-3.1 \left( \frac{A'}{A} \right)$ | (1.3) | $5.1 \pm 0.5$ | f |
| $Q_1 \rightarrow (\pi\pi)_\epsilon K$ | $+ ( \frac{1}{9})^{1/2} A' \vec{q}$ | $+0.04 \left( \frac{A'}{A} \right)$ | | $1.6 \pm 0.4$ | f, j |

| $Q_2 \rightarrow [(K\pi)_{K^*} \pi]_S$ | $+ ( \frac{1}{3})^{1/2} S $ | $+16$ |  | $13 \pm 3$ | f |
| $Q_2 \rightarrow [(K\pi)_{K^*} \pi]_D$ | $- ( \frac{1}{24})^{1/2} A \vec{q}^2$ | $-0.1$ | (1.4) | $2.6 \pm 0.9$ | b,f,i,j |
| $Q_2 \rightarrow [(\pi\pi)_\rho K]_S$ | $+ ( \frac{1}{3})^{1/2} S$ | $+4.2$ | | $2.4 \pm 2.3$ | b,f,i |
| $Q_2 \rightarrow [(\pi\pi)_\rho K]_D$ | $- ( \frac{1}{24})^{1/2} A \vec{q}^2$ | $-1.8$ | (1.4) | not seen | b,f,i |
| $Q_2 \rightarrow [\omega K]_S$ | $- ( \frac{1}{9})^{1/2} S$ | $-3.0$ | | $1.1 \pm 1.0$ | b,i |
| $Q_2 \rightarrow [\omega K]_D $ | $+ ( \frac{1}{72})^{1/2}A \vec{q}^2$ | $+0.8$ | (1.4) | not seen | b,i |
| $Q_2 \rightarrow (K\pi)_\kappa \pi$ | $- (\frac{2}{3})^{1/2} A \vec{q}$ | $-1.4$ | | ~0 | b,f,i |
| $Q_2 \rightarrow (\pi\pi)_\epsilon K$ | $- (\frac{2}{9})^{1/2} A\vec{q}$ | $-2.0$ | | $1.7 \pm 1.4$ | b,f,i |

$^1P_1$ (nonstrange)
| $B \rightarrow [\omega\pi]_S $ | $- ( \frac{2}{9})^{1/2}S$ | $-11$ | 11 | fit |
| $B \rightarrow [\omega\pi]_D $ | $- ( \frac{1}{9})^{1/2} A \vec{q}^2$ | $-3.0$ | (1.5) | $3.2 \pm 0.5$ | e,f |
| $B \rightarrow (\eta\pi)_{\delta_2} \pi $ | $- (\frac{8}{9})^{1/2} A'\vec{q}$ | $-4.4 \left( \frac{A'}{A} \right)$ | | $18 \pm 2$ | f |

| $H \rightarrow [\rho\pi]_S $ | $+ ( \frac{2}{3})^{1/2} S$ | $+19$ | |  | f |
| $H \rightarrow [\rho\pi]_D $ | $+ ( \frac{1}{3})^{1/2} A \vec{q}^2$ | $+5.2$ | |  | f |
| $H' \rightarrow [\rho\pi]_S $ | $0$ | $\sim 0$ | |  | f |
| $H' \rightarrow [\rho\pi]_D $ | $0$ | $\sim 0$ | |  | f |
| $H' \rightarrow [(K\pi)_{K^*} \overline{K}]_S$ | $+ ( \frac{2}{9})^{1/2} S$ | $+8.0$ | |  | f |
| $H' \rightarrow [(K\pi)_{K^*} \overline{K}]_D$ | $+ ( \frac{1}{9})^{1/2} A \vec{q}^2$ | $+1.2$ | |  | f |

## PDF Page 16

***

TABLE V. (Continued).

| Decay              | Harmonic-oscillator amplitude |   (MeV$^{1/2}$)  |   "Realistic" factor  |   Experiment$^{a}$ (MeV$^{1/2}$) | References, footnotes |
|:-------------------------|:----------------------------:|:-------------------:|:----------------------:|:---------------------:|:---------------------:|
|                         $\mathbf{1^3P_0}$ |

$\delta_2 \rightarrow \eta\pi$       $+\left( \frac{1}{3} \right)^{1/2} S$    +$14$  (1.3) see Sec. VD e  
$\delta_2 \rightarrow KK$       $-\left( \frac{1}{3} \right)^{1/2} S$    $-11$  (1.3) see Sec. VD e  
$\delta_2 \rightarrow \eta'\pi$      $+\left( \frac{1}{3} \right)^{1/2} S$    below threshold (1.3) see Sec. VD e  
$\epsilon \rightarrow \pi\pi$       $-S$    $-27$  (1.3) see Sec. VD  
$\epsilon \rightarrow KK$       $-\left( \frac{1}{3} \right)^{1/2} S$    $-11$  (1.3) see Sec. VD  
$\epsilon \rightarrow \eta\eta$     $+\left( \frac{1}{12} \right)^{1/2} S$   below threshold (1.3) see Sec. VD  
$\epsilon' \rightarrow \pi\pi$      0    $\sim$0  (1.3) see Sec. VD  
$\epsilon' \rightarrow KK$       $-\left( \frac{2}{3} \right)^{1/2} S$    $-21$  (1.2) see Sec. VD  
$\epsilon' \rightarrow \eta\eta$    $+\left( \frac{1}{6} \right)^{1/2} S$    $+10$  (1.2) see Sec. VD  
$\kappa \rightarrow K\pi$       $+\left(\frac{1}{2}\right)^{1/2} S$      $+19$  (1.2) see Sec. VD  
$\kappa \rightarrow K\eta$       $-\left[ \frac{\sqrt{2}-1}{\sqrt{12}} \right] S$    $-2.8$  (1.2) see Sec. VD  
$\kappa \rightarrow K\eta'$      $+\left[ \frac{\sqrt{2}+1}{\sqrt{12}} \right] S$    below threshold (1.2) see Sec. VD  

  $0^{++}$ cryptoexotics see Sec. VD  
  see Sec. VD  
  see Sec. VD  
  see Sec. VD  
  see Sec. VD  

***

$\mathbf{1^3D_3}$

$\delta \rightarrow \eta\pi$       ---  
$\delta \rightarrow KK$         ---  
$S^* \rightarrow \pi\pi$        ---  
$S^* \rightarrow KK$         ---  

$g \rightarrow \pi\pi$       $-\left(\frac{1}{560}\right)^{1/2} A \tilde{q}^3$    $-5.6$  (1.8)  $6.9 \pm 0.6$  
$g \rightarrow \omega\pi$     $+\left(\frac{1}{420}\right)^{1/2} A \tilde{q}^3$   $+3.1$  (1.8)  $5 \pm 2$  
$g \rightarrow KK$        $-\left(\frac{1}{1120}\right)^{1/2} A \tilde{q}^3$  $-2.2$  $1.7 \pm 0.3$  
$g \rightarrow K^* \bar{K}$   $-\left(\frac{1}{840}\right)^{1/2} A \tilde{q}^3$  $-0.7$  (1.8)  $<2.6$  
$\omega \rightarrow \rho\pi$    $-\left(\frac{1}{140}\right)^{1/2} A \tilde{q}^3$  $-5.2$  (1.8)  $8 \pm 2$  
$\omega \rightarrow \omega\eta$  $+\left(\frac{1}{840}\right)^{1/2} A \tilde{q}^3$  $+0.9$  (1.8)    
$\omega \rightarrow [B\pi]_D$  $-\left(\frac{2}{5}\right)^{1/2} A'' \tilde{q}^2$  $-6.1$  $\left| \frac{A''}{A} \right|$ $7 \pm 2$  
$\omega \rightarrow [B\pi]_G$  $-\left(\frac{3}{7840}\right)^{1/2} A\tilde{q}^4$  $-0.2$    
$\omega \rightarrow KK$      $-\left(\frac{1}{1120}\right)^{1/2} A\tilde{q}^3$  $-2.1$    
$\omega \rightarrow K^* \bar{K}$ $-\left(\frac{1}{840}\right)^{1/2} A\tilde{q}^3$  $-0.6$  $7 \pm 2$  
$\phi \rightarrow \rho\pi$        0    0      
$\phi \rightarrow \omega\eta$      0    0       
$\phi \rightarrow [B\pi]_D$      0    0       
$\phi \rightarrow [B\pi]_G$      0    0       
$\phi \rightarrow K K$      $+\left(\frac{1}{560}\right)^{1/2} A\tilde{q}^3$  $+5.2$  
$\phi \rightarrow K^* \bar{K}$  $-\left(\frac{1}{420}\right)^{1/2} A\tilde{q}^3$  $-2.8$  (1.6)  $5 \pm 2$  
$K^* \rightarrow K\pi$     $-\left(\frac{3}{2240}\right)^{1/2} A\tilde{q}^3$  $-4.6$  $4.9 \pm 1.0$  
$K^* \rightarrow K^*\pi$    $+\left(\frac{1}{560}\right)^{1/2} A\tilde{q}^3$  $+2.8$  (1.6)  seen  
$K^* \rightarrow \rho K$    $-\left(\frac{1}{560}\right)^{1/2} A\tilde{q}^3$  $-2.3$  (1.6)  seen  
$K^* \rightarrow \omega K$   $+\left(\frac{1}{1680}\right)^{1/2} A\tilde{q}^3$  $+1.3$  (1.6)  seen

## PDF Page 17

TABLE V. (Continued).

| Decay                                   | Harmonic-oscillator amplitude $(\vec{q} = q / \beta)$                                         | (MeV$^{1/2}$) | "Realistic" factor | Experiment$^{a}$ (MeV$^{1/2}$) | References, footnotes |
|------------------------------------------|----------------------------------------------------------------------------------------|--------------|--------------------|-------------------------------|----------------------|
| $\rho \rightarrow [\omega \pi]_P$        | $-\left(\frac{1}{36}\right)^{1/2} D \bar{q}$                                            | $-7.8$       | $\left[ \frac{D}{S} \right]$ | (1.7)                        |                      |
| $\rho \rightarrow [\omega \pi]_F$        | $+\left(\frac{1}{600}\right)^{1/2} A \bar{q}^3$                                         | $+2.7$       |                      | (1.7)                        |                      |
| $\rho \rightarrow [\rho \eta]_P$         | $-\left(\frac{1}{72}\right)^{1/2} D \bar{q}$                                            | $-4.2$       | $\left[ \frac{D}{S} \right]$ | (1.7)                        |                      |
| $\rho \rightarrow [\rho \eta]_F$         | $+\left(\frac{1}{1200}\right)^{1/2} A \bar{q}^3$                                        | $+0.9$       |                      | (1.7)                        |                      |
| $\rho \rightarrow [K^* K]_P$             | $+\left(\frac{1}{72}\right)^{1/2} D \bar{q}$                                            | $+3.7$       | $\left[ \frac{D}{S} \right]$ | (1.7)                        |                      |
| $\rho \rightarrow [K^* \overline{K}]_F$  | $-\left(\frac{1}{1200}\right)^{1/2} A \bar{q}^3$                                        | $-0.7$       |                      | (1.7)                        |                      |
| $\omega \rightarrow [\rho \pi]_P$        | $+\left(\frac{1}{12}\right)^{1/2} D \bar{q}$                                            | $+14$        | $\left[ \frac{D}{S} \right]$ | (1.7)                        |                      |
| $\omega \rightarrow [\rho \pi]_F$        | $-\left(\frac{1}{200}\right)^{1/2} A \bar{q}^3$                                         | $-4.7$       |                      | (1.7)                        |                      |
| $\omega \rightarrow [\omega \eta]_P$     | $-\left(\frac{1}{72}\right)^{1/2} D \bar{q}$                                            | $-4.1$       | $\left[ \frac{D}{S} \right]$ | (1.7)                        |                      |
| $\omega \rightarrow [\omega \eta]_F$     | $+\left(\frac{1}{1200}\right)^{1/2} A \bar{q}^3$                                        | $+0.9$       |                      | (1.7)                        |                      |
| $\omega \rightarrow [K^* K]_P$           | $+\left(\frac{1}{72}\right)^{1/2} D \bar{q}$                                            | $+3.7$       | $\left[ \frac{D}{S} \right]$ | (1.7)                        |                      |
| $\omega \rightarrow [K^* K]_F$           | $-\left(\frac{1}{1200}\right)^{1/2} A \bar{q}^3$                                        | $-0.7$       |                      | (1.7)                        |                      |
| $\phi \rightarrow [\rho \pi]_P$          | $0$                                                                                    | $0$          |                      | (1.7)                        |                      |
| $\phi \rightarrow [\rho \pi]_F$          | $0$                                                                                    | $0$          |                      | (1.7)                        |                      |
| $\phi \rightarrow [\omega \eta]_P$       | $0$                                                                                    | $0$          |                      | (1.7)                        |                      |
| $\phi \rightarrow [\omega \eta]_F$       | $0$                                                                                    | $0$          |                      | (1.7)                        |                      |
| $\phi \rightarrow [K^* K]_P$             | $+\left(\frac{1}{36}\right)^{1/2} D \bar{q}$                                            | $+7.5$       | $\left[ \frac{D}{S} \right]$ | (1.6)                        |                      |
| $\phi \rightarrow [K^* K]_F$             | $-\left(\frac{1}{600}\right)^{1/2} A \bar{q}^3$                                        | $-2.4$       |                      | (1.6)                        |                      |

$1^3D_2$ and $1^1D_2$ (nonstrange)

| $Q_1 \rightarrow [\rho K]_P$             | $+\left(\frac{1}{72}\right)^{1/2} D \bar{q}$                                            | $+7.7$       | $\left[ \frac{D}{S} \right]$ | (1.6)                        | c                    |
| $Q_1 \rightarrow [\rho K]_F$             | $+\left(\frac{3}{1600}\right)^{1/2} A \bar{q}^3$                                        | $+0.9$       |                      | (1.6)                        | c                    |
| $Q_1 \rightarrow [\omega K]_P$           | $-\left(\frac{1}{216}\right)^{1/2} D \bar{q}$                                           | $-4.4$       | $\left[ \frac{D}{S} \right]$ | (1.6)                        | c                    |
| $Q_1 \rightarrow [\omega K]_F$           | $-\left(\frac{1}{1600}\right)^{1/2} A \bar{q}^3$                                       | $-0.5$       |                      | (1.6)                        | c                    |
| $Q_1 \rightarrow [K^* \pi]_P$            | $+\left(\frac{1}{72}\right)^{1/2} D \bar{q}$                                            | $+0.9$       | $\left[ \frac{D}{S} \right]$ | (1.6)                        | c                    |
| $Q_1 \rightarrow [K^* \pi]_F$            | $+\left(\frac{3}{1600}\right)^{1/2} A \bar{q}^3$                                       | $+3.6$       |                      | (1.6)                        | c                    |
| $Q_1 \rightarrow [K^* \eta]_P$           | $+\left[ \frac{\sqrt{2}+1}{\sqrt{432}} \right] D\bar{q}$                                | $+3.8$       | $\left[ \frac{D}{S} \right]$ | (1.6)                        | c                    |
| $Q_1 \rightarrow [K^* \eta]_F$           | $+\left[ \frac{\sqrt{2}+1}{\sqrt{3200}} \right] A \bar{q}^3$                            | $+0.9$       |                      | (1.6)                        | c                    |
| $Q_2 \rightarrow [\rho K]_P$             | $+\left(\frac{1}{48}\right)^{1/2} D \bar{q}$                                            | $+2.6$       | $\left[ \frac{D}{S} \right]$ | (1.6)                        | c                    |
| $Q_2 \rightarrow [\rho K]_F$             | $-\left(\frac{1}{800}\right)^{1/2} A \bar{q}^3$                                         | $-3.1$       |                      | (1.6)                        | c                    |
| $Q_2 \rightarrow [\omega K]_P$           | $-\left(\frac{1}{144}\right)^{1/2} D \bar{q}$                                           | $-1.5$       | $\left[ \frac{D}{S} \right]$ | (1.6)                        | c                    |
| $Q_2 \rightarrow [\omega K]_F$           | $+\left(\frac{1}{2400}\right)^{1/2} A \bar{q}^3$                                        | $+1.7$       |                      | (1.6)                        | c                    |

## PDF Page 18

### TABLE V. (Continued).

| Decay | Harmonic-oscillator amplitude $(\vec{q} \equiv q/\beta)$ | (MeV$^{1/2}$) |  | "Realistic" factor | Experiment$^{a}$ (MeV$^{1/2}$) | References, footnotes |
|-------|-------------------------------------------------------|--------------|---|----------------------|--------------------------|-------------------------|

#### $Q_2 \rightarrow [K^* \pi]_P$
- $-(\frac{1}{48})^{1/2} D \tilde{q}$
- $-8.9 \left[ \frac{D}{S} \right]$ 
- $7 \pm 2 $  c

#### $Q_2 \rightarrow [K^* \pi]_F$
- $+(\frac{1}{800})^{1/2}A \tilde{q}^3$
- $+0.1$  c

#### $Q_2 \rightarrow [K^* \eta]_P$
- $+ \left[ \frac{\sqrt{2} - 1}{\sqrt{288}} \right] D \tilde{q}$
- $-1.6 \left[ \frac{D}{S} \right]$ 
- $(1.6)$  c

#### $Q_2 \rightarrow [K^* \eta]_F$
- $- \frac{\sqrt{2} - 1}{\sqrt{4800}} A \tilde{q}^3$
- $-1.0$  c

#### $A_3 \rightarrow [\rho \pi]_P$
- $+(\frac{1}{27})^{1/2} D \tilde{q}$
- $+8.9 \left[ \frac{D}{S} \right]$ 
- $9.5 \pm 1 $  (nonstrange)

#### $A_3 \rightarrow [\rho \pi]_F$
- $+(\frac{1}{200})^{1/2}A \tilde{q}^3$
- $+4.5$  (1.7)

#### $A_3 \rightarrow [K^* \bar{K}]_P$
- $+(\frac{1}{108})^{1/2} D \tilde{q}$
- $+2.8 \left[ \frac{D}{S} \right]$ 
- $3.1 \pm 0.9$  (1.7)

#### $A_3 \rightarrow [K^* \bar{K}]_F$
- $+(\frac{1}{800})^{1/2} A \tilde{q}^3$
- $+0.7$  (1.7)

#### $\omega \rightarrow [K^* \bar{K}]_P$
- $+(\frac{1}{108})^{1/2} D \tilde{q}$
- $+2.8 \left[ \frac{D}{S} \right]$ 
- (1.6)

#### $\omega \rightarrow [K^* \bar{K}]_F$
- $+(\frac{1}{800})^{1/2}A \tilde{q}^3$
- $+0.7$  (1.6)

#### $\phi \rightarrow [K^* K]_P$
- $-(\frac{1}{54})^{1/2} D \tilde{q}$
- $-6.0 \left[ \frac{D}{S} \right]$ 
- 

#### $\phi \rightarrow [K^* K]_F$
- $-(\frac{1}{400})^{1/2}A\tilde{q}^3$
- $-2.8$

***

#### $1^3 D_1$

#### $\rho_D \rightarrow \pi \pi$
- $-(\frac{5}{162})^{1/2} D \tilde{q}$
- $-10 \left[ \frac{D}{S} \right]$
- see Sec. VA

#### $\rho_D \rightarrow \omega \pi$
- $-(\frac{5}{324})^{1/2} D \tilde{q}$
- $-5.5 \left[ \frac{D}{S} \right]$
- see Sec. VA

#### $\rho_D \rightarrow K \bar{K}$
- $-(\frac{5}{324})^{1/2} D \tilde{q}$
- $-5.9 \left[ \frac{D}{S} \right]$
- see Sec. VA

#### $\rho_D \rightarrow \rho \eta$
- $-(\frac{5}{648})^{1/2} D \tilde{q}$
- $-2.9 \left[ \frac{D}{S} \right]$
- see Sec. VA

#### $\rho_D \rightarrow K^* \bar{K}$
- $+(\frac{5}{648})^{1/2} D \tilde{q}$
- $+2.5 \left[ \frac{D}{S} \right]$
- see Sec. VA

#### $\omega_D \rightarrow \rho \pi$
- $+(\frac{5}{108})^{1/2} D \tilde{q}$
- $+9.7 \left[ \frac{D}{S} \right]$
- see Sec. VA

#### $\omega_D \rightarrow \omega \eta$
- $-(\frac{5}{648})^{1/2} D \tilde{q}$
- $-2.8 \left[ \frac{D}{S} \right]$
- see Sec. VA

#### $\omega_D \rightarrow \omega \eta'$
- $-(\frac{5}{648})^{1/2} D \tilde{q}$
- "below threshold"
- see Sec. VA

#### $\omega_D \rightarrow K K$
- $-(\frac{5}{324})^{1/2} D \tilde{q}$
- $-5.9 \left[ \frac{D}{S} \right]$
- see Sec. VA

## PDF Page 19

TABLE V. (Continued).

| Decay | Harmonic-oscillator amplitude $(\tilde{q} \equiv q / \beta)$ | (MeV$^{1/2}$) | "Realistic" factor | Experiment$^a$ (MeV$^{1/2}$) | References, footnotes |
|--------------------------|------------------------------------------|---------------|----------------------|----------------------------|-----------------------|
| $\omega_D \to K^* \bar{K}$ | $+(\frac{5}{648})^{1/2} D \tilde{q}$ | $+2.5$ | $\left[\frac{D}{S}\right]$ | see Sec. V A | |
| $\phi_D \to \phi \eta$    | $+(\frac{5}{324})^{1/2} D \tilde{q}$ | $+4.1$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $\phi_D \to K \bar{K}$    | $+(\frac{5}{162})^{1/2} D \tilde{q}$ | $+10$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $\phi_D \to K^* \bar{K}$  | $+(\frac{5}{324})^{1/2} D \tilde{q}$ | $+5.4$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $K_D^* \to K \pi$         | $-(\frac{5}{216})^{1/2} D \tilde{q}$  | $-8.9$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $K_D^* \to K \eta$        | $-\left[ \frac{\sqrt{10}+\sqrt{5}}{\sqrt{1296}} \right] D \tilde{q}$ | $-7.7$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $K_D^* \to K \eta'$       | $+\left[ \frac{\sqrt{10}-\sqrt{5}}{\sqrt{1296}} \right] D \tilde{q}$ | $+0.8$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $K_D^* \to K^* \pi$       | $-(\frac{5}{432})^{1/2} D \tilde{q}$  | $-5.0$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $K_D^* \to K^* \eta$      | $+\left[ \frac{\sqrt{10}-\sqrt{5}}{\sqrt{2592}} \right] D \tilde{q}$ | $+0.6$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $K_D^* \to \rho K$        | $+(\frac{5}{432})^{1/2} D \tilde{q}$  | $+4.6$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $K_D^* \to \omega K$      | $-(\frac{5}{1296})^{1/2} D \tilde{q}$ | $-2.6$ | $\left[\frac{D}{S}\right]$ | see Sec. V A |  |
| $K_D^* \to \phi K$        | $-(\frac{5}{648})^{1/2} D \tilde{q}$ | $-2.5$ | $\left[\frac{D}{S}\right]$ | see Sec. V A | $2^1 S_0$ |
| $\pi' \to \rho \pi$       | $-(\frac{1}{27})^{1/2} P \tilde{q}$ | $-4.8$ | $\left[\frac{P}{S}\right]$ | $\sim 14$ |  |
| $\pi' \to K^* \bar{K}$    | $-(\frac{1}{108})^{1/2} P \tilde{q}$ | below threshold |  |  |  |
| $\eta_r \to (K \pi)_{K^*}$ | $-(\frac{1}{108})^{1/2} P \tilde{q}$ | $-0.11 - 1.1$ | $\left[\frac{P}{S}\right]$ |  | $f$ |
| $\eta'_r \to K^* \bar{K}$ | $+(\frac{1}{54})^{1/2} P \tilde{q}$  | $-0.46 + 2.4$ | $\left[\frac{P}{S}\right]$ |  |  |
| $K' \to K^* \pi$          | $-(\frac{1}{72})^{1/2} P \tilde{q}$  | $-3.2$ | $\left[\frac{P}{S}\right]$ | $10 \pm 2$ |  |
| $K' \to K^* \eta$         | $- \left[ \frac{\sqrt{2}+1}{\sqrt{432}} \right] P \tilde{q}$ | $-0.5$ | $\left[\frac{P}{S}\right]$ |  |  |

## PDF Page 20

### TABLE V. (Continued)

| Decay | Harmonic-oscillator amplitude $(\vec{q} \equiv q/\beta)$ | (MeV$^{1/2}$) | "Realistic" factor | Experiment$^a$ (MeV$^{1/2}$) | References, footnotes |
|-------|-----------------------------------------|------------------|--------------------|-------------------------------|----------------------|
| $K' \rightarrow \rho K$ | $-\left(\frac{1}{72}\right)^{1/2} P \tilde{q}$ | $-2.4 \left(\frac{P}{S}\right)$ | | $6 \pm 1$ |  |
| $K' \rightarrow \omega K$ | $+\left(\frac{1}{216}\right)^{1/2} P \tilde{q}$ | $+1.3 \left(\frac{P}{S}\right)$ | |  |  |
| $K' \rightarrow \phi K$ | $-\left(\frac{1}{108}\right)^{1/2} P \tilde{q}$ | below threshold $2^3 S_1$ | |  |  |
| $\rho_S \rightarrow \pi \pi$ | $+\left(\frac{1}{162}\right)^{1/2} P \tilde{q}$ | $+4.0 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |
| $\rho_S \rightarrow \omega \pi$ | $-\left(\frac{1}{81}\right)^{1/2} P \tilde{q}$ | $-3.7 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |
| $\rho_S \rightarrow (\pi \pi)_\rho \eta$ | $-\left(\frac{1}{162}\right)^{1/2} P \tilde{q}$ | $-1.3 \left(\frac{P}{S}\right)$ | | see Sec. VA | f |
| $\rho_S \rightarrow K \bar{K}$ | $+\left(\frac{1}{324}\right)^{1/2} P \tilde{q}$ | $+2.0 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |
| $\rho_S \rightarrow (K \pi)_{K^*} \bar{K}$ | $+\left(\frac{1}{162}\right)^{1/2} P \tilde{q}$ | $+0.5 \left(\frac{P}{S}\right)$ | | see Sec. VA | f |
| $\omega_S \rightarrow (\pi \pi)_\rho \pi$ | $+\left(\frac{1}{27}\right)^{1/2} P \tilde{q}$ | $-5.8 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |
| $\omega_S \rightarrow \omega \eta$ | $-\left(\frac{1}{162}\right)^{1/2} P \tilde{q}$ | $-1.3 \left(\frac{P}{S}\right)$ | | see Sec. VA | f |
| $\omega_S \rightarrow \omega \eta'$ | $-\left(\frac{1}{162}\right)^{1/2} P \tilde{q}$ | below threshold | | see Sec. VA |  |
| $\omega_S \rightarrow K \bar{K}$ | $+\left(\frac{1}{324}\right)^{1/2} P \tilde{q}$ | $+2.0 \left(\frac{P}{S}\right)$ | | see Sec. VA | f |
| $\omega_S \rightarrow (K \pi)_{K^*} \bar{K}$ | $+\left(\frac{1}{162}\right)^{1/2} P \tilde{q}$ | $+0.6 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |
| $\phi_S \rightarrow \rho \pi$ | $0$ | $0$ | | see Sec. VA |  |
| $\phi_S \rightarrow \phi \eta$ | $+\left(\frac{1}{81}\right)^{1/2} P \tilde{q}$ | $+1.8 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |
| $\phi_S \rightarrow K \bar{K}$ | $-\left(\frac{1}{162}\right)^{1/2} P \tilde{q}$ | $-3.8 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |
| $\phi_S \rightarrow K^* \bar{K}$ | $+\left(\frac{1}{81}\right)^{1/2} P \tilde{q}$ | $+3.4 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |
| $K_S^* \rightarrow K \pi$ | $+\left(\frac{1}{216}\right)^{1/2} P \tilde{q}$ | $+3.4 \left(\frac{P}{S}\right)$ | | see Sec. VA |  |

## PDF Page 21

TABLE V. (Continued).

| Decay                           | Harmonic-oscillator amplitude | (MeV$^{1/2}$) | "Realistic" factor | Experiment$^{a}$ (MeV$^{1/2}$) | References, footnotes |
|----------------------------------|-------------------------------|-----------------|-------------------|---------------------|----------------------|
|                                 | ($\vec{q}\equiv q/\beta$)   |                 |                   |                     |                      |
| $K_S^* \rightarrow K \eta$      | $+ \left[ \frac{\sqrt{2}+1}{\sqrt{1296}} \right] P\vec{q} $      | $+2.7$            | $\left[ \frac{P}{S} \right]$ | see Sec. VA           |                      |
| $K_S^* \rightarrow K \eta'$     | $- \left[ \frac{\sqrt{2}-1}{\sqrt{1296}} \right] P\vec{q} $      | $-0.2$            | $\left[ \frac{P}{S} \right]$ | see Sec. VA           |                      |
| $K_S^* \rightarrow \rho K$      | $+\left(-\frac{1}{108}\right)^{1/2}P\vec{q} $             | $+2.9$            | $\left[ \frac{P}{S} \right]$ | see Sec. VA           |                      |
| $K_S^* \rightarrow \omega K$    | $-\left(\frac{1}{324}\right)^{1/2}P\vec{q} $              | $-1.6$            | $\left[ \frac{P}{S} \right]$ | see Sec. VA           |                      |
| $K_S^* \rightarrow \phi K$      | $-\left( \frac{1}{162} \right)^{1/2}P\vec{q} $            | $-0.8$            | $\left[ \frac{P}{S} \right]$ | see Sec. VA           |                      |
| $K_S^* \rightarrow K^* \pi$     | $-\left( \frac{1}{108} \right)^{1/2}P\vec{q} $            | $-3.4$            | $\left[ \frac{P}{S} \right]$ | see Sec. VA           |                      |
| $K_S^* \rightarrow (K \pi) K^* \eta $ | $+ \left[ \frac{\sqrt{2}-1}{\sqrt{648}} \right] P\vec{q} $   | $+0.3$            | $\left[ \frac{P}{S} \right]$ | see Sec. VA           | f                    |
| $\delta \rightarrow \rho \pi$   | $+ \left( \frac{1}{12\,096} \right)^{1/2} A \vec{q}^4 $   | $+2.7$            | $1^3F_4$           | (2.1)                  |                      |
| $\delta \rightarrow \eta \pi$   | $- \left( \frac{1}{60\,480} \right)^{1/2} A \vec{q}^4 $   | $-1.7$            |                       |                         |                      |
| $\delta \rightarrow K K$        | $+ \left( \frac{1}{60\,480} \right)^{1/2} A \vec{q}^4 $   | $+1.4$            |                       |                         |                      |
| $\delta \rightarrow \eta' \pi$  | $- \left( \frac{1}{60\,480} \right)^{1/2} A \vec{q}^4 $   | $-0.8$            |                       |                         |                      |
| $\delta \rightarrow K^* \bar{K}$| $+ \left( \frac{1}{48\,384} \right)^{1/2} A \vec{q}^4 $   | $+0.7$            |                       |                         |                      |
| $h \rightarrow \pi \pi$         | $+ \left( \frac{1}{20\,160} \right)^{1/2} A \vec{q}^4 $   | $+3.9$            |                       |                         |                      |
| $h \rightarrow K K$             | $+ \left( \frac{1}{60\,480} \right)^{1/2} A \vec{q}^4 $   | $+1.4$            |                       |                         |                      |
| $h \rightarrow \eta \eta$       | $+ \left( \frac{1}{241\,920} \right)^{1/2} A \vec{q}^4 $  | $+1.0$            |                       |                         |                      |
| $h \rightarrow \eta \eta'$      | $- \left( \frac{1}{120\,960} \right)^{1/2} A \vec{q}^4 $  | $-0.3$            |                       |                         |                      |
| $h \rightarrow \eta' \eta'$     | $+ \left( \frac{1}{241\,920} \right)^{1/2} A \vec{q}^4 $  | $+0.6$            |                       |                         |                      |
| $h \rightarrow K^* \bar{K}$     | $+ \left( \frac{1}{48\,384} \right)^{1/2} A \vec{q}^4 $   | $+0.7$            |                       |                         |                      |
| $h' \rightarrow \pi \pi$        | $0$                                                 | $0$              |                       | (2.1)                  |                      |
| $h' \rightarrow K K$            | $+ \left( \frac{1}{30\,240} \right)^{1/2} A \vec{q}^4 $   | $+3.1$            |                       |                         |                      |
| $h' \rightarrow \eta \eta$      | $+ \left( \frac{1}{120\,960} \right)^{1/2} A \vec{q}^4 $  | $+1.4$            |                       |                         |                      |
| $h' \rightarrow \eta \eta'$     | $+ \left( \frac{1}{60\,480} \right)^{1/2} A \vec{q}^4 $   | $+1.0$            |                       |                         |                      |
| $h' \rightarrow \eta' \eta'$    | $- \left( \frac{1}{120\,960} \right)^{1/2} A \vec{q}^4 $  | $-0.2$            |                       |                         |                      |
| $h' \rightarrow K^* \bar{K}$    | $- \left( \frac{1}{24\,192} \right)^{1/2} A \vec{q}^4 $   | $-1.9$            |                       |                         |                      |
| $K^* \rightarrow K \pi$         | $+ \left( \frac{1}{40\,320} \right)^{1/2} A \vec{q}^4 $   | $-2.7$            |                       |                         |                      |
| $K^* \rightarrow K^* \pi$       | $+ \left( \frac{1}{32\,256} \right)^{1/2} A \vec{q}^4 $   | $+1.7$            |                       | $5 \pm 1$            |                      |
| $K^* \rightarrow \rho K$        | $+ \left( \frac{1}{32\,256} \right)^{1/2} A \vec{q}^4 $   | $+1.6$            |                       | $1.2 \pm 0.4$        |                      |
| $K^* \rightarrow \omega K$      | $- \left( \frac{1}{96\,768} \right)^{1/2} A \vec{q}^4 $   | $-0.9$            |                       |                         |                      |
| $K^* \rightarrow K \eta$        | $+ \left[ \frac{\sqrt{2}-1}{\sqrt{241\,920}} \right] A \vec{q}^4 $ | $+0.3$  |                       |                         |                      |
| $K^* \rightarrow K \eta'$       | $- \left[ \frac{\sqrt{2}+1}{\sqrt{241\,920}} \right]A \vec{q}^4 $ | $-0.9$  |                       | $4 \pm 1$            |                      |

## PDF Page 22

210

STEPHEN GODFREY AND NATHAN ISGUR

32

TABLE V. (Continued).

-------------------------------------------------------------------------------
| Decay                                    | Harmonic-oscillator amplitude  | "Realistic" factor | Experiment\textsuperscript{a} | References, |
|                                          | $(\tilde{q} \equiv q/\beta)$   | (MeV$^{1/2}$)      | (MeV$^{1/2}$)                  | footnotes   |
-------------------------------------------------------------------------------

Charmed mesons $1^3S_1$

$D^{*+} \rightarrow D^0 \pi^+$ 
$-\left( \frac{2}{3} \right)^{1/2} A_c \tilde{q}$
\hspace{2.5cm} $-0.34 \left[ \frac{A_c}{A} \right]$

$D^{*+} \rightarrow D^+ \pi^0$
$+\left( \frac{1}{3} \right)^{1/2} A_c \tilde{q}$
\hspace{2.5cm} $+0.24 \left[ \frac{A_c}{A} \right]$

$D^{*0} \rightarrow D^0 \pi^0$
$-\left( \frac{1}{3} \right)^{1/2} A_c \tilde{q}$
\hspace{2.5cm} $-0.27 \left[ \frac{A_c}{A} \right]$

Charmed mesons $1^3P_2$

$K_c^* \rightarrow D\pi$
$-\left( \frac{1}{5} \right)^{1/2} A_c \tilde{q}^2 \frac{m_c \beta}{(m_c + m_d) \beta_c}$
\hspace{2.5cm} $-7.3 \left[ \frac{A_c}{A} \right]$ \hspace{1cm} d

$K_c^* \rightarrow D^* \pi$
$-\left( \frac{3}{10} \right)^{1/2} A_c \tilde{q}^2 \frac{m_c \beta}{(m_c + m_d) \beta_c}$
\hspace{2.5cm} $-5.1 \left[ \frac{A_c}{A} \right]$ \hspace{1cm} d

$1^3P_1$ and $1^1P_1$

$Q_{1c} \rightarrow [D^* \pi]_S$
$+\left(\frac{1}{6}\right)^{1/2} S_c$
\hspace{2.5cm} $-1.5 \left[ \frac{S_c}{S} \right]$ \hspace{1cm} b,d

$Q_{1c} \rightarrow [D^* \pi]_D$
$+\left(\frac{1}{3}\right)^{1/2} A_c \tilde{q}^2 \frac{m_c \beta}{(m_c + m_d) \beta_c}$
\hspace{2.5cm} $+5.3 \left[ \frac{A_c}{A} \right]$ \hspace{1cm} b,d

$Q_{2c} \rightarrow [D^* \pi]_S$
$+\left( \frac{1}{3} \right)^{1/2} S_c$
\hspace{2.5cm} $+15 \left[ \frac{S_c}{S} \right]$ \hspace{1cm} b,d

$Q_{2c} \rightarrow [D^* \pi]_D$
$-\left( \frac{1}{6} \right)^{1/2} A_c \tilde{q}^2 \frac{m_c \beta}{(m_c + m_d) \beta_c}$
\hspace{2.5cm} $+0.7 \left[ \frac{A_c}{A} \right]$ \hspace{1cm} b,d

Charmed mesons $1^3P_0$

$\kappa_c \rightarrow D\pi$
$-\left( \frac{1}{2} \right)^{1/2} S_c$
\hspace{2.5cm} $-15 \left[ \frac{S_c}{S} \right]$ \hspace{1cm} d

***

\textsuperscript{a}From the Particle Data Group, Rev. Mod. Phys. \textbf{56}, S1 (1984), unless otherwise noted.  
\textsuperscript{b}We quote the pure $Q_{B(1)}$ (i.e., pure singlet (triplet)) amplitude formulas under $Q_{1(2)}$, where $Q_{1(2)}$ is the lower (higher) state in mass.  
\textsuperscript{c}We quote the pure $^1D_2$ ($^3D_2$) amplitude formulas under $Q_{1(2)}$, where $Q_{1(2)}$ is the lower (higher) state in mass.  
\textsuperscript{d}We denote a charmed $P$-wave meson by the name of its $I = \frac{1}{2}$ $u$-$d$-$s$ analog with a subscript $c$; $\beta_c$ denotes the harmonic-oscillator parameter appropriate to charmed mesons which in our numerical results we have taken equal to $\beta$. We also note that the form factor for these states has become
$$
\exp\left[ -\frac{1}{4} \left( \frac{m_c}{m_c + m_d} \right)^2 q^2/\beta^2 \right].
$$

\textsuperscript{e}Here, as discussed in Sec. V D, $\delta_2$ is the $q\bar{q}$ state at about 1100 MeV that decays to $\eta \pi$ and $KK$ with widths of about 200 MeV and 125 MeV, respectively. We have assumed that the observed $KK\pi$ and $\eta\pi\pi$ decays arise from this channel.  
\textsuperscript{f}We have allowed for final-state-particle width in calculating our width. Our notation is $M^* \rightarrow (AB)_M P$ for the process $M^* \rightarrow MP \rightarrow ABP$.  
\textsuperscript{g}Reference 14.  
\textsuperscript{h}Reference 9.  
\textsuperscript{i}Reference 15.  
\textsuperscript{j}This result is extremely sensitive to the $Q_A \cdot Q_B$ mixing angle of Fig. 4.

## PDF Page 23

are $A = 1.67$ and $S = 3.27$. (The structure-independent amplitude $A$ is remarkably close to the value 1.76 one would have predicted from the baryon analysis of Ref. 11.) The harmonic-oscillator parameter $\beta$ is fitted to be $0.40$ GeV, a reasonable reflection of the values one would deduce by fitting to our exact wave functions (fits to $\langle r^{2} \rangle_{\pi}$, $\langle p^{2} \rangle_{\pi}$, $\langle r^{2} \rangle_{\rho}$, $\langle p^{2} \rangle_{\rho}$, $\langle r^{2} \rangle_{A_{2}}$, $\langle p^{2} \rangle_{A_{2}}$, $\langle r^{2} \rangle_{g}$, and $\langle p^{2} \rangle_{g}$ give 0.62 GeV, 0.76 GeV, 0.43 GeV, 0.45 GeV, 0.36 GeV, 0.37 GeV, 0.33 GeV, and 0.34 GeV, respectively). The results of the analysis are given in Table V.

Since we are in possession of more realistic wave functions for the mesons than those provided by the harmonic-oscillator model, it is possible to elevate this calculation to a second (but still low) level of sophistication. As an illustration, consider the two decays $A_{2} \rightarrow \rho \pi$ and $f \rightarrow \pi \pi \pi$. In the above harmonic-oscillator approximation the $\rho$ and $\pi$ have the same wave function while the actual wave functions make the $\pi$ considerably smaller than the $\rho$ as a result of hyperfine interactions; a more realistic calculation of decay amplitudes would certainly reflect this difference. Unfortunately, it is difficult to take these effects into account with a spectator model of the sort we are using here: the real $M^{*} \rightarrow MP$ amplitude necessarily involves an integral that is linear in each of the three meson wave functions and which only collapses to spectator form in the limit that $P$ becomes pointlike. A completely satisfactory treatment must therefore await a more dynamical calculation involving $q\bar{q}$ pair creation. In the meantime we can semiquantitatively explore these wave-function effects by considering the ratio of spatial matrix elements,
$$
\frac{\langle {}^{3}S_{1}| r^{L-1}| M^{*} \rangle}{\langle {}^{1}S_{0} | r^{L-1} | M^{*} \rangle}
\tag{20}
$$
for type-$A$ amplitudes ($L$ is the orbital angular momentum of the decay products) and
$$
\frac{\langle {}^{1}S_{0} | p | M^{*} \rangle}{\langle {}^{3}S_{1} | p | M^{*} \rangle}
\tag{21}
$$
for type-$S$ amplitudes. We show the first ratio in parentheses to the right of the naive harmonic-oscillator results for the $M^{*} \rightarrow P \, {}^{3}S_{1}$ decays (which were fit to $\rho \rightarrow \pi \pi$) and the second to the right of the $M^{*} \rightarrow P \, {}^{1}S_{0}$ decays [which were fit to $B \rightarrow (\omega \pi)_{S}$] as a rough indication of the probable (and apparently beneficial) effect that the use of more realistic wave functions would have on our results.

This last observation indicates that a more sophisticated treatment of meson decays could be worthwhile. Such alternative treatments are certainly required in any event if one is to be convinced that the above extremely crude model is not in some cases at least qualitatively misleading. The preliminary results of such a treatment$^{16}$ indicate that the above analysis is actually extraordinarily similar to the results of a more sophisticated analysis.

### C. Photon emission

While the data on strong decay amplitudes is richer, photon emission by excited mesons is much better founded theoretically: in the latter case we know the interaction Hamiltonian and we also know that a spectator model should be approximately valid. We can as a result come closer to a complete calculation of the operative current matrix elements $\langle M | j_{\text{em}}^{\mu} | M^{*} \rangle$. Our calculation of these radiative decay amplitudes (and of many other amplitudes discussed in the next section) is based on the results of Ref. 8 on relativistic modifications to decay amplitudes. In analogy to Sec. II, our prescription for taking these effects qualitatively into account is to modify the leading nonrelativistic decay amplitudes by smearing factors$^{17}$ of the type $(m/E)^{f}$ where the exponents $f$ are given in the tables following; the bag-model result that magnetic moments of relativistic quarks are proportional to $1/E$ instead of $1/m$ is probably the best known example of this phenomenon. Of course, these relativistic modifications to the quark matrix elements we encounter will not make our calculations fully relativistic: other less tractable but comparable effects are necessarily being neglected. Here we also follow Ref. 8 and deal with the resulting ambiguities by assuming a correspondence of real mesons with weakly bound ones ("mock mesons"). For details see Appendix D.

In the nonrelativistic limit (which is a good approximation for the heavy-quark systems) the photon helicity amplitudes are obtained by taking $\mathbf{q} = q \hat{z}$ and the corresponding photon polarization vectors $\mathbf{e}_{\gamma}$ in
$$
A_{s's} = \langle M(s') | \frac{e_{Q}}{m_{Q}} \left( -i \mathbf{p}' \cdot \mathbf{e}_{\gamma}^{*} + \frac{1}{2} \boldsymbol{\sigma} \cdot \mathbf{q} \times \mathbf{e}_{\gamma}^{*} \right) e^{-i \mathbf{q} \cdot \mathbf{r}_{q}}
+ (Q \rightarrow \bar{Q}) | M^{*}(s) \rangle
\tag{22}
$$
with
$$
\Gamma(M^{*} \rightarrow M\gamma) = \left( \frac{1}{2j^{*} + 1} \right) \frac{q}{2\pi} \sum_{s's} |A_{s's}|^{2} .
$$

[See the related discussion for meson emission in Appendix C and compare this formula to the meson-emission amplitudes of Eq. (19)]. In the nonrelativistic limit the mock-meson method and Eq. (22) of course agree.$^{18,19}$ In this subsection we actually use a hybrid of the two approaches which allows us to take into account the $m \leftrightarrow E$ ambiguity (see the end of Table VI). In neither this section nor the next do we consider QCD corrections to decays although in many cases they have been calculated. Aside from the fact that the systematic treatment of such effects is beyond the scope of this paper, we note that our decay results already contain some QCD effects: our wave functions have been affected not only by the Coulomb-type $-4\alpha_{s}/3r$ term but also by various transverse-exchange effects like the hyperfine interaction. The extraction of those effects we have not already thereby included from the available calculations would definitely be nontrivial. In any event one must expect that for light-quark states residual effects of this type could lead to corrections of order $\alpha_{s}^{\text{critical}} \simeq \frac{1}{2}$, although appropriate ratios should be more reliable than this.

## PDF Page 24

TABLE VI. Photon-decay amplitudes. The factors $I_m(x,y)$ and $E_n^i(x,y)$ appearing in the table are (see Appendix D)

$$
I_i(x,y) = \frac{(4\tilde{m}_x\tilde{m}_y)^{1/2}}{\tilde{m}_x + \tilde{m}_y}
\int dp\, p^2 \Phi_x^*(p) \Phi_y(p) \frac{1}{m_i} \left| \frac{m_i}{E_i} \right|^{0.7},
$$

$$
E_n^i(x,y) = \left| \frac{m_i}{(\langle E_i \rangle_x \langle E_i \rangle_y)^{1/2}} \right|^{0.5} \int dr\, r^2 R_x^*(r) R_y(r) r^n,
$$

where $E = (p^2 + m^2)^{1/2}$, $\Phi(p) = \Phi(p) Y_{lm}(\Omega_p)$, and $\psi(r) = R(r) Y_{lm}(\Omega_r)$, where $\Phi(p)$ and $\psi(r)$ are full normalized momentum and space wave functions. To proceed from the tabulated magnetic moments to rates, use

$$
\Gamma(V \rightarrow P \gamma) = \frac{\alpha}{3M_N^2} (\mu/\mu_N)^2 \omega^3
$$

and

$$
\Gamma(P \rightarrow V \gamma) = \frac{\alpha}{M_N^2} (\mu/\mu_N)^2 \omega^3
$$

The squares of other amplitudes are the corresponding rates. Note that the exponents of $m/E$ in the expressions for $I_m$ and $E_n^i$ were chosen to fit $\rho \rightarrow \pi \gamma$ and $A_2 \rightarrow \pi \gamma$, respectively; of course such factors only affect light mesons. We have allowed for isoscalar mixing in our numerical results. The results in the formula column are for ideal mixing [as in (B10) and (B11)] in every nonet except $^1S_0$ where we show perfect-mixing formulas from (B14) and (B15). The isoscalar mixings are taken from Table III (using P1 for pseudoscalars). A mixing angle not given explicitly in Table III is assumed, for now, to be zero.

***

### Magnetic-dipole decays

| Decay                  | Magnetic moment ($\mu$)    |
|------------------------|----------------------------|
| $\rho \rightarrow \pi\gamma$          | $+\frac{1}{3}I_d(\pi,\rho)$                 |
| $\rho \rightarrow \eta\gamma$         | $+\frac{1}{\sqrt{2}} I_d(\eta_{ns},\rho)$    |
| $\eta' \rightarrow \rho\gamma$        | $+\frac{1}{\sqrt{2}} I_d(\eta_{ns},\rho)$    |
| $\omega \rightarrow \pi\gamma$        | $+I_d(\pi,\omega)$                           |
| $\omega \rightarrow \eta\gamma$       | $+\frac{1}{3\sqrt{2}} I_d(\eta_{ns},\omega)$ |
| $\eta' \rightarrow \omega\gamma$      | $-\frac{1}{3\sqrt{2}} I_d(\eta_{ns},\omega)$ |
| $\phi \rightarrow \pi\gamma$          | $0$                                          |
| $\phi \rightarrow \eta\gamma$         | $+\frac{\sqrt{2}}{3} I_s(\eta_s, \phi)$      |
| $\phi \rightarrow \eta'\gamma$        | $-\frac{\sqrt{2}}{3} I_s(\eta_s, \phi)$      |
| $\phi \rightarrow \eta'\gamma$        | $+\frac{1}{3} I_d(\eta_{ns}, \omega)$        |
| $\eta_{s} \rightarrow \omega\gamma$   | $+I_d(\eta_{ns},\omega)$                     |
| $\eta_{s} \rightarrow \phi\gamma$     | $0$                                          |
| $\eta_{s}' \rightarrow \omega\gamma$  | $0$                                          |
| $\eta_{s}' \rightarrow \rho\gamma$    | $0$                                          |
| $\eta_{s}' \rightarrow \phi\gamma$    | $-\frac{2}{3} I_s(\eta_s, \phi)$             |
| $K^{*+} \rightarrow K^+\gamma$        | $+\frac{2}{3} I_d(K, K^*) - \frac{1}{3} I_s(K, K^*)$  |
| $K^{*0} \rightarrow K^0\gamma$        | $-\frac{1}{3} I_d(K, K^*) - \frac{1}{3} I_s(K, K^*)$  |
| $D^{*+} \rightarrow D^+\gamma$        | $+\frac{2}{3} I_c(D,D^*) - \frac{1}{3} I_d(D, D^*)$   |
| $D^{*0} \rightarrow D^0\gamma$        | $+\frac{2}{3} I_c(D,D^*) + \frac{1}{3} I_d(D, D^*)$   |
| $F^{*} \rightarrow F\gamma$           | $+\frac{2}{3} I_c(F, F^*) - \frac{1}{3} I_s(F, F^*)$  |
| $B^{*+} \rightarrow B^+ \gamma$       | $-\frac{1}{3} I_b(B,B^*) - \frac{1}{3} I_d(B,B^*)$    |
| $B^{*0} \rightarrow B^0\gamma$        | $-\frac{1}{3} I_b(B,B^*) - \frac{1}{3} I_d(B,B^*)$    |
| $F_b^* \rightarrow F_b \gamma$        | $-\frac{1}{3} I_b(F_b,F_b^*) - \frac{1}{3} I_c(F_b,F_b^*)$  |
| $\psi \rightarrow \eta_c \gamma$      | $+\frac{4}{3} I_c(\eta_c, \psi)$              |


***

| Magnetic-dipole decays | Predicted $\mu$ (in units of $\mu_N$) | Experimental $\mu$ (in units of $\mu_N$) | References, footnotes |
|-----------------------|----------------------------------------|------------------------------------------|-----------------------|
| $\rho \rightarrow \pi\gamma$     | $+0.69$    | $0.69 \pm 0.04$     | fit              |
| $\rho \rightarrow \eta\gamma$    | $+1.53$    | $1.62 \pm 0.21$     |                  |
| $\eta' \rightarrow \rho\gamma$   | $+1.85$    | $1.43 \pm 0.31$     |                  |
| $\omega \rightarrow \pi\gamma$   | $+2.07$    | $2.38 \pm 0.13$     |                  |
| $\omega \rightarrow \eta\gamma$  | $+0.50$    | $0.37 \pm 0.15$     |                  |
| $\eta' \rightarrow \omega\gamma$ | $+0.63$    | $0.44 \pm 0.12$     |                  |
| $\phi \rightarrow \pi\gamma$     | $(+0.06)$  | $0.13 \pm 0.02$     | a,d              |
| $\phi \rightarrow \eta\gamma$    | $+0.71$    | $0.69 \pm 0.07$     |                  |
| $\phi \rightarrow \eta'\gamma$   | $-0.66$    |                      |                  |
| $\eta_{s} \rightarrow \omega\gamma$  | $-0.18$ |                      |                  |
| $\eta_{s} \rightarrow \rho\gamma$    | $+0.57$ |                      |                  |
| $\eta_{s}' \rightarrow \omega\gamma$ | $+0.37$ |                      |                  |
| $\eta_{s}' \rightarrow \rho\gamma$   | $+0.01$ |                      |                  |
| $\eta_{s}' \rightarrow \phi\gamma$   | $+0.008$|                      |                  |
| $K^{*+} \rightarrow K^+\gamma$       | $-0.29$ |                      |                  |
| $K^{*0} \rightarrow K^0\gamma$       | $+0.91$ | $0.86 \pm 0.11$      |                  |
| $D^{*+} \rightarrow D^+\gamma$       | $-1.20$ | $0.95 \pm 0.22$      |                  |
| $D^{*0} \rightarrow D^0\gamma$       | $-0.35$ |                      | b                |
| $F^{*} \rightarrow F\gamma$          | $+1.78$ |                      | b                |
| $B^{*+} \rightarrow B^+ \gamma$      | $-0.13$ |                      |                  |
| $B^{*0} \rightarrow B^0\gamma$       | $+1.37$ |                      |                  |
| $F_b^* \rightarrow F_b \gamma$       | $-0.78$ |                      |                  |
| $\psi \rightarrow \eta_c \gamma$     | $-0.55$ |                      | f                |
| $\psi \rightarrow \eta_c \gamma$     | $+0.69$ | $0.44 \pm 0.10$      | i                |

## PDF Page 25

TABLE VI. (Continued).

Magnetic-dipole decays

| Decay | Magnetic moment ($\mu$) (in units of $e/2$) | Predicted $\mu$ (in units of $\mu_N$) | Experimental $\mu$ (in units of $\mu_N$) | References, footnotes |
|-------|---------------------------------------------|--------------------------------------|-----------------------------------------|----------------------|
| $\psi' \rightarrow \eta_c' \gamma$ | $+\frac{4}{3} I_c(\eta_c',\psi')$ | $+0.68$ | $0.8 \pm 0.4$ |  |
| $\psi' \rightarrow \eta_c \gamma$ | $+\frac{4}{3} \left[ I_c(\eta_c, \psi') - \frac{q^2}{24m_c} E_2^{c}(\eta_c, \psi) \right]$ | $-0.056$ | $0.029 \pm 0.007$ | c |
| $\Upsilon \rightarrow \eta_b \gamma$ | $-\frac{2}{3} I_b(\eta_b, \Upsilon)$ | $-0.13$ |  |  |
| $\Upsilon' \rightarrow \eta_b' \gamma$ | $-\frac{2}{3} I_b(\eta_b', \Upsilon')$ | $-0.12$ |  |  |
| $\Upsilon' \rightarrow \eta_b \gamma$ | $-\frac{2}{3} \left[ I_b(\eta_b, \Upsilon') - \frac{q^2}{24m_b} E_2^{b}(\eta_b, \Upsilon') \right]$ | $+0.007$ |  | c |
| $\Upsilon'' \rightarrow \eta_b'' \gamma$ | $-\frac{2}{3} I_b(\eta_b'', \Upsilon'')$ | $-0.12$ |  |  |
| $\Upsilon'' \rightarrow \eta_b \gamma$ | $-\frac{2}{3} \left[ I_b(\eta_b, \Upsilon'') - \frac{q^2}{24m_b} E_2^{b}(\eta_b, \Upsilon'') \right]$ | $+0.007$ |  | c |
| $\Upsilon''' \rightarrow \eta_b \gamma$ | $-\frac{2}{3} \left[ I_b(\eta_b, \Upsilon''') - \frac{q^2}{24 m_b} E_2^{b}(\eta_b, \Upsilon''') \right]$ | $-0.004$ |  | c |
| $\psi \rightarrow \eta \gamma$ | $0$ | $(+0.001)$ | $0.003 \pm 0.001$ | d |
| $\psi \rightarrow \eta' \gamma$ | $0$ | $(+0.003)$ | $0.005 \pm 0.002$ | d |
| $\psi \rightarrow \pi^0 \gamma$ | $0$ | $(+0.0007)$ | $0.0007 \pm 0.0003$ | d,g |
| $\psi \rightarrow \eta' \gamma$ | $0$ | $(-0.002)$ | $0.007 \pm 0.002$ | d,e |
| $\psi \rightarrow \pi \gamma$ | $0$ | $(-0.002)$ |  | d |
| $\psi \rightarrow \pi \gamma$ | $0$ | $(+8 \times 10^{-5})$ |  | d |
| $\Upsilon \rightarrow \eta' \gamma$ | $0$ | $(+3 \times 10^{-6})$ |  | d |
| $\Upsilon \rightarrow \pi^0 \gamma$ | $0$ | $(\sim 0)$ |  | d,g |
| $\Upsilon \rightarrow \eta \gamma$ | $0$ | $(-9 \times 10^{-5})$ |  | d |
| $\Upsilon \rightarrow \eta' \gamma$ | $0$ | $(+3 \times 10^{-4})$ |  | d |

***

Other electric and magnetic multipole decays

| Decay | Amplitude formula (in units of $\sqrt{\alpha q}$) | Theory (MeV$^{1/2}$) | Experiment (MeV$^{1/2}$) | References, footnotes |
|-------|----------------------------------------------------|------------------------|----------------------------|----------------------|
| $A_2^{\pm} \rightarrow \pi^{\pm} \gamma$ | $\frac{q^2}{\sqrt{60} m_u} E_1^{u}(\pi, A_2)$ | $+0.55$ | $0.55 \pm 0.09$ | h, fit |
| $K^*+(1420) \rightarrow K^+ \gamma$ | $\frac{q^2}{\sqrt{60}} \left[ \frac{2}{3m_u} E_1^{u}(K, K^*) + \frac{1}{3m_s} E_1^{s}(K, K^*) \right]$ | $+0.48$ | $0.50 \pm 0.05$ | h |
| $A_2^{+,0,-} \rightarrow \rho^{+,0,-} \gamma$ | $\frac{q}{9} E_1^{u}(\rho, A_2)$ | $+0.15$ |  |  |
| $f' \rightarrow \phi \gamma$ | $-\frac{2q}{9} E_1^{s}(\phi, f')$ | $-0.31$ |  |  |
| $A_2^0 \rightarrow \omega \gamma$ | $\frac{q}{3} E_1^{u}(\omega, A_2)$ | $+0.44$ |  |  |
| $A_1^{\pm} \rightarrow \pi^{\pm} \gamma$ | $\frac{q^2}{6m_u} E_1^{u}(\pi, A_1)$ | $+0.56$ | $1.0 \pm 0.3$ | h |
| $B^{\pm} \rightarrow \pi^{\pm} \gamma$ | $\frac{\sqrt{2}q}{3} E_1^{u}(\pi, B)$ | $+0.63$ | $0.45 \pm 0.04$ | h |
| $\chi_{2c} \rightarrow \psi \gamma$ | $\frac{4q}{9} E_1^{c}(\chi_{2c}, \psi)$ | $+0.50$ | $0.58 \pm 0.18$ |  |
| $\chi_{1c} \rightarrow \psi \gamma$ | $\frac{4q}{9} E_1^{c}(\psi, \chi_{1c})$ | $+0.44$ | $<0.86$ |  |

## PDF Page 26

TABLE VI. (Continued).

Other electric and magnetic multipole decays

| Decay                                | Amplitude formula (in units of $\sqrt{\alpha q}$) | Theory (MeV$^{1/2}$) | Experiment (MeV$^{1/2}$) | References, footnotes |
|--------------------------------------|----------------------------------------------------|----------------------|--------------------------|----------------------|
| $\chi_{0c} \to \psi\gamma$           | $\frac{4q}{9} E_1^c(\psi, \chi_{0c})$              | $+0.30$              | $0.33 \pm 0.08$           |                      |
| $\psi' \to \chi_{2c}\gamma$          | $\frac{4q}{9} \left[\frac{5}{3}\right]^{1/2} E_1^c(\chi_{2c}, \psi')$ | $+0.14$ | $0.13 \pm 0.02$ | |
| $\psi' \to \chi_{1c}\gamma$          | $\frac{4q}{9} E_1^c(\chi_{1c}, \psi')$            | $+0.15$              | $0.13 \pm 0.02$           |                      |
| $\psi' \to \chi_{0c}\gamma$          | $\frac{4q}{9} \left[\frac{1}{3}\right]^{1/2} E_1^c(\chi_{0c}, \psi')$ | $+0.14$ | $0.13 \pm 0.02$ | |
| $\chi_{2b} \to \Upsilon\gamma$       | $-\frac{2q}{9} E_1^b(\Upsilon, \chi_{2b})$         | $-0.18$              |                          |                      |
| $\chi_{1b} \to \Upsilon\gamma$       | $-\frac{2q}{9} E_1^b(\Upsilon, \chi_{1b})$         | $-0.17$              |                          |                      |
| $\chi_{0b} \to \Upsilon\gamma$       | $-\frac{2q}{9} E_1^b(\Upsilon, \chi_{0b})$         | $-0.16$              |                          |                      |
| $\Upsilon' \to \chi_{2b}\gamma$      | $-\frac{2q}{9} \left[\frac{5}{3}\right]^{1/2} E_1^b(\chi_{2b}, \Upsilon')$ | $-0.040$ | $0.043 \pm 0.008$ | |
| $\Upsilon' \to \chi_{1b}\gamma$      | $-\frac{2q}{9} E_1^b(\chi_{1b}, \Upsilon')$        | $-0.038$             | $0.042 \pm 0.008$         |                      |
| $\Upsilon' \to \chi_{0b}\gamma$      | $-\frac{2q}{9} \left[\frac{1}{3}\right]^{1/2} E_1^b(\chi_{0b}, \Upsilon')$ | $-0.025$ | $0.032 \pm 0.009$ | |
| $\chi'_{2b} \to \Upsilon'\gamma$     | $-\frac{2q}{9} E_1^b(\Upsilon', \chi'_{2b})$       | $-0.12$              |                          |                      |
| $\chi'_{1b} \to \Upsilon'\gamma$     | $-\frac{2q}{9} E_1^b(\Upsilon', \chi'_{1b})$       | $+0.09$              |                          |                      |
| $\chi'_{1b} \to \Upsilon'\gamma$     | $-\frac{2q}{9} E_1^b(\Upsilon', \chi'_{1b})$       | $-0.11$              |                          |                      |
| $\chi'_{0b} \to \Upsilon'\gamma$     | $-\frac{2q}{9} E_1^b(\Upsilon', \chi'_{0b})$       | $+0.069$             |                          |                      |
| $\chi_{0b} \to \Upsilon'\gamma$      | $-\frac{2q}{9} E_1^b(\Upsilon', \chi_{0b})$        | $+0.10$              |                          |                      |
| $\Upsilon'' \to \chi'_{2b}\gamma$    | $-\frac{2q}{9} \left[\frac{5}{3}\right]^{1/2} E_1^b(\chi'_{2b}, \Upsilon'')$ | $+0.047$ | | |
| $\Upsilon'' \to \chi_{2b}\gamma$     | $-\frac{2q}{9} \left[\frac{5}{3}\right]^{1/2} E_1^b(\chi_{2b}, \Upsilon'')$ | $-0.044$ | $0.048 \pm 0.015$ | |
| $\Upsilon'' \to \chi'_{1b}\gamma$    | $-\frac{2q}{9} E_1^b(\chi'_{1b}, \Upsilon'')$      | $-0.020$             |                          |                      |
| $\Upsilon'' \to \chi_{1b}\gamma$     | $-\frac{2q}{9} E_1^b(\chi_{1b}, \Upsilon'')$       | $-0.043$             | $0.053 \pm 0.015$         |                      |
| $\Upsilon'' \to \chi'_{0b}\gamma$    | $-\frac{2q}{9} \left[\frac{1}{3}\right]^{1/2} E_1^b(\chi'_{0b}, \Upsilon'')$ | $-0.008$ | | |
| $\Upsilon'' \to \chi_{0b}\gamma$     | $-\frac{2q}{9} \left[\frac{1}{3}\right]^{1/2} E_1^b(\chi_{0b}, \Upsilon'')$ | $-0.030$ | $0.037 \pm 0.014$ | |
| $\psi \to A_2\gamma$                 | $0$                                             | $-0.002$             |                          |                      |
| $\psi \to f \gamma$                  | $0$                                             | $(+0.002)$           | $<7 \times 10^{-3}$       | d                    |
| $\psi \to f' \gamma$                 | $0$                                             | $(+0.024)$           | $0.010 \pm 0.002$         | d                    |
| $\Upsilon \to A_2\gamma$             | $0$                                             | $(+0.010)$           | $<4 \times 10^{-3}$       | d                    |
| $\Upsilon \to f \gamma$              | $0$                                             | $(+0.002)$           |                          | d                    |
| $\Upsilon \to f' \gamma$             | $0$                                             | $(-0.006)$           |                          | d                    |
|                                      |                                                 | $(-0.002)$           |                          | d                    |

***

a This moment includes a contribution of $+0.01$ which we calculate will arise from $\pi^0$-$\eta$ mixing.  
b While absolute widths are not known, the predicted branching ratios can be favorably compared to experiment by using Table V.  
c We have included the recoil term in this decay since the direct term is so small.

## PDF Page 27

32            MESONS IN A RELATIVISTIC QUARK MODEL ...                215

***

TABLE VI. (Continued).

***

d These are the amplitudes which arise from emission of the photon either from the non-$s \bar{s}$ ($c\bar{c}$, $b\bar{b}$) components of the $\phi$ ($\psi$, $\Upsilon$) or to the $s\bar{s}$ ($c\bar{c}$, $b\bar{b}$) component of the final meson as listed in Table III; note that the predicted amplitude should, but does not, include contributions for photon emission during the Zweig-rule-violating $s\bar{s}$ ($c\bar{c}$, $b\bar{b}$) to $q\bar{q}$ transition which could be comparable to the contributions for photon emission before or after Zweig violation. The "predicted" transition amplitudes, shown in parentheses, can therefore only be considered order-of-magnitude estimates.

e We tentatively identify the $\eta_r$ with the $\iota(1440)$; see Sec. V A.

f We use $F_b$ to denote a $b\bar{s}$ meson.

g Since the recoil in this decay is being absorbed by a light-quark system, we include a form factor $\exp(-q^2 / 16\beta^2)$ with $\beta = 0.40$ GeV as in Table V.

h Reference 20.

i Reference 23.

***

### D. Miscellaneous, weak, electromagnetic, and strong couplings

In addition to the very large body of data on $M^* \rightarrow MP$ and $M^* \rightarrow M \gamma$, considerable information on meson structure resides in weak-pseudoscalar decays ($P \rightarrow l \bar{\nu}$), leptonic-pair decays of vector mesons ($V \rightarrow l^+ l^-$), two-photon decays ($P \rightarrow \gamma\gamma$), gluonic decays of heavy quarkonia ($Q\bar{Q} \rightarrow$ gluons), and meson charge radii. We once again take into account relativistic corrections as described in Appendix D. That this is appropriate, in view of our approach to relativistic corrections to spectroscopy, is especially clear in some of these cases. For example,\super{21} the usual result\super{22} that $V \rightarrow l^+ l^-$ proceeds via $\Psi(0)$, the spatial wave function at zero relative coordinate, is modified by the same smearing effect that regulates the hyperfine interaction $\delta$ function as described in Sec. II. This both affects “allowed” decays and allows “forbidden” ones such as $\tau \rightarrow A_1 \nu_\tau$ and $^3D_1 \rightarrow e^+ e^-$. For convenience we have reproduced several relevant definitions in Appendix D. Our results are given in Table VII. The comments in subsection C above regarding QCD corrections apply with equal force in these cases.

***

## V. DISCUSSION

For the most part we believe that the figures and tables of the preceding sections can speak for themselves, and so we will comment rather briefly on most of our results. The major exceptions to this rule occur in our discussions of the $2^3S_1$–$1^3D_1$ complex, of the pseudoscalar self-conjugate isoscalar mesons, and of the scalar mesons.

***

### A. Pseudoscalars and vectors

Mesons with $0^{-+}$ and $1^{--}$ quantum numbers ($^1S_0$ and $^3S_1$ or $^3D_1$) are well described in the model. The spectroscopy of the ground-state $S$-wave hyperfine pairs $\rho$-$\pi$, $K^*$-$K$, $D^*$-$D$, $F^*$-$F$, and $\psi$-$\eta_c$ (along with splittings like $A_2$-$B$ and $g$-$A_3$ of subsections B and C below) provide strong evidence for the Fermi contact term with the properties predicted from one-gluon exchange. As already noted, it is very significant that the strength of this term is that expected from the Coulombic term $-4 \alpha_s(r) / 3r$ required to fit the gross spectroscopy of mesons. The success of this term depends on its short range, its $(m_i m_j)^{-1}$ mass dependence, and on the fact that $\alpha_s$ depends on $r$. It is also from these states that we have our best evidence that the model simultaneously describes orbital and radial excitations via the heavy-quark sequences ($\psi$, $\psi'$, $\psi''$, $\psi'''$, ...) and ($\Upsilon$, $\Upsilon'$, $\Upsilon''$, $\Upsilon'''$, ...). In addition, the decay analysis indicates that the internal structure of these states is well understood. Perhaps the simplest example is that the ratios of the leptonic widths of the $n\ ^3S_1$ radial excitations is in accord with the expected behavior of the ratios of the quantities $|\Psi_{n\ ^3S_1}(0)|^2$.

Along with these successes, however, we encounter two puzzles in the light-quark mesons with these quantum numbers. The first puzzle relates to the identification of the $2^3S_1$ and $1^3D_1$ states with observed meson resonances; the second concerns the self-conjugate isoscalar-pseudoscalar mesons.

In the $c\bar{c}$ system it is well known that the $1^{3}D_1$ $\psi'(3770)$ level lies just above the $2^3S_1$ $\psi'(3685)$ and we expect that for heavier systems this splitting will asymptotically approach (from above) the hydrogenic value of $+\frac{5}{27}$ of the $2S$-$1S$ splitting. Conversely, according to this picture this splitting will increase as we move toward light-quark systems so that we would expect the isovector $2^3S_1$ state to now lie about twice as far below $1^3D_1$. We know experimentally that the $1^3D_3$ state is at about 1700 MeV where we predict it to be; it follows that we have good grounds for supposing that the $2^3S_1$ state should be at about 1450 MeV, where it is predicted to be. This does not, however, accord very well with either the position of the $\rho'(1600)$, which is usually taken to be the radial excitation of the $\rho$, nor with various other usual suppositions about the radially excited vector-meson nonet. We propose an alternative interpretation of the situation. Whether or not the details of what follows are correct, we would stress that it seems to us certain that there are two sets of $1^{--}$ resonances in the 1400–1900-MeV region and we would argue that analysis of this region should always allow for this possibility.

The inhabitants of our alternative scenario are, accordingly, the $2^3S_1$ states $\rho_S(1.45)$, $\omega_S(1.46)$, $K^*_S(1.58)$, $\phi_S(1.69)$, and the $1^3D_1$ states $\rho_D(1.66)$, $\omega_D(1.66)$, $K^*_D(1.78)$, $\phi_D(1.88)$. From the decay analysis of Sec. IV we know many of their characteristics: from Table VII we see that, via relativistic corrections, the $\rho_D$, $\omega_D$, and $\phi_D$ are coupled to $e^+e^-$ with strengths comparable to those of $\rho_S$, $\omega_S$, and $\phi_S$, and from Table V (using the $A_3$ to tell us that $D/S \simeq 1$ and the $\pi'$ and $K'$ to tell us that $P/S \simeq 3$) we can roughly estimate that

## PDF Page 28

216

STEPHEN GODFREY AND NATHAN ISGUR

***

TABLE VII. Leptonic, $\gamma\gamma$, and gluonic decays; charge radii.

(a) Leptonic decays. The factors $P_P$, $P'_{A_1}$, $V_V$, and $V'_V$ are

$$
P_P = M_P \tilde{M}_P^{-1}(2\pi)^{-3/2} \int d^3 p (4\pi)^{-1/2} \Phi_P(p) \left( \frac{m_1 m_2}{E_1 E_2} \right)^{1/2},
$$

$$
P'_{A_1} = M_{A_1} \tilde{M}_{A_1}^{-2}(2\pi)^{-3/2} \int d^3 p (4\pi)^{-1/2} \Phi_{A_1}(p) \left( \frac{mp}{E^2} \right),
$$

$$
V_V = M_V^{-2} \tilde{M}_V^{1/2}(2\pi)^{-3/2} \int d^3 p (4\pi)^{-1/2} \Phi_V(p) \left( \frac{m_1 m_2}{E_1 E_2} \right)^{1/2},
$$

$$
V'_V = M_V^{-2} \tilde{M}_V^{1/2}(2\pi)^{-3/2} \int d^3 p (4\pi)^{-1/2} \Phi_V(p) \frac{m}{E} \left( 1 - \frac{m}{E} \right).
$$

See Appendix D for further details. Note that $\phi(p) \equiv \Phi(p) Y_{lm}(\Omega_p)$. As usual, the formulas shown are for unmixed states; the effects of mixing are taken into account in the numerical results. Note that an implicit exponent for $m/E$ has been chosen to be unity in these formulas.

(b) $\gamma\gamma$ decays. These predictions are based on the formula

$$
A(P \rightarrow \gamma\gamma) = \sqrt{6} e_q^2 \frac{\alpha}{m} \left( \frac{M_P}{\tilde{M}_P} \right)^{3/2} \frac{1}{2\pi} \int d^3 p \, \phi_P(p) \left( \frac{m}{E} \right)
$$

for an $S$-wave pseudoscalar made of $q\bar{q}$ of mass $m$ in a wave function $\phi_P(p)$; and

$$
A({}^3P_2 \rightarrow \gamma\gamma) = - \left( \frac{4}{5} \right)^{1/2} e_q^2 \frac{\alpha}{m} \left( \frac{M_{^3P_2}}{\tilde{M}_{^3P_2}} \right)^{3/2} \left( \frac{2}{\pi} \right)^{1/2} \int dp\, p^2 \Phi(p) \left( \frac{mp}{E^2} \right)
$$

for a $P$-wave meson with radial wave function $\Phi(p)$; see Appendix D and Ref. 21 for details. Note that an implicit exponent for $m/E$ has once again been chosen to be unity in these formulas.

(c) Gluonic decays. We use here the lowest-order QCD formulas

$$
\Gamma({}^1S_0 \rightarrow 2g) = \frac{8\pi \alpha_s^2}{3m_Q^2} |S_0(\Psi)|^2, \quad \Gamma({}^3S_1 \rightarrow 3g) = \frac{40(\pi^2 - 9)}{81 m_Q^2} \alpha_s^3 |S_0(\Psi)|^2,
$$

and

$$
\Gamma({}^3P_2 \rightarrow 2g) = \frac{32\pi \alpha_s^2}{45 m_Q^2} |S_1(\Psi)|^2, \quad \Gamma({}^3P_0 \rightarrow 2g) = \frac{8\pi \alpha_s^2}{3 m_Q^2} |S_1(\Psi)|^2,
$$

where $\alpha_s = \alpha_s(\mu^2)$ with $\mu$ the mass of the decaying meson and where $S_L(\Psi)$ is defined in Eq. (17).

(d) Charge radii. These predictions are based on the result (from Ref. 8 modified in the usual way)

$$
r_E^2 = \sum_i e_i \left[ \left< r_i^2 \right> + \frac{3}{4 m_i^2} \int d^3 p\, |\phi(p)|^2 \left( \frac{m_i}{E_i} \right)^{2f} \right],
$$

where the first term is the ordinary expectation value of the quark radius vector and the second an approximate formula for the relativistic smearing of the quark-position operator. We stress that the second term (with the value $f = 0.2$ we obtain by fitting $r_E^2$ of the $\pi^+$) is very important for light-quark systems and accounts for the usual discrepancy between the size of the wave function of such systems and their measured charge radius.

***

Leptonic decays

\begin{tabular}{lccccc}
Typical & Amplitude & Amplitude & Predicted & Experimental & References,\\
decay & Amplitude & formula & amplitude & amplitude & footnotes \\
\hline
$\pi \rightarrow \mu \nu$ & $f_\pi$ & $2\sqrt{3}P_\pi$ & $+1.3$ & $0.95$ \\
$K \rightarrow \mu \nu$ & $f_K$ & $2\sqrt{3}P_K$ & $+0.47$ & $0.32$ \\
$D \rightarrow \mu \nu$ & $f_D$ & $2\sqrt{3}P_D$ & $+0.13$ \\
\end{tabular}

## PDF Page 29

TABLE VII. (Continued).

-------------------------------------------------------------------------------

Typical decay | Amplitude | Amplitude formula | Leptonic decays | Predicted amplitude | Experimental amplitude | References, footnotes
--- | --- | --- | --- | --- | --- | ---
$F \rightarrow \mu \nu$ | $f_F$ | $2\sqrt{3} P_F$ | $+0.17$
$B \rightarrow \mu \nu$ | $f_B$ | $2\sqrt{3} P_B$ | $+0.033$
$1 \, ^1 S_0 (b\bar{c}) \rightarrow \mu \bar{\nu}$ | $f_{b\bar{c}}$ | $2\sqrt{3} P_{b\bar{c}}$ | $+0.091$
$\tau \rightarrow A_1 \nu$ | $f_{A_1}$ | $2\sqrt{2} P'_{A_1}$ | $+0.11$ | $0.10 \pm 0.04$
$\tau \rightarrow K^*(892)\nu$ | $f_{K^*}$ | $2\sqrt{3} V^{K^*}$ | $+0.30$ | $0.28 \pm 0.07$
$\rho \rightarrow e^+ e^-$ | $f_\rho$ | $\sqrt{6} V_\rho$ | $+0.20$ | $0.20 \pm 0.02$
$\omega \rightarrow e^+ e^-$ | $f_\omega$ | $(\frac{2}{3})^{1/2} V_\omega$ | $+0.07$ | $0.07 \pm 0.01$
$\phi \rightarrow e^+ e^-$ | $f_\phi$ | $-(\frac{4}{3})^{1/2} V_\phi$ | $-0.11$ | $0.07 \pm 0.01$
$\rho_S \rightarrow e^+ e^-$ | $f_{\rho_S}$ | $\sqrt{6} V_{\rho_S}$ | $-0.037$
$\rho_D \rightarrow e^+ e^-$ | $f_{\rho_D}$ | $(\frac{4}{3})^{1/2} V'_{\rho_D}$ | $+0.019$
$\omega_S \rightarrow e^+ e^-$ | $f_{\omega_S}$ | $(\frac{2}{3})^{1/2} V_{\omega_S}$ | $-0.012$
$\omega_D \rightarrow e^+ e^-$ | $f_{\omega_D}$ | $(\frac{4}{27})^{1/2} V'_{\omega_D}$ | $+0.006$
$\phi_S \rightarrow e^+ e^-$ | $f_{\phi_S}$ | $-(\frac{4}{3})^{1/2} V_{\phi_S}$ | $+0.027$
$\phi_D \rightarrow e^+ e^-$ | $f_{\phi_D}$ | $-(\frac{8}{27})^{1/2} V'_{\phi_D}$ | $-0.012$
$\psi \rightarrow e^+ e^-$ | $f_\psi$ | $(\frac{16}{3})^{1/2} V_\psi$ | $+0.12$ | $0.08 \pm 0.01$
$\psi' \rightarrow e^+ e^-$ | $f_{\psi'}$ | $(\frac{16}{3})^{1/2} V_{\psi'}$ | $-0.063$ | $0.049 \pm 0.007$
$\psi'' \rightarrow e^+ e^-$ | $f_{\psi''}$ | $(\frac{16}{3})^{1/2} V_{\psi''}$ | $+0.011$ | $0.020 \pm 0.003$
$\psi''' \rightarrow e^+ e^-$ | $f_{\psi'''}$ | $(\frac{16}{3})^{1/2} V_{\psi'''}$ | $+0.046$ | $0.028 \pm 0.007$
$\Upsilon \rightarrow e^+ e^-$ | $f_{\Upsilon}$ | $(\frac{4}{3})^{1/2} V_{\Upsilon}$ | $-0.026$ | $0.023 \pm 0.002$
$\Upsilon' \rightarrow e^+ e^-$ | $f_{\Upsilon'}$ | $-(\frac{4}{3})^{1/2} V_{\Upsilon'}$ | $+0.017$ | $0.015 \pm 0.002$
$\Upsilon'' \rightarrow e^+ e^-$ | $f_{\Upsilon''}$ | $-(\frac{4}{3})^{1/2} V_{\Upsilon''}$ | $-0.014$ | $0.012 \pm 0.002$
$\Upsilon''' \rightarrow e^+ e^-$ | $f_{\Upsilon'''}$ | $+(\frac{16}{3})^{1/2} V_{\Upsilon'''}$ | $+0.012$ | $0.010 \pm 0.002$
$1\, ^3D_1(b\bar{b}) \rightarrow e^+ e^-$ | $f_{1 ^3D_1(b\bar{b})}$ | $-(\frac{8}{27})^{1/2} V'_{1 ^3D_1(b\bar{b})}$ | $-0.0008$
$\xi \rightarrow e^+ e^-$ | $f_\xi$ | $+(\frac{16}{3})^{1/2} V_\xi$ | $+0.01$ | hypothetical $\bar{t}t$ $1\,^3S_1$ meson

***

$\gamma\gamma$ decays

| Decay | Predicted amplitude | Experimental amplitude | References, footnotes |
|-------|---------------------|------------------------|----------------------|
$\pi \rightarrow \gamma\gamma$ | $2.6$ eV$^{1/2}$ | $2.8 \pm 0.1$ eV$^{1/2}$
$\pi' \rightarrow \gamma\gamma$ | $-1.0$ keV$^{1/2}$
$\eta \rightarrow \gamma\gamma$ | $0.50$ keV$^{1/2}$ | $0.57 \pm 0.05$ keV$^{1/2}$
$\eta' \rightarrow \gamma\gamma$ | $1.3$ keV$^{1/2}$ | $2.3 \pm 0.6$ keV$^{1/2}$
$\eta_r \rightarrow \gamma\gamma$ | $-2.7$ keV$^{1/2}$
$\eta'_r \rightarrow \gamma\gamma$ | $-2.2$ keV$^{1/2}$
$\eta_c \rightarrow \gamma\gamma$ | $2.6$ keV$^{1/2}$ | $< 5$ keV$^{1/2}$ | a
$\eta'_c \rightarrow \gamma\gamma$ | $-2.2$ keV$^{1/2}$
$\eta_b \rightarrow \gamma\gamma$ | $0.62$ keV$^{1/2}$
$\eta_t \rightarrow \gamma\gamma$ | $1.5$ keV$^{1/2}$ | hypothetical $\bar{t}t$ $1\,^1S_0$ meson
$A_2 \rightarrow \gamma\gamma$ | $-1.2$ keV$^{1/2}$ | $0.9 \pm 0.1$ keV$^{1/2}$
$f \rightarrow \gamma\gamma$ | $-1.9$ keV$^{1/2}$ | $1.7 \pm 0.1$ keV$^{1/2}$
$f'' \rightarrow \gamma\gamma$ | $-0.25$ keV$^{1/2}$ | $0.3 \pm 0.1$ keV$^{1/2}$ | b, c

***

Gluonic decays

| Decay | Predicted amplitude (MeV$^{1/2}$) | Experimental amplitude (MeV$^{1/2}$) | References, footnotes |
|-------|-----------------------|------------------------|----------------------|
$\eta_c \rightarrow 2g$ | $4.7$ | $3.5 \pm 0.6$
$\psi \rightarrow 3g$ | $0.42$ | $0.21 \pm 0.01$
$\eta'_c \rightarrow 2g$ | $-2.7$ | $< 3$

## PDF Page 30

TABLE VII. (Continued).

***

#### Gluonic decays

| Decay                | Predicted amplitude (MeV$^{1/2}$) | Experimental amplitude (MeV$^{1/2}$) | References, footnotes |
|----------------------|-------------------------------------|---------------------------------------|----------------------|
| $\psi' \to 3g$         | $-0.28$       | $0.22 \pm 0.08$  |        |
| $\eta_b \to 2g$        | $2.5$         |                   |        |
| $\Upsilon \to 3g$      | $0.21$        | $0.21 \pm 0.04$   |        |
| $\eta_b' \to 2g$       | $-1.7$        |                   |        |
| $\Upsilon' \to 3g$     | $-0.15$       |                   |        |
| $\Upsilon'' \to 3g$    | $0.13$        |                   |        |
| $\Upsilon''' \to 3g$   | $-0.11$       |                   |        |
| $\eta_t \to 2g$        | $1.1$         | hypothetical $t\bar{t}$ $^1S_0$ meson | |
| $\zeta \to 3g$         | $0.10$        | hypothetical $t\bar{t}$ $^3S_1$ meson | |
| $\chi_{2c} \to 2g$     | $0.88$        | $1.6 \pm 0.3$     |      |
| $\chi_{0c} \to 2g$     | $2.5$         | $4 \pm 0.5$       |      |
| $\chi_{2b} \to 2g$     | $0.35$        |                   |      |
| $\chi_{0b} \to 2g$     | $0.82$        |                   |      |
| $\chi_{2b}' \to 2g$    | $-0.37$       |                   |      |
| $\chi_{0b}' \to 2g$    | $-0.82$       |                   |      |

***
#### Charge radii

| Meson     | Predicted $r_E^2$ (fm$^2$)      | Measured $r_E^2$ (fm$^2$)     | References, footnotes |
|-----------|----------------------------------|-------------------------------|----------------------|
| $\pi^+$   | $+ (0.66)^2$                    | $+ (0.66 \pm 0.02)^2$         | d, fit               |
| $K^+$     | $+ (0.59)^2$                    | $+ (0.53 \pm 0.07)^2$         | e                    |
| $K^0$     | $-(0.30)^2$                     | $- (0.23 \pm 0.05)^2$         | f                    |

***

aReference 23.  
bThis amplitude is very sensitive to $f$-$f'$ mixing.  
cReference 24.  
dReference 25.  
eReference 26.  
fReference 27.

***

$\Gamma(\rho_S \to \pi\pi) \simeq \Gamma(\rho_S \to \omega\pi) \simeq 125\ \mathrm{MeV}$,

$\Gamma(\omega_S \to \rho\pi) \simeq 300\ \mathrm{MeV}$,

$\Gamma(K_S^* \to K\pi) \simeq \Gamma(K_S^* \to K \eta) \simeq \Gamma(K_S^* \to \rho K) \simeq \Gamma(K_S^* \to K^* \pi)$

$\qquad \simeq 100\ \mathrm{MeV}$,

$\Gamma(\phi_S \to K^* \bar{K} + cc) \simeq 2\Gamma(\phi_S \to K\bar{K}) \simeq 200\ \mathrm{MeV}$,

and

$\Gamma(\rho_D \to \pi\pi) \simeq 4\Gamma(\rho_D \to \omega\pi) \simeq 4\Gamma(\rho_D \to K\bar{K}) \simeq 100\ \mathrm{MeV}$,

$\Gamma(\omega_D \to \rho\pi) \simeq 4\Gamma(\omega_D \to K\bar{K}) \simeq 100\ \mathrm{MeV}$,

$\Gamma(K_D^* \to K\pi) \simeq \Gamma(K_D^* \to K\eta) \simeq 2\Gamma(K_D^* \to \rho K)$

$\qquad \simeq 2\Gamma(K_D^* \to \rho K) \simeq 75\ \mathrm{MeV}$,

and

$\Gamma(\phi_D \to K\bar{K}) \simeq 2\Gamma(\phi_D \to K^* \bar{K} + cc) \simeq 100\ \mathrm{MeV}$.

***

We now argue that the properties of these states are not only consistent with what is known empirically but also that they also could explain some puzzling features of the data. Let us consider first the isovector channel. From the simple decay modes of Table V above we have $\Gamma_{\rho_S} \gtrsim 300\ \mathrm{MeV}$ and $\Gamma_{\rho_D} \gtrsim 200\ \mathrm{MeV}$ so that after taking into account other modes it would not be surprising to find $\Gamma_{\rho_S} \approx 500\ \mathrm{MeV}$ and $\Gamma_{\rho_D} \approx 300 \ \mathrm{MeV}$. Recalling that their mass difference is only 200 MeV, it might thus be that the reported $\rho'$ is a mixture of the $\rho_D$ and $\rho_S$ weighted in favor of the narrower $\rho_D$; the fact that these states will interfere and will appear with different relative strengths in different experiments may account for the large range of masses reported for the $\rho'$ (from around 1450 to around 1700 MeV). For the $\omega$-like states we find from the simple decay modes $\Gamma_{\omega_S} \gtrsim 300\ \mathrm{MeV}$ while $\Gamma_{\omega_D} \gtrsim 150\ \mathrm{MeV}$, so taking into account possible widths to other modes it seems plausible that states reported in this channel (see, e.g., the report in Ref. 28 of an $\omega$-like state at 1.67 GeV) will be the $\omega_D$ with $\omega_S$ appearing, if at all, as ...

## PDF Page 31

32

a very broad background. When we consider the $\phi$-like states, however, this pattern is broken: here the modes of Table V should approximately exhaust the Zweig-allowed decays so we estimate $\Gamma_{\phi_S} \simeq 300$ MeV and $\Gamma_{\phi_D} \simeq 150$ MeV. We would therefore expect in this case to see $\phi_S$, and we would indeed argue for identifying this state with the $\phi'(1680)$: it has close to the predicted mass, and its width and branching ratios are consistent with expectation [e.g., $\phi'(1680)$ decays preferentially to $K^* \overline{K} + cc$ rather than $K\overline{K}$ as $\phi_S$ is supposed to do and contrary to the predictions for $\phi_D$]. Moreover, we note that this identification avoids the otherwise embarrassing interpretation of $\phi'(1680)$ as a partner of $\rho'(1600)$: empirically, an $ss$ state in this region should be about $200$ MeV heavier than its isovector partner. Finally, turning to the $K^*$-like states, we find $\Gamma_{K^*_S} \gtrsim 350$ MeV and $\Gamma_{K^*_D} \gtrsim 200$ MeV. The sightings of such a state range in mass from about 1500 to 1800 MeV with widths from about 200 to 400 MeV. This certainly seems consistent with two such broad overlapping resonances; the measurement of a $K\pi$ width of $\simeq 75$ MeV for this effect is also consistent with expectations. As a final bit of evidence in favor of this interpretation, we mention the fact that the pseudoscalar partners of the radially excited vector mesons are seen at approximately their expected masses. To summarize: we argue that there is already some indication that there are two $1\,^{--}$ nonets close together and often overlapping in the region normally associated with the radial excitation of the ground-state nonet, and that analysis of the data in this region should take this possibility into account.

The second puzzle in this sector which we would like to discuss concerns the pseudoscalar mesons $\eta,\,\eta'$, and their radial excitations. As discussed in Sec. II, there are reasons to believe that the annihilation amplitude in this sector receives large positive nonperturbative contributions for light mesons and that it changes sign for very heavy mesons. We have little faith in our ability to describe this complex situation, but, as an illustration of two possibilities, discuss the consequences of our models P1 (18a) and P2 (18b) here.

We begin by discussing the simple model P1 in which there are very strong positive annihilation amplitudes for masses below $m_\eta$ which die out in the $1-2$ GeV region. The result is a fairly standard picture of the $\eta$ and $\eta'$ as nearly "perfectly mixed" (see the first of Refs. 5) states although the $\eta'$ has in this model a significant admixture of the radial excitation of the $2^{-1/2}(u\overline{u} + d\overline{d})$ state. With a weaker but still very significant positive annihilation amplitude in the mass range of the $\pi'$, its isoscalar partner is heavier ($m_{\eta'}-m_{\pi'}\simeq 100$ MeV) and it is natural to identify this third isoscalar state with the $\iota(1440)$ seen in $\psi$ radiative decays.\textsuperscript{29} Since the mechanism of Cohen, Lipkin, and Isgur\textsuperscript{29} (in which the wave function at the origin is automatically shifted to higher states by a positive-annihilation amplitude) is in operation here, there is no problem with the branching ratios for $\psi$ radiative decays into these states: see Table III and also footnote d to Table VI which explains why the usual pole-model calculations of radiative decay rates based on mixing coefficients like those of Table III can only be used as a qualita-

tive guide in any event. This model predicts a further state should be seen in $\psi$ radiative decays at around 1.6 GeV. P1 is, however, incorrect if there really is a $q\overline{q}$ isoscalar in this channel at around 1.28 GeV as indicated by one experiment.\textsuperscript{30} Before discussing the alternative P2, which is a direct response to this possible difficulty with P1, we note that there are a number of conceivable problems with the interpretation of the $\pi^- p \to n\eta\pi\pi$ data of Ref. 30. These problems are mainly related to the fact that the observed $\eta\pi\pi$ peak, with which it was not possible to associate resonant phase motion, was seen in an isobar analysis in the $\delta\pi$ channel where it both lies directly under the $D(1285)$ and is superimposed on a reasonably large $0^-$ background. Another possible problem is associated with the peculiar nature of the $\delta$ (see Sec. V D below) which could distort the attempted $\delta\pi$ isobar analysis.

In the model P2 the annihilation mass matrix depends on the annihilation-channel mass so that its $i$th eigenvector must be found by iterating until the $i$th eigenvalue $M_i$ coincides with the input mass $M$. At low masses (18b) creates strong positive mixing between the $2^{-1/2}(u\overline{u} + d\overline{d})$ and $ss$ $1\,^1S_0$ states to produce once again a nearly "perfectly mixed" $\eta$. As $M$ approaches 1 GeV the annihilation weakens so that the $\eta'$ is half $1\,^1S_0$ $ss$ with components of both $1\,^1S_0$ and $2\,^1S_0$ nonstrange states. The third isoscalar state appears at 1.27 GeV, just above the annihilation zero $M_0$, and can now be identified with the $\eta(1.28)$,\textsuperscript{30} while the fourth isoscalar at 1.55 GeV might now be associated with the $\iota(1440)$. In P2, $\iota$ is a dominantly $2\,^1S_0$ $ss$ state which has been shifted down in mass by annihilation. (Note that this predicted state is split from its isovector companion by 0.25 GeV, compared to the $\iota$-$\pi'$ splitting of $0.14 \pm 0.11$ GeV; in view of our experience with the spectroscopy of the nonisoscalar $n\,^3S_1$—$n\,^1S_0$ splittings, we would consider our absolute prediction in this case to be acceptable). Table III shows that the $cc$ content of these states is, in view of the uncertainties just mentioned, once again acceptable.

These two pseudoscalar-mixing models are most sharply distinguished by the presence or absence of the state at around 1.28 GeV, but they can also be distinguished by the fact that in P1 the $\iota$ is mainly $n=2\ (u\overline{u}+d\overline{d})/\sqrt{2}$, while in P2 it is mainly $n=2\ ss$. Thus, for example, $\iota \to \rho\gamma$ would be expected to be large in P1 and small in P2 as indicated in Table VI.\textsuperscript{31}

The $\iota(1440)$ is widely considered as a candidate for a glueball.\textsuperscript{29} Our models, while certainly not compelling, illustrate that with our ignorance of annihilation forces in this channel it is possible to accommodate this state as a radial excitation. So long as such an interpretation of the $\iota(1440)$ is possible, its credentials as a glueball will probably remain questionable.

***

**B. The $P$ waves**

The spectroscopy of $P$-wave mesons is most clearly viewed in the $cc$ system. In perturbation theory (which is a good approximation in $cc$) we would find that

$$
E(2^{++}) = E_0 + \frac{1}{4}S - \frac{1}{5}T + L\ , \tag{23}
$$
$$
E(1^{++}) = E_0 + \frac{1}{4}S + T - L\ , \tag{24}
$$

## PDF Page 32

$$
E(0^{++}) = E_0 + \frac{1}{4}S - 2T - 2L \tag{25}
$$
$$
E(1^{+-}) = E_0 - \frac{3}{4}S , \tag{26}
$$

where $S$ arises from the contact term (and is small but nonzero due to relativistic smearing), $T$ arises from the tensor term, and $L$ from the spin-orbit interaction. Since the $1^{+-}$ state has not yet been found, we concentrate on the $C = +$ states from which we can deduce that $T = +20 \pm 2$ MeV and $L = +34 \pm 3$ MeV, which compare reasonably well with the calculated values of $+13$ MeV and $+28$ MeV, respectively. We stress that these values are not fit, but rather follow once the parameters of the spin-independent parts of the potential are determined. It is particularly worth noting that the small value of $L$ is a result of a strong cancellation between the color-magnetic spin-orbit interaction (6) and the Thomas-precession term (7).

The situation in light-quark $P$-wave mesons is less satisfactory because the relativistic modification factors $(m/E)^{1+2e_\alpha}$ play an important role. This means that although a strong cancellation between $L_{cm}$ and $L_{tp}$ certainly occurs here as well, the residual is less certain. We believe it would be extremely valuable to study these effects in a fully relativistic model to determine less phenomenologically the nature of the required modifications. In the meantime we can note that our values for the $e$'s are reasonable and have roughly the characteristics we anticipated in the discussion of Appendix A. Furthermore, our method of implementing these effects (which, after all, are certainly present in some form) can be tested in several ways. It would, for example, be encouraging to find that the $s\bar{s}$ spectrum with its considerably smaller values of $E/m$, is indeed correspondingly closer to having unmodified spin-orbit effects. It will also be interesting to see if the inversion of spin-orbit multiplets at high $l$ predicted$^{32}$ as a consequence of the color-magnetic–Thomas-precession cancellation, which our results confirm, actually occurs.

Despite these uncertainties, there seems to be little doubt that the main characteristics of the $P$-wave states are correctly given by the model after examining their decays:

(1) The successful description of the $E1$ transitions in charmonium and $b$-quarkonium indicates that the spatial wave functions are basically correct; they even provide reasonable evidence that the distortions expected from spin-orbit effects are present.$^{19}$

(2) The successful description of the decays of the light-quark $P$ waves indicates that the expected relations between these multiplets exist.

(3) The mixing in the $Q_1$-$Q_2$ system is automatically explained in sign and magnitude by the same forces that generate the spectroscopy,$^{32}$ although it appears that there may be a significant decay-channel mixing in operation here as well.$^{33}$

### C. Higher-$L$ excitations

Apart from the $\psi''$ (3770), which is mainly $^3D_1$, our knowledge of mesons with $L > 1$ comes from the SU(3) meson families, where we have the sequences $(\pi, B, A_3, \ldots)$, $(\rho, A_2, g, A_2^*, \ldots)$, and $(K, Q_B, L, \ldots)$,
$(K^*, K^*(1780), K^*(2060), \ldots)$. The overall quality of the model’s spectroscopic description of these states is good. Some of the details of this spectroscopy bear further comment however: (1) For a long-range spin-spin interaction, the $S=0$ state of any $L$ would be split from the center of gravity of its $S=1$ partners by (approximately) the $^1S_0$–$^3S_1$ ground-state splitting. Clearly, as mentioned already, the facts that the $B(1235)$ and $A_3(1680)$ are *not* 640 MeV below their $S=1$ counterparts indicates once again that the $S_i S_j$ piece of the hyperfine interaction is short range.
(2) The spacings of the light-quark excitations, falling as they do in the linear part of the potential, provide a strong constraint on the slope parameter $b$. It is impressive that this slope agrees with the one required to explain the spectroscopy of the radial excitations of the $\psi$ and $\Upsilon$ systems (see Fig. 12).

Once again, any questions about the model’s ability to describe the internal structure of orbital excitations would seem to be dispelled by the successful description the model affords of the couplings of these states. Especially impressive is the fact that the pionic decays of the whole natural-parity sequence can be described by the single decay partner $A$.

***

### D. The scalar mesons

As one would expect from the charmonium $P$ waves, our model predicts that the $^3P_0$ $q\bar{q}$ states lie (in each flavor sector) roughly 200 MeV below the $^3P_2$ states $A_2(1310)$, $f(1270)$, $K^*(1430)$, and $f'(1515)$. It is at first sight tempting to associate the two lowest such states with the nearly degenerate states $\delta(980)$ ($IJ^{PC} = 10^{++}$) and $S^*(980)$ ($IJ^{PC} = 00^{++}$), but such an assignment is very difficult to implement. Putting aside for the moment any discrepancies with mass predictions, we can probably rule out this assignment by noting that these states have predicted couplings which are inconsistent with their observed properties: the width of the $S^*$ would, if $S^*$ were a $2^{-1/2}(u\bar{u}+d\bar{d})\,^3P_0$ state, be nearly 1000 MeV (versus its measured total width of about 30 MeV) while if the $\delta$ were a $2^{-1/2}(u\bar{u}-d\bar{d})\,^3P_0$ state its width would be about 400 MeV (versus a measured total width of approximately 50 MeV). Clearly there is something wrong.

It seems to us very unlikely that the decay model can be blamed since it correctly predicts the closely related $S$-wave decays $B \to [\omega\pi]_{S}$ and $A_1 \to [\rho\pi]_{S}$. (Of course we cannot rule out the possibility of anomalous behavior in the channel with vacuum quantum numbers, but preliminary results from a study of more realistic decay models$^{16}$ confirm these conclusions.)

Especially in view of this glaring discrepancy, but also considering the spectroscopy and branching ratios of these states, the most attractive alternative seems to us to be the one given by Ref. 34, that the $\delta$ and $S^*$ are $K\bar{K}$ bound states. This identification, which was based on a variational calculation for the $qq\bar{q}\bar{q}$ system using a simplified form of our Hamiltonian, is closely related to the earlier suggestion that the $\delta$ and $S^*$ could be identified as two members of the low-lying "cryptoexotic" $qq\bar{q}\bar{q}$ $0^{++}$ nonet of the bag model.$^{35}$ According to Ref. 34, while basically confirming the bag-model proposal for these two

## PDF Page 33

states, the potential-model approach helps to strengthen the $qq\bar{q}\bar{q}$ interpretation of these states in several ways:

(1) The $qq\bar{q}\bar{q}$ wave function in the potential model naturally clusters into two weakly bound color-singlet mesons. This explains why the $\delta$ and $S^*$ are found just below $K\bar{K}$ threshold.

(2) The color hyperfine interaction plays a crucial role in the binding process. If the quarks are too heavy, this interaction is too weak to cause binding; if they are too light, then the hyperfine interactions create a strong intermeson attraction but the mesons become too light to be bound. Thus a complete cryptoexotic nonet is not expected and in Ref. 34 it was argued that probably only the $K\bar{K}$ system binds. The potential model thus avoids the prediction not only of unseen nonet partners to the $S^*$ and $\delta$ but also, perhaps, of $qq\bar{q}\bar{q}$ resonances in all other channels (for which there is no apparent place in the spectrum).

(3) The correct absolute widths of the $\delta$ and $S^*$ emerge from the $K\bar{K}$ bound-state picture. The decay of the $S^*$ to $\pi\pi$, for example, occurs in this view as a consequence of the inelastic collision $K\bar{K} \rightarrow \pi\pi$ of the two weakly bound kaons. The predicted narrow widths $\Gamma(S^*\rightarrow\pi\pi)\simeq 15$ MeV and $\Gamma(\delta\rightarrow\eta\pi)\simeq 40$ MeV are in fact a consequence of the weak binding since the $K\bar{K}$ wave function is so large as to make such collisions improbable. These predicted narrow widths are in marked contrast to the $q\bar{q}$ interpretation, as already mentioned, but they may also be contrasted with those of the bag-model interpretation where $\delta\rightarrow\eta\pi$ would be a fast "fall-apart" mode.

We propose to accept the $qq\bar{q}\bar{q}$ potential-model calculation which identifies the $\delta(980)$ and $S^*(980)$ as weakly bound $K\bar{K}$ states and to argue that together with our predicted $q\bar{q}$ $^3P_0$ nonet (and their radial excitations) these states may provide an adequate description of the complicated experimental situation in $0^+$ channels.$^{36-38}$

Let us begin by considering the $IJ^{PC}=00^{++}$ channel of $\pi\pi$ and $K\bar{K}$ scattering.$^{37}$ The states which may be expected to contribute to this channel are listed in Table VIII along with some of their properties. If we neglect for the moment any channels except $\pi\pi$ and $K\bar{K}$, then the $S$ matrix must be of the form

***

TABLE VIII. States contributing to $I=0$, S-wave, $\pi\pi$ and $K\bar{K}$ scattering. Note that we have not considered here the possibility of mixing between these states via annihilation.

| States | (mass in MeV) | Widths on resonance | Comments |
|--------|---------------|--------------------|----------|
| $S^*(980)$ |  | $\Gamma \simeq 15$ MeV <br> $x_{\pi\pi} \simeq 1$ | $K\bar{K}$ bound state |
| $\epsilon(1090)$ |  | $\Gamma \simeq 850$ MeV <br> $x_{\pi\pi} \simeq 0.86$ <br> $x_{K\bar{K}} \simeq 0.14$ | $\left(\frac{1}{2}\right)^{1/2}(u\bar{u}+d\bar{d}) \ ^3P_0$ |
| $\epsilon'(1360)$ |  | $\Gamma \simeq 600$ MeV <br> $x_{\pi\pi} \simeq 0$ <br> $x_{K\bar{K}} \simeq 0.8$ <br> $x_{\eta\eta} \simeq 0.2$ | $s\bar{s} \ ^3P_0$ |
| $\epsilon_r(1780)$ |  | $\Gamma \simeq 1000$ MeV? <br> $x_{\pi\pi} \simeq 3x_{K\bar{K}} \ll 1$ | radial excitation of $\epsilon(1090)$ |

***

$$
S = 
\begin{bmatrix}
\eta e^{2i\delta_{\pi\pi}} & i(1-\eta^2)^{1/2} e^{i(\delta_{\pi\pi} + \delta_{K\bar{K}})} \\
i(1-\eta^2)^{1/2} e^{i(\delta_{\pi\pi} + \delta_{K\bar{K}})} & \eta e^{2i\delta_{K\bar{K}}}
\end{bmatrix}
\tag{27}
$$

and $T = (1/2i)(S-1)$. For a single Breit-Wigner resonance we would have

$$
T_{ij} = (x_i x_j)^{1/2} \left[ \frac{\Gamma/2}{M-E-i\Gamma/2} \right] = (x_i x_j)^{1/2} e^{i\delta_{BW}(E)} \sin \delta_{BW}(E) = a_{ij} e^{i\phi_{ij}},
\tag{28}
$$

where $x_i = \Gamma_i / \Gamma$ is the branching ratio to channel $i$. Using these results and the values of Table VIII as input, we show in Fig. 16 the effect that each of the states in the table would have on the coupled-channel $T$ matrix *if they acted alone*. Note that in this limit

$$
\eta = [1 - 4x_{\pi\pi}(E)x_{K\bar{K}}(E)\sin^2\delta_{BW}(E)]^{1/2}.
\tag{29}
$$

We have also shown in the figure the effect on the scattering of the potential $V$ which binds the $S^*$.

We believe that these curves of individual contributions make plausible the eventual explanation of the $IJ^{PC} = 00^{++}$ channel in terms of these elements. Of

![Figure 16](attachment:image3)

**FIG. 16.** Some effects in $IJ^{PC} = 00^{++}$ scattering.

## PDF Page 34

course to demonstrate that such an explanation is correct, one must simultaneously consider the effects of $V$, $S^*$, $\epsilon$, $\epsilon'$, and $\epsilon_r$ (as well as other possible inelastic modes like $\eta \eta$ and multibody final states). There is, _a priori_, no unique way of doing this: there are an infinity of ways in which unitarity can be implemented in such a system, and only by understanding the physics of the interplay of the contributing effects can the correct method be chosen. Let us try to illustrate this point by going below $K\bar{K}$ threshold and considering $\pi\pi$ scattering in the $S^*, \epsilon$ system. From the physical picture of the $S^*$ as a $K\bar{K}$ bound state, it is most natural to suppose that these two states interact (in first order) only very indirectly via the chain $\pi\pi \rightarrow \epsilon \rightarrow K\bar{K} \rightarrow S^* \rightarrow K\bar{K} \rightarrow \epsilon \rightarrow \pi\pi$. Obviously such an interaction can lead to different results from one where, for example, both contributing states couple directly to $\pi\pi$ so that chains like $\pi\pi \rightarrow R_1 \rightarrow \pi\pi \rightarrow R_2 \rightarrow \pi\pi$ can occur. We hope that by carefully considering dynamics of this sort it may be possible to substantiate an explanation of the $IJ^{PC}=00^{++}$ channel along the lines we have outlined here.

We next turn to the $IJ^{PC}=10^{++}$ channel which may be expected to receive contributions from the states listed in Table IX. States in this channel have for the most part been seen in production experiments where the very broad $\delta_2(1090)$ might be lost (see Fig. 17). Since these states can only be produced in $\pi N$ reactions by $\eta$, $B, \ldots$ exchange, they are considerably more difficult to study than their isoscalar companions. There are some indications, however, that the $\delta(980)$ behaves much like the $S^*(980)$ in the $K\bar{K}$ channel, lending further support to an interpretation of both as $K\bar{K}$ bound states.$^{37,38}$

Finally, we consider the $IJ^P = \frac{1}{2} 0^{+}$ channel with the states of Table X contributing. Figure 18 shows the contributions of these states along with the potential term $V$ analogous to that present in $K\bar{K}$ scattering, but too weak to produce a bound state. To emphasize our earlier comments, especially in view of our prediction that the $^3\!P_0$ $\kappa$ lies considerably below the apparent resonance in this channel, we make a further simple comment here on the question of whether the three effects noted in the figure might reproduce the observed phase shifts. In the region around 1250 MeV where $\phi \simeq \delta_\kappa \simeq 90^\circ$, unitarity could be preserved either by adding to $T_\kappa$ some small T-matrix element $\Delta T$ (from $V$ or $\kappa'$) or by adding to the $\kappa$ phase shift

***

TABLE IX. States contributing to $I=1$, S-wave, $\eta\pi$ and $K\bar{K}$ scattering.

| State (mass in MeV) | Widths on resonance | Comments |
|---------------------|--------------------|----------|
| $\delta(980)$       | $\Gamma \simeq 40$ MeV, $x_{\eta\pi} \simeq 1$    | $K\bar{K}$ bound state |
| $\delta_2(1090)$    | $\Gamma \simeq 400$ MeV, $x_{\eta\pi} \simeq 0.6$, $x_{K\bar{K}}\simeq 0.4$ | $q\bar{q}$ $^3P_0$ |
| $\delta_2'(1780)$   | $\Gamma \simeq 500$ MeV?, $x_{\eta\pi} \simeq x_{K\bar{K}} \ll 1$ | radial excitation of $\delta_2(1090)$ |

***

![FIG. 17](fig17.png)

**FIG. 17.** Incoherent production rate (in arbitrary units) of $\delta(980)$ and $\delta_2(1090)$ in the ratio 1:5.

***

$\delta_\kappa$ the corresponding $\Delta \delta$ (from $V$ or $\kappa'$). These two methods lead to $\delta = \delta_\kappa + \Delta\delta$, respectively. The correct result, as remarked earlier, can only be decided on the basis of the physical mechanisms by which the various effects interact.

Of course, in all three of these $IJ^P$ channels we are dealing with resonances which would be among the widest known hadrons, so there is also the possibility that they have suffered larger than normal shifts in their masses from the narrow-resonance approximation we have assumed. Even without this uncertainty, however, we believe that our interpretation of these states is, for the moment at least, viable.

***

### E. The missing resonances

As was the case with the related analysis of baryon decays,$^{11}$ an examination of the predicted couplings of missing meson resonances usually provides a reasonable understanding of why they have not yet been found, as well as in most cases providing guidance in how a search for them might most profitably be made. While such an analysis is beyond the scope of the present work, we shall mention, by way of illustration, two examples:

(1) The $A_3$ is very prominent because it has a predicted (and observed) width of about 100 MeV to the simple

***

![FIG. 18](fig18.png)

**FIG. 18.** Some effects in $IJ^P = \frac{1}{2} 0^{+}$ scattering.

## PDF Page 35

TABLE X. States contributing to $I = \frac{1}{2}$, S-wave, $K\pi$ scattering.

| State (mass in MeV) | Widths on resonance | Comments |
|---------------------|--------------------|----------|
| $\kappa(1240)$      | $\Gamma \simeq 400$ MeV  | $qq\, ^3P_0$ |
|                     | $x_{K\pi}\simeq 1$       |             |
| $\kappa'(1890)$     | $\Gamma \sim 500$ MeV?   | radial excitation |
|                     | $x_{K\pi} \ll 1$         | of $\kappa(1240)$ |

channel $\rho\pi$. Its isoscalar partner, on the other hand, has as its only available simple decay channel $K^* \bar{K} + \bar{K}^* K$ for which it has a width of only about 20 MeV. Clearly it would have been difficult to find such an inelastic state by bump-hunting, and a high statistics isobar-type analysis will be required to see this state in this mode. (It is also possible that this state may be prominent in a mode like $A_2\pi$ but a more sophisticated decay analysis will be required to delineate such possibilities.)

(2) The $I=1\ ^{13}D_2$ state along with all its nonet partners (with the possible exception of the $u\bar{s}$ state), is missing. In this case the state is predicted to have a prominent decay mode to $[\omega\pi]_P$ with a width of about 50 MeV and so could be seen by an experiment sensitive to this channel. The main uncertainty in examining these sorts of questions is the total width of a given high-mass state: even a state with a healthy width to some simple mode will probably be lost if its total width is much above 400 MeV. Since highly excited states would in fact be expected to have strong decay modes to other excited states, a more sophisticated decay model than the one we have used here will once again be required to make any study of missing resonances complete.

***

$^{3}P_2$ $\quad$—$\quad x_b\quad$—$\quad x_c\quad$—$\quad f'$ $\quad$—$\quad K^{**}$—$\quad A_2$
             
$^{3}S_1$ — T — $\psi$ — $\phi$ — $K^*$ — $\rho$  
$^{1}S_0$ - - - - - $\eta_b$ — $\eta_c$
         - - - - - $\eta, \eta'$
                         — K
                                             - π

FIG. 19. A graphic illustration of the universality of meson dynamics from the $\pi$ to the $\Upsilon$, showing the splittings of $^{3}P_2$ and $^{1}S_0$ from $^{3}S_1$ in the $b\bar{b}$, $c\bar{c}$, $s\bar{s}$, $u\bar{s}$, and $u\bar{d}$ families.

***

VI. CONCLUSIONS AND COMMENTS

In our view, this attempt to construct a quark model with chromodynamics for mesons has been reasonably successful. Although if one examines the ability of the model to describe any particular aspect of meson structure it is usually less accurate than more *ad hoc* schemes, we believe its inaccuracy is more than offset by its unity and breadth of application (see Fig. 19).

On the negative side, perhaps the most disappointing feature of our model is that to cope with relativistic effects we have been forced to introduce a new set of parameters which we should in principle be able to calculate. While the concomitant uncertainties have little effect in heavy-quark sectors, we believe that future progress in light-quark spectroscopy within this framework lies in the direction of handling these effects more rigorously and more economically. Nevertheless, we believe that we have managed to take into account semiquantitatively the main relativistic modifications to the nonrelativistic quark model and to demonstrate their importance to a unified picture of meson structure. A second shortcoming of the model, which it shares with other attempts, is its weakness in isoscalar channels where annihilation effects can be important. As discussed in the text, this problem is especially severe in the pseudoscalar channel. Of course, one can concentrate on channels with net flavor, but in practical terms this is an unpleasant limitation: over a third of all known light mesons are isoscalars. It remains to be seen, on the other hand, whether the model's inability to provide an immediate explanation for some states, like the possible resonance $\theta(1660)$ seen in $\psi$ radiative decays, is a shortcoming. One exciting possibility is, of course, that states which cannot be accommodated in this scheme will prove to be new kinds of hadronic objects, like pure glue states or hybrids.

Finally, despite our reservations, which are more extensive than the main ones we have just mentioned, we find it difficult to doubt the basic validity of the model. From their masses and quantum numbers to their static properties, electromagnetic, weak, and strong decays, it seems to us to consistently provide a recognizable (though at times very coarse) portrait of the whole meson family.$^{39}$

***

**ACKNOWLEDGMENTS**

This work owes a great deal to many people. Without the computing power made available to us by George Luste and Anna Pezacki it could not have been done at all. Bruce Campbell, Simon Capstick, Cameron Hayne, Rick Kokoski, Kim Maltman, and John Weinstein at Toronto all engaged actively in discussions of the ongoing work. R. L. Jaffe, G. Karl, H. J. Lipkin, R. Longacre, and J. Paton all made important comments and suggestions. We would also like to acknowledge L. Montanet for suggesting the problem in his rapporteur's talk at Madison 1980. This research was supported in part by a grant from the Natural Sciences and Engineering Research Council of Canada.

***

**APPENDIX A: BUILDING THE EFFECTIVE $q\bar{q}$ POTENTIAL**

We arrive at the effective $q\bar{q}$ potential of the text by first constructing a potential that reproduces on-shell $q\bar{q}$ scattering amplitudes in the center-of-mass frame, motivated by simple arguments based on QCD.$^{40}$ The two key ingredients of the assumed interaction are a short-range $\gamma^{\mu} \otimes \gamma_{\mu}$ interaction

## PDF Page 36

224                                                                                                        32

STEPHEN GODFREY AND NATHAN ISGUR

$$
G(Q^2) = \frac{-4\alpha_s(Q^2)}{3} \frac{4\pi}{Q^2}
$$

and a long-range $1\otimes 1$ linear confining interaction $S(Q^2)$ suggested by lattice QCD calculations.$^{41}$ We include the effects of asymptotic freedom by using the Ansatz (12) for $\alpha_s(Q^2)$.

By this method we find for $|\mathbf{p}| = |\mathbf{p}'| \equiv p$ that

$$
\chi_s^{\dagger} \chi_{\bar{s}}^{\dagger} V_{\text{eff}}(\mathbf{P}, \mathbf{r}) \chi_s \chi_{\bar{s}}
$$
$$
= \frac{1}{(2\pi)^3} \int d^3Q \, e^{i\mathbf{Q} \cdot \mathbf{r}} \overline{U} (\mathbf{p}', s') \overline{V}(-\mathbf{p}, \bar{s}) I(Q^2)
\times U(\mathbf{p}, s) V(-\mathbf{p}', \bar{s}'), 
\tag{A1}
$$
where

$$
I(Q^2) = G(Q^2)(\gamma^{\mu})_q (\gamma_{\mu})_{\bar{q}} - S(Q^2)(1)_q (1)_{\bar{q}}
\tag{A2}
$$

with $\mathbf{Q} = \mathbf{p}' - \mathbf{p}$ and $\mathbf{P} = (\mathbf{p} + \mathbf{p}')/2$. With

$$
f(\mathbf{P}, \mathbf{r}) \equiv \frac{1}{(2\pi)^3} \int d^3 Q \, e^{i \mathbf{Q}\cdot\mathbf{r}} f(\mathbf{P}, \mathbf{Q})
\tag{A3}
$$

we find

$$
V_{\text{eff}}(\mathbf{P}, \mathbf{r}) = G_{\text{eff}}(\mathbf{P}, \mathbf{r}) + S_{\text{eff}}(\mathbf{P}, \mathbf{r}),
\tag{A4}
$$

where

$$
\begin{aligned}
G_{\text{eff}}(\mathbf{P}, \mathbf{r}) = {}& \frac{1}{(2\pi)^3} \int d^3 Q \, e^{i \mathbf{Q} \cdot \mathbf{r}} G(Q^2)
\left\{
\left[
1 - \frac{Q^2}{4E(E+m)} + \frac{i\mathbf{Q} \times \mathbf{P} \cdot \mathbf{S}_q}{E(E+m)}
\right]
\left[
1 - \frac{Q^2}{4\overline{E}(\overline{E}+\overline{m})} + \frac{i\mathbf{Q} \times \mathbf{P} \cdot \mathbf{S}_{\bar{q}}}{\overline{E}(\overline{E}+\overline{m})}
\right]
\right. \\
&\left.
+ \left[
\frac{\mathbf{P} - 2i\mathbf{Q} \times \mathbf{S}_q}{2E}
\right]
\left[
\frac{\mathbf{P} - 2i\mathbf{Q} \times \mathbf{S}_{\bar{q}}}{2\overline{E}}
\right]
\right\}
\tag{A5}
\end{aligned}
$$

and

$$
S_{\text{eff}}(\mathbf{P},\mathbf{r}) = \frac{1}{(2\pi)^3} \int d^3 Q \, e^{i \mathbf{Q} \cdot \mathbf{r}} S(Q^2)
\left[
\frac{m}{E} + \frac{Q^2}{4E(E+m)} - \frac{i\mathbf{Q} \times \mathbf{P} \cdot \mathbf{S}_q}{E(E+m)}
\right]
\left[
\frac{\overline{m}}{\overline{E}} + \frac{Q^2}{4\overline{E}(\overline{E}+\overline{m})} - \frac{i\mathbf{Q} \times \mathbf{P} \cdot \mathbf{S}_{\bar{q}}}{\overline{E}(\overline{E}+\overline{m})}
\right].
\tag{A6}
$$

where $m$ and $\overline{m}$ are the quark and antiquark masses, $\mathbf{S}_q$ and $\mathbf{S}_{\bar{q}}$ are their spins, and $E = (p^2 + m^2)^{1/2}$, $\overline{E} = (p^2 + \overline{m}^2)^{1/2}$.

For on-shell $q\bar{q}$ scattering at c.m. momentum $p$ these results are exact, i.e., these potentials will reproduce exactly the scattering amplitudes we have assumed. There are several reasons, however, for not using these potentials directly in (1). The first is that as they stand they do not adequately reflect the full expected momentum dependence of the potential $V(\mathbf{P}, \mathbf{r})$ in (1) since it will in general have off-energy-shell behavior not considered in (A1). In field theory the Schrödinger equation (1) arises in the $q\bar{q}$ sector of Fock space by integrating over more complex components of Fock space such as $| q\bar{q}g \rangle$, and this integration will introduce additional $\mathbf{P}$ dependence in $V$ not seen in (A1). Related to this deficiency is the fact that potentials like (A5) and (A6) have the usual ambiguity in the ordering of the classical quantities $E$ and $\overline{E}$ into a quantum operator: a matrix element involving $E$ will in general depend on $\mathbf{p}$ and $\mathbf{p}'$. Another reason for not using the above potentials without modification has specifically to do with (A6). It seems to us very unlikely that the confinement potential is a simple $1 \times 1$ interaction: the picture that emerges from studies of lattice QCD indicates that it is spin-independent, but that it arises from a distortion of a Coulomb interaction. If in fact we use one-dimensional QED as a guide, then we would expect the confinement potential to not only be spin-independent, but also $p$-independent.

We respond to this situation by treating (A5) and (A6) as a framework on which to build a semiquantitative model of relativistic effects. We roughly classify these effects into three categories: (a) the strengths of the various interactions will depend on the c.m. momentum of the interacting quarks, (b) the interactions will, since they depend on both $\mathbf{P}$ and $\mathbf{Q}$, be nonlocal, and (c) the interactions will, through $\mathbf{Q}$ dependence, take on new $\mathbf{r}$ dependences. Based on (b) and (c) we introduce a smearing function for a meson $q_i \bar{q}_j$

$$
\rho_{ij}(\mathbf{r}-\mathbf{r}') = \frac{\sigma_{ij}^3}{\pi^{3/2}} e^{-\sigma_{ij}^2(\mathbf{r}-\mathbf{r}')^2}
\tag{A7}
$$

which we apply to our basic potentials $G(r)$ and $S(r)$ to obtain smeared potentials $\widetilde{G}(r)$ and $\widetilde{S}(r)$ via

$$
\widetilde{f}_{ij}(r) \equiv \int d^3 r' \, \rho_{ij}(\mathbf{r}-\mathbf{r}') f(r')
\tag{A8}
$$
with the prescription

$$
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
$$

## PDF Page 37

where $\sigma_0$ and $s$ are the universal parameters given in Table II. The parameter $\sigma_0$ reflects the fact that in a confined system the smearing must be limited, while $s$ is the coefficient of the expected linear relation $\sigma_{Q\bar{Q}} = sm_Q$ for a heavy $Q\bar{Q}$ system. The complicated $m_i, m_j$ dependence of the $\sigma_0$ term in (A9), designed to reflect the fact that in a $Q\bar{q}$ system the light quark is more relativistic than in $q\bar{q}$, is significant mainly for pseudoscalar mesons. (If, however, the coefficient of $\sigma_0^2$ were replaced by unity these states would shift by only of order 30 MeV.) Using (12) for the running coupling constant and Fourier transforming leads to

$$
G(r) = -\sum_k \frac{4\alpha_k}{3r} \left[ \frac{2}{\sqrt{\pi}} \int_0^{\gamma_k r} e^{-x^2} dx \right] ,
\tag{A10}
$$

and

$$
S(r) = br + c .
\tag{A11}
$$

The result (13) follows from the definition

$$
G(r) = -\frac{4\alpha_s(r)}{3r} .
$$

The resulting smeared potentials are

$$
\tilde{G}(r) = -\sum_k \frac{4\alpha_k}{3r} \left[ \frac{2}{\sqrt{\pi}} \int_0^{\tau_{kij} r} e^{-x^2} dx \right]
\tag{A12}
$$

and

$$
\tilde{S}(r) = br
\left[
\frac{e^{-\sigma_{ij}^2 r^2}}{\sqrt{\pi}\sigma_{ij} r}
+ \left( 1 + \frac{1}{2\sigma_{ij}^2 r^2} \right) \frac{2}{\sqrt{\pi}} \int_0^{\sigma_{ij} r} e^{-x^2} dx
\right] + c ,
\tag{A13}
$$

where

$$
\frac{1}{\tau_{kij}^2} = \frac{1}{\gamma_k^2} + \frac{1}{\sigma_{ij}^2} .
\tag{A14}
$$

On the other hand, we take into account the effect (a) by introducing momentum-dependent factors in the various interactions which go to unity in the nonrelativistic limit to give back the potentials (3)-(7) of the text. Since (A12) and (A13) should already contain the $Q^2$-dependence-induced modifications of the form of the potentials, we examine (A5) and (A6) as $Q^2 \rightarrow 0$ and conclude the following.

(1) The Coulomb term should be modified according to

$$
\tilde{G}(r) \rightarrow \left[1 + \frac{p^2}{EE}\right]^{1/2} \tilde{G}(r) \left[1 + \frac{p^2}{EE}\right]^{1/2} .
$$

(2) The contact, tensor, vector spin-orbit, and scalar spin-orbit potentials should be modified according to

$$
\frac{\tilde{V}_i(r)}{m_1 m_2} \rightarrow \left(\frac{m_1 m_2}{E_1 E_2}\right)^{1/2 + \epsilon_i} \frac{\tilde{V}_i(r)}{m_1 m_2} \left(\frac{m_1 m_2}{E_1 E_2}\right)^{1/2 + \epsilon_i},
$$

where $i=$ contact $(c)$, tensor $(t)$, vector spin-orbit [so$(v)$], scalar spin-orbit [so$(s)$]. If $\epsilon_i=0$ then these modifications have the effect of replacing the nonrelativistic mass dependences $1/m_\alpha m_\beta$ $(=1/m^2,~ 1/\bar{m}^2,~ {\rm or}~ 1/m\bar{m})$ of these potentials by $1/E_\alpha E_\beta$. The parameters $\epsilon_i$, which are therefore expected to be small, are given in Table II.

(3) On the basis of the solutions to the (at least superficially) similar example of QED in one dimension, which confines with a linear potential, we assume that $\tilde{S}(r)$ is unmodified by relativistic corrections.$^{41}$ We remark that if we were to introduce a parameter $\epsilon_{\text{linear}}$ for this potential, we would be forced phenomenologically to set $\epsilon_{\text{linear}} \simeq0$; this observation offers some indirect support for the one-dimensional flux-tube model of confinement in QCD.

Since both experimentally and, partly as a consequence of $m/E$ suppressions, theoretically the spin-orbit interactions are relatively weak, we ignore the "second-order" spin-orbit terms of the form $(\mathbf{Q} \cdot \mathbf{p} \times \mathbf{S}_q)(\mathbf{Q} \cdot \mathbf{p} \times \mathbf{S}_q)$ in both $G_{\rm eff}$ and $S_{\rm eff}$. Reverting (as allowed by the resulting symmetry) to the 1,2 labeling of the text, this leaves us with the approximate forms of the potentials which we use in our calculations. Defining for compactness

$$
f^{i}_{\alpha \beta}(r) = \left(\frac{m_\alpha m_\beta}{E_\alpha E_\beta}\right)^{1/2 + \epsilon_i} f(r) \left(\frac{m_\alpha m_\beta}{E_\alpha E_\beta}\right)^{1/2 + \epsilon_i}
$$

we have

$$
\begin{aligned}
G_{\rm eff}(r) = {}& \left[ 1 + \frac{p^2}{E_1 E_2} \right]^{1/2} \tilde{G}(r) \left[ 1 + \frac{p^2}{E_1 E_2} \right]^{1/2} \\
&+ \left[
\frac{\mathbf{S}_1 \cdot \mathbf{L}}{2 m_1^2} \frac{1}{r} \frac{\partial \tilde{G}_{11}^{\rm so(v)}}{\partial r}
+ \frac{\mathbf{S}_2 \cdot \mathbf{L}}{2 m_2^2} \frac{1}{r} \frac{\partial \tilde{G}_{22}^{\rm so(v)}}{\partial r}
+ \frac{(\mathbf{S}_1 + \mathbf{S}_2) \cdot \mathbf{L}}{m_1 m_2} \frac{1}{r} \frac{\partial \tilde{G}_{12}^{\rm so(v)}}{\partial r}
\right] \\
&+ \frac{2\mathbf{S}_1 \cdot \mathbf{S}_2}{3 m_1 m_2} \nabla^2 \tilde{G}_{12}^{c}
- \left[
\frac{\mathbf{S}_1 \cdot \hat{\mathbf{r}}\,\mathbf{S}_2 \cdot \hat{\mathbf{r}} - \frac{1}{3} \mathbf{S}_1 \cdot \mathbf{S}_2}{m_1 m_2}
\right]
\left[
\frac{\partial^2}{\partial r^2} - \frac{1}{r} \frac{\partial}{\partial r}
\right]
\tilde{G}_{12}^{t}
\tag{A15}
\end{aligned}
$$

## PDF Page 38

and
$$
S_{\rm eff}(r) = \tilde{S}(r) - \frac{\mathbf{S}_1 \cdot \mathbf{L}}{2 m_1^2} \frac{1}{r} \frac{\partial \tilde{S}_{11}^{\rm so(s)}}{\partial r} - \frac{\mathbf{S}_2 \cdot \mathbf{L}}{2 m_2^2} \frac{1}{r} \frac{\partial \tilde{S}_{22}^{\rm so(s)}}{\partial r} .
\tag{A16}
$$

To actually perform calculations with these potentials we diagonalize the Hamiltonian matrix obtained from (1) in a (large) harmonic-oscillator basis. The harmonic basis is particularly useful in this context. For example, to take a matrix element of an operator of the form $f(p)g(r)$ we use

$$
\langle i | f(p)g(r) | j \rangle = \sum_n \langle i | f(p) | n \rangle \langle n | g(r) | j \rangle .
\tag{A17}
$$

The matrix elements $\langle i | f(p) | n \rangle$ are then calculated in momentum space, while the matrix elements $\langle n | g(r) | j \rangle$ are calculated in configuration space. Since for a harmonic-oscillator basis these two sets of wave functions are simple polynomials of the same form, the calculations are especially simple. Once the elements of the Hamiltonian matrix are calculated this way, the choice of the Gaussian parameter $\beta$ which characterizes the harmonic-oscillator basis is optimized in accordance with the variational principle. (Note that the best energy for a given state can always be obtained by minimizing it with respect to $\beta$. Of course to generate an orthogonal set of wave functions in a given sector, a single value of $\beta$ must be used, and in practice we take the value that minimizes the energy of the last state of the set for this purpose. Our basis is so large that very little error is introduced in this approximation.)

**APPENDIX B: WAVE-FUNCTION CONVENTIONS**

Here we make all of our wave-function conventions explicit so that our results may be more readily used. First, since we use the “natural” SU(3) conventions for quark and antiquark transfer operators ($q_i \rightarrow q_j$ is always $+1$ and $\bar{q}_i \rightarrow \bar{q}_j$ is always $-1$) our flavor wave functions, referred to the names of the pseudoscalar octet but generally applicable, are

$$
\pi^+ = -u \bar{d}, \tag{B1}
$$

$$
\pi^0 = \frac{1}{\sqrt{2}} (u\bar{u} - d\bar{d}), \tag{B2}
$$

$$
\pi^- = d\bar{u}, \tag{B3}
$$

$$
K^+ = -u\bar{s}, \tag{B4}
$$

$$
K^0 = -d\bar{s}, \tag{B5}
$$

$$
\bar{K}^0 = -s\bar{d}, \tag{B6}
$$

$$
K^- = s\bar{u}, \tag{B7}
$$

$$
\eta_8 = \frac{1}{\sqrt{6}} (u\bar{u} + d\bar{d} - 2s\bar{s}), \tag{B8}
$$

$$
\eta_1 = \frac{1}{\sqrt{3}} (u\bar{u} + d\bar{d} + s\bar{s}) \tag{B9}
$$

which satisfy the de Swart conventions on phases. For ideally mixed isoscalar mesons we take

$$
M_{ns} = \frac{1}{\sqrt{2}} (u\bar{u} + d\bar{d}), \tag{B10}
$$

$$
M_s = s\bar{s}, \tag{B11}
$$

and normally define mixing angles $\phi$ relative to this basis via

$$
M = M_{ns} \cos\phi - M_s \sin\phi, \tag{B12}
$$

$$
M' = M_s \cos\phi + M_{ns} \sin\phi, \tag{B13}
$$

so that $\phi = \theta_{{SU(3)}} - \theta_{{\rm ideal}}$, where $\theta_{\rm ideal} \simeq 35.3^\circ$. (Note that for $\theta_{SU(3)} \rightarrow 0$ we get $M \rightarrow + “\eta_1”, M' \rightarrow - “\eta_8”$ as a consequence of our conventions.) In the special case of the pseudoscalar mesons, we will often use the “perfect-mixing” states (see the first of Refs. 5),

$$
\eta = \frac{1}{\sqrt{2}} (M_{ns} - M_s), \tag{B14}
$$

$$
\eta' = \frac{1}{\sqrt{2}} (M_{ns} + M_s), \tag{B15}
$$

which correspond to an SU(3) mixing angle of $\theta_{\rm ideal} - 45^\circ \simeq -10^\circ$. (Note that these states follow from (B12) and (B13) by taking $\phi = -45^\circ$ so that $M \rightarrow \eta' \simeq \eta_1$, $M' \rightarrow -\eta \simeq -\eta_8$.) For heavy-quark mesons we ignore all symmetries except isospin and simply take our state vectors to be $+|Q\bar{q}\rangle$, where $Q$ is the heavy quark and $q$ any other quark except the $d$, including $Q$ itself. We use $-|Q\bar{d}\rangle$ so that $( -|Q\bar{d}\rangle, |Q\bar{u}\rangle )$ form an isospin multiplet analogous to $(\bar{K}^0, K^-)$.

The SU(3) flavor wave functions lead to coupling operators $X_q^i$ in (19) given by

$$
X_q^{\pi^+} = -\sqrt{2}\left[ \frac{\lambda_1 - i \lambda_2}{2} \right] (u \rightarrow -\sqrt{2}d ), \tag{B16}
$$

$$
X_q^{\pi^0} = +\lambda_3~(u \rightarrow u, d \rightarrow -d ), \tag{B17}
$$

$$
X_q^{\pi^-} = +\sqrt{2} \left[ \frac{\lambda_1 + i \lambda_2}{2} \right] (d \rightarrow \sqrt{2}u ), \tag{B18}
$$

$$
X_q^{K^+} = -\sqrt{2} \left[ \frac{\lambda_4 - i \lambda_5}{2} \right] (u \rightarrow -\sqrt{2}s ), \tag{B19}
$$

$$
X_q^{K^0} = -\sqrt{2} \left[ \frac{\lambda_6 - i \lambda_7}{2} \right] (d \rightarrow -\sqrt{2}s ), \tag{B20}
$$

$$
X_q^{\bar{K}^0} = -\sqrt{2} \left[ \frac{\lambda_6 + i \lambda_7}{2} \right] (s \rightarrow -\sqrt{2}d ), \tag{B21}
$$

## PDF Page 39

$$
X_q^{K^-} = + \sqrt{2} \left[ \frac{\lambda_4 + i \lambda_5}{2} \right] (s \rightarrow \sqrt{2}u) , \tag{B22}
$$

$$
X_q^{\eta} = \left[ \frac{\sqrt{2}+1}{\sqrt{6}} \right] \lambda_8 + \left[ \frac{\sqrt{2}-1}{\sqrt{6}} \right] \left[\frac{2}{3}\right]^{1/2} 1
$$

$$
\left[ u \rightarrow \frac{u}{\sqrt{2}}, d \rightarrow \frac{d}{\sqrt{2}}, s \rightarrow -s \right] , \tag{B23}
$$

$$
X_q^{\eta'} = \left[ \frac{\sqrt{2}+1}{\sqrt{6}} \right] \left[ \frac{2}{3} \right]^{1/2} 1 - \left[ \frac{\sqrt{2}-1}{\sqrt{6}} \right] \lambda_8
$$

$$
\left[ u \rightarrow \frac{u}{\sqrt{2}}, d \rightarrow \frac{d}{\sqrt{2}}, s \rightarrow s \right] , \tag{B24}
$$

where $\eta$ and $\eta'$ are defined in (B14) and (B15) and the $\lambda$'s are the usual Gell-Mann matrices. The operators $X_{\bar{q}}^{i}$ are given by the same formulas but with
$$
\lambda_i \rightarrow - \lambda_i^*
\tag{B25}
$$
with our conventions.

Our spin wave functions are the obvious ones:

$$
\chi_0 = 
\left[
\begin{array}{c}
1 \\
2
\end{array}
\right]^{1/2}
(\uparrow\downarrow - \downarrow\uparrow),
\tag{B26}
$$
and
$$
\chi_{11} = \uparrow\uparrow , \tag{B27}
$$
$$
\chi_{10} = 
\left[
\begin{array}{c}
1 \\
2
\end{array}
\right]^{1/2}
(\uparrow\downarrow+\downarrow\uparrow) , \tag{B28}
$$
$$
\chi_{1-1} = \downarrow\downarrow ,
\tag{B29}
$$
and our color wave function,
$$
\phi_{\mathrm{color}}=
\left[\frac{1}{3}\right]^{1/2}
(\overline{R}R + \overline{B}B + \overline{Y}Y)
\tag{B30}
$$
is even more obvious.

Our conventions for spatial wave functions are more involved. We choose them so that in the harmonic limit they go over into the following harmonic-oscillator wave functions $\Psi_{nLM}$ (where $n$ is the number of radial nodes minus one):

$$
\Psi_{000}=+\frac{\beta^{3/2}}{\pi^{3/4}}e^{-\beta^{2}r^{2}/2}, \tag{B31}
$$
$$
\Psi_{011} = -\frac{\beta^{5/2}}{\pi^{3/4}} r_+ e^{-\beta^{2}r^{2}/2}, \tag{B32}
$$
$$
\Psi_{022} = +\frac{\beta^{7/2}}{\pi^{3/4}}\frac{1}{\sqrt{2}} r_+^2 e^{-\beta^{2}r^{2}/2} , \tag{B33}
$$
$$
\Psi_{033} = -\frac{\beta^{9/2}}{\pi^{3/4}}\frac{1}{\sqrt{6}} r_+^3 e^{-\beta^{2}r^{2}/2}, \tag{B34}
$$
etc., where $r_+ = x + iy$ and the alternating sign follows the sign of $Y_{LL}(\theta, \phi)$, and where the lower states of each $L$ multiplet follow from the Condon-Shortley convention. For radial excitations we follow the convention that as $r \rightarrow \infty$ the ratio of a radial excitation to its ground state should be positive. Thus

$$
\Psi_{100} = +\frac{\beta^{7/2}}{\pi^{3/4}}\left[\frac{2}{3}\right]^{1/2}(r^2 - \frac{3}{2}\beta^{-2})e^{-\beta^{2}r^{2}/2}, \tag{B35}
$$
$$
\Psi_{111} = -\frac{\beta^{9/2}}{\pi^{3/4}}\left[\frac{5}{2}\right]^{1/2}r_+(r^2 - \frac{5}{2}\beta^{-2})e^{-\beta^{2}r^{2}/2}. \tag{B36}
$$

We complete these conventions by noting that when constructing the states $|^{2S+1}L_{JM}\rangle$ we combine angular momenta in the $L\cdot S$ order with Wigner’s conventional Clebsch-Gordon coefficients. Thus,
$$
|A_1^{+}(+)\rangle = |I=I_3=1, ^{3}P_1\rangle
$$
$$
=\phi_{\mathrm{color}}(-|u\bar{d}\rangle) \left[\frac{1}{2}\right]^{1/2} (\Psi_{011}\chi_{10}-\Psi_{010}\chi_{11}). \tag{B37}
$$

***

APPENDIX C: CONVERSION FROM HELICITY TO PARTIAL-WAVE AMPLITUDES

In terms of the helicity amplitude
$$
H_m = A[M^{*}_{j^{*}m} \rightarrow M_{jm} + P(q\hat{z})] \tag{C1}
$$
defined by Eq. (19), the decay rate is (see below)
$$
\Gamma(M^{*} \rightarrow M + P) = \frac{1}{2j^{*}+1} \left[ \frac{q}{2\pi} \right] \sum_{m} | \tilde{H}_m |^2 , \tag{C2}
$$
where
$$
\tilde{H}_m = \frac{(2\pi)^{9/2}}{(2M^{*}2E)^{1/2}} H_m .
$$

It is convenient to define $h_m = \sqrt{2} H_m$ for $m > 0$ and $h_0=H_0$ so that we can write
$$
\sum_m |H_m|^2 = \sum_{m \geq 0} |h_m|^2 . \tag{C3}
$$

We then find that if $M$ is a spin-zero particle, the partial-wave amplitudes $A_L = h_0$ and that if $M$ is a vector particle we have the results of Table XI. For photon amplitudes we quote our results directly as helicity amplitudes so that these conversions are unnecessary.

In proceeding to the above formula for the rate, we have set a factor of $E/M^{*} = 1$, consistent with our nonrelativistic calculation. There is, however, some further motivation for doing this: in the mock-hadron method of Ref. 8 and Appendix D, this factor is automatically absent in the simple cases like magnetic-dipole decays which it can treat.

***

APPENDIX D: THE MOCK-MESON METHOD

The mock-meson method$^8$ assumes a correspondence between Lorentz-invariant amplitudes of real hadrons and those of free quarks. The prescription is to (1) express the physical matrix element $\mathcal{M}$ in terms of Lorentz covariants with scalar coefficients $A$, (2) with each hadron $H$ of mass $M_H$ associate a mock hadron $\bar{H}$ (consisting of free quarks with the wave function of the bound quarks in $H$)

## PDF Page 40

***

TABLE XI. Conversion from helicity to partial-wave amplitudes in $M^*_j \rightarrow V+P$.

| $j^*$ | Positive parity $M^*$                      | Negative parity $M^*$                        |
|-------|--------------------------------------------|----------------------------------------------|
| 0     |                                            | $A_P = -h_0$                                 |
| 1     | $A_S=\sqrt{2/3} h_1 + \sqrt{1/3} h_0$      | $A_P = -h_1$                                 |
|       | $A_D=\sqrt{1/3} h_1 - \sqrt{2/3} h_0$      |                                              |
| 2     | $A_D = - h_1$                              | $A_P = \sqrt{3/5} h_1 + \sqrt{2/5} h_0$      |
|       |                                            | $A_F = \sqrt{2/5} h_1 - \sqrt{3/5} h_0$      |
| 3     | $A_D = \sqrt{4/7} h_1 + \sqrt{3/7} h_0$    | $A_F = -h_1$                                 |
|       | $A_G = \sqrt{3/7} h_1 - \sqrt{4/7} h_0$    |                                              |
| 4     | $A_G = - h_1$                              | $A_F = \sqrt{5/9} h_1 + \sqrt{4/9} h_0$      |
|       |                                            | $A_H = \sqrt{4/9} h_1 - \sqrt{5/9} h_0$      |
| 5     | $A_G = \sqrt{6/11} h_1 + \sqrt{5/11} h_0$  | $A_H = - h_1$                                |
|       | $A_I = \sqrt{5/11} h_1 - \sqrt{6/11} h_0$  |                                              |

***

and a mock mass $\widetilde{M}_H$ equal to the mean total energy of the (free) quarks in $H$, (3) calculate $\widetilde{\mathcal{M}}$, the mock matrix element, in terms of free-quark amplitudes, and (4) if (as is the case in many simple circumstances) $\widetilde{\mathcal{M}}$ has the same form as $\mathcal{M}$ take $A = \widetilde{A}$. Such a procedure is obviously not completely satisfactory, but it at least amounts to a partial relativization of the quark model.

A simple example will perhaps help to clarify the prescription: consider $\omega \rightarrow \pi \gamma$. The relevant hadronic matrix element is

$$
\langle \pi(k') | j^{\mu}_{\text{em}}(0) | \omega(e, k) \rangle
$$
$$
= \frac{1}{(2\pi)^3} \mu_{\pi\omega} \epsilon^{\mu\nu\rho\gamma} e_\nu (k' - k)_\rho (k' + k)_\gamma . \tag{D1}
$$

We wish to know $\mu_{\pi\omega}$ so we calculate instead $\widetilde{\mu}_{\pi\omega}$ by taking the matrix element of $j^\mu$ in mock mesons. For example, as $k \rightarrow 0$ (we always must work with states nearly at rest for which $\widetilde{E} \simeq \widetilde{M}$)

$$
| \widetilde{\omega}(+, \mathbf{k}) \rangle
$$
$$
= (2 \widetilde{M}_\omega)^{1/2} \int \sum_i d^3p \, \phi_\omega(p)
$$
$$
\qquad \times a_i \left| q_i \left[ \frac{\mathbf{k}}{2} + \mathbf{p}, \uparrow \right] \bar{q}_i \left[ \frac{\mathbf{k}}{2} - \mathbf{p}, \uparrow \right] \right\rangle , \tag{D2}
$$

where $a_u = a_d = (1/2)^{1/2}$ are flavor factors, $\phi_\omega(p)$ is the normalized momentum-space wave function, and

***

$\widetilde{M}_\omega = 2 \int d^3p \, E \, | \phi_\omega(p) |^2 ;$

we then find

$$
\mu_{\pi\omega} = \widetilde{\mu}_{\pi\omega} = \frac{e}{2} \cos(\theta_V - \theta_{\text{ideal}})
\left(
\frac{2\widetilde{M}_\pi^{1/2} \widetilde{M}_\omega^{1/2}}
     {\widetilde{M}_\pi + \widetilde{M}_\omega}
\right)
$$

$$
\qquad\qquad  \times \int d^3p \, \phi^*_{\pi} \phi_{\omega} 
\left[
\frac{m + 2E}{3E^2}
\right]
\tag{D3}
$$

which of course reduces to the usual nonrelativistic result in the limit that $\langle p^2 \rangle \rightarrow 0$.

We do not take the precise forms of the relativistic modifications to amplitudes, e.g., the factor $(m + 2E / 3E^2)$ in (D3) too seriously, but use it as a guide: our prescription is to modify the leading behavior of all such matrix elements by inserting factors of $(m_i m_j / E_i E_j)^f/2$ with $f$ a fitted parameter.

For convenience we now list various other definitions and results we use in the text. First the definitions: for $P \rightarrow l \nu$ and related decays we use

$$
\langle 0 | A^\mu_1(0) | P(k) \rangle = \frac{1}{(2\pi)^{3/2}} i f_P M_P k^\mu \tag{D4}
$$

with all matrix elements defined in terms of the appropriate axial-vector current with unit strength; for $V \rightarrow l^+ l^-$ and related decays we take

$$
\langle 0 | j^\mu_{\text{em}}(0) | V(e, k) \rangle = -\frac{1}{(2\pi)^{3/2}} e f_V M_V^2 e^{\mu} \tag{D5}
$$

while for $\tau \rightarrow A_1 \nu_\tau$ the analogous definition

$$
\langle 0 | A^\mu_{1+i2}(0) | A_1(e, k) \rangle = -\frac{1}{(2\pi)^{3/2}} f_{A_1} M_{A_1}^2 e^{\mu} \tag{D6}
$$

is used. In terms of these couplings it follows that

$$
\Gamma(P \rightarrow l \nu) = \frac{G^2 f_P^2 m_l^2}{8 M_P \pi} (M_P^2 - M_l^2)^2 , \tag{D7}
$$

$$
\Gamma(V \rightarrow l^+ l^-) = \frac{4\pi}{3} \alpha^2 M_V f_V^2 , \tag{D8}
$$

$$
\Gamma(\tau \rightarrow A_1 \nu_\tau) = \frac{G^2 f_{A_1}^2 m_\tau^3 M_{A_1}^2}{16\pi}
\left[
1 - \frac{M_{A_1}^2}{m_\tau^2}
\right]
\left[
1 + \frac{2 M_{A_1}^2}{m_\tau^2}
\right] , \tag{D9}
$$

$$
\Gamma(V \rightarrow P \gamma) = \frac{4}{3} \alpha \left( \frac{\mu_{PV}}{e} \right)^2 \omega_\gamma^3 , \tag{D10}
$$

and

$$
\Gamma(P \rightarrow V \gamma) = 4 \alpha \left( \frac{\mu_{PV}}{e} \right)^2 \omega_\gamma^3 . \tag{D11}
$$

## PDF Page 41

***

1. T. Appelquist and H. D. Politzer, Phys. Rev. Lett. **34**, 43 (1975); Phys. Rev. D. **12**, 1404 (1975); E. Eichten _et al._, Phys. Rev. Lett. **34**, 369 (1975); T. Appelquist _et al._, _ibid._ 34, 365 (1975).

2. E. Eichten _et al._, Phys. Rev. D **17**, 3090 (1978); **21**, 203 (1980); W. Celmaster, H. Georgi, and M. Machacek, _ibid._ **17**, 879 (1978); **17**, 886 (1978); H. Krasemann and S. Ono, Nucl. Phys. **B154**, 283 (1979); W. Buchmüller and S.-H. H. Tye, Phys. Rev. D **24**, 132 (1981); T. Appelquist, R. M. Barnett, and K. D. Lane, Ann. Rev. Nucl. Sci. **28**, 387 (1978); M. Krammer and H. Krasemann, Acta Phys. Austriaca, Suppl. XXI, 259 (1979); E. Eichten in _Experimental Meson Spectroscopy—1980_, proceedings of the Sixth International Conference, Brookhaven, edited by S. U. Chung and S. J. Lindenbaum (AIP, New York, 1981), p. 387; R. D. Viollier and J. Rafelski, Helv. Phys. Acta **53**, 352 (1981).

3. A. De Rújula, H. Georgi, and S. L. Glashow, Phys. Rev. D **12**, 147 (1975) was a seminal paper for the application of potential-model ideas to light quarks. For reviews of this general area, see N. Isgur, in _The New Aspects of Subnuclear Physics_, proceedings of the XVI International School of Subnuclear Physics, Erice, 1978, edited by A. Zichichi (Plenum, New York, 1980), p. 107; O. W. Greenberg, Ann. Rev. of Nucl. Part. Phys. **28**, 327 (1978); A. J. G. Hey, in _Proceedings of the European Physical Society International Conference on High Energy Physics, Geneva, 1979_, edited by A. Zichichi (CERN, Geneva, 1980); J. Rosner, in _Techniques and Concepts of High Energy Physics_, proceedings of the NATO Advanced Study Institute, St. Croix, 1980, edited by T. Ferbel (Plenum, New York, 1981); and R. H. Dalitz, Nucl. Phys. **A353**, 251 (1981).

4. While our model is unique, it shares many features with the (vast) literature on this subject. It is especially close in spirit to the excellent work of D. P. Stanley and D. Robson, Phys. Rev. D **21**, 3180 (1980), although quite different in detail, and also to the recent work of J. Carlson, J. Kogut, and V. R. Pandaripande, _ibid._ **27**, 233 (1983); **28**, 2807 (1983). See also H. J. Schnitzer, Phys. Lett. **65B**, 239 (1976); **69B**, 477 (1977); Phys. Rev. D **18**, 3482 (1978); B. R. Martin and L. J. Reinders, Nucl. Phys. **B143**, 309 (1978); Phys. Lett. **78B**, 144 (1978); A. B. Henriques, B. Kellet and R. G. Moorhouse, Phys. Lett. **64B**, 85 (1976); R. Barbieri _et al._, Nucl. Phys. **B105**, 125 (1976); E. Eichten _et al._, Phys. Rev. D **17**, 3090 (1978); L. J. Reinders, Z. Phys. C **4**, 95 (1980); I. Cohen and H. J. Lipkin, Nucl. Phys. **B112**, 213 (1976); J. Arafune, M. Fukugita, and Y. Oyanagi, Phys. Rev. D **16**, 772 (1977); T. Barnes, Z. Phys. C **11**, 135 (1981); N. Barik and S. N. Jena, Phys. Rev. D **21**, 2647 (1980); **22**, 1704 (1980); **24**, 680 (1981); Phys. Lett. **97B**, 261 (1980); **97B**, 265 (1980); S. Mahmood and M. M. Illyas, Can. J. Phys. **59**, 387 (1981); A. Bradley, J. Phys. G **4**, 1517 (1978); Z. Phys. C **5**, 239 (1980); D. Beavis _et al._, Phys. Rev. D **20**, 743 (1979); R. Levine and Y. Tomozawa, _ibid._ **19**, 1572 (1979); **21**, 840 (1980); W. Celmaster and F. S. Henyey, _ibid._ **17**, 3268 (1978); **18**, 1688 (1978); M. Kauragi _et al._, Phys. Lett. **97B**, 143 (1980); R. K. Bhaduri, L. E. Cohler, and Y. Nogami, Phys. Rev. Lett. **44**, 1369 (1980); D. B. Lichtenberg and J. G. Wills, _ibid._ **35**, 1055 (1975); J. F. Gunion and R. S. Wiley, Phys. Rev. D **12**, 174 (1975); A. B. Henriques, Z. Phys. C **11**, 31 (1981); A. Martin and J. M. Richard, Phys. Lett. **115B**, 323 (1982); R. W. Childers, _ibid. 126B_, 485 (1983); S. Ono and F. Schöberl, _ibid. 118B_, 419 (1982); A. T. Aerts and L. Heller, Phys. Rev. D **29**, 513 (1984); A. B. Henriques, Z. Phys. C **18**, 213 (1983); J. Morishitz _et al._, _ibid._ **19**, 167 (1983); D. B. Lichtenberg, W. Namgung, and J. G. Wills, _ibid._ **19**, 19 (1983); K. Sink, J. Plum. Phys. **2**, 212 (1984); see also Ref. 10.

5. N. Isgur, Phys. Rev. D **12**, 3770 (1975); **13**, 122 (1976), and in Ref. 3; A. De Rújula, H. Georgi, and S. L. Glashow, in Ref. 3; H. Fritzsch and P. Minkowski, Nuovo Cimento **30A**, 393 (1975); H. Fritzsch and J. D. Jackson, Phys. Lett. **66B**, 365 (1977); C. De La Vaissiere, Nuovo Cimento **41A**, 419 (1977); R. H. Capps, Phys. Rev. D **18**, 848 (1978); H. F. Jones and M. D. Scadron, Nucl. Phys. **B155**, 409 (1979); I. Cohen and H. J. Lipkin, _ibid._ **B151**, 16 (1979); P. J. O'Donnell and R. H. Graham, Phys. Rev. D **19**, 284 (1979); A. T. Filippov, Yad. Fiz. **29**, 1035 (1979) [Sov. J. Nucl. Phys. **29**, 534 (1979)]; K. Hirata, T. Kobayashi, and Y. Takaiwa, Phys. Rev. D **18**, 236 (1978); F. D. Gault and A. B. Rimmer, Z. Phys. C **8**, 353 (1981); A. B. Rimmer, Phys. Lett. **103B**, 459 (1981).

6. E. Witten, Nucl. Phys. **B156**, 269 (1979).

7. S. Godfrey and N. Isgur (work in progress).

8. C. Hayne and N. Isgur, Phys. Rev. D **25**, 1944 (1982).

9. The status of the $E$ is in a state of flux as the old $E$ was partly the newly discovered pseudoscalar state at 1.44 GeV (see Sec. V A). Ph. Gavillet _et al_. [Z. Phys. C **16**, 119 (1982)] have recently argued that the residual $1^{++}$ signal in the $1^{+} S(K^* \overline{K}^*_{cc})$ wave is not resonant and that the true $E$ is an object they call $D'$ at $1.53 \pm 0.01$ GeV with a width of $0.11 \pm 0.02$ GeV. In view of these developments, we do not show an experimental value for the $E$ mass in Fig. 5, nor do we take too seriously the experimental $E$ widths of Table V.

10. Our potential turns out to be quite similar to others which have appeared previously, mainly in the context of heavy-quark spectroscopy. See, e.g., J. L. Richardson, Phys. Lett. **82B**, 272 (1979); G. Bhanot and S. Rudaz, _ibid. 78B_, 119 (1978); A. Martin, _ibid. 90B_, 338 (1980); _ibid. 100B_, 511 (1981); _High Energy Physics—1980_, proceedings of the XX International Conference, Madison, Wisconsin, edited by Loy-al Durand and Lee G. Pondrom (AIP, New York, 1981), p. 715; D. B. Lichtenberg and J. G. Wills, Nuovo Cimento **47A**, 483 (1978); A. Billiore and A. Mord, Nucl. Phys. **B135**, 131 (1978); C. Quigg and J. L. Rosner, Phys. Lett. **71B**, 153 (1977); AIP Conf. Proc. **68**, 719 (1980); Phys. Rev. D **23**, 2625 (1981); D. Pignon and C. A. Piketty, Phys. Lett. **74B**, 108 (1978); see also Ref. 4.

11. R. Koniuk and N. Isgur, Phys. Rev. Lett. **44**, 845 (1980); Phys. Rev. D **21**, 1868 (1980); **23**, 818 (E) (1981); R. Koniuk in _Baryon 1980_, proceedings of the IVth International Conference on Baryon Resonances, Toronto, edited by N. Isgur (University of Toronto, Toronto, 1981), p. 217.

12. C. Becchi and G. Morpurgo, Phys. Rev. **149**, 1284 (1966); **140B**, 687 (1965); Phys. Lett. **17**, 352 (1965); A. N. Mitra and M. Ross, Phys. Rev. **158**, 1630 (1967); D. Faiman and A. W. Hendry, _ibid._ **173**, 1720 (1968); H. J. Lipkin, Phys. Rep. **8C**, 173 (1973); J. L. Rosner, Phys. Rep. **11C**, 189 (1974); R. Horgan, in _Proceedings of the Topical Conference on Baryon Resonances, Oxford, 1976_, edited by R. T. Ross and D. H. Saxon (Rutherford Laboratory, Chilton, Didcot, England, 1977), p. 435; A. Le Yaouanc _et al._, Phys. Rev. D **11**, 1272 (1975); L. A. Copley, G. Karl, and E. Obryk, Phys. Lett. **29B**, 177 (1969); Nucl. Phys. **B13**, 303 (1969); D. Faiman, A. W. Hendry, Phys. Rev. **180**, 1572 (1969); Hoichi Ohta, Phys. Rev. Lett. **43**, 1201 (1979); R. G. Moorhouse, _ibid. 16_, 771 (1966); R. P. Feynman, M. Kislinger, and F. Ravndal, Phys. Rev. D **3**, 2706 (1971); R. G. Moorhouse and N. H. Parsons, Nucl. Phys. B**62**, 109 (1973).

## PDF Page 42

***

1.  H. J. Lipkin and S. Meshkov, Phys. Rev. Lett. 14, 670 (1965); D. Faiman and A. W. Hendry, Phys. Rev. **173**, 1720 (1968); **180**, 1609 (1969); E. W. Colglazier and J. L. Rosner, Nucl. Phys. **B27**, 349 (1971); W. Petersen and J. Rosner, Phys. Rev. D **6**, 820 (1972); A. J. G. Hey, P. J. Litchfield, and R. J. Cashmore, Nucl. Phys. **B95**, 516 (1975); F. Gilman and I. Karliner, Phys. Rev. D **10**, 2194 (1974); J. Babcock and J. Rosner, Ann. Phys. (N.Y.) **96**, 191 (1976); J. Babcock _et al._, Nucl. Phys. **B126**, 87 (1977); D. Faiman and D. E. Plane, _ibid._ **B50**, 379 (1972).

2.  F. Binon _et al._, Nuovo Cimento **78A**, 313 (1983).

3.  W. Hoogland, in _New Flavors and Hadron Spectroscopy_, proceedings of the XVI Rencontre de Moriond, Les Arcs, France, 1981, edited by J. Trân Thanh Vân (Editions Frontières, Dreux, France, 1981), p. 209; C. Daum _et al._, Phys. Lett. **99B**, 270 (1980); **99B**, 281 (1980); **99B**, 285 (1980); Nucl. Phys. **B182**, 269 (1981); **B187**, 1 (1981).

4.  R. Kokoski, Ph.D. thesis, University of Toronto, 1983; R. Kokoski and N. Isgur, University of Toronto Report No. UTPT-85-05, 1984 (unpublished).

5.  Relativistic smearing effects have often been considered: See J. F. Gunion and L. F. Li, Phys. Rev. D **13**, 82 (1976); A. Bradley, Phys. Lett. **77B**, 422 (1978); J. S. Kang and J. Sucher, Phys. Rev. D **18**, 2698 (1978); E. C. Poggio and H. J. Schnitzer, _ibid._ **20**, 1175 (1979); P. Ditsas, N. A. McDougall, and R. G. Moorhouse, Nucl. Phys. **B146**, 191 (1978); J. E. Paschalis and G. J. Gounaris, _ibid._ **B222**, 473 (1983).

6.  For other work on the magnetic-dipole decays of mesons, see G. Feinberg and J. Sucher, Phys. Rev. Lett. **35**, 1740 (1975); N. Isgur, _ibid._ **36**, 1262 (1976); A. B. Govorkov and S. B. Drenszky, Yad. Fiz. **26**, 851 (1977) [Sov. J. Nucl. Phys. **26**, 446 (1977)]; L. Maharana and S. P. Misra, Phys. Rev. D **18**, 2530 (1978); L. P. Singh, Phys. Rev. D **19**, 2812 (1979); J. L. Rosner, in _High Energy Physics––1980_ (Ref. 10), p. 540; H. Grötch and K. J. Sebastian, Phys. Rev. D **25**, 2944 (1982); G. Cocho, M. Fortes, and H. Vucetich, _ibid._ **16**, 3339 (1977); and for a recent review see P. J. O’Donnell, Rev. Mod. Phys. **53**, 673 (1981).

7.  For E1 decays we use a hybrid form of Eq. (22) after application of Seigert’s theorem [see R. G. Sachs and N. Austern, Phys. Rev. **81**, 705 (1951); **81**, 710 (1951)] to transform Eq. (22) into ordinary dipole form. We owe our appreciation of the importance of this theorem to R. McClary and N. Byers, Phys. Rev. D **28**, 1692 (1983), and private communication. Our results are quite similar to theirs. See also H. Krasemann, Phys. Lett. **101B**, 259 (1981); P. Moxhay and J. L. Rosner, Phys. Rev. D **28**, 1132 (1983); and also the more recent work by R. McClary and N. Byers, UCLA Report No. UCLA/TEP/83/15 (unpublished).

8.  D. Berg _et al._, in _Experimental Meson Spectroscopy—1983_, proceedings of the Seventh International Conference, Brookhaven, edited by S. J. Lindenbaum (AIP, New York, 1984), p. 157.

9.  In addition to Ref. 8, see R. Barbieri, R. Gatto, R. Kögerler, and Z. Kunszt, Phys. Lett. **57B**, 455 (1975); Y. Abe _et al._, Prog. Theor. Phys. **60**, 639 (1978); **61**, 1566 (1979); **63**, 1078 (1980); L. Bergström, H. Snellman, and G. Tengstrand, Phys. Lett. **80B**, 242 (1979); P. M. Fishbane, D. Horn, and S. Meshkov, Phys. Rev. D **19**, 288 (1979); W. Kummer, Nucl. Phys. **B179**, 365 (1981); **B185**, 41 (1981); H. Kraseman, Phys. Lett. **96B**, 397 (1980); G. R. Goldstein and J. Maharana, Z. Phys. C **12**, 23 (1982).

10. R. Van Royen and V. F. Weisskopf, Nuovo Cimento **50**, 617 (1967); **51**, 583 (1967).

***

23. E. D. Bloom and C. W. Peck, Ann. Rev. Nucl. Part. Phys. **33**, 143 (1983).

24. M. Althoff _et al._, Phys. Lett. **121B**, 216 (1983). Note that the results of this experiment are misquoted by the Particle Data Group, Rev. Mod. Phys. **56**, S1(1984), we presume because a handwritten European 11 was transcribed into 77.

25. E. B. Dally _et al._, Phys. Rev. Lett. **48**, 375 (1982).

26. E. B. Dally _et al._, Phys. Rev. Lett. **45**, 232 (1980).

27. W. R. Molzon _et al._, Phys. Rev. Lett. **41**, 1213 (1978).

28. A. Cordier _et al._, Phys. Lett. **106B**, 155 (1981); M. Atkinson _et al._, _ibid._ **127B**, 132 (1983).

29. The possibility that ι(1440) is a glueball has been widely discussed. H. J. Lipkin in Phys. Lett. **106B**, 114 (1981); I. Cohen, N. Isgur, and H. J. Lipkin in Phys. Rev. Lett. **48**, 1074 (1982); and H. J. Lipkin and I. Cohen, Phys. Lett. **135B**, 215 (1984) argue for the interpretation of ι(1440) as a $q \bar q$ radial excitation, although not exactly as outlined in the text. See Gault and Rimmer and also Rimmer in Ref. 5 for a similar discussion. The glueball interpretation has been strongly advocated by M. Chanowitz, Phys. Rev. Lett. **46**, 981 (1981); in _Strong Interactions_, proceedings of the 9th SLAC Summer Institute on Particle Physics, 1981, edited by Anne Mosher (SLAC Report No. 245, 1982); and in _Particles and Fields—1981: Testing the Standard Model_, proceedings of the Annual Meeting of the Division of Particles and Fields of the APS, Santa Cruz, California, edited by C. A. Heusch and W. T. Kirk (AIP, New York, 1982), p. 85. See also J. Donoghue, K. Johnson, and B. A. Li, Phys. Lett. **99B**, 416 (1981); K. Ishikawa, Phys. Rev. Lett. **46**, 978 (1981); J. F. Donoghue and H. Gomm, Phys. Lett. **112B**, 409 (1982); T. Barnes in _New Flavors and Hadron Spectroscopy_ (Ref. 15), p. 175; M. Frank and P. J. O’Donnell, Phys. Lett. **133B**, 253 (1983); Phys. Rev. D **29**, 921 (1984).

30. N. Stanton _et al._, Phys. Rev. Lett. **42**, 346 (1979).

31. We are grateful to Gabriel Karl for pointing this out to us.

32. See H. J. Schnitzer, Phys. Lett. **76B**, 461 (1978); Nucl. Phys. **B207**, 131 (1982).

33. For another possibility (which could be in effect in parallel to ours) see H. J. Lipkin, Phys. Lett. **72B**, 249 (1977).

34. J. Weinstein and N. Isgur, Phys. Rev. Lett. **48**, 659 (1982); Phys. Rev. D **27**, 588 (1983).

35. R. L. Jaffe, Phys. Rev. D **15**, 267 (1977); **15**, 281 (1977); R. L. Jaffe and K. Johnson, Phys. Lett..**60B**, 201 (1976); R. L. Jaffe, Phys. Rev. D **17**, 1444 (1978).

36. For analyses of experiments on the $0^+$ channels see P. Estabrooks, Phys. Rev. D **19**, 2678 (1979); A. C. Irving, A. D. Martin, and P. J. Done, Z. Phys., C **10**, 45 (1981); N. N. Achasov, S. A. Devyanin, G. N. Shestakov, Phys. Lett. **96B**, 168 (1980); Yad. Fiz. 32, 1098 (1980) [Sov. J. Nucl. Phys. **32**, 566 (1980)]; see also Refs. 37 and 38.

37. Our suggestion here has much in common with the analysis of A. B. Wicklund in _New Flavors and Hadron Spectroscopy_ (Ref. 15), p. 339; and A. B. Wicklund _et al._, Phys. Rev. Lett. **45**, 1469 (1980). An interpretation quite different from ours has been given in N. A. Törnquist, Phys. Rev. Lett. **49**, 624 (1982); see also, however, N. N. Achasov, S. A. Devyanin, and G. N. Shestakov, Z. Phys. C **22**, 53 (1984).

38. D. Cohen _et al._, Phys. Rev. D **22**, 2595 (1980).

39. Preliminary results from this work have been quoted by one of us (N.I.) in _New Flavors and Hadron Spectroscopy_ (Ref. 15), p. 247, and in _Particles and Fields—1981: Testing the Standard Model_ (Ref. 29), p. 1. A more complete account of this work will appear shortly in the University of Toronto Ph.D. thesis of S.G.

## PDF Page 43

32

MESONS IN A RELATIVISTIC QUARK MODEL…                                                                 231

***

40 Variants of this approach (e.g., use of the Bethe-Salpeter equation) have been used by many authors to derive relativistic corrections: D. Gromes, Nucl. Phys. B131, 80 (1977); W. Celmaster and F. Henyey in Ref. 4; H. J. Schnitzer, Phys. Rev. D 13, 74 (1976); Phys. Rev. Lett. 35, 1540 (1975); and in Refs. 4 and 32; J. Pumplin, W. Repko, and A. Sato, Phys. Rev. Lett. 35, 1538 (1975); E. Eichten and F. Feinberg, ibid. 43, 1205 (1979); Phys. Rev. 23, 2724 (1981); A. B. Henriques, B. H. Kellet, and R. G. Moorhouse in Ref. 4; Lai-Him Chan, Phys. Lett. 71B, 422 (1977); L. J. Reinders, J. Phys. G 4, 1241 (1978); and in Baryon 1980 (Ref. 11), p. 203; B. Durand and L. Durand, Phys. Rev. D 25, 2312 (1982); M. G. Olsson and K. J. Miller, Phys. Rev. D 28, 674 (1983).

41 Our picture of this interaction is that of the flux-tube model where the long-range force arises continuously out of the color-electric potential but is, unlike the short-range force, spin independent. Although such a potential is simple nonrelativistically, it cannot easily be fully represented relativistically since we are asking that it have the $J^{PC}$ properties of a $\gamma^0 \gamma^0$ interaction, namely, $0^{++}$ overall but $0^{+-}$ (an exotic combination) at each vertex. We can therefore use the interaction to derive the Thomas-precession interaction but we must impose the vertex properties under $q \leftrightarrow \bar{q}$ (relevant for meson-baryon connections) by hand.

