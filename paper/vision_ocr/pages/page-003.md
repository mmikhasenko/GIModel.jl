<!-- PDF page 3; model gpt-4.1 -->

and $\alpha_s(r)$ is the running coupling constant of QCD which we will discuss below in more detail. Since

$$
\langle \mathbf{F}_i \cdot \mathbf{F}_j \rangle =
\begin{cases}
- \frac{4}{3} & \text{in a meson} \,, \\
- \frac{2}{3} & \text{in a baryon} \,,
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

---

**FIG. 1.** The origin of the annihilation term $H_A$: a typical graph.
