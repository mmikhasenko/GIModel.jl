<!-- PDF page 40; model gpt-4.1 -->

228

STEPHEN GODFREY AND NATHAN ISGUR

---

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

---

and a mock mass $\widetilde{M}_H$ equal to the mean total energy of the (free) quarks in $H$, (3) calculate $\widetilde{\mathcal{M}}$, the mock matrix element, in terms of free-quark amplitudes, and (4) if (as is the case in many simple circumstances) $\widetilde{\mathcal{M}}$ has the same form as $\mathcal{M}$ take $A = \widetilde{A}$. Such a procedure is obviously not completely satisfactory, but it at least amounts to a partial relativization of the quark model.

A simple example will perhaps help to clarify the prescription: consider $\omega \rightarrow \pi \gamma$. The relevant hadronic matrix element is

\[
\langle \pi(k') | j^{\mu}_{\text{em}}(0) | \omega(e, k) \rangle
\]
\[
= \frac{1}{(2\pi)^3} \mu_{\pi\omega} \epsilon^{\mu\nu\rho\gamma} e_\nu (k' - k)_\rho (k' + k)_\gamma . \tag{D1}
\]

We wish to know $\mu_{\pi\omega}$ so we calculate instead $\widetilde{\mu}_{\pi\omega}$ by taking the matrix element of $j^\mu$ in mock mesons. For example, as $k \rightarrow 0$ (we always must work with states nearly at rest for which $\widetilde{E} \simeq \widetilde{M}$)

\[
| \widetilde{\omega}(+, \mathbf{k}) \rangle
\]
\[
= (2 \widetilde{M}_\omega)^{1/2} \int \sum_i d^3p \, \phi_\omega(p)
\]
\[
\qquad \times a_i \left| q_i \left[ \frac{\mathbf{k}}{2} + \mathbf{p}, \uparrow \right] \bar{q}_i \left[ \frac{\mathbf{k}}{2} - \mathbf{p}, \uparrow \right] \right\rangle , \tag{D2}
\]

where $a_u = a_d = (1/2)^{1/2}$ are flavor factors, $\phi_\omega(p)$ is the normalized momentum-space wave function, and

---

$\widetilde{M}_\omega = 2 \int d^3p \, E \, | \phi_\omega(p) |^2 ;$

we then find

\[
\mu_{\pi\omega} = \widetilde{\mu}_{\pi\omega} = \frac{e}{2} \cos(\theta_V - \theta_{\text{ideal}})
\left(
\frac{2\widetilde{M}_\pi^{1/2} \widetilde{M}_\omega^{1/2}}
     {\widetilde{M}_\pi + \widetilde{M}_\omega}
\right)
\]

\[
\qquad\qquad  \times \int d^3p \, \phi^*_{\pi} \phi_{\omega} 
\left[
\frac{m + 2E}{3E^2}
\right]
\tag{D3}
\]

which of course reduces to the usual nonrelativistic result in the limit that $\langle p^2 \rangle \rightarrow 0$.

We do not take the precise forms of the relativistic modifications to amplitudes, e.g., the factor $(m + 2E / 3E^2)$ in (D3) too seriously, but use it as a guide: our prescription is to modify the leading behavior of all such matrix elements by inserting factors of $(m_i m_j / E_i E_j)^f/2$ with $f$ a fitted parameter.

For convenience we now list various other definitions and results we use in the text. First the definitions: for $P \rightarrow l \nu$ and related decays we use

\[
\langle 0 | A^\mu_1(0) | P(k) \rangle = \frac{1}{(2\pi)^{3/2}} i f_P M_P k^\mu \tag{D4}
\]

with all matrix elements defined in terms of the appropriate axial-vector current with unit strength; for $V \rightarrow l^+ l^-$ and related decays we take

\[
\langle 0 | j^\mu_{\text{em}}(0) | V(e, k) \rangle = -\frac{1}{(2\pi)^{3/2}} e f_V M_V^2 e^{\mu} \tag{D5}
\]

while for $\tau \rightarrow A_1 \nu_\tau$ the analogous definition

\[
\langle 0 | A^\mu_{1+i2}(0) | A_1(e, k) \rangle = -\frac{1}{(2\pi)^{3/2}} f_{A_1} M_{A_1}^2 e^{\mu} \tag{D6}
\]

is used. In terms of these couplings it follows that

\[
\Gamma(P \rightarrow l \nu) = \frac{G^2 f_P^2 m_l^2}{8 M_P \pi} (M_P^2 - M_l^2)^2 , \tag{D7}
\]

\[
\Gamma(V \rightarrow l^+ l^-) = \frac{4\pi}{3} \alpha^2 M_V f_V^2 , \tag{D8}
\]

\[
\Gamma(\tau \rightarrow A_1 \nu_\tau) = \frac{G^2 f_{A_1}^2 m_\tau^3 M_{A_1}^2}{16\pi}
\left[
1 - \frac{M_{A_1}^2}{m_\tau^2}
\right]
\left[
1 + \frac{2 M_{A_1}^2}{m_\tau^2}
\right] , \tag{D9}
\]

\[
\Gamma(V \rightarrow P \gamma) = \frac{4}{3} \alpha \left( \frac{\mu_{PV}}{e} \right)^2 \omega_\gamma^3 , \tag{D10}
\]

and

\[
\Gamma(P \rightarrow V \gamma) = 4 \alpha \left( \frac{\mu_{PV}}{e} \right)^2 \omega_\gamma^3 . \tag{D11}
\]
