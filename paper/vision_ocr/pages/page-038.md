<!-- PDF page 38; model gpt-4.1 -->

226                                           STEPHEN GODFREY AND NATHAN ISGUR            32

and
\[
S_{\rm eff}(r) = \tilde{S}(r) - \frac{\mathbf{S}_1 \cdot \mathbf{L}}{2 m_1^2}\frac{1}{r} \frac{\partial \tilde{S}_{11}^{so(s)}}{\partial r} - \frac{\mathbf{S}_2 \cdot \mathbf{L}}{2 m_2^2}\frac{1}{r} \frac{\partial S_{22}^{so(s)}}{\partial r} .
\tag{A16}
\]

To actually perform calculations with these potentials we diagonalize the Hamiltonian matrix obtained from (1) in a (large) harmonic-oscillator basis. The harmonic basis is particularly useful in this context. For example, to take a matrix element of an operator of the form $f(p)g(r)$ we use

\[
\langle i | f(p)g(r) | j \rangle = \sum_n \langle i | f(p) | n \rangle \langle n | g(r) | j \rangle .
\tag{A17}
\]

The matrix elements $\langle i | f(p) | n \rangle$ are then calculated in momentum space, while the matrix elements $\langle n | g(r) | j \rangle$ are calculated in configuration space. Since for a harmonic-oscillator basis these two sets of wave functions are simple polynomials of the same form, the calculations are especially simple. Once the elements of the Hamiltonian matrix are calculated this way, the choice of the Gaussian parameter $\beta$ which characterizes the harmonic-oscillator basis is optimized in accordance with the variational principle. (Note that the best energy for a given state can always be obtained by minimizing it with respect to $\beta$. Of course to generate an orthogonal set of wave functions in a given sector, a single value of $\beta$ must be used, and in practice we take the value that minimizes the energy of the last state of the set for this purpose. Our basis is so large that very little error is introduced in this approximation.)

**APPENDIX B: WAVE-FUNCTION CONVENTIONS**

Here we make all of our wave-function conventions explicit so that our results may be more readily used. First, since we use the “natural” SU(3) conventions for quark and antiquark transfer operators ($q_i \rightarrow q_j$ is always $+1$ and $\bar{q}_i \rightarrow \bar{q}_j$ is always $-1$) our flavor wave functions, referred to the names of the pseudoscalar octet but generally applicable, are

\[
\pi^+ = -u \bar{d}, \tag{B1}
\]

\[
\pi^0 = \frac{1}{\sqrt{2}} (u\bar{u} - d\bar{d}), \tag{B2}
\]

\[
\pi^- = d\bar{u}, \tag{B3}
\]

\[
K^+ = -u\bar{s}, \tag{B4}
\]

\[
K^0 = -d\bar{s}, \tag{B5}
\]

\[
\bar{K}^0 = -s\bar{d}, \tag{B6}
\]

\[
K^- = s\bar{u}, \tag{B7}
\]

\[
\eta_8 = \frac{1}{\sqrt{6}} (u\bar{u} + d\bar{d} - 2s\bar{s}), \tag{B8}
\]

\[
\eta_1 = \frac{1}{\sqrt{3}} (u\bar{u} + d\bar{d} + s\bar{s}) \tag{B9}
\]

which satisfy the de Swart conventions on phases. For ideally mixed isoscalar mesons we take

\[
M_{ns} = \frac{1}{\sqrt{2}} (u\bar{u} + d\bar{d}), \tag{B10}
\]

\[
M_s = s\bar{s}, \tag{B11}
\]

and normally define mixing angles $\phi$ relative to this basis via

\[
M = M_{ns} \cos\phi - M_s \sin\phi, \tag{B12}
\]

\[
M' = M_s \cos\phi + M_{ns} \sin\phi, \tag{B13}
\]

so that $\phi = \theta_{{SU(3)}} - \theta_{{\rm ideal}}$, where $\theta_{\rm ideal} \simeq 35.3^\circ$. (Note that for $\theta_{SU(3)} \rightarrow 0$ we get $M \rightarrow + “\eta_1”, M' \rightarrow - “\eta_8”$ as a consequence of our conventions.) In the special case of the pseudoscalar mesons, we will often use the “perfect-mixing” states (see the first of Refs. 5),

\[
\eta = \frac{1}{\sqrt{2}} (M_{ns} - M_s), \tag{B14}
\]

\[
\eta' = \frac{1}{\sqrt{2}} (M_{ns} + M_s), \tag{B15}
\]

which correspond to an SU(3) mixing angle of $\theta_{\rm ideal} - 45^\circ \simeq -10^\circ$. (Note that these states follow from (B12) and (B13) by taking $\phi = -45^\circ$ so that $M \rightarrow \eta' \simeq \eta_1$, $M' \rightarrow -\eta \simeq -\eta_8$.) For heavy-quark mesons we ignore all symmetries except isospin and simply take our state vectors to be $+|Q\bar{q}\rangle$, where $Q$ is the heavy quark and $q$ any other quark except the $d$, including $Q$ itself. We use $-|Q\bar{d}\rangle$ so that $( -|Q\bar{d}\rangle, |Q\bar{u}\rangle )$ form an isospin multiplet analogous to $(\bar{K}^0, K^-)$.

The SU(3) flavor wave functions lead to coupling operators $X_q^i$ in (19) given by

\[
X_q^{\pi^+} = -\sqrt{2}\left[ \frac{\lambda_1 - i \lambda_2}{2} \right] (u \rightarrow -\sqrt{2}d ), \tag{B16}
\]

\[
X_q^{\pi^0} = +\lambda_3~(u \rightarrow u, d \rightarrow -d ), \tag{B17}
\]

\[
X_q^{\pi^-} = +\sqrt{2} \left[ \frac{\lambda_1 + i \lambda_2}{2} \right] (d \rightarrow \sqrt{2}u ), \tag{B18}
\]

\[
X_q^{K^+} = -\sqrt{2} \left[ \frac{\lambda_4 - i \lambda_5}{2} \right] (u \rightarrow -\sqrt{2}s ), \tag{B19}
\]

\[
X_q^{K^0} = -\sqrt{2} \left[ \frac{\lambda_6 - i \lambda_7}{2} \right] (d \rightarrow -\sqrt{2}s ), \tag{B20}
\]

\[
X_q^{\bar{K}^0} = -\sqrt{2} \left[ \frac{\lambda_6 + i \lambda_7}{2} \right] (s \rightarrow -\sqrt{2}d ), \tag{B21}
\]
