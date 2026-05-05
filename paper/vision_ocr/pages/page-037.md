<!-- PDF page 37; model gpt-4.1 -->

32 MESONS IN A RELATIVISTIC QUARK MODEL...

where $\sigma_0$ and $s$ are the universal parameters given in Table II. The parameter $\sigma_0$ reflects the fact that in a confined system the smearing must be limited, while $s$ is the coefficient of the expected linear relation $\sigma_{Q\bar{Q}} = sm_Q$ for a heavy $Q\bar{Q}$ system. The complicated $m_i, m_j$ dependence of the $\sigma_0$ term in (A9), designed to reflect the fact that in a $Q\bar{q}$ system the light quark is more relativistic than in $q\bar{q}$, is significant mainly for pseudoscalar mesons. (If, however, the coefficient of $\sigma_0^2$ were replaced by unity these states would shift by only of order 30 MeV.) Using (12) for the running coupling constant and Fourier transforming leads to

\[
G(r) = -\sum_k \frac{4\alpha_k}{3r} \left[ \frac{2}{\sqrt{\pi}} \int_0^{\gamma_k r} e^{-x^2} dx \right] ,
\tag{A10}
\]

and

\[
S(r) = br + c .
\tag{A11}
\]

The result (13) follows from the definition

\[
G(r) = -\frac{4\alpha_s(r)}{3r} .
\]

The resulting smeared potentials are

\[
\tilde{G}(r) = -\sum_k \frac{4\alpha_k}{3r} \left[ \frac{2}{\sqrt{\pi}} \int_0^{\tau_{kij} r} e^{-x^2} dx \right]
\tag{A12}
\]

and

\[
\tilde{S}(r) = br
\left[
\frac{e^{-\sigma_{ij}^2 r^2}}{\sqrt{\pi}\sigma_{ij} r}
+ \left( 1 + \frac{1}{2\sigma_{ij}^2 r^2} \right) \frac{2}{\sqrt{\pi}} \int_0^{\sigma_{ij} r} e^{-x^2} dx
\right] + c ,
\tag{A13}
\]

where

\[
\frac{1}{\tau_{kij}^2} = \frac{1}{\gamma_k^2} + \frac{1}{\sigma_{ij}^2} .
\tag{A14}
\]

On the other hand, we take into account the effect (a) by introducing momentum-dependent factors in the various interactions which go to unity in the nonrelativistic limit to give back the potentials (3)-(7) of the text. Since (A12) and (A13) should already contain the $Q^2$-dependence-induced modifications of the form of the potentials, we examine (A5) and (A6) as $Q^2 \rightarrow 0$ and conclude the following.

(1) The Coulomb term should be modified according to

\[
\tilde{G}(r) \rightarrow \left[1 + \frac{p^2}{EE}\right]^{1/2} \tilde{G}(r) \left[1 + \frac{p^2}{EE}\right]^{1/2} .
\]

(2) The contact, tensor, vector spin-orbit, and scalar spin-orbit potentials should be modified according to

\[
\frac{\tilde{V}_i(r)}{m_1 m_2} \rightarrow \left(\frac{m_1 m_2}{E_1 E_2}\right)^{1/2 + \epsilon_i} \frac{\tilde{V}_i(r)}{m_1 m_2} \left(\frac{m_1 m_2}{E_1 E_2}\right)^{1/2 + \epsilon_i},
\]

where $i=$ contact $(c)$, tensor $(t)$, vector spin-orbit [so$(v)$], scalar spin-orbit [so$(s)$]. If $\epsilon_i=0$ then these modifications have the effect of replacing the nonrelativistic mass dependences $1/m_\alpha m_\beta$ $(=1/m^2,~ 1/\bar{m}^2,~ {\rm or}~ 1/m\bar{m})$ of these potentials by $1/E_\alpha E_\beta$. The parameters $\epsilon_i$, which are therefore expected to be small, are given in Table II.

(3) On the basis of the solutions to the (at least superficially) similar example of QED in one dimension, which confines with a linear potential, we assume that $\tilde{S}(r)$ is unmodified by relativistic corrections.$^{41}$ We remark that if we were to introduce a parameter $\epsilon_{\text{linear}}$ for this potential, we would be forced phenomenologically to set $\epsilon_{\text{linear}} \simeq0$; this observation offers some indirect support for the one-dimensional flux-tube model of confinement in QCD.

Since both experimentally and, partly as a consequence of $m/E$ suppressions, theoretically the spin-orbit interactions are relatively weak, we ignore the "second-order" spin-orbit terms of the form $(\mathbf{Q} \cdot \mathbf{p} \times \mathbf{S}_q)(\mathbf{Q} \cdot \mathbf{p} \times \mathbf{S}_q)$ in both $G_{\rm eff}$ and $S_{\rm eff}$. Reverting (as allowed by the resulting symmetry) to the 1,2 labeling of the text, this leaves us with the approximate forms of the potentials which we use in our calculations. Defining for compactness

\[
f^{i}_{\alpha \beta}(r) = \left(\frac{m_\alpha m_\beta}{E_\alpha E_\beta}\right)^{1/2 + \epsilon_i} f(r) \left(\frac{m_\alpha m_\beta}{E_\alpha E_\beta}\right)^{1/2 + \epsilon_i}
\]

we have

\[
G_{\rm eff}(r) = \left[ 1 + \frac{p^2}{E_1 E_2} \right]^{1/2} \tilde{G}(r) \left[ 1 + \frac{p^2}{E_1 E_2} \right]^{1/2}
\]

\[
+ \left[ \frac{\mathbf{S}_1 \cdot \mathbf{L}}{2 m_1^2} \frac{1}{r} \frac{\partial \tilde{G}_{11}^{\rm so(v)}}{\partial r}
+ \frac{\mathbf{S}_2 \cdot \mathbf{L}}{2 m_2^2} \frac{1}{r} \frac{\partial \tilde{G}_{22}^{\rm so(v)}}{\partial r}
+ \frac{(\mathbf{S}_1 + \mathbf{S}_2) \cdot \mathbf{L}}{m_1 m_2} \frac{1}{r} \frac{\partial \tilde{G}_{12}^{\rm so(v)}}{\partial r}
\right]
\]

\[
+ \frac{2\mathbf{S}_1 \cdot \mathbf{S}_2}{3 m_1 m_2} \nabla^2 \tilde{G}_{12}^{c}
- \left[
\frac{\mathbf{S}_1 \cdot \hat{\mathbf{r}}~ \mathbf{S}_2 \cdot \hat{\mathbf{r}} - \frac{1}{3} \mathbf{S}_1 \cdot \mathbf{S}_2}{m_1 m_2}
\right]
\left[
\frac{\partial^2}{\partial r^2} - \frac{1}{r} \frac{\partial}{\partial r}
\right]
\tilde{G}_{12}^{t}
\tag{A15}
\]
