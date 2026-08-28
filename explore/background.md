# Takagi-Taupin Equation

$$
\begin{aligned}
\frac{dD_0}{ds} &= \frac{ik}{2}\Big[\chi_0 D_0 + C\,\chi_{\bar h}\,e^{+i\phi} D_h\Big]\\
\frac{dD_h}{ds} &= \frac{ik}{2}\Big[C\,\chi_h\,e^{-i\phi} D_0 + (\chi_0 + 2\delta) D_h\Big]
\end{aligned}
$$

## Eliminate $\chi_0$ (Gauge Transformation)
Redefine wavefunction
### $\tilde{D}_0(s)=D_0(s)e^{-i\frac{k}{2}\chi_0s},\ \tilde{D}_h=D_h(s)e^{-i\frac{k}{2}\chi_0s}$

Differentiate $D_0$ and $D_h$, and insert them to Takagi-Taupin Equation.
The result is

$$
\begin{aligned}
\frac{d\tilde{D}_0}{ds}&=\frac{ik}{2}C\chi_{\bar{h}}e^{+i\phi}\tilde{D}_h\\
\frac{d\tilde{D}_h}{ds}&=\frac{ik}{2}C\chi_he^{-i\phi}\tilde{D_0}+\frac{ik}{2}2\delta\tilde{D}_h
\end{aligned}
$$

The real paths of x-ray ($s_0, s_h$) have a relation $s=\frac{z}{cos\theta_B}$

## Projection to z
$z=s\cos\theta_B$

$\frac{d}{ds}=\cos\theta_B\frac{d}{dz}$

$$
\begin{aligned}
\frac{d\tilde{D}_0}{dz}&=\frac{ik}{2\cos\theta_B}C\chi_{\bar{h}}e^{+i\phi}\tilde{D}_h\\
\frac{d\tilde{D}_h}{dz}&=\frac{ik}{2\cos\theta_B}C\chi_he^{-i\phi}\tilde{D_0}+\frac{ik}{2}2\delta\tilde{D}_h
\end{aligned}
$$

## Nondimensionalization
Substitute $\zeta$ 