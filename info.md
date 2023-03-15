
# Basic model
The CMB lensing problem in this repo, returned by `load_cmb_lensing_problem`, is defined by:

$$
\begin{align}
r &\sim \mathcal{P}(r) \propto 1/r \\
A_\phi &\sim \mathcal{P}(A_\phi) \propto 1/A_\phi \\
f &\sim \mathcal{N}(0, \mathbb{C}_f(r)) \\
\phi &\sim \mathcal{N}(0, \mathbb{C}_\phi(A_\phi)) \\
\tilde f &= \mathbb{L}(\phi) f \\
d &\sim \mathcal{N}(\mathbb{A} \tilde f, \mathbb{C}_n)
\end{align} 
$$

where

* $r$ is the tensor-to-scalar ratio (scales the tensor contribution) to $B$-mode maps
* $A_\phi$ is the lensing amplitudes (scales total lensing covariance, $\mathbb{C}_\phi \rightarrow A_\phi \mathbb{C}_\phi$)
* $f$ is the unlensed CMB polarization
* $\phi$ is the lensing potential 
* $\mathbb{C}_f$, $\mathbb{C}_\phi$, and $\mathbb{C}_n$, and the unlensed CMB, lensing, and noise covariances
* $\mathbb{L}$ is the lensing operator (implemented with LenseFlow)
* $\mathbb{A}$ includes beams and an $\ell_{\rm max}=3000$ highpass. 

The simulated data has noise levels similar to CMB-S4, including:

* $1\,\mu{\rm K\,arcmin}$ noise in temperature ($\sqrt 2$ higher in polarization)
* $1\,{\rm arcmin}$ FWHM Gaussian beams
* $3\,{\rm arcmin}$ pixels
* The flat-sky approximation
* true values of $r=0.2$ and $A_\phi=1$. 

Using the shorthand $\theta\equiv(r,A_\phi)$, this fully describes the posterior in terms of the physical variables:

$$\log\mathcal{P}(f,\phi,\theta\,|\,d)$$


# Change-of-variables

 We further perform two changes-of-variables to help Gaussianize the posterior.

First, we do the standard lensing change-of-variables ([Millea et al, 2022](https://arxiv.org/abs/2002.00965)):

$$
\begin{align}
f^\prime &= \mathbb{L}(\phi) \, \mathbb{D}(r) \, f \\
\phi^\prime &= \mathbb{G}(A_\phi) \, \phi
\end{align}
$$

where 

$$
\begin{align}
\mathbb{G}(A_\phi) &= \frac{\big(\mathbb{I} + 2 \, \mathbb{N}_\phi \, \mathbb{C}_\phi(A_\phi)^{-1}\big)^{\frac{1}{2}}}{\big(\mathbb{I} + 2 \, \mathbb{N}_\phi \, \mathbb{C}_\phi(A_\phi=1)^{-1}\big)^\frac{1}{2}} \\
\mathbb{D}(r) &= \frac{\big( \mathbb{C}_f(r) + \mathbb{N}_{\rm len} + 2\,\mathbb{C}_n \big)^{\frac{1}{2}}}{ \mathbb{C}_f(r)^{\frac{1}{2}} }
\end{align}
$$

and $\mathbb{N}_\phi$ is an estimate of the effective reconstruction noise which we take to be the $N_0$ quadratic estimator noise and $\mathbb{N}_{\rm len}$ is an estimate of lensed $B$-mode noise, which we take to be $5\,\mu{\rm K\,arcmin}$ white noise. 

This leads to an $f^\prime$ and $\phi^\prime$ which are much less correlated than $f$ and $\phi$, and also removes some (but not all) of the funnel structure in $(r,f)$ and $(A_\phi, \phi)$.

Second, we take $\theta^\prime = \log \theta$. This gives $\theta^\prime$ support on $(-\infty,\infty)$, rather than just $(0,\infty)$. Given this change-of-variables and the original prior on $\theta$, note that this corresponds to a uniform prior on $\log \theta^\prime$.

The necessary jacobian determinant factors which arise from these two changes-of-variables are included.

This gives us the final posterior which we work with:

$$
\mathcal{P}(f^\prime, \phi^\prime, \theta^\prime\,|\,d)
$$
