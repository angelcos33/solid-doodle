#### Brief Theoretical Introduction

Chance, Prock and Silbey developed a model to describe dipole emissions and energy transfer near interfaces 
(CPS). In CPS theory excitons are modeled as point dipoles, which allows for the description of their electric fields using a particular class of Green's functions called dyadic Green's functions.

Since LED and OLEDs are multilayer devices, CPS theory is, by itself, unable to describe the energy transfer behavior of such devices due to the need to model the interlayer energy transfer for each of the stack's layers. To remedy this, Celebi, Heidel and Baldo introduced an extension of CPS theory which allows for the analytical determination of the Poynting vector and thus for the calculation of interlayer energy transfer. 

The library presented here applies this extension of CPS theory and uses its main results to provide simple and computationally efficient tools to model several aspects of the energy emission behavior of LEDs and OLEDS, such as angular emission profile, and energy transfer efficiency.

As previously mentioned, the electric fields of the point dipoles can be described by using dyadic Green's functions, this is done by integrating the product of the Green's functions and the current as follows:

\f[
	\mathbf{E}(\mathbf{R}) = \omega \mu_0 \int  \mathbf{J}(\mathbf{R'}) \cdot \mathbf{G}(\mathbf{R}|\mathbf{R'}) \, d\mathbf{R'}

\f]

The dyadic Green's functions for two-dimensionally symmetric multilayer stack is given by two independent sets of eigenfunctions in cylindrical coordinates, as reported by Chance, Prock and Silbey:

\f[
\mathbf{M}_{n\kappa}^e(h) = e^{ihz} \left[ \pm \frac{\eta J_n(\kappa r)}{r} 
\begin{array}{cc}
\sin n\phi \\
\cos n\phi 
\end{array} 
\hat{r} - \frac{\partial J_n(\kappa r)}{\partial r} 
\begin{array}{cc}
\cos n\phi \\
\sin n\phi 
\end{array} 
\hat{\phi} \right]
\f]


\f[
\mathbf{N}_{n\kappa}^e(h) = \frac{e^{ihz}}{k_j} \left[ ih \frac{\partial J_n(\kappa r)}{\partial r} 
\begin{array}{cc}
\cos n\phi \\
\sin n\phi 
\end{array}
\hat{r} \pm inh \frac{J_n(\kappa r)}{r} 
\begin{array}{cc}
\sin n\phi \\
\cos n\phi 
\end{array} 
\hat{\phi} + \kappa^2 J_n(\kappa r) 
\begin{array}{cc}
\cos n\phi \\
\sin n\phi 
\end{array} 
\hat{z} \right]

\f]

Here kappa and h are the amplitudes of the parallel and perpendicular components of the propagation vector k, and Jn refers to a Bessel function of the first type. Moreover, j refers to the layer index, and even and odd eigenfunctions are represented by e and o. For the dipole source and the corresponding scattering, we then obtain:

\f[
\mathbf G_0(\mathbf R \mid \mathbf R \prime ) = \frac{i}{4\pi }\int\limits_0^{\infty}d \kappa \sum_{\begin{array}{cc} n = 0 \\ t = e,o \\  \end{array}}^{\infty}\frac{2- \delta_{n0}}{\kappa h_s}\left[\begin{array}{cc} \mathbf M_{\mathrm{tn\kappa}}(h_s)\mathbf M_{\mathrm{tn\kappa}}^{\prime}(- h_s) + \mathbf N_{\mathrm{tn\kappa}}(h_s)\mathbf N_{\mathrm{tn\kappa}}^{\prime}(- h_s) \\ \mathbf M_{\mathrm{tn\kappa}}(- h_s)\mathbf M_{\mathrm{tn\kappa}}^{\prime}(h_s) + \mathbf N_{\mathrm{tn\kappa}}(- h_s)\mathbf N_{\mathrm{tn\kappa}}^{\prime}(h_s) \\  \end{array}\right]\begin{array}{cc} z \geq 0 \\ z \leq 0 \\  \end{array} 
		
\f]

\f[ \begin{array}{c}
	\mathbf G_j(\mathbf R \mid \mathbf R \prime ) = \frac{i}{4\pi }\int\limits_0^{\infty}d\kappa \sum_{\begin{array}{cc} n = 0 \\ t = e,o \\  \end{array}}^{\infty}\frac{2- \delta_{n0}}{\kappa h_s}[c_j\mathbf M_{\mathrm{tn}\kappa}(- h_j)\mathbf M_{\mathrm{tn}\kappa}^{\prime}(h_s) + f_j\mathbf N_{\mathrm{tn}\kappa}(- h_j)\mathbf N_{\mathrm{tn}\kappa}^{\prime}(h_s) + \\
	 c_j \prime \mathbf M_{\mathrm{tn}\kappa}(h_j)\mathbf M_{\mathrm{tn}\kappa}^{\prime}(h_s) + f_j \prime \mathbf N_{tn\kappa}(h_j)\mathbf N_{\mathrm{tn}\kappa}^{\prime}(h_s)]	\end{array}
\f]

Solving the Maxwell's equation at the interfaces gives rise to a number of relations between the coefficients. These relations allow us to solve for the corresponding coefficients iteratively as implemented in this library's BaseSolver class. 
\f[
	c_je^{- ih_jz_j} + c\prime_je^{ih_jz_j} = c_{j + 1}e^{- ih_{j + 1}z_j} + c\prime_{j + 1}e^{ih_{j + 1}z_j}
\f]

\f[
	\frac{h_j}{k_j}(- f_je^{- ih_jz_j} + f\prime_je^{ih_jz_j}) = \frac{h_{j + 1}}{k_{j + 1}}(- f_{j + 1}e^{- ih_{j + 1}z_j} + f\prime_{j + 1}e^{ih_{j + 1}z_j})
\f]

\f[
	- c_jh_je^{- ih_jz_j} + c\prime_jh_je^{ih_jz_j} = - h_{j + 1}c_{j + 1}e^{- ih_{j + 1}z_j} + h_{j + 1}c\prime_{j + 1}e^{ih_{j + 1}z_j}
\f]

\f[
	k_jf_je^{- ih_jz_j} + k_jf\prime_je^{ih_jz_j} = k_{j + 1}f_{j + 1}e^{- ih_{j + 1}z_j} + k_{j + 1}f\prime_{j + 1}e^{ih_{j + 1}z_j}	
\f]

Since no external sources of radiation are considered. the initial c, c' and f, f' coefficients are zero. Next, we need to calculate the interlayer energy transmission efficiency for both perpendicular and parallel components. The respective formulas for both components are given below:
\f[
	\frac{b^{\bot}}{b_0} = 1- q + q\{ 1 + \frac{3}{2}\mathrm{Re}[\int\limits_0^{\infty}d\kappa \frac{\kappa^3}{h_sk_s^3}(f_s + f_s \prime )]\} 
\f]

\f[
	\frac{b^{\parallel}}{b_0} = 1- q + q\{ 1 + \frac{3}{4}\mathrm{Re}[\int\limits_0^{\infty}d\kappa \frac{\kappa }{h_sk_s}(c_s + c_s \prime  + \frac{h_s^2}{k_s^2}(f_s- f_s \prime ))]\} 
\f]

Finally, by simplifying the integral of the Poynting's vector as follows:

\f[
	\int \nabla \ldotp \mathbf SdV = \oint \mathbf S\ldotp d\mathbf A \approx \int S_zdA\ldotp 
\f]

\f[
S_z = \frac{i}{2\mu_0\omega }[E_r\left(\frac{\partial E_r}{\partial z}- \frac{\partial E_z}{\partial r}\right)^* + E_{\phi}\left(\frac{\partial E_{\phi}}{\partial z}- \frac{1}{r}\frac{\partial E_z}{\partial \phi }\right)^*]
\f]

Finally, the power flow can be determined from the real part of the integral of the Poynting vector with respect to the surface area of the interface:


\f[
	\mathrm{Re}(\int S_{z,j}^{\bot^*}\mathrm{dA}) = \frac{3q}{4}\mathrm{Re}[\int\limits_0^{\infty}\mathrm{du}\frac{u^3(\sqrt{\varepsilon_j})^*}{ \mid 1- u^2 \mid \sqrt{\varepsilon_j}}\left(\frac{\varepsilon_j}{\varepsilon_s- u^2}\right)^{\frac{1}{2}}(f_j \prime e^{ih_jz}- f_je^{- ih_jz})(f_j \prime e^{ih_jz} + f_je^{- ih_jz})^*]
\f]


\f[\begin{array}{c}
	\mathrm{Re}(\int S_{z,j}^{\parallel}\mathop{dA}\limits^*) = \frac{3q}{8}\mathrm{Re}[\int\limits_0^{\infty}du\frac{u(\sqrt{\varepsilon_j})^*}{\sqrt{\varepsilon_j}}\left(\frac{\varepsilon_j}{\varepsilon_s- u^2}\right)^{\frac{1}{2}}(f_j \prime e^{ih_jz}- f_je^{- ih_jz})(f_j \prime e^{ih_jz} + f_je^{- ih_jz})^* + \\
	 \int\limits_0^{\infty}du\frac{u((\frac{\varepsilon_j}{\varepsilon_s- u^2})^{\frac{1}{2}})^*}{ \mid 1- u^2 \mid }(c_j \prime e^{ih_jz} + c_je^{- ih_jz})(c_j \prime e^{ih_jz}- c_je^{- ih_jz})^*] \end{array}
\f]

Where u = kappa/ks. These calculations build the basis for both the simulation and fitting functionalities. Therefore, as previously mentioned, they are implemented in the BaseSolver class, which is inherited by both the Simulation and Fitting classes. Notably, this strong relation between these three classes allowed us to implement the BaseSolver class as an abstract class, such that objects from child classes can be manipulated via pointers to the base class. 
