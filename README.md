# Documentation for Latent Alpha Model Estimation Code

*Author: Philippe Casgrain* 

Email: <img src="http://safemail.justlikeed.net/e/f5105fb08d53f22bc44ef60f94c8f474.png" border="0" align="center" title="Email image created with safemail.justlikeed.net">



This collection of MATLAB and C code can be used for the estimation of paramters in a pure-jump latent alpha model via an EM algorithm. For more information on these models and on the estimation algorithm, as well as their use in algorithmic trading, see *[Casgrain, Jaimungal (2016)](https://arxiv.org/abs/1806.04472)*[^fn1].


*Note:* Much of the C code found in this repository is based on a C/mex implementaton of the Forward-Backward algorithm for classical HMMs by [Aurélien Garivier](https://www.math.univ-toulouse.fr/~agarivie/?q=node/34). The original code for this algorithm and comparisons of differnt implementations [here](https://www.math.univ-toulouse.fr/~agarivie/Telecom/code/index.php).


## The Price Process Model
We consider a continuous-time model for an asset price process $S_t$, which is driven by a latent hidden markov chain $\Theta$. We assume that the dynamics of this particular model are expressed as
$$
dS_t = \delta \left( dN_t^+ - dN_t^- \right)
\;,
$$
where $\delta>0$ represents the tick size and  $N_t^\pm$ are Poisson processes with respective stochastic intensities $\lambda_t^{\pm}$. We assume that the intensity processes takes the form
$$
\lambda_t^{\pm} = \sigma + 
\kappa(\Theta_t - S_t)_{\pm}
\;,
$$
for $\kappa,\sigma \geq 0$ and where $(x)_{+}$ indicates the positive part of $x\in\mathbb{R}$ and vice-versa. 

The parameter $\sigma$ can be interpreted as controlling the base amount of activity in the process, whereas the parameter $\kappa$ represents the strength of the mean-reversion of the process, which can be seen by computing the expression
$$
\partial_h \mathbb{E}\left[ S_{t+h} \lvert \Theta_t, S_t \right] = \kappa (\Theta_t - S_t)
\;,
$$
which indicates a linear rate of mean reversion of the process $S_t$ towards $\Theta_t$ controlled by the magnitude of the parameter $\kappa$.


The process $\Theta_t$ is assumed to be an $M$-state continuous time Markov chain, taking values in the set $\{\theta_i\}_{i=1}^M$. $\Theta_t$ has a generator matrix $Q$, which defines its transition probabilities through the relation $\mathbb{P}\left(\Theta_{t+h} = \theta_j \lvert \Theta_t = \theta_i \right) = \left( e^{Qh}\right)_{i,j}$. We also assume here that at time $t=0$, the process $\Theta_t$ begins in state $\theta_i$ with probability $\pi_0^i$.

## The Estimation Problem

Given independent observations of discrete paths of $S_t$, we would like to obtain frequentist maximum likelihood estimates of the model parameters $\boldsymbol{\Psi} = (\boldsymbol{\pi}_0,Q,\boldsymbol{\theta},\kappa,\sigma)$. To do this, we must maximize the likelihood of the observed data with respect to this set of parameters. Since the paths of $\Theta_t$ are not known, the likelihood will correspond to the marginal likelihood of $S$, which is not tractable in general. To get around this issue, we instead resort to the use of the EM algorithm, which achieves this by maximizing a computable lower bound to the marginal likelihood.

## Information on the Code

This repo includes some example price path data in `./Example_Data.mat`.
The files which run the EM algorithm on example data are contained in the root of `./PureJumpEstimation/`.

Parts of the code are written in C/mex. Before you run the code for the first time, you should compile these files through mex by `cd`-ing to the appropriate folder and running `HMM_makeMex.m` and `HMM_makeMex.m`.


[^fn1]: Casgrain, P. and S. Jaimungal (2016, Nov). Trading algorithms with learning in latent alpha models. Mathematical Finance, Forthcoming. [arXiv:1806.04472](https://arxiv.org/abs/1806.04472)
