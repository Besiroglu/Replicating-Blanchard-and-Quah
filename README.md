# Replicating-Blanchard-and-Quah

I replicate the decompositions of output into permanent and transitory shocks in [Blanchard and Quah, 1988](https://uh.edu/~bsorense/BlanchardQuah1989.pdf). I replicate their quantitative results, and find that Blanchard and Quah, 1998 likely over-estimate the unemployment that results from positive supply disturbances but under-estimate the long-run effect on GNP growth of positive supply disturbances. Moreover, by evaluating Bayesian information criteria for a variety of lag-lengths, I find  evidence that Blanchard and Quah might have over-fit their their vector autoregressive (VAR) models.

## Repo Overview

In [Replication.m](https://github.com/Besiroglu/Replicating-Blanchard-and-Quah/blob/main/Replication.m) I replicate the decompositions of output into permanent and transitory shocks in [Blanchard and Quah, 1988](https://uh.edu/~bsorense/BlanchardQuah1989.pdf). Along with Blanchard and Quah, we assume that demand-shocks have no long-run effect on output, and impose this long-run restriction through a Cholesky Decomposition. [Replication.m](https://github.com/Besiroglu/Replicating-Blanchard-and-Quah/blob/main/Replication.m) graphs the asccoiated impulse response response functions for a unit supply and demand disturbances estimated using the base-case  VAR(8)  model (figure reproduced below). Shaded  regions  represent  95%  confidence  interval  generated  by  bootstrapping  with 1,000 replications.

![](https://github.com/Besiroglu/Replicating-Blanchard-and-Quah/blob/main/Images/1.jpg)



We consider lags ranging from 2 (as is suggested by theBayesian information criteria) to 8 (as is adopted in Blanchard and Quah, 1988) and the intermediate numberof lags, spaced out by 2 lags.


![](https://github.com/Besiroglu/Replicating-Blanchard-and-Quah/blob/main/Images/2.jpg)


![](https://github.com/Besiroglu/Replicating-Blanchard-and-Quah/blob/main/Images/3.jpg)


![](https://github.com/Besiroglu/Replicating-Blanchard-and-Quah/blob/main/Images/4.jpg)


![](https://github.com/Besiroglu/Replicating-Blanchard-and-Quah/blob/main/Images/5.jpg)
