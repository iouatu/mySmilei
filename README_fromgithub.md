# mySmilei
includes TunnelBSIrate in its correct (standard_rate * au_to_wc) form

######################################################


* 3 rates implemented, as in I. Yu. Kostyukov, A. A. Golovanov (2019): https://arxiv.org/abs/1906.01358
1. Tunnel rate (PPT): already implemented in SMILEI
2. Quadratic Bauer-Mulser rate, Phys. Rev. A 59, 569 (1999) https://journals.aps.org/pra/abstract/10.1103/PhysRevA.59.569
3. Linear I. I. Artmenko, I. Yu. Kostyukov Phys. Rev. A 96, 032106 (2017) https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.032106

* int continuity_tool( , ... , ) helps choosing from the 3 rates above, at each timestep

* Multiple Ionizations during 1 timestep can happen for any of the rates above, 
* The exact procedure is described in R.Nuter et al. PoP 19, 033107 (2011) and was already implemented in SMILEI.
