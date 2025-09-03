

# File canopy\_bioemi\_mod.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_bioemi\_mod.F90**](canopy__bioemi__mod_8F90.md)

[Go to the source code of this file](canopy__bioemi__mod_8F90_source.md)

_Biogenic Emissions Module._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**canopy\_bioemi\_mod**](namespacecanopy__bioemi__mod.md) <br> |








## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**above**](#variable-above)  <br> |
|  real(rk) | [**agro**](#variable-agro)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**air**](#variable-air)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), dimension(default=0/viirs), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**al**](#variable-al)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**all**](#variable-all)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**alpha\_p\_shade**](#variable-alpha_p_shade)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**alpha\_p\_sun**](#variable-alpha_p_sun)  <br> |
|  real(rk) | [**amat**](#variable-amat)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**and**](#variable-and)  <br> |
|  real(rk) | [**anew**](#variable-anew)  <br> |
|  real(rk) | [**aold**](#variable-aold)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**applied**](#variable-applied)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**applying**](#variable-applying)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**aq**](#variable-aq)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**aqopt**](#variable-aqopt)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**area**](#variable-area)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**atmospheric**](#variable-atmospheric)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ave**](#variable-ave)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**average**](#variable-average)  <br> |
|  real(rk) | [**beta**](#variable-beta)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**between**](#variable-between)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**biogenic**](#variable-biogenic)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**biogenics**](#variable-biogenics)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**calculated**](#variable-calculated)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**calculation**](#variable-calculation)  <br> |
|  real(rk) | [**canloss\_fac**](#variable-canloss_fac)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**canopy**](#variable-canopy)  <br> |
|  real(rk) | [**caq**](#variable-caq)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cce**](#variable-cce)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cell**](#variable-cell)  <br> |
|  real(rk) | [**ceo**](#variable-ceo)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**chemical**](#variable-chemical)  <br> |
|  real(rk) | [**cht**](#variable-cht)  <br> |
|  real(rk) | [**chw**](#variable-chw)  <br> |
|  real(rk) | [**clt**](#variable-clt)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cm**](#variable-cm)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cm2**](#variable-cm2)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**co2**](#variable-co2)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**co2opt**](#variable-co2opt)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**co2set**](#variable-co2set)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**coefficient**](#variable-coefficient)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**conc**](#variable-conc)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**constant**](#variable-constant)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**correction**](#variable-correction)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**cp\_shade**](#variable-cp_shade)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**cp\_sun**](#variable-cp_sun)  <br> |
|  real(rk) | [**ct1**](#variable-ct1)  <br> |
|  real(rk), parameter | [**ct2**](#variable-ct2)   = `230.0\_rk`<br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**current**](#variable-current)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**currentlai**](#variable-currentlai)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**daily**](#variable-daily)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**daily\_maxt2**](#variable-daily_maxt2)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**daily\_maxws10**](#variable-daily_maxws10)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**daily\_mint2**](#variable-daily_mint2)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**days**](#variable-days)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**depth**](#variable-depth)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**dominant**](#variable-dominant)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**downward**](#variable-downward)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**dswrf**](#variable-dswrf)  <br> |
|  real(rk) | [**dtaq**](#variable-dtaq)  <br> |
|  real(rk) | [**dtht**](#variable-dtht)  <br> |
|  real(rk) | [**dthw**](#variable-dthw)  <br> |
|  real(rk) | [**dtlt**](#variable-dtlt)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**e\_opt**](#variable-e_opt)  <br> |
|  real(rk) | [**ef**](#variable-ef)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**emi\_ind**](#variable-emi_ind)  <br> |
|  real(rk), dimension(:), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**emi\_out**](#variable-emi_out)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), dimension(kg [**m**](canopy__bioemi__mod_8F90.md#variable-m)-3 s-1), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**emissions**](#variable-emissions)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**environment**](#variable-environment)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**et**](#variable-et)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**factor**](#variable-factor)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fch**](#variable-fch)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fclai**](#variable-fclai)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**flai**](#variable-flai)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**for**](#variable-for)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fraction**](#variable-fraction)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**friction**](#variable-friction)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**from**](#variable-from)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fsun**](#variable-fsun)  <br> |
|  real(rk) | [**gammaaq**](#variable-gammaaq)  <br> |
|  real(rk) | [**gammaco2**](#variable-gammaco2)  <br> |
|  real(rk) | [**gammaht**](#variable-gammaht)  <br> |
|  real(rk) | [**gammahw**](#variable-gammahw)  <br> |
|  real(rk) | [**gammaleafage**](#variable-gammaleafage)  <br> |
|  real(rk) | [**gammalt**](#variable-gammalt)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammappfd\_shade\_ldf**](#variable-gammappfd_shade_ldf)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammappfd\_sun\_ldf**](#variable-gammappfd_sun_ldf)  <br> |
|  real(rk) | [**gammasoim**](#variable-gammasoim)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_ppfd\_ave**](#variable-gammatleaf_ppfd_ave)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_ppfd\_ldf**](#variable-gammatleaf_ppfd_ldf)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_ppfd\_lif**](#variable-gammatleaf_ppfd_lif)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_shade\_ldf**](#variable-gammatleaf_shade_ldf)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_shade\_ldf\_den**](#variable-gammatleaf_shade_ldf_den)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_shade\_ldf\_num**](#variable-gammatleaf_shade_ldf_num)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_shade\_lif**](#variable-gammatleaf_shade_lif)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_sun\_ldf**](#variable-gammatleaf_sun_ldf)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_sun\_ldf\_den**](#variable-gammatleaf_sun_ldf_den)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_sun\_ldf\_num**](#variable-gammatleaf_sun_ldf_num)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gammatleaf\_sun\_lif**](#variable-gammatleaf_sun_lif)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**gauss**](#variable-gauss)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**gfs**](#variable-gfs)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**grid**](#variable-grid)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**height**](#variable-height)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**heights**](#variable-heights)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**high**](#variable-high)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**hours**](#variable-hours)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**hr**](#variable-hr)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**htopt**](#variable-htopt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**hwopt**](#variable-hwopt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**i**](#variable-i)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**index**](#variable-index)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), dimension(&gt; 0), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**indices**](#variable-indices)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**inhibition**](#variable-inhibition)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**input**](#variable-input)  <br> |
|  integer, intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**integer**](#variable-integer)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**integration**](#variable-integration)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**interpolated**](#variable-interpolated)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lai**](#variable-lai)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**layer**](#variable-layer)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**layers**](#variable-layers)  <br> |
|  real(rk) | [**ldf**](#variable-ldf)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**leaf**](#variable-leaf)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), dimension(0=on, 1=off i.e. [**gammaleafage**](canopy__bioemi__mod_8F90.md#variable-gammaleafage)=1, [**in**](canopy__calcs_8F90.md#variable-in) canopy\_readnml.f90), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**leafage\_opt**](#variable-leafage_opt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**leafageopt**](#variable-leafageopt)  <br> |
|  real(rk), dimension(umol phot/m2 s), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**leaves**](#variable-leaves)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lifetime**](#variable-lifetime)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**loss**](#variable-loss)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**loss\_opt**](#variable-loss_opt)   = `2 (Default = 0.96)`<br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lossind**](#variable-lossind)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lossopt**](#variable-lossopt)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lossset**](#variable-lossset)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**low**](#variable-low)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ltopt**](#variable-ltopt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lu**](#variable-lu)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lu\_opt**](#variable-lu_opt)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**m**](#variable-m)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**m3**](#variable-m3)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**mapped**](#variable-mapped)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**massman**](#variable-massman)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**maximum**](#variable-maximum)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**megan**](#variable-megan)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**minimum**](#variable-minimum)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**model**](#variable-model)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**modlays**](#variable-modlays)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**modres**](#variable-modres)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**moisture**](#variable-moisture)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**number**](#variable-number)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**of**](#variable-of)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**only**](#variable-only)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), dimension(default=0/no [**integration**](canopy__bioemi__mod_8F90.md#variable-integration)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**option**](#variable-option)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**or**](#variable-or)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**output**](#variable-output)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ozone**](#variable-ozone)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**past**](#variable-past)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pastlai**](#variable-pastlai)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**photolysis**](#variable-photolysis)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**point**](#variable-point)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd**](#variable-ppfd)  <br> |
|  real(rk), parameter | [**ppfd0\_shade**](#variable-ppfd0_shade)   = `50.0`<br> |
|  real(rk), parameter | [**ppfd0\_sun**](#variable-ppfd0_sun)   = `200.0`<br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd240\_shade**](#variable-ppfd240_shade)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd240\_sun**](#variable-ppfd240_sun)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd24\_shade**](#variable-ppfd24_shade)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd24\_sun**](#variable-ppfd24_sun)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd\_shade**](#variable-ppfd_shade)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd\_sun**](#variable-ppfd_sun)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppm**](#variable-ppm)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppmv**](#variable-ppmv)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**proportion**](#variable-proportion)  <br> |
|  real(rk), dimension(w/m2), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**radiation**](#variable-radiation)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**resolution**](#variable-resolution)  <br> |
|  real(rk) | [**roota**](#variable-roota)  <br> |
|  real(rk) | [**rootb**](#variable-rootb)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**set**](#variable-set)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**shaded**](#variable-shaded)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**shortwave**](#variable-shortwave)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soid1**](#variable-soid1)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soid2**](#variable-soid2)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soid3**](#variable-soid3)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soid4**](#variable-soid4)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soil**](#variable-soil)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soim1**](#variable-soim1)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soim2**](#variable-soim2)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soim3**](#variable-soim3)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soim4**](#variable-soim4)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**soimopt**](#variable-soimopt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**specie**](#variable-specie)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), dimension(=0), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**species**](#variable-species)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**specific**](#variable-specific)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**speed**](#variable-speed)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**stress**](#variable-stress)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**summing**](#variable-summing)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**sunlit**](#variable-sunlit)  <br> |
|  real(rk) | [**tabovecanopy**](#variable-tabovecanopy)  <br> |
|  real(rk) | [**taq**](#variable-taq)  <br> |
|  real(rk), dimension([**k**](canopy__calcs_8F90.md#variable-k)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**temp**](#variable-temp)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**temp2**](#variable-temp2)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**temperature**](#variable-temperature)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**the**](#variable-the)  <br> |
|  real(rk) | [**tht**](#variable-tht)  <br> |
|  real(rk) | [**thw**](#variable-thw)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**tka**](#variable-tka)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**tleaf240\_ave**](#variable-tleaf240_ave)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**tleaf24\_ave**](#variable-tleaf24_ave)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**tleaf\_opt**](#variable-tleaf_opt)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**tleaf\_shade**](#variable-tleaf_shade)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**tleaf\_sun**](#variable-tleaf_sun)  <br> |
|  real(rk) | [**tlt**](#variable-tlt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**to**](#variable-to)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**top**](#variable-top)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**total**](#variable-total)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**tsteplai**](#variable-tsteplai)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**type**](#variable-type)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**used**](#variable-used)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**user**](#variable-user)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ustar**](#variable-ustar)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**value**](#variable-value)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vegetation**](#variable-vegetation)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**velocity**](#variable-velocity)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vert**](#variable-vert)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vertical**](#variable-vertical)  <br> |
|  real(rk), dimension(s), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**voc**](#variable-voc)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**volume**](#variable-volume)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**volumetric**](#variable-volumetric)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**vpgwt**](#variable-vpgwt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vtype**](#variable-vtype)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**w126**](#variable-w126)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**w126\_ref**](#variable-w126_ref)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**w126\_set**](#variable-w126_set)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**when**](#variable-when)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**wilt**](#variable-wilt)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**wilting**](#variable-wilting)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**wind**](#variable-wind)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**with**](#variable-with)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zk**](#variable-zk)  <br> |












































## Detailed Description


This module contains subroutines for calculating parameterized canopy biogenic emissions based on the algorithms described in Clifton et al. (2022). The module handles various biogenic volatile organic compounds (BVOCs) including isoprene, monoterpenes, and other biogenic species.




**Author:**

Patrick C. Campbell 




**Date:**

January 2023




**Version:**


* Jan 2023 P.C. Campbell: Initial canopy isoprene only version
* Feb 2023 P.C. Campbell: Modified for multiple biogenic species
* Jul 2023 P.C. Campbell: Restructured to use FSUN, TLEAF, and PPFD as inputs
* Sept 2023 QZ Rasool: Modifications for LeafAge Response for multiple BVOCs




\references Clifton, O. E. et al. (2022). Large eddy simulation for investigating coupled forest canopy and turbulence influences on atmospheric chemistry. Journal of Advances in Modeling Earth Systems, 14, e2022MS003078. [https://doi.org/10.1029/2022MS003078](https://doi.org/10.1029/2022MS003078) 


    
## Public Attributes Documentation




### variable above 

```Fortran
real(rk), intent(in) above;
```




<hr>



### variable agro 

```Fortran
real(rk) agro;
```




<hr>



### variable air 

```Fortran
real(rk), intent(in) air;
```




<hr>



### variable al 

```Fortran
integer, dimension (default = 0/viirs), intent(in) al;
```




<hr>



### variable all 

```Fortran
integer, intent(in) all;
```




<hr>



### variable alpha\_p\_shade 

```Fortran
real(rk), dimension(size(zk)) alpha_p_shade;
```




<hr>



### variable alpha\_p\_sun 

```Fortran
real(rk), dimension(size(zk)) alpha_p_sun;
```




<hr>



### variable amat 

```Fortran
real(rk) amat;
```




<hr>



### variable and 

```Fortran
real(rk), intent(in) and;
```




<hr>



### variable anew 

```Fortran
real(rk) anew;
```




<hr>



### variable aold 

```Fortran
real(rk) aold;
```




<hr>



### variable applied 

```Fortran
real(rk), intent(in) applied;
```




<hr>



### variable applying 

```Fortran
integer, intent(in) applying;
```




<hr>



### variable aq 

```Fortran
integer, intent(in) aq;
```




<hr>



### variable aqopt 

```Fortran
integer, intent(in) aqopt;
```




<hr>



### variable area 

```Fortran
real(rk), intent(in) area;
```




<hr>



### variable atmospheric 

```Fortran
real(rk), intent(in) atmospheric;
```




<hr>



### variable ave 

```Fortran
real(rk), intent(in) ave;
```




<hr>



### variable average 

```Fortran
real(rk), intent(in) average;
```




<hr>



### variable beta 

```Fortran
real(rk) beta;
```




<hr>



### variable between 

```Fortran
real(rk), intent(in) between;
```




<hr>



### variable biogenic 

```Fortran
integer, intent(in) biogenic;
```




<hr>



### variable biogenics 

```Fortran
integer, intent(in) biogenics;
```




<hr>



### variable calculated 

```Fortran
real(rk), intent(in) calculated;
```




<hr>



### variable calculation 

```Fortran
integer, intent(in) calculation;
```




<hr>



### variable canloss\_fac 

```Fortran
real(rk) canloss_fac;
```




<hr>



### variable canopy 

```Fortran
real(rk), intent(out) canopy;
```




<hr>



### variable caq 

```Fortran
real(rk) caq;
```




<hr>



### variable cce 

```Fortran
real(rk), intent(in) cce;
```




<hr>



### variable cell 

```Fortran
integer, intent(in) cell;
```




<hr>



### variable ceo 

```Fortran
real(rk) ceo;
```




<hr>



### variable chemical 

```Fortran
real(rk), intent(in) chemical;
```




<hr>



### variable cht 

```Fortran
real(rk) cht;
```




<hr>



### variable chw 

```Fortran
real(rk) chw;
```




<hr>



### variable clt 

```Fortran
real(rk) clt;
```




<hr>



### variable cm 

```Fortran
real(rk), intent(in) cm;
```




<hr>



### variable cm2 

```Fortran
real(rk), intent(in) cm2;
```




<hr>



### variable co2 

```Fortran
real(rk), intent(in) co2;
```




<hr>



### variable co2opt 

```Fortran
integer, intent(in) co2opt;
```




<hr>



### variable co2set 

```Fortran
real(rk), intent(in) co2set;
```




<hr>



### variable coefficient 

```Fortran
real(rk), intent(in) coefficient;
```




<hr>



### variable conc 

```Fortran
real(rk), intent(in) conc;
```




<hr>



### variable constant 

```Fortran
real(rk), intent(in) constant;
```




<hr>



### variable correction 

```Fortran
real(rk), intent(in) correction;
```




<hr>



### variable cp\_shade 

```Fortran
real(rk), dimension(size(zk)) cp_shade;
```




<hr>



### variable cp\_sun 

```Fortran
real(rk), dimension(size(zk)) cp_sun;
```




<hr>



### variable ct1 

```Fortran
real(rk) ct1;
```




<hr>



### variable ct2 

```Fortran
real(rk), parameter ct2;
```




<hr>



### variable current 

```Fortran
real(rk), intent(in) current;
```




<hr>



### variable currentlai 

```Fortran
real(rk), intent(in) currentlai;
```




<hr>



### variable daily 

```Fortran
real(rk), intent(in) daily;
```




<hr>



### variable daily\_maxt2 

```Fortran
real(rk), intent(in) daily_maxt2;
```




<hr>



### variable daily\_maxws10 

```Fortran
real(rk), intent(in) daily_maxws10;
```




<hr>



### variable daily\_mint2 

```Fortran
real(rk), intent(in) daily_mint2;
```




<hr>



### variable days 

```Fortran
real(rk), intent(in) days;
```




<hr>



### variable depth 

```Fortran
real(rk), intent(in) depth;
```




<hr>



### variable dominant 

```Fortran
integer, intent(in) dominant;
```




<hr>



### variable downward 

```Fortran
real(rk), intent(in) downward;
```




<hr>



### variable dswrf 

```Fortran
real(rk), intent(in) dswrf;
```




<hr>



### variable dtaq 

```Fortran
real(rk) dtaq;
```




<hr>



### variable dtht 

```Fortran
real(rk) dtht;
```




<hr>



### variable dthw 

```Fortran
real(rk) dthw;
```




<hr>



### variable dtlt 

```Fortran
real(rk) dtlt;
```




<hr>



### variable e\_opt 

```Fortran
real(rk), dimension(size(zk)) e_opt;
```




<hr>



### variable ef 

```Fortran
real(rk) ef;
```




<hr>



### variable emi\_ind 

```Fortran
integer, intent(in) emi_ind;
```




<hr>



### variable emi\_out 

```Fortran
real(rk), dimension(:), intent(out) emi_out;
```




<hr>



### variable emissions 

```Fortran
real(rk), dimension (kg m-3 s-1), intent(out) emissions;
```




<hr>



### variable environment 

```Fortran
real(rk), intent(in) environment;
```




<hr>



### variable et 

```Fortran
integer, intent(in) et;
```




<hr>



### variable factor 

```Fortran
integer, intent(in) factor;
```




<hr>



### variable fch 

```Fortran
real(rk), intent(in) fch;
```




<hr>



### variable fclai 

```Fortran
real(rk), dimension(:), intent(in) fclai;
```




<hr>



### variable flai 

```Fortran
real(rk), dimension(size(zk)) flai;
```




<hr>



### variable for 

```Fortran
integer, intent(in) for;
```




<hr>



### variable fraction 

```Fortran
real(rk), intent(in) fraction;
```




<hr>



### variable friction 

```Fortran
real(rk), intent(in) friction;
```




<hr>



### variable from 

```Fortran
integer, intent(in) from;
```




<hr>



### variable fsun 

```Fortran
real(rk), dimension(:), intent(in) fsun;
```




<hr>



### variable gammaaq 

```Fortran
real(rk) gammaaq;
```




<hr>



### variable gammaco2 

```Fortran
real(rk) gammaco2;
```




<hr>



### variable gammaht 

```Fortran
real(rk) gammaht;
```




<hr>



### variable gammahw 

```Fortran
real(rk) gammahw;
```




<hr>



### variable gammaleafage 

```Fortran
real(rk) gammaleafage;
```




<hr>



### variable gammalt 

```Fortran
real(rk) gammalt;
```




<hr>



### variable gammappfd\_shade\_ldf 

```Fortran
real(rk), dimension(size(zk)) gammappfd_shade_ldf;
```




<hr>



### variable gammappfd\_sun\_ldf 

```Fortran
real(rk), dimension(size(zk)) gammappfd_sun_ldf;
```




<hr>



### variable gammasoim 

```Fortran
real(rk) gammasoim;
```




<hr>



### variable gammatleaf\_ppfd\_ave 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_ppfd_ave;
```




<hr>



### variable gammatleaf\_ppfd\_ldf 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_ppfd_ldf;
```




<hr>



### variable gammatleaf\_ppfd\_lif 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_ppfd_lif;
```




<hr>



### variable gammatleaf\_shade\_ldf 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_shade_ldf;
```




<hr>



### variable gammatleaf\_shade\_ldf\_den 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_shade_ldf_den;
```




<hr>



### variable gammatleaf\_shade\_ldf\_num 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_shade_ldf_num;
```




<hr>



### variable gammatleaf\_shade\_lif 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_shade_lif;
```




<hr>



### variable gammatleaf\_sun\_ldf 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_sun_ldf;
```




<hr>



### variable gammatleaf\_sun\_ldf\_den 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_sun_ldf_den;
```




<hr>



### variable gammatleaf\_sun\_ldf\_num 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_sun_ldf_num;
```




<hr>



### variable gammatleaf\_sun\_lif 

```Fortran
real(rk), dimension(size(zk)) gammatleaf_sun_lif;
```




<hr>



### variable gauss 

```Fortran
real(rk), dimension(size(zk)) gauss;
```




<hr>



### variable gfs 

```Fortran
real(rk), intent(in) gfs;
```




<hr>



### variable grid 

```Fortran
integer, intent(in) grid;
```




<hr>



### variable height 

```Fortran
real(rk), dimension (m), intent(in) height;
```




<hr>



### variable heights 

```Fortran
real(rk), dimension (m), intent(in) heights;
```




<hr>



### variable high 

```Fortran
integer, intent(in) high;
```




<hr>



### variable hours 

```Fortran
real(rk), intent(in) hours;
```




<hr>



### variable hr 

```Fortran
real(rk), intent(in) hr;
```




<hr>



### variable htopt 

```Fortran
integer, intent(in) htopt;
```




<hr>



### variable hwopt 

```Fortran
integer, intent(in) hwopt;
```




<hr>



### variable i 

```Fortran
integer i;
```




<hr>



### variable index 

```Fortran
integer, intent(in) index;
```




<hr>



### variable indices 

```Fortran
integer, dimension (> 0), intent(in) indices;
```




<hr>



### variable inhibition 

```Fortran
integer, intent(in) inhibition;
```




<hr>



### variable input 

```Fortran
integer, intent(in) input;
```




<hr>



### variable integer 

```Fortran
integer, intent(in) integer;
```




<hr>



### variable integration 

```Fortran
integer, intent(in) integration;
```




<hr>



### variable interpolated 

```Fortran
real(rk), intent(in) interpolated;
```




<hr>



### variable lai 

```Fortran
integer save lai;
```




<hr>



### variable layer 

```Fortran
real(rk), intent(out) layer;
```




<hr>



### variable layers 

```Fortran
integer, intent(in) layers;
```




<hr>



### variable ldf 

```Fortran
real(rk) ldf;
```




<hr>



### variable leaf 

```Fortran
real(rk), intent(in) leaf;
```




<hr>



### variable leafage\_opt 

```Fortran
integer, dimension (0= on, 1= off i.e. gammaleafage =1, in canopy_readnml.f90), intent(in) leafage_opt;
```




<hr>



### variable leafageopt 

```Fortran
integer, intent(in) leafageopt;
```




<hr>



### variable leaves 

```Fortran
real(rk), dimension (umol phot/m2 s), intent(in) leaves;
```




<hr>



### variable lifetime 

```Fortran
real(rk), intent(in) lifetime;
```




<hr>



### variable loss 

```Fortran
integer, intent(in) loss;
```




<hr>



### variable loss\_opt 

```Fortran
real(rk), intent(in) loss_opt;
```




<hr>



### variable lossind 

```Fortran
integer, intent(in) lossind;
```




<hr>



### variable lossopt 

```Fortran
integer, intent(in) lossopt;
```




<hr>



### variable lossset 

```Fortran
real(rk), intent(in) lossset;
```




<hr>



### variable low 

```Fortran
integer, intent(in) low;
```




<hr>



### variable ltopt 

```Fortran
integer, intent(in) ltopt;
```




<hr>



### variable lu 

```Fortran
integer, intent(in) lu;
```




<hr>



### variable lu\_opt 

```Fortran
integer, intent(in) lu_opt;
```




<hr>



### variable m 

```Fortran
real(rk), intent(in) m;
```




<hr>



### variable m3 

```Fortran
real(rk), intent(in) m3;
```




<hr>



### variable mapped 

```Fortran
integer, intent(in) mapped;
```




<hr>



### variable massman 

```Fortran
integer, intent(in) massman;
```




<hr>



### variable maximum 

```Fortran
real(rk), intent(in) maximum;
```




<hr>



### variable megan 

```Fortran
integer, intent(in) megan;
```




<hr>



### variable minimum 

```Fortran
real(rk), intent(in) minimum;
```




<hr>



### variable model 

```Fortran
integer save model;
```




<hr>



### variable modlays 

```Fortran
integer, intent(in) modlays;
```




<hr>



### variable modres 

```Fortran
real(rk), intent(in) modres;
```




<hr>



### variable moisture 

```Fortran
real(rk), intent(in) moisture;
```




<hr>



### variable number 

```Fortran
integer save number;
```




<hr>



### variable of 

```Fortran
integer save of;
```




<hr>



### variable only 

```Fortran
integer, intent(in) only;
```




<hr>



### variable option 

```Fortran
integer, dimension (default = 0/no integration), intent(in) option;
```




<hr>



### variable or 

```Fortran
integer, intent(in) or;
```




<hr>



### variable output 

```Fortran
real(rk), intent(out) output;
```




<hr>



### variable ozone 

```Fortran
real(rk), intent(in) ozone;
```




<hr>



### variable past 

```Fortran
real(rk), intent(in) past;
```




<hr>



### variable pastlai 

```Fortran
real(rk), intent(in) pastlai;
```




<hr>



### variable photolysis 

```Fortran
real(rk), intent(in) photolysis;
```




<hr>



### variable point 

```Fortran
real(rk), intent(in) point;
```




<hr>



### variable ppfd 

```Fortran
real(rk), intent(in) ppfd;
```




<hr>



### variable ppfd0\_shade 

```Fortran
real(rk), parameter ppfd0_shade;
```




<hr>



### variable ppfd0\_sun 

```Fortran
real(rk), parameter ppfd0_sun;
```




<hr>



### variable ppfd240\_shade 

```Fortran
real(rk), dimension(:), intent(in) ppfd240_shade;
```




<hr>



### variable ppfd240\_sun 

```Fortran
real(rk), dimension(:), intent(in) ppfd240_sun;
```




<hr>



### variable ppfd24\_shade 

```Fortran
real(rk), dimension(:), intent(in) ppfd24_shade;
```




<hr>



### variable ppfd24\_sun 

```Fortran
real(rk), dimension(:), intent(in) ppfd24_sun;
```




<hr>



### variable ppfd\_shade 

```Fortran
real(rk), dimension(:), intent(in) ppfd_shade;
```




<hr>



### variable ppfd\_sun 

```Fortran
real(rk), dimension(:), intent(in) ppfd_sun;
```




<hr>



### variable ppm 

```Fortran
real(rk), intent(in) ppm;
```




<hr>



### variable ppmv 

```Fortran
real(rk), intent(in) ppmv;
```




<hr>



### variable proportion 

```Fortran
real(rk), intent(in) proportion;
```




<hr>



### variable radiation 

```Fortran
real(rk), dimension (w/m2), intent(in) radiation;
```




<hr>



### variable resolution 

```Fortran
real(rk), dimension (m), intent(in) resolution;
```




<hr>



### variable roota 

```Fortran
real(rk) roota;
```




<hr>



### variable rootb 

```Fortran
real(rk) rootb;
```




<hr>



### variable set 

```Fortran
real(rk), intent(in) set;
```




<hr>



### variable shaded 

```Fortran
real(rk), intent(in) shaded;
```




<hr>



### variable shortwave 

```Fortran
real(rk), intent(in) shortwave;
```




<hr>



### variable soid1 

```Fortran
real(rk), intent(in) soid1;
```




<hr>



### variable soid2 

```Fortran
real(rk), intent(in) soid2;
```




<hr>



### variable soid3 

```Fortran
real(rk), intent(in) soid3;
```




<hr>



### variable soid4 

```Fortran
real(rk), intent(in) soid4;
```




<hr>



### variable soil 

```Fortran
real(rk), intent(in) soil;
```




<hr>



### variable soim1 

```Fortran
real(rk), intent(in) soim1;
```




<hr>



### variable soim2 

```Fortran
real(rk), intent(in) soim2;
```




<hr>



### variable soim3 

```Fortran
real(rk), intent(in) soim3;
```




<hr>



### variable soim4 

```Fortran
real(rk), intent(in) soim4;
```




<hr>



### variable soimopt 

```Fortran
integer, intent(in) soimopt;
```




<hr>



### variable specie 

```Fortran
integer, intent(in) specie;
```




<hr>



### variable species 

```Fortran
integer, dimension (=0), intent(in) species;
```




<hr>



### variable specific 

```Fortran
integer, intent(in) specific;
```




<hr>



### variable speed 

```Fortran
real(rk), intent(in) speed;
```




<hr>



### variable stress 

```Fortran
integer, intent(in) stress;
```




<hr>



### variable summing 

```Fortran
integer, intent(in) summing;
```




<hr>



### variable sunlit 

```Fortran
real(rk), intent(in) sunlit;
```




<hr>



### variable tabovecanopy 

```Fortran
real(rk) tabovecanopy;
```




<hr>



### variable taq 

```Fortran
real(rk) taq;
```




<hr>



### variable temp 

```Fortran
real(rk), dimension (k), intent(in) temp;
```




<hr>



### variable temp2 

```Fortran
real(rk), intent(in) temp2;
```




<hr>



### variable temperature 

```Fortran
real(rk), dimension (m/s), intent(in) temperature;
```




<hr>



### variable the 

```Fortran
integer save the;
```




<hr>



### variable tht 

```Fortran
real(rk) tht;
```




<hr>



### variable thw 

```Fortran
real(rk) thw;
```




<hr>



### variable tka 

```Fortran
real(rk), dimension(:), intent(in) tka;
```




<hr>



### variable tleaf240\_ave 

```Fortran
real(rk), dimension(:), intent(in) tleaf240_ave;
```




<hr>



### variable tleaf24\_ave 

```Fortran
real(rk), dimension(:), intent(in) tleaf24_ave;
```




<hr>



### variable tleaf\_opt 

```Fortran
real(rk), dimension(size(zk)) tleaf_opt;
```




<hr>



### variable tleaf\_shade 

```Fortran
real(rk), dimension(:), intent(in) tleaf_shade;
```




<hr>



### variable tleaf\_sun 

```Fortran
real(rk), dimension(:), intent(in) tleaf_sun;
```




<hr>



### variable tlt 

```Fortran
real(rk) tlt;
```




<hr>



### variable to 

```Fortran
integer, intent(in) to;
```




<hr>



### variable top 

```Fortran
integer, intent(in) top;
```




<hr>



### variable total 

```Fortran
integer, intent(in) total;
```




<hr>



### variable tsteplai 

```Fortran
real(rk), intent(in) tsteplai;
```




<hr>



### variable type 

```Fortran
integer, intent(in) type;
```




<hr>



### variable used 

```Fortran
real(rk), intent(in) used;
```




<hr>



### variable user 

```Fortran
real(rk), intent(in) user;
```




<hr>



### variable ustar 

```Fortran
real(rk), intent(in) ustar;
```




<hr>



### variable value 

```Fortran
real(rk), intent(in) value;
```




<hr>



### variable vegetation 

```Fortran
integer, intent(in) vegetation;
```




<hr>



### variable velocity 

```Fortran
real(rk), dimension (m/s), intent(in) velocity;
```




<hr>



### variable vert 

```Fortran
integer, intent(in) vert;
```




<hr>



### variable vertical 

```Fortran
integer, intent(in) vertical;
```




<hr>



### variable voc 

```Fortran
real(rk), dimension (s), intent(in) voc;
```




<hr>



### variable volume 

```Fortran
real(rk), intent(out) volume;
```




<hr>



### variable volumetric 

```Fortran
real(rk), intent(in) volumetric;
```




<hr>



### variable vpgwt 

```Fortran
real(rk), dimension(size(zk)) vpgwt;
```




<hr>



### variable vtype 

```Fortran
integer, intent(in) vtype;
```




<hr>



### variable w126 

```Fortran
real(rk), intent(in) w126;
```




<hr>



### variable w126\_ref 

```Fortran
real(rk), intent(in) w126_ref;
```




<hr>



### variable w126\_set 

```Fortran
real(rk), intent(in) w126_set;
```




<hr>



### variable when 

```Fortran
integer, intent(in) when;
```




<hr>



### variable wilt 

```Fortran
real(rk), intent(in) wilt;
```




<hr>



### variable wilting 

```Fortran
real(rk), intent(in) wilting;
```




<hr>



### variable wind 

```Fortran
real(rk), intent(in) wind;
```




<hr>



### variable with 

```Fortran
real(rk), intent(in) with;
```




<hr>



### variable zk 

```Fortran
real(rk), dimension(:), intent(in) zk;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_bioemi_mod.F90`

