

# Namespace canopy\_utils\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_utils\_mod**](namespacecanopy__utils__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calccair**](#function-calccair) (pmbi pmbi, tki tki) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calcdx**](#function-calcdx) (lat lat, dlon dlon) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calcflameh**](#function-calcflameh) (frp frp, dx dx, [**lu\_opt**](canopy__profile__mod_8F90.md#variable-lu_opt) lu\_opt, [**vtype**](canopy__profile__mod_8F90.md#variable-vtype) vtype, flameh\_cal flameh\_cal) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calcpai**](#function-calcpai) (ch ch, canfrac canfrac) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calcpressure**](#function-calcpressure) ([**zk**](canopy__var3din__mod_8F90.md#variable-zk) zk, zref zref, pmbzref pmbzref, tref tref, t0 t0) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calcrelhum**](#function-calcrelhum) (tki tki, pmbi pmbi, qhi qhi) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calcrib**](#function-calcrib) (tak tak, tsk tsk, ubari ubari, [**d**](canopy__calcs_8F90.md#variable-d) d, zref zref) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calcspechum**](#function-calcspechum) (rhi rhi, tki tki, pmbi pmbi) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**calctemp**](#function-calctemp) ([**zk**](canopy__var3din__mod_8F90.md#variable-zk) zk, zref zref, taref taref, tsref tsref) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**convert\_qh\_to\_h2o**](#function-convert_qh_to_h2o) (qhi qhi, cairi cairi) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**effhenryslawcoeff**](#function-effhenryslawcoeff) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, ispec ispec) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**esat**](#function-esat) (tki tki) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**get\_canloss\_bio**](#function-get_canloss_bio) ([**loss\_opt**](canopy__bioemi__mod_8F90.md#variable-loss_opt) loss\_opt, [**lifetime**](canopy__bioemi__mod_8F90.md#variable-lifetime) lifetime, [**ustar**](canopy__drydep__mod_8F90.md#variable-ustar) ustar, ch ch) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**get\_gamma\_aq**](#function-get_gamma_aq) (aq\_opt aq\_opt, w126\_ozone w126\_ozone, [**w126\_set**](canopy__bioemi__mod_8F90.md#variable-w126_set) w126\_set, [**caq**](canopy__bioemi__mod_8F90.md#variable-caq) caq, [**taq**](canopy__bioemi__mod_8F90.md#variable-taq) taq, [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq) dtaq) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**get\_gamma\_co2**](#function-get_gamma_co2) (co2\_opt co2\_opt, co2\_set co2\_set) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**get\_gamma\_ht**](#function-get_gamma_ht) (ht\_opt ht\_opt, maxt2m maxt2m, [**cht**](canopy__bioemi__mod_8F90.md#variable-cht) cht, [**tht**](canopy__bioemi__mod_8F90.md#variable-tht) tht, [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht) dtht) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**get\_gamma\_hw**](#function-get_gamma_hw) (hw\_opt hw\_opt, maxws10m maxws10m, [**chw**](canopy__bioemi__mod_8F90.md#variable-chw) chw, [**thw**](canopy__bioemi__mod_8F90.md#variable-thw) thw, [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw) dthw) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**get\_gamma\_leafage**](#function-get_gamma_leafage) ([**leafage\_opt**](canopy__bioemi__mod_8F90.md#variable-leafage_opt) leafage\_opt, LAIpast LAIpast, LAIcurrent LAIcurrent, [**tsteplai**](canopy__calcs_8F90.md#variable-tsteplai) tsteplai, TABOVE TABOVE, Anew Anew, Agro Agro, Amat Amat, Aold Aold) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**get\_gamma\_lt**](#function-get_gamma_lt) (lt\_opt lt\_opt, mint2m mint2m, [**clt**](canopy__bioemi__mod_8F90.md#variable-clt) clt, [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt) tlt, [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt) dtlt) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**get\_gamma\_soim**](#function-get_gamma_soim) (soim\_opt soim\_opt, [**soim1**](canopy__bioemi__mod_8F90.md#variable-soim1) soim1, [**soim2**](canopy__bioemi__mod_8F90.md#variable-soim2) soim2, [**soim3**](canopy__bioemi__mod_8F90.md#variable-soim3) soim3, [**soim4**](canopy__bioemi__mod_8F90.md#variable-soim4) soim4, [**soid1**](canopy__bioemi__mod_8F90.md#variable-soid1) soid1, [**soid2**](canopy__bioemi__mod_8F90.md#variable-soid2) soid2, [**soid3**](canopy__bioemi__mod_8F90.md#variable-soid3) soid3, [**soid4**](canopy__bioemi__mod_8F90.md#variable-soid4) soid4, [**wilt**](canopy__bioemi__mod_8F90.md#variable-wilt) wilt, [**roota**](canopy__bioemi__mod_8F90.md#variable-roota) roota, [**rootb**](canopy__bioemi__mod_8F90.md#variable-rootb) rootb) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**integratetrapezoid**](#function-integratetrapezoid) (x x, y y) <br>_Numerical integration using trapezoidal rule._  |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**interp\_linear1\_internal**](#function-interp_linear1_internal) (x x, y y, xout xout) <br>_Linear interpolation function._  |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**lebasmvgas**](#function-lebasmvgas) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, ispec ispec) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**molarmassgas**](#function-molarmassgas) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, ispec ispec) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**molecdiff**](#function-molecdiff) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, ispec ispec, tkx tkx, pmbx pmbx) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**rav**](#function-rav) ([**ubar**](canopy__drydep__mod_8F90.md#variable-ubar) ubar, zref zref, [**d**](canopy__calcs_8F90.md#variable-d) d, hc hc, [**rib**](canopy__calcs_8F90.md#variable-rib) rib) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**rbl**](#function-rbl) ([**mdiffl**](canopy__drydep__mod_8F90.md#variable-mdiffl) mdiffl, ubari ubari) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**rcl**](#function-rcl) ([**hstarl**](canopy__drydep__mod_8F90.md#variable-hstarl) hstarl, [**f01**](canopy__drydep__mod_8F90.md#variable-f01) f01) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**reactivityparam**](#function-reactivityparam) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, ispec ispec) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**reactivityparamhno3**](#function-reactivityparamhno3) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, ispec ispec) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**rml**](#function-rml) ([**hstarl**](canopy__drydep__mod_8F90.md#variable-hstarl) hstarl, [**f01**](canopy__drydep__mod_8F90.md#variable-f01) f01) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**rs\_zhang\_gas**](#function-rs_zhang_gas) ([**mdiffl**](canopy__drydep__mod_8F90.md#variable-mdiffl) mdiffl, tki tki, pmbi pmbi, ppfdi ppfdi, [**srad**](canopy__drydep__mod_8F90.md#variable-srad) srad, relhumi relhumi) <br> |
|  subroutine, public | [**seteffhenryslawcoeffs**](#function-seteffhenryslawcoeffs) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, hstar hstar) <br> |
|  subroutine, public | [**setlebasmvgas**](#function-setlebasmvgas) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, mvg mvg) <br> |
|  subroutine, public | [**setmolarmassgas**](#function-setmolarmassgas) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, mmg mmg) <br> |
|  subroutine, public | [**setmolecdiffstp**](#function-setmolecdiffstp) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, mdiffstp mdiffstp) <br> |
|  subroutine, public | [**setreactivityhno3**](#function-setreactivityhno3) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, ar ar) <br> |
|  subroutine, public | [**setreactivityparams**](#function-setreactivityparams) ([**chemmechgas\_opt**](canopy__drydep__mod_8F90.md#variable-chemmechgas_opt) chemmechgas\_opt, [**chemmechgas\_tot**](canopy__drydep__mod_8F90.md#variable-chemmechgas_tot) chemmechgas\_tot, f0 f0) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**soilrbg**](#function-soilrbg) ([**ubar**](canopy__drydep__mod_8F90.md#variable-ubar) ubar) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**soilresist**](#function-soilresist) ([**mdiffl**](canopy__drydep__mod_8F90.md#variable-mdiffl) mdiffl, [**socat**](canopy__drydep__mod_8F90.md#variable-socat) socat, [**sotyp**](canopy__drydep__mod_8F90.md#variable-sotyp) sotyp, [**dsoil**](canopy__drydep__mod_8F90.md#variable-dsoil) dsoil, [**stheta**](canopy__drydep__mod_8F90.md#variable-stheta) stheta) <br> |
|  real(rk) [**function**](canopy__var3din__mod_8F90.md#variable-function), public | [**waterrbw**](#function-waterrbw) ([**ustar**](canopy__drydep__mod_8F90.md#variable-ustar) ustar) <br> |




























## Public Functions Documentation




### function calccair 

```Fortran
real(rk) function , public canopy_utils_mod::calccair (
    pmbi pmbi,
    tki tki
) 
```




<hr>



### function calcdx 

```Fortran
real(rk) function , public canopy_utils_mod::calcdx (
    lat lat,
    dlon dlon
) 
```




<hr>



### function calcflameh 

```Fortran
real(rk) function , public canopy_utils_mod::calcflameh (
    frp frp,
    dx dx,
    lu_opt lu_opt,
    vtype vtype,
    flameh_cal flameh_cal
) 
```




<hr>



### function calcpai 

```Fortran
real(rk) function , public canopy_utils_mod::calcpai (
    ch ch,
    canfrac canfrac
) 
```





**Parameters:**


* `ch` Input Grid cell canopy height (m)
* `canfrac` Input Grid cell canopy fraction 




        

<hr>



### function calcpressure 

```Fortran
real(rk) function , public canopy_utils_mod::calcpressure (
    zk zk,
    zref zref,
    pmbzref pmbzref,
    tref tref,
    t0 t0
) 
```




<hr>



### function calcrelhum 

```Fortran
real(rk) function , public canopy_utils_mod::calcrelhum (
    tki tki,
    pmbi pmbi,
    qhi qhi
) 
```




<hr>



### function calcrib 

```Fortran
real(rk) function , public canopy_utils_mod::calcrib (
    tak tak,
    tsk tsk,
    ubari ubari,
    d d,
    zref zref
) 
```




<hr>



### function calcspechum 

```Fortran
real(rk) function , public canopy_utils_mod::calcspechum (
    rhi rhi,
    tki tki,
    pmbi pmbi
) 
```




<hr>



### function calctemp 

```Fortran
real(rk) function , public canopy_utils_mod::calctemp (
    zk zk,
    zref zref,
    taref taref,
    tsref tsref
) 
```




<hr>



### function convert\_qh\_to\_h2o 

```Fortran
real(rk) function , public canopy_utils_mod::convert_qh_to_h2o (
    qhi qhi,
    cairi cairi
) 
```




<hr>



### function effhenryslawcoeff 

```Fortran
real(rk) function , public canopy_utils_mod::effhenryslawcoeff (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    ispec ispec
) 
```




<hr>



### function esat 

```Fortran
real(rk) function , public canopy_utils_mod::esat (
    tki tki
) 
```




<hr>



### function get\_canloss\_bio 

```Fortran
real(rk) function , public canopy_utils_mod::get_canloss_bio (
    loss_opt loss_opt,
    lifetime lifetime,
    ustar ustar,
    ch ch
) 
```




<hr>



### function get\_gamma\_aq 

```Fortran
real(rk) function , public canopy_utils_mod::get_gamma_aq (
    aq_opt aq_opt,
    w126_ozone w126_ozone,
    w126_set w126_set,
    caq caq,
    taq taq,
    dtaq dtaq
) 
```




<hr>



### function get\_gamma\_co2 

```Fortran
real(rk) function , public canopy_utils_mod::get_gamma_co2 (
    co2_opt co2_opt,
    co2_set co2_set
) 
```




<hr>



### function get\_gamma\_ht 

```Fortran
real(rk) function , public canopy_utils_mod::get_gamma_ht (
    ht_opt ht_opt,
    maxt2m maxt2m,
    cht cht,
    tht tht,
    dtht dtht
) 
```




<hr>



### function get\_gamma\_hw 

```Fortran
real(rk) function , public canopy_utils_mod::get_gamma_hw (
    hw_opt hw_opt,
    maxws10m maxws10m,
    chw chw,
    thw thw,
    dthw dthw
) 
```




<hr>



### function get\_gamma\_leafage 

```Fortran
real(rk) function , public canopy_utils_mod::get_gamma_leafage (
    leafage_opt leafage_opt,
    LAIpast LAIpast,
    LAIcurrent LAIcurrent,
    tsteplai tsteplai,
    TABOVE TABOVE,
    Anew Anew,
    Agro Agro,
    Amat Amat,
    Aold Aold
) 
```




<hr>



### function get\_gamma\_lt 

```Fortran
real(rk) function , public canopy_utils_mod::get_gamma_lt (
    lt_opt lt_opt,
    mint2m mint2m,
    clt clt,
    tlt tlt,
    dtlt dtlt
) 
```




<hr>



### function get\_gamma\_soim 

```Fortran
real(rk) function , public canopy_utils_mod::get_gamma_soim (
    soim_opt soim_opt,
    soim1 soim1,
    soim2 soim2,
    soim3 soim3,
    soim4 soim4,
    soid1 soid1,
    soid2 soid2,
    soid3 soid3,
    soid4 soid4,
    wilt wilt,
    roota roota,
    rootb rootb
) 
```




<hr>



### function integratetrapezoid 

_Numerical integration using trapezoidal rule._ 
```Fortran
real(rk) function , public canopy_utils_mod::integratetrapezoid (
    x x,
    y y
) 
```



Calculates the integral of an array y with respect to x using the trapezoid approximation. Note that the mesh spacing of x does not have to be uniform. 

**Parameters:**


* `x` Variable x array 
* `y` Function y(x) array 



**Returns:**

IntegrateTrapezoid Integral ∫y(x)·dx 




**Parameters:**


* `x` Variable x 
* `y` Function y(x) 



**Returns:**

Integral ∫y(x)·dx 





        

<hr>



### function interp\_linear1\_internal 

_Linear interpolation function._ 
```Fortran
real(rk) function , public canopy_utils_mod::interp_linear1_internal (
    x x,
    y y,
    xout xout
) 
```



Interpolates for the y value at the desired x value, given x and y values around the desired point. 

**Parameters:**


* `x` Two x-values surrounding the interpolation point 
* `y` Two y-values corresponding to the x-values 
* `xout` Desired x-value for interpolation 



**Returns:**

yout Interpolated y-value at xout 




**Parameters:**


* `xout` Input arrays and target value 



**Returns:**

Interpolated result 





        

<hr>



### function lebasmvgas 

```Fortran
real(rk) function , public canopy_utils_mod::lebasmvgas (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    ispec ispec
) 
```




<hr>



### function molarmassgas 

```Fortran
real(rk) function , public canopy_utils_mod::molarmassgas (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    ispec ispec
) 
```




<hr>



### function molecdiff 

```Fortran
real(rk) function , public canopy_utils_mod::molecdiff (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    ispec ispec,
    tkx tkx,
    pmbx pmbx
) 
```




<hr>



### function rav 

```Fortran
real(rk) function , public canopy_utils_mod::rav (
    ubar ubar,
    zref zref,
    d d,
    hc hc,
    rib rib
) 
```




<hr>



### function rbl 

```Fortran
real(rk) function , public canopy_utils_mod::rbl (
    mdiffl mdiffl,
    ubari ubari
) 
```




<hr>



### function rcl 

```Fortran
real(rk) function , public canopy_utils_mod::rcl (
    hstarl hstarl,
    f01 f01
) 
```




<hr>



### function reactivityparam 

```Fortran
real(rk) function , public canopy_utils_mod::reactivityparam (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    ispec ispec
) 
```




<hr>



### function reactivityparamhno3 

```Fortran
real(rk) function , public canopy_utils_mod::reactivityparamhno3 (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    ispec ispec
) 
```




<hr>



### function rml 

```Fortran
real(rk) function , public canopy_utils_mod::rml (
    hstarl hstarl,
    f01 f01
) 
```




<hr>



### function rs\_zhang\_gas 

```Fortran
real(rk) function , public canopy_utils_mod::rs_zhang_gas (
    mdiffl mdiffl,
    tki tki,
    pmbi pmbi,
    ppfdi ppfdi,
    srad srad,
    relhumi relhumi
) 
```




<hr>



### function seteffhenryslawcoeffs 

```Fortran
subroutine, public canopy_utils_mod::seteffhenryslawcoeffs (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    hstar hstar
) 
```




<hr>



### function setlebasmvgas 

```Fortran
subroutine, public canopy_utils_mod::setlebasmvgas (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    mvg mvg
) 
```




<hr>



### function setmolarmassgas 

```Fortran
subroutine, public canopy_utils_mod::setmolarmassgas (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    mmg mmg
) 
```




<hr>



### function setmolecdiffstp 

```Fortran
subroutine, public canopy_utils_mod::setmolecdiffstp (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    mdiffstp mdiffstp
) 
```




<hr>



### function setreactivityhno3 

```Fortran
subroutine, public canopy_utils_mod::setreactivityhno3 (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    ar ar
) 
```




<hr>



### function setreactivityparams 

```Fortran
subroutine, public canopy_utils_mod::setreactivityparams (
    chemmechgas_opt chemmechgas_opt,
    chemmechgas_tot chemmechgas_tot,
    f0 f0
) 
```




<hr>



### function soilrbg 

```Fortran
real(rk) function , public canopy_utils_mod::soilrbg (
    ubar ubar
) 
```




<hr>



### function soilresist 

```Fortran
real(rk) function , public canopy_utils_mod::soilresist (
    mdiffl mdiffl,
    socat socat,
    sotyp sotyp,
    dsoil dsoil,
    stheta stheta
) 
```




<hr>



### function waterrbw 

```Fortran
real(rk) function , public canopy_utils_mod::waterrbw (
    ustar ustar
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_utils_mod.F90`

