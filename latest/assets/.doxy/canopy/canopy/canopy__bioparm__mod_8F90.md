

# File canopy\_bioparm\_mod.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_bioparm\_mod.F90**](canopy__bioparm__mod_8F90.md)

[Go to the source code of this file](canopy__bioparm__mod_8F90_source.md)

_Biogenic Parameters Module._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**canopy\_bioparm\_mod**](namespacecanopy__bioparm__mod.md) <br> |








## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), dimension(default=0/viirs), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**al**](#variable-al)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**biogenic**](#variable-biogenic)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cell**](#variable-cell)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**dominant**](#variable-dominant)  <br> |
|  real(rk), dimension((ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf) !&gt; light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**beta**](canopy__bioemi__mod_8F90.md#variable-beta) !&gt; empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1) !&gt; [**out**](canopy__bioparm__mod_8F90.md#variable-out) activation energy(kj/mol) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo) !&gt; [**out**](canopy__bioparm__mod_8F90.md#variable-out) empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**anew**](canopy__bioemi__mod_8F90.md#variable-anew), [**agro**](canopy__bioemi__mod_8F90.md#variable-agro), [**amat**](canopy__bioemi__mod_8F90.md#variable-amat), [**aold**](canopy__bioemi__mod_8F90.md#variable-aold) !&gt; empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**roota**](canopy__bioemi__mod_8F90.md#variable-roota), [**rootb**](canopy__bioemi__mod_8F90.md#variable-rootb) !&gt; coefficients a [**and**](canopy__calcs_8F90.md#variable-and) b [**used**](canopy__bioemi__mod_8F90.md#variable-used) [**for**](canopy__calcs_8F90.md#variable-for) pft dependent cumulative root [**depth**](canopy__bioemi__mod_8F90.md#variable-depth) [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction)[[**m**](canopy__bioemi__mod_8F90.md#variable-m)-1] real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**caq**](canopy__bioemi__mod_8F90.md#variable-caq) !&gt; [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**taq**](canopy__bioemi__mod_8F90.md#variable-taq) !&gt; threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq) !&gt; delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**cht**](canopy__bioemi__mod_8F90.md#variable-cht) !&gt; [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**tht**](canopy__bioemi__mod_8F90.md#variable-tht) !&gt; threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht) !&gt; delta threshold [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**clt**](canopy__bioemi__mod_8F90.md#variable-clt) !&gt; [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt) !&gt; threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt) !&gt; delta threshold [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**chw**](canopy__bioemi__mod_8F90.md#variable-chw) !&gt; [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**thw**](canopy__bioemi__mod_8F90.md#variable-thw) !&gt; threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) ::[**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw) !&gt; delta threshold [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) !&gt; \} !&gt; \defgroup bioparm\_local\_vars local [**variables**](canopy__calcs_8F90.md#variable-variables) !! \brief local [**variables**](canopy__calcs_8F90.md#variable-variables) [**for**](canopy__calcs_8F90.md#variable-for) parameter assignment !! \{ real(rk) ::ef1, ef2, ef3, ef4, ef5, ef6, ef7 !&gt; [**plant**](canopy__var3din__mod_8F90.md#variable-plant) emission factors(ef)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk) ::ef8, ef9, ef10, ef11, ef12, ef13 !&gt; [**plant**](canopy__var3din__mod_8F90.md#variable-plant) emission factors(ef)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk) ::ef14, ef15 !&gt; [**plant**](canopy__var3din__mod_8F90.md#variable-plant) emission factors(ef)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) !&gt; \} !&gt; \defgroup bioparm\_isop\_params isoprene parameters !! \brief [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent emission capacity factors [**for**](canopy__calcs_8F90.md#variable-for) isoprene [**from**](canopy__var3din__mod_8F90.md#variable-from) tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al.(2012) !! \{ !&gt; \brief needleleaf evergreen temperate tree isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_isop= 600.0\_rk !&gt; \brief needleleaf evergreen boreal tree isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef2\_isop= 3000.0\_rk !&gt; \brief needleleaf deciduous boreal tree isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef3\_isop= 1.0\_rk !&gt; \brief broadleaf evergreen tropical tree isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef4\_isop= 7000.0\_rk !&gt; \brief broadleaf evergreen temperate tree isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef5\_isop= 10000.0\_rk !&gt; \brief broadleaf deciduous tropical tree isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef6\_isop= 7000.0\_rk !&gt; \brief broadleaf deciduous temperate tree isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef7\_isop= 10000.0\_rk !&gt; \brief broadleaf deciduous boreal tree isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef8\_isop= 11000.0\_rk !&gt; \brief broadleaf evergreen temperate shrub isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef9\_isop= 2000.0\_rk !&gt; \brief broadleaf deciduous temperate shrub isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef10\_isop= 4000.0\_rk !&gt; \brief broadleaf deciduous boreal shrub isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef11\_isop= 4000.0\_rk !&gt; \brief arctic c3 grass isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef12\_isop= 1600.0\_rk !&gt; \brief cool c3 grass isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef13\_isop= 800.0\_rk !&gt; \brief warm c4 grass isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef14\_isop= 200.0\_rk !&gt; \brief crop1 isoprene ef(μg/m²/[**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef15\_isop= 1.0\_rk !&gt; \brief isoprene [**leaf**](canopy__phot__mod_8F90.md#variable-leaf) age [**factor**](canopy__phot__mod_8F90.md#variable-factor) [**for**](canopy__calcs_8F90.md#variable-for) new [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage)(table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012) real(rk), parameter ::anew\_isop=0.05\_rk !&gt; \brief isoprene [**leaf**](canopy__phot__mod_8F90.md#variable-leaf) age [**factor**](canopy__phot__mod_8F90.md#variable-factor) [**for**](canopy__calcs_8F90.md#variable-for) growing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage)(table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012) real(rk), parameter ::agro\_isop=0.6\_rk !&gt; \brief isoprene [**leaf**](canopy__phot__mod_8F90.md#variable-leaf) age [**factor**](canopy__phot__mod_8F90.md#variable-factor) [**for**](canopy__calcs_8F90.md#variable-for) mature [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage)(table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012) real(rk), parameter ::amat\_isop=1.0\_rk !&gt; \brief isoprene [**leaf**](canopy__phot__mod_8F90.md#variable-leaf) age [**factor**](canopy__phot__mod_8F90.md#variable-factor) [**for**](canopy__calcs_8F90.md#variable-for) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage)(table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012) real(rk), parameter ::aold\_isop=0.9\_rk !&gt; \} ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) myrcene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_myrc= 70.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_myrc= 70.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_myrc= 60.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_myrc= 80.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_myrc= 30.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_myrc= 80.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_myrc= 30.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_myrc= 30.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_myrc= 30.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_myrc= 50.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_myrc= 30.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_myrc= 0.3\_rk ! arctic c3 grass real(rk), parameter ::ef13\_myrc= 0.3\_rk ! cool c3 grass real(rk), parameter ::ef14\_myrc= 0.3\_rk ! warm c4 grass real(rk), parameter ::ef15\_myrc= 0.3\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) myrcene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_myrc=2.0\_rk real(rk), parameter ::agro\_myrc=1.8\_rk real(rk), parameter ::amat\_myrc=1.0\_rk real(rk), parameter ::aold\_myrc=1.05\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) sabinene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_sabi= 70.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_sabi= 70.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_sabi= 40.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_sabi= 80.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_sabi= 50.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_sabi= 80.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_sabi= 50.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_sabi= 50.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_sabi= 50.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_sabi= 70.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_sabi= 50.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_sabi= 0.7\_rk ! arctic c3 grass real(rk), parameter ::ef13\_sabi= 0.7\_rk ! cool c3 grass real(rk), parameter ::ef14\_sabi= 0.7\_rk ! warm c4 grass real(rk), parameter ::ef15\_sabi= 0.7\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) sabinene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_sabi=2.0\_rk real(rk), parameter ::agro\_sabi=1.8\_rk real(rk), parameter ::amat\_sabi=1.0\_rk real(rk), parameter ::aold\_sabi=1.05\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) limonene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_limo= 100.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_limo= 100.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_limo= 130.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_limo= 80.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_limo= 80.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_limo= 80.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_limo= 80.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_limo= 80.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_limo= 60.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_limo= 100.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_limo= 60.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_limo= 0.7\_rk ! arctic c3 grass real(rk), parameter ::ef13\_limo= 0.7\_rk ! cool c3 grass real(rk), parameter ::ef14\_limo= 0.7\_rk ! warm c4 grass real(rk), parameter ::ef15\_limo= 0.7\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) limonene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_limo=2.0\_rk real(rk), parameter ::agro\_limo=1.8\_rk real(rk), parameter ::amat\_limo=1.0\_rk real(rk), parameter ::aold\_limo=1.05\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) 3-carene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_care= 160.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_care= 160.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_care= 80.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_care= 40.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_care= 30.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_care= 40.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_care= 30.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_care= 30.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_care= 30.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_care= 100.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_care= 30.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_care= 0.3\_rk ! arctic c3 grass real(rk), parameter ::ef13\_care= 0.3\_rk ! cool c3 grass real(rk), parameter ::ef14\_care= 0.3\_rk ! warm c4 grass real(rk), parameter ::ef15\_care= 0.3\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) 3-carene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_care=2.0\_rk real(rk), parameter ::agro\_care=1.8\_rk real(rk), parameter ::amat\_care=1.0\_rk real(rk), parameter ::aold\_care=1.05\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) t-[**beta**](canopy__bioemi__mod_8F90.md#variable-beta)-ocimene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_ocim= 70.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_ocim= 70.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_ocim= 60.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_ocim= 150.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_ocim= 120.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_ocim= 150.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_ocim= 120.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_ocim= 120.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_ocim= 90.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_ocim= 150.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_ocim= 90.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_ocim= 2.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_ocim= 2.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_ocim= 2.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_ocim= 2.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) t-[**beta**](canopy__bioemi__mod_8F90.md#variable-beta)-ocimene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_ocim=2.0\_rk real(rk), parameter ::agro\_ocim=1.8\_rk real(rk), parameter ::amat\_ocim=1.0\_rk real(rk), parameter ::aold\_ocim=1.05\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)-pinene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_bpin= 300.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_bpin= 300.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_bpin= 200.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_bpin= 120.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_bpin= 130.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_bpin= 120.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_bpin= 130.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_bpin= 130.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_bpin= 100.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_bpin= 150.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_bpin= 100.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_bpin= 1.5\_rk ! arctic c3 grass real(rk), parameter ::ef13\_bpin= 1.5\_rk ! cool c3 grass real(rk), parameter ::ef14\_bpin= 1.5\_rk ! warm c4 grass real(rk), parameter ::ef15\_bpin= 1.5\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)-pinene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_bpin=2.0\_rk real(rk), parameter ::agro\_bpin=1.8\_rk real(rk), parameter ::amat\_bpin=1.0\_rk real(rk), parameter ::aold\_bpin=1.05\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) alpha-pinene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_apin= 500.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_apin= 500.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_apin= 510.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_apin= 600.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_apin= 400.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_apin= 600.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_apin= 400.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_apin= 400.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_apin= 200.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_apin= 300.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_apin= 200.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_apin= 2.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_apin= 2.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_apin= 2.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_apin= 2.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) alpha-pinene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_apin=2.0\_rk real(rk), parameter ::agro\_apin=1.8\_rk real(rk), parameter ::amat\_apin=1.0\_rk real(rk), parameter ::aold\_apin=1.05\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) other monoterpenes(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) ! ! other monoterpenes category(34 compounds): see table 1 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al.(2012) real(rk), parameter ::ef1\_mono= 180.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_mono= 180.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_mono= 170.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_mono= 150.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_mono= 150.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_mono= 150.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_mono= 150.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_mono= 150.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_mono= 110.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_mono= 200.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_mono= 110.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_mono= 5.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_mono= 5.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_mono= 5.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_mono= 5.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) other monoterpenes as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_mono=2.0\_rk real(rk), parameter ::agro\_mono=1.8\_rk real(rk), parameter ::amat\_mono=1.0\_rk real(rk), parameter ::aold\_mono=1.05\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) alpha-farnesene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_farn= 40.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_farn= 40.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_farn= 40.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_farn= 60.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_farn= 40.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_farn= 60.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_farn= 40.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_farn= 40.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_farn= 40.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_farn= 40.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_farn= 40.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_farn= 3.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_farn= 3.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_farn= 3.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_farn= 4.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) alpha-farnesene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_farn=0.4\_rk real(rk), parameter ::agro\_farn=0.6\_rk real(rk), parameter ::amat\_farn=1.0\_rk real(rk), parameter ::aold\_farn=0.95\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)-caryophyllene(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_cary= 80.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_cary= 80.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_cary= 80.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_cary= 60.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_cary= 40.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_cary= 60.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_cary= 40.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_cary= 40.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_cary= 50.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_cary= 50.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_cary= 50.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_cary= 1.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_cary= 1.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_cary= 1.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_cary= 4.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)-caryophyllene as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_cary=0.4\_rk real(rk), parameter ::agro\_cary=0.6\_rk real(rk), parameter ::amat\_cary=1.0\_rk real(rk), parameter ::aold\_cary=0.95\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) other sesquieterpenes(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) ! other sesquiterpenes category(30 compounds): see table 1 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al.(2012) real(rk), parameter ::ef1\_sesq= 120.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_sesq= 120.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_sesq= 120.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_sesq= 120.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_sesq= 100.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_sesq= 120.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_sesq= 100.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_sesq= 100.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_sesq= 100.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_sesq= 100.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_sesq= 100.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_sesq= 1.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_sesq= 1.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_sesq= 1.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_sesq= 1.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) other sesquieterpenes as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_sesq=0.4\_rk real(rk), parameter ::agro\_sesq=0.6\_rk real(rk), parameter ::amat\_sesq=1.0\_rk real(rk), parameter ::aold\_sesq=0.95\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) 232-mbo(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_mbol= 700.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_mbol= 60.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_mbol= 0.01\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_mbol= 0.01\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_mbol= 0.01\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_mbol= 0.01\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_mbol= 0.01\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_mbol= 2.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_mbol= 0.01\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_mbol= 0.01\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_mbol= 0.01\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_mbol= 0.01\_rk ! arctic c3 grass real(rk), parameter ::ef13\_mbol= 0.01\_rk ! cool c3 grass real(rk), parameter ::ef14\_mbol= 0.01\_rk ! warm c4 grass real(rk), parameter ::ef15\_mbol= 0.01\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) 232-mbo as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_mbol=0.05\_rk real(rk), parameter ::agro\_mbol=0.6\_rk real(rk), parameter ::amat\_mbol=1.0\_rk real(rk), parameter ::aold\_mbol=0.9\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) methanol(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_meth= 900.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_meth= 900.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_meth= 900.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_meth= 500.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_meth= 900.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_meth= 500.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_meth= 900.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_meth= 900.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_meth= 900.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_meth= 900.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_meth= 900.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_meth= 500.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_meth= 500.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_meth= 500.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_meth= 900.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) methanol as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_meth=3.5\_rk real(rk), parameter ::agro\_meth=3.0\_rk real(rk), parameter ::amat\_meth=1.0\_rk real(rk), parameter ::aold\_meth=1.2\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) acetone(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_acet= 240.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_acet= 240.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_acet= 240.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_acet= 240.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_acet= 240.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_acet= 240.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_acet= 240.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_acet= 240.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_acet= 240.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_acet= 240.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_acet= 240.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_acet= 80.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_acet= 80.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_acet= 80.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_acet= 80.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) acetone as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_acet=1.0\_rk real(rk), parameter ::agro\_acet=1.0\_rk real(rk), parameter ::amat\_acet=1.0\_rk real(rk), parameter ::aold\_acet=1.0\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) carbon monoxide(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) real(rk), parameter ::ef1\_co= 600.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_co= 600.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_co= 600.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_co= 600.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_co= 600.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_co= 600.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_co= 600.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_co= 600.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_co= 600.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_co= 600.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_co= 600.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_co= 600.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_co= 600.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_co= 600.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_co= 600.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) carbon monoxide as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_co=1.0\_rk real(rk), parameter ::agro\_co=1.0\_rk real(rk), parameter ::amat\_co=1.0\_rk real(rk), parameter ::aold\_co=1.0\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) bidi [**voc**](canopy__bioemi__mod_8F90.md#variable-voc) [**species**](canopy__bioemi__mod_8F90.md#variable-species)(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) ! bidirectional [**voc**](canopy__bioemi__mod_8F90.md#variable-voc)(5 compounds):see table 1 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al.(2012) real(rk), parameter ::ef1\_bvoc= 500.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_bvoc= 500.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_bvoc= 500.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_bvoc= 500.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_bvoc= 500.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_bvoc= 500.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_bvoc= 500.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_bvoc= 500.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_bvoc= 500.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_bvoc= 500.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_bvoc= 500.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_bvoc= 80.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_bvoc= 80.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_bvoc= 80.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_bvoc= 80.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) bidi [**voc**](canopy__bioemi__mod_8F90.md#variable-voc) [**species**](canopy__bioemi__mod_8F90.md#variable-species) as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_bvoc=1.0\_rk real(rk), parameter ::agro\_bvoc=1.0\_rk real(rk), parameter ::amat\_bvoc=1.0\_rk real(rk), parameter ::aold\_bvoc=1.0\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) vocs(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) ! [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) [**voc**](canopy__bioemi__mod_8F90.md#variable-voc)(15 compounds):see table 1 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al.(2012) real(rk), parameter ::ef1\_svoc= 300.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_svoc= 300.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_svoc= 300.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_svoc= 300.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_svoc= 300.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_svoc= 300.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_svoc= 300.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_svoc= 300.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_svoc= 300.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_svoc= 300.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_svoc= 300.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_svoc= 300.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_svoc= 300.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_svoc= 300.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_svoc= 300.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) vocs as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_svoc=1.0\_rk real(rk), parameter ::agro\_svoc=1.0\_rk real(rk), parameter ::amat\_svoc=1.0\_rk real(rk), parameter ::aold\_svoc=1.0\_rk ! [**plant**](canopy__var3din__mod_8F90.md#variable-plant)-dependent [**emissions**](canopy__bioparm__mod_8F90.md#variable-emissions) capacity/factors(efs) [**for**](canopy__calcs_8F90.md#variable-for) other vocs(tables 2-3 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012)(ug/m2 [**hr**](canopy__bioemi__mod_8F90.md#variable-hr)) ! other [**voc**](canopy__bioemi__mod_8F90.md#variable-voc)(49 compounds):see table 1 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al.(2012) real(rk), parameter ::ef1\_ovoc= 140.0\_rk ! needleleaf evergreen temperate tree real(rk), parameter ::ef2\_ovoc= 140.0\_rk ! needleleaf evergreen boreal tree real(rk), parameter ::ef3\_ovoc= 140.0\_rk ! needleleaf deciduous boreal tree real(rk), parameter ::ef4\_ovoc= 140.0\_rk ! broadleaf evergreen tropical tree real(rk), parameter ::ef5\_ovoc= 140.0\_rk ! broadleaf evergreen temperate tree real(rk), parameter ::ef6\_ovoc= 140.0\_rk ! broadleaf deciduous tropical tree real(rk), parameter ::ef7\_ovoc= 140.0\_rk ! broadleaf deciduous temperate tree real(rk), parameter ::ef8\_ovoc= 140.0\_rk ! broadleaf deciduous boreal tree real(rk), parameter ::ef9\_ovoc= 140.0\_rk ! broadleaf evergreen temperate shrub real(rk), parameter ::ef10\_ovoc= 140.0\_rk ! broadleaf deciduous temperate shrub real(rk), parameter ::ef11\_ovoc= 140.0\_rk ! broadleaf deciduous boreal shrub real(rk), parameter ::ef12\_ovoc= 140.0\_rk ! arctic c3 grass real(rk), parameter ::ef13\_ovoc= 140.0\_rk ! cool c3 grass real(rk), parameter ::ef14\_ovoc= 140.0\_rk ! warm c4 grass real(rk), parameter ::ef15\_ovoc= 140.0\_rk ! crop1 !empirical factors [**or**](canopy__bioemi__mod_8F90.md#variable-or) coefficients for:growing, mature, [**and**](canopy__calcs_8F90.md#variable-and) old/senescing [**foliage**](canopy__var3din__mod_8F90.md#variable-foliage), [**for**](canopy__calcs_8F90.md#variable-for) other vocs as per table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012 real(rk), parameter ::anew\_ovoc=1.0\_rk real(rk), parameter ::agro\_ovoc=1.0\_rk real(rk), parameter ::amat\_ovoc=1.0\_rk real(rk), parameter ::aold\_ovoc=1.0\_rk ! [**species**](canopy__bioemi__mod_8F90.md#variable-species)-dependent parameterized [**canopy**](canopy__var3din__mod_8F90.md#variable-canopy) [**model**](canopy__phot__mod_8F90.md#variable-model) parameters(table 4 [**of**](canopy__var3din__mod_8F90.md#variable-of) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012) real(rk), parameter ::ldf\_isop= 1.0\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_isop= 0.13\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_isop= 95.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_isop= 2.0\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_isop= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_isop= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_isop= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_isop= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_isop= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_isop= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_isop= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_isop= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_isop= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_isop= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_isop= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_isop= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_myrc= 0.6\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_myrc= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_myrc= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_myrc= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_myrc= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_myrc= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_myrc= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_myrc= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_myrc= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_myrc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_myrc= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_myrc= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_myrc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_myrc= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_myrc= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_myrc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_sabi= 0.6\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_sabi= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_sabi= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_sabi= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_sabi= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_sabi= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_sabi= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_sabi= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_sabi= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_sabi= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_sabi= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_sabi= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_sabi= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_sabi= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_sabi= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_sabi= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_limo= 0.2\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_limo= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_limo= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_limo= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_limo= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_limo= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_limo= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_limo= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_limo= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_limo= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_limo= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_limo= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_limo= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_limo= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_limo= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_limo= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_care= 0.2\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_care= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_care= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_care= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_care= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_care= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_care= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_care= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_care= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_care= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_care= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_care= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_care= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_care= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_care= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_care= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_ocim= 0.8\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_ocim= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_ocim= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_ocim= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_ocim= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_ocim= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_ocim= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_ocim= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_ocim= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_ocim= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_ocim= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_ocim= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_ocim= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_ocim= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_ocim= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_ocim= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_bpin= 0.2\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_bpin= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_bpin= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_bpin= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_bpin= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_bpin= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_bpin= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_bpin= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_bpin= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_bpin= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_bpin= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_bpin= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_bpin= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_bpin= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_bpin= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_bpin= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_apin= 0.6\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_apin= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_apin= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_apin= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_apin= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_apin= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_apin= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_apin= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_apin= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_apin= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_apin= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_apin= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_apin= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_apin= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_apin= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_apin= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_mono= 0.4\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_mono= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_mono= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_mono= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_mono= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_mono= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_mono= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_mono= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_mono= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_mono= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_mono= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_mono= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_mono= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_mono= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_mono= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_mono= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_farn= 0.5\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_farn= 0.17\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_farn= 130.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_farn= 2.37\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_farn= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_farn= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_farn= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_farn= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_farn= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_farn= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_farn= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_farn= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_farn= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_farn= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_farn= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_farn= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_cary= 0.5\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_cary= 0.17\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_cary= 130.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_cary= 2.37\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_cary= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_cary= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_cary= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_cary= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_cary= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_cary= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_cary= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_cary= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_cary= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_cary= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_cary= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_cary= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_sesq= 0.5\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_sesq= 0.17\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_sesq= 130.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_sesq= 2.37\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_sesq= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_sesq= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_sesq= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_sesq= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_sesq= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_sesq= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_sesq= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_sesq= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_sesq= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_sesq= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_sesq= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_sesq= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_mbol= 1.0\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_mbol= 0.13\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_mbol= 95.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_mbol= 2.0\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_mbol= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_mbol= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_mbol= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_mbol= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_mbol= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_mbol= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_mbol= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_mbol= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_mbol= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_mbol= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_mbol= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_mbol= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_meth= 0.8\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_meth= 0.08\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_meth= 60.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_meth= 1.6\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_meth= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_meth= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_meth= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_meth= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_meth= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_meth= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_meth= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_meth= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_meth= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_meth= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_meth= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_meth= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_acet= 0.2\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_acet= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_acet= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_acet= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_acet= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_acet= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_acet= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_acet= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_acet= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_acet= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_acet= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_acet= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_acet= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_acet= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_acet= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_acet= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_co= 1.0\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_co= 0.08\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_co= 60.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_co= 1.6\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_co= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_co= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_co= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_co= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_co= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_co= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_co= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_co= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_co= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_co= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_co= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_co= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_bvoc= 0.8\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_bvoc= 0.13\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_bvoc= 95.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_bvoc= 2.0\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_bvoc= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_bvoc= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_bvoc= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_bvoc= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_bvoc= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_bvoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_bvoc= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_bvoc= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_bvoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_bvoc= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_bvoc= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_bvoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_svoc= 0.8\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_svoc= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_svoc= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_svoc= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_svoc= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_svoc= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_svoc= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_svoc= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_svoc= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_svoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_svoc= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_svoc= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_svoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_svoc= 5.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_svoc= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_svoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::ldf\_ovoc= 0.2\_rk !light-dependent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::beta\_ovoc= 0.1\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) dependence [**of**](canopy__var3din__mod_8F90.md#variable-of) light-independent [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction) real(rk), parameter ::ct1\_ovoc= 80.0\_rk !activation energy(kj/mol) real(rk), parameter ::ceo\_ovoc= 1.83\_rk !empirical [**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) real(rk), parameter ::caq\_ovoc= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::taq\_ovoc= 20.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::dtaq\_ovoc= 30.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) poor [**air**](canopy__bioemi__mod_8F90.md#variable-air) quality [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**ppm**](canopy__bioemi__mod_8F90.md#variable-ppm)-[**hours**](canopy__bioemi__mod_8F90.md#variable-hours)) real(rk), parameter ::cht\_ovoc= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tht\_ovoc= 313.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtht\_ovoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::clt\_ovoc= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::tlt\_ovoc= 283.15\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::dtlt\_ovoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**low**](canopy__bioemi__mod_8F90.md#variable-low) [**temperature**](canopy__bioemi__mod_8F90.md#variable-temperature) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**k**](canopy__calcs_8F90.md#variable-k)) real(rk), parameter ::chw\_ovoc= 1.0\_rk ![**coefficient**](canopy__bioemi__mod_8F90.md#variable-coefficient) [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress) real(rk), parameter ::thw\_ovoc= 12.0\_rk !threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) real(rk), parameter ::dthw\_ovoc= 8.0\_rk !delta threshold [**for**](canopy__calcs_8F90.md#variable-for) [**high**](canopy__bioemi__mod_8F90.md#variable-high) [**wind**](canopy__bioemi__mod_8F90.md#variable-wind) [**stress**](canopy__bioemi__mod_8F90.md#variable-stress)([**m**](canopy__bioemi__mod_8F90.md#variable-m)/s) ! [**set**](canopy__bioemi__mod_8F90.md#variable-set) tree [**and**](canopy__calcs_8F90.md#variable-and) [**species**](canopy__bioemi__mod_8F90.md#variable-species) dependent coefficients if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 1) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_isop [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_isop [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_isop [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_isop ef1=ef1\_isop ef2=ef2\_isop ef3=ef3\_isop ef4=ef4\_isop ef5=ef5\_isop ef6=ef6\_isop ef7=ef7\_isop ef8=ef8\_isop ef9=ef9\_isop ef10=ef10\_isop ef11=ef11\_isop ef12=ef12\_isop ef13=ef13\_isop ef14=ef14\_isop ef15=ef15\_isop [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_isop [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_isop [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_isop [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_isop [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_isop [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_isop [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_isop [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_isop [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_isop [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_isop [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_isop [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_isop [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_isop [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_isop [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_isop [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_isop else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 2) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_myrc [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_myrc [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_myrc [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_myrc ef1=ef1\_myrc ef2=ef2\_myrc ef3=ef3\_myrc ef4=ef4\_myrc ef5=ef5\_myrc ef6=ef6\_myrc ef7=ef7\_myrc ef8=ef8\_myrc ef9=ef9\_myrc ef10=ef10\_myrc ef11=ef11\_myrc ef12=ef12\_myrc ef13=ef13\_myrc ef14=ef14\_myrc ef15=ef15\_myrc [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_myrc [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_myrc [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_myrc [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_myrc [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_myrc [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_myrc [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_myrc [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_myrc [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_myrc [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_myrc [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_myrc [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_myrc [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_myrc [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_myrc [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_myrc [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_myrc else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 3) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_sabi [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_sabi [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_sabi [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_sabi ef1=ef1\_sabi ef2=ef2\_sabi ef3=ef3\_sabi ef4=ef4\_sabi ef5=ef5\_sabi ef6=ef6\_sabi ef7=ef7\_sabi ef8=ef8\_sabi ef9=ef9\_sabi ef10=ef10\_sabi ef11=ef11\_sabi ef12=ef12\_sabi ef13=ef13\_sabi ef14=ef14\_sabi ef15=ef15\_sabi [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_sabi [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_sabi [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_sabi [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_sabi [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_sabi [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_sabi [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_sabi [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_sabi [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_sabi [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_sabi [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_sabi [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_sabi [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_sabi [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_sabi [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_sabi [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_sabi else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 4) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_limo [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_limo [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_limo [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_limo ef1=ef1\_limo ef2=ef2\_limo ef3=ef3\_limo ef4=ef4\_limo ef5=ef5\_limo ef6=ef6\_limo ef7=ef7\_limo ef8=ef8\_limo ef9=ef9\_limo ef10=ef10\_limo ef11=ef11\_limo ef12=ef12\_limo ef13=ef13\_limo ef14=ef14\_limo ef15=ef15\_limo [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_limo [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_limo [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_limo [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_limo [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_limo [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_limo [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_limo [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_limo [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_limo [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_limo [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_limo [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_limo [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_limo [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_limo [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_limo [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_limo else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 5) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_care [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_care [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_care [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_care ef1=ef1\_care ef2=ef2\_care ef3=ef3\_care ef4=ef4\_care ef5=ef5\_care ef6=ef6\_care ef7=ef7\_care ef8=ef8\_care ef9=ef9\_care ef10=ef10\_care ef11=ef11\_care ef12=ef12\_care ef13=ef13\_care ef14=ef14\_care ef15=ef15\_care [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_care [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_care [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_care [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_care [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_care [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_care [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_care [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_care [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_care [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_care [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_care [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_care [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_care [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_care [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_care [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_care else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 6) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_ocim [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_ocim [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_ocim [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_ocim ef1=ef1\_ocim ef2=ef2\_ocim ef3=ef3\_ocim ef4=ef4\_ocim ef5=ef5\_ocim ef6=ef6\_ocim ef7=ef7\_ocim ef8=ef8\_ocim ef9=ef9\_ocim ef10=ef10\_ocim ef11=ef11\_ocim ef12=ef12\_ocim ef13=ef13\_ocim ef14=ef14\_ocim ef15=ef15\_ocim [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_ocim [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_ocim [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_ocim [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_ocim [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_ocim [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_ocim [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_ocim [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_ocim [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_ocim [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_ocim [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_ocim [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_ocim [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_ocim [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_ocim [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_ocim [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_ocim else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 7) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_bpin [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_bpin [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_bpin [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_bpin ef1=ef1\_bpin ef2=ef2\_bpin ef3=ef3\_bpin ef4=ef4\_bpin ef5=ef5\_bpin ef6=ef6\_bpin ef7=ef7\_bpin ef8=ef8\_bpin ef9=ef9\_bpin ef10=ef10\_bpin ef11=ef11\_bpin ef12=ef12\_bpin ef13=ef13\_bpin ef14=ef14\_bpin ef15=ef15\_bpin [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_bpin [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_bpin [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_bpin [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_bpin [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_bpin [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_bpin [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_bpin [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_bpin [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_bpin [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_bpin [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_bpin [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_bpin [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_bpin [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_bpin [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_bpin [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_bpin else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 8) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_apin [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_apin [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_apin [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_apin ef1=ef1\_apin ef2=ef2\_apin ef3=ef3\_apin ef4=ef4\_apin ef5=ef5\_apin ef6=ef6\_apin ef7=ef7\_apin ef8=ef8\_apin ef9=ef9\_apin ef10=ef10\_apin ef11=ef11\_apin ef12=ef12\_apin ef13=ef13\_apin ef14=ef14\_apin ef15=ef15\_apin [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_apin [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_apin [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_apin [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_apin [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_apin [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_apin [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_apin [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_apin [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_apin [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_apin [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_apin [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_apin [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_apin [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_apin [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_apin [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_apin else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 9) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_mono [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_mono [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_mono [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_mono ef1=ef1\_mono ef2=ef2\_mono ef3=ef3\_mono ef4=ef4\_mono ef5=ef5\_mono ef6=ef6\_mono ef7=ef7\_mono ef8=ef8\_mono ef9=ef9\_mono ef10=ef10\_mono ef11=ef11\_mono ef12=ef12\_mono ef13=ef13\_mono ef14=ef14\_mono ef15=ef15\_mono [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_mono [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_mono [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_mono [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_mono [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_mono [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_mono [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_mono [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_mono [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_mono [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_mono [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_mono [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_mono [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_mono [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_mono [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_mono [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_mono else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 10) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_farn [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_farn [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_farn [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_farn ef1=ef1\_farn ef2=ef2\_farn ef3=ef3\_farn ef4=ef4\_farn ef5=ef5\_farn ef6=ef6\_farn ef7=ef7\_farn ef8=ef8\_farn ef9=ef9\_farn ef10=ef10\_farn ef11=ef11\_farn ef12=ef12\_farn ef13=ef13\_farn ef14=ef14\_farn ef15=ef15\_farn [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_farn [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_farn [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_farn [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_farn [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_farn [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_farn [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_farn [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_farn [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_farn [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_farn [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_farn [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_farn [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_farn [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_farn [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_farn [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_farn else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 11) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_cary [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_cary [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_cary [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_cary ef1=ef1\_cary ef2=ef2\_cary ef3=ef3\_cary ef4=ef4\_cary ef5=ef5\_cary ef6=ef6\_cary ef7=ef7\_cary ef8=ef8\_cary ef9=ef9\_cary ef10=ef10\_cary ef11=ef11\_cary ef12=ef12\_cary ef13=ef13\_cary ef14=ef14\_cary ef15=ef15\_cary [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_cary [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_cary [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_cary [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_cary [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_cary [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_cary [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_cary [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_cary [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_cary [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_cary [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_cary [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_cary [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_cary [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_cary [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_cary [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_cary else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 12) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_sesq [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_sesq [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_sesq [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_sesq ef1=ef1\_sesq ef2=ef2\_sesq ef3=ef3\_sesq ef4=ef4\_sesq ef5=ef5\_sesq ef6=ef6\_sesq ef7=ef7\_sesq ef8=ef8\_sesq ef9=ef9\_sesq ef10=ef10\_sesq ef11=ef11\_sesq ef12=ef12\_sesq ef13=ef13\_sesq ef14=ef14\_sesq ef15=ef15\_sesq [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_sesq [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_sesq [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_sesq [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_sesq [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_sesq [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_sesq [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_sesq [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_sesq [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_sesq [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_sesq [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_sesq [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_sesq [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_sesq [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_sesq [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_sesq [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_sesq else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 13) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_mbol [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_mbol [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_mbol [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_mbol ef1=ef1\_mbol ef2=ef2\_mbol ef3=ef3\_mbol ef4=ef4\_mbol ef5=ef5\_mbol ef6=ef6\_mbol ef7=ef7\_mbol ef8=ef8\_mbol ef9=ef9\_mbol ef10=ef10\_mbol ef11=ef11\_mbol ef12=ef12\_mbol ef13=ef13\_mbol ef14=ef14\_mbol ef15=ef15\_mbol [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_mbol [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_mbol [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_mbol [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_mbol [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_mbol [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_mbol [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_mbol [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_mbol [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_mbol [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_mbol [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_mbol [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_mbol [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_mbol [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_mbol [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_mbol [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_mbol else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 14) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_meth [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_meth [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_meth [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_meth ef1=ef1\_meth ef2=ef2\_meth ef3=ef3\_meth ef4=ef4\_meth ef5=ef5\_meth ef6=ef6\_meth ef7=ef7\_meth ef8=ef8\_meth ef9=ef9\_meth ef10=ef10\_meth ef11=ef11\_meth ef12=ef12\_meth ef13=ef13\_meth ef14=ef14\_meth ef15=ef15\_meth [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_meth [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_meth [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_meth [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_meth [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_meth [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_meth [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_meth [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_meth [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_meth [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_meth [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_meth [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_meth [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_meth [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_meth [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_meth [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_meth else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 15) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_acet [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_acet [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_acet [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_acet ef1=ef1\_acet ef2=ef2\_acet ef3=ef3\_acet ef4=ef4\_acet ef5=ef5\_acet ef6=ef6\_acet ef7=ef7\_acet ef8=ef8\_acet ef9=ef9\_acet ef10=ef10\_acet ef11=ef11\_acet ef12=ef12\_acet ef13=ef13\_acet ef14=ef14\_acet ef15=ef15\_acet [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_acet [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_acet [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_acet [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_acet [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_acet [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_acet [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_acet [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_acet [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_acet [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_acet [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_acet [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_acet [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_acet [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_acet [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_acet [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_acet else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 16) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_co [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_co [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_co [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_co ef1=ef1\_co ef2=ef2\_co ef3=ef3\_co ef4=ef4\_co ef5=ef5\_co ef6=ef6\_co ef7=ef7\_co ef8=ef8\_co ef9=ef9\_co ef10=ef10\_co ef11=ef11\_co ef12=ef12\_co ef13=ef13\_co ef14=ef14\_co ef15=ef15\_co [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_co [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_co [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_co [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_co [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_co [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_co [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_co [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_co [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_co [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_co [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_co [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_co [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_co [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_co [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_co [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_co else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 17) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_bvoc [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_bvoc [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_bvoc [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_bvoc ef1=ef1\_bvoc ef2=ef2\_bvoc ef3=ef3\_bvoc ef4=ef4\_bvoc ef5=ef5\_bvoc ef6=ef6\_bvoc ef7=ef7\_bvoc ef8=ef8\_bvoc ef9=ef9\_bvoc ef10=ef10\_bvoc ef11=ef11\_bvoc ef12=ef12\_bvoc ef13=ef13\_bvoc ef14=ef14\_bvoc ef15=ef15\_bvoc [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_bvoc [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_bvoc [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_bvoc [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_bvoc [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_bvoc [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_bvoc [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_bvoc [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_bvoc [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_bvoc [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_bvoc [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_bvoc [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_bvoc [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_bvoc [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_bvoc [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_bvoc [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_bvoc else if([**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind) .eq. 18) then [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_svoc [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_svoc [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_svoc [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_svoc ef1=ef1\_svoc ef2=ef2\_svoc ef3=ef3\_svoc ef4=ef4\_svoc ef5=ef5\_svoc ef6=ef6\_svoc ef7=ef7\_svoc ef8=ef8\_svoc ef9=ef9\_svoc ef10=ef10\_svoc ef11=ef11\_svoc ef12=ef12\_svoc ef13=ef13\_svoc ef14=ef14\_svoc ef15=ef15\_svoc [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_svoc [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_svoc [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_svoc [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_svoc [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_svoc [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_svoc [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_svoc [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_svoc [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_svoc [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_svoc [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_svoc [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_svoc [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_svoc [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_svoc [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_svoc [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_svoc else ! [**emi\_ind**](canopy__bioparm__mod_8F90.md#variable-emi_ind)=19 [**ldf**](canopy__bioemi__mod_8F90.md#variable-ldf)=ldf\_ovoc [**beta**](canopy__bioemi__mod_8F90.md#variable-beta)=beta\_ovoc [**ct1**](canopy__bioemi__mod_8F90.md#variable-ct1)=ct1\_ovoc [**ceo**](canopy__bioemi__mod_8F90.md#variable-ceo)=ceo\_ovoc ef1=ef1\_ovoc ef2=ef2\_ovoc ef3=ef3\_ovoc ef4=ef4\_ovoc ef5=ef5\_ovoc ef6=ef6\_ovoc ef7=ef7\_ovoc ef8=ef8\_ovoc ef9=ef9\_ovoc ef10=ef10\_ovoc ef11=ef11\_ovoc ef12=ef12\_ovoc ef13=ef13\_ovoc ef14=ef14\_ovoc ef15=ef15\_ovoc [**anew**](canopy__bioemi__mod_8F90.md#variable-anew)=anew\_ovoc [**agro**](canopy__bioemi__mod_8F90.md#variable-agro)=agro\_ovoc [**amat**](canopy__bioemi__mod_8F90.md#variable-amat)=amat\_ovoc [**aold**](canopy__bioemi__mod_8F90.md#variable-aold)=aold\_ovoc [**caq**](canopy__bioemi__mod_8F90.md#variable-caq)=caq\_ovoc [**taq**](canopy__bioemi__mod_8F90.md#variable-taq)=taq\_ovoc [**dtaq**](canopy__bioemi__mod_8F90.md#variable-dtaq)=dtaq\_ovoc [**cht**](canopy__bioemi__mod_8F90.md#variable-cht)=cht\_ovoc [**tht**](canopy__bioemi__mod_8F90.md#variable-tht)=tht\_ovoc [**dtht**](canopy__bioemi__mod_8F90.md#variable-dtht)=dtht\_ovoc [**clt**](canopy__bioemi__mod_8F90.md#variable-clt)=clt\_ovoc [**tlt**](canopy__bioemi__mod_8F90.md#variable-tlt)=tlt\_ovoc [**dtlt**](canopy__bioemi__mod_8F90.md#variable-dtlt)=dtlt\_ovoc [**chw**](canopy__bioemi__mod_8F90.md#variable-chw)=chw\_ovoc [**thw**](canopy__bioemi__mod_8F90.md#variable-thw)=thw\_ovoc [**dthw**](canopy__bioemi__mod_8F90.md#variable-dthw)=dthw\_ovoc end if if([**lu\_opt**](canopy__profile__mod_8F90.md#variable-lu_opt) .eq. 0 .or. [**lu\_opt**](canopy__profile__mod_8F90.md#variable-lu_opt) .eq. 1) then !viirs [**or**](canopy__bioemi__mod_8F90.md#variable-or) modis [**lu**](canopy__bioparm__mod_8F90.md#variable-lu) types ! simple [**megan**](canopy__bioemi__mod_8F90.md#variable-megan)(table 3 [**in**](canopy__calcs_8F90.md#variable-in) guenther [**et**](canopy__bioparm__mod_8F90.md#variable-et) al., 2012) pft [**to**](canopy__bioparm__mod_8F90.md#variable-to) viirs/modis [**vtype**](canopy__profile__mod_8F90.md#variable-vtype) mapping if([**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .eq. 1) then !viirs cat 1 evergreen needleleaf !--&gt; [**average**](canopy__bioemi__mod_8F90.md#variable-average) needleleaf evergreen temperate tree [**and**](canopy__calcs_8F90.md#variable-and) needleleaf evergreen boreal tree ef=(ef1+ef2)/2.0\_rk ![**set**](canopy__bioemi__mod_8F90.md#variable-set) pft dependent a [**and**](canopy__calcs_8F90.md#variable-and) b coefficients [**for**](canopy__calcs_8F90.md#variable-for) cumulative root [**depth**](canopy__bioemi__mod_8F90.md#variable-depth) [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction)(zeng 2001) !see table 2 [**for**](canopy__calcs_8F90.md#variable-for) igpb classification at: https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541\_2001\_002\_0525\_gvrdfl\_2\_0\_co\_2.xml [**roota**](canopy__bioemi__mod_8F90.md#variable-roota)=6.706\_rk [**rootb**](canopy__bioemi__mod_8F90.md#variable-rootb)=2.175\_rk else if([**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .eq. 2) then !viirs/modis cat 2 evergreen broadleaf !--&gt; [**average**](canopy__bioemi__mod_8F90.md#variable-average) broadleaf evergreen tropical tree [**and**](canopy__calcs_8F90.md#variable-and) broadleaf evergreen temperate tree ef=(ef4+ef5)/2.0\_rk ![**set**](canopy__bioemi__mod_8F90.md#variable-set) pft dependent a [**and**](canopy__calcs_8F90.md#variable-and) b coefficients [**for**](canopy__calcs_8F90.md#variable-for) cumulative root [**depth**](canopy__bioemi__mod_8F90.md#variable-depth) [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction)(zeng 2001) !see table 2 [**for**](canopy__calcs_8F90.md#variable-for) igpb classification at: https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541\_2001\_002\_0525\_gvrdfl\_2\_0\_co\_2.xml [**roota**](canopy__bioemi__mod_8F90.md#variable-roota)=7.344\_rk [**rootb**](canopy__bioemi__mod_8F90.md#variable-rootb)=1.303\_rk else if([**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .eq. 3) then !viirs/modis cat 3 deciduous needleaf !--&gt; [**average**](canopy__bioemi__mod_8F90.md#variable-average) needleleaf deciduous boreal tree ef=ef3 ![**set**](canopy__bioemi__mod_8F90.md#variable-set) pft dependent a [**and**](canopy__calcs_8F90.md#variable-and) b coefficients [**for**](canopy__calcs_8F90.md#variable-for) cumulative root [**depth**](canopy__bioemi__mod_8F90.md#variable-depth) [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction)(zeng 2001) !see table 2 [**for**](canopy__calcs_8F90.md#variable-for) igpb classification at: https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541\_2001\_002\_0525\_gvrdfl\_2\_0\_co\_2.xml [**roota**](canopy__bioemi__mod_8F90.md#variable-roota)=7.066\_rk [**rootb**](canopy__bioemi__mod_8F90.md#variable-rootb)=1.953\_rk else if([**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .eq. 4) then !viirs/modis cat 4 deciduous broadleaf !--&gt; [**average**](canopy__bioemi__mod_8F90.md#variable-average) broadleaf deciduous tropical tree, broadleaf deciduous temperate tree, ! [**and**](canopy__calcs_8F90.md#variable-and) broadleaf deciduous boreal tree ef=(ef6+ef7+ef8)/3.0\_rk ![**set**](canopy__bioemi__mod_8F90.md#variable-set) pft dependent a [**and**](canopy__calcs_8F90.md#variable-and) b coefficients [**for**](canopy__calcs_8F90.md#variable-for) cumulative root [**depth**](canopy__bioemi__mod_8F90.md#variable-depth) [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction)(zeng 2001) !see table 2 [**for**](canopy__calcs_8F90.md#variable-for) igpb classification at: https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541\_2001\_002\_0525\_gvrdfl\_2\_0\_co\_2.xml [**roota**](canopy__bioemi__mod_8F90.md#variable-roota)=5.990\_rk [**rootb**](canopy__bioemi__mod_8F90.md#variable-rootb)=1.955\_rk else if([**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .eq. 5) then !viirs/modis cat 5 mixed forests !--&gt; avearge [**of**](canopy__var3din__mod_8F90.md#variable-of) [**all**](canopy__bioemi__mod_8F90.md#variable-all) [**above**](canopy__var3din__mod_8F90.md#variable-above) ef1-ef8 pfts. ef=(ef1+ef2+ef3+ef4+ef5+ef6+ef7+ef8)/8.0\_rk ![**set**](canopy__bioemi__mod_8F90.md#variable-set) pft dependent a [**and**](canopy__calcs_8F90.md#variable-and) b coefficients [**for**](canopy__calcs_8F90.md#variable-for) cumulative root [**depth**](canopy__bioemi__mod_8F90.md#variable-depth) [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction)(zeng 2001) !see table 2 [**for**](canopy__calcs_8F90.md#variable-for) igpb classification at: https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541\_2001\_002\_0525\_gvrdfl\_2\_0\_co\_2.xml [**roota**](canopy__bioemi__mod_8F90.md#variable-roota)=4.453\_rk [**rootb**](canopy__bioemi__mod_8F90.md#variable-rootb)=1.631\_rk else if([**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .ge. 6 .and. [**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .le. 7) then !viirs/modis cat 6-7 closed/open shrublands !--&gt; avearge broadleaf evergreen temperate shrub, broadleaf deciduous temperate shrub, ! [**and**](canopy__calcs_8F90.md#variable-and) broadleaf deciduous boreal shrub ef=(ef9+ef10+ef11)/3.0\_rk ![**set**](canopy__bioemi__mod_8F90.md#variable-set) pft dependent a [**and**](canopy__calcs_8F90.md#variable-and) b coefficients [**for**](canopy__calcs_8F90.md#variable-for) cumulative root [**depth**](canopy__bioemi__mod_8F90.md#variable-depth) [**fraction**](canopy__bioemi__mod_8F90.md#variable-fraction)(zeng 2001) !see table 2 [**for**](canopy__calcs_8F90.md#variable-for) igpb classification at: https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541\_2001\_002\_0525\_gvrdfl\_2\_0\_co\_2.xml [**roota**](canopy__bioemi__mod_8F90.md#variable-roota)=(6.326\_rk+7.718\_rk)/2.0\_rk [**rootb**](canopy__bioemi__mod_8F90.md#variable-rootb)=(1.567\_rk+1.262\_rk)/2.0\_rk else if([**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .ge. 8 .and. [**vtype**](canopy__profile__mod_8F90.md#variable-vtype) .le. 11) then !viirs/modis cat 8-10 savannas [**and**](canopy__calcs_8F90.md#variable-and) grasslands !--&gt; avearge arctic c3 grass, cool c3 grass, warm c4 grass), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**ef**](#variable-ef)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**emi\_ind**](#variable-emi_ind)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**emissions**](#variable-emissions)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**et**](#variable-et)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**for**](#variable-for)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**from**](#variable-from)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**grid**](#variable-grid)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**index**](#variable-index)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**input**](#variable-input)  <br> |
|  integer, intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**integer**](#variable-integer)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lu**](#variable-lu)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lu\_opt**](#variable-lu_opt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**mapped**](#variable-mapped)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**massman**](#variable-massman)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**model**](#variable-model)  <br> |
|  real(rk), intent(out) | [**out**](#variable-out)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**to**](#variable-to)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**type**](#variable-type)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vegetation**](#variable-vegetation)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vtype**](#variable-vtype)  <br> |












































## Detailed Description


This module contains the CANOPY\_BIOP subroutine which provides biogenic emission factors and parameters from MEGAN2.1 (Model of Emissions of Gases and Aerosols from Nature). The module contains extensive parameter tables for different biogenic volatile organic compounds (BVOCs) and vegetation types.




**Author:**

Patrick C. Campbell 




**Date:**

February 2023


\references Guenther, A. B., et al.: The Model of Emissions of Gases and Aerosols from Nature version 2.1 (MEGAN2.1): an extended and updated framework for modeling biogenic emissions, Geosci. Model Dev., 5, 1471–1492, [https://doi.org/10.5194/gmd-5-1471-2012](https://doi.org/10.5194/gmd-5-1471-2012), 2012. 


    
## Public Attributes Documentation




### variable al 

```Fortran
integer, dimension (default = 0/viirs), intent(in) al;
```




<hr>



### variable biogenic 

```Fortran
integer, intent(in) biogenic;
```




<hr>



### variable cell 

```Fortran
integer, intent(in) cell;
```




<hr>



### variable dominant 

```Fortran
integer, intent(in) dominant;
```




<hr>



### variable ef 

```Fortran
real(rk), dimension ((ug/m2 hr)
        real(rk),    intent( out )      :: ldf             !> light-dependent fraction
        real(rk),    intent( out )      :: beta            !> empirical coefficient for temperature dependence of light-independent fraction
        real(rk),    intent( out )      :: ct1             !> out activation energy (kj/mol)
        real(rk),    intent( out )      :: ceo             !> out empirical coefficient
        real(rk),    intent( out )      :: anew, agro, amat, aold   !> empirical factors or coefficients for: growing, mature, and old/senescing foliage, as per table 4 of guenther et al., 2012
        real(rk),    intent( out )      :: roota, rootb    !> coefficients a and b used for pft dependent cumulative root depth fraction [m-1]
        real(rk),    intent( out )      :: caq             !> coefficient for poor air quality stress
        real(rk),    intent( out )      :: taq             !> threshold for poor air quality stress (ppm-hours)
        real(rk),    intent( out )      :: dtaq            !> delta threshold for poor air quality stress (ppm-hours)
        real(rk),    intent( out )      :: cht             !> coefficient for high temperature stress
        real(rk),    intent( out )      :: tht             !> threshold for high temperature stress (k)
        real(rk),    intent( out )      :: dtht            !> delta threshold high temperature stress (k)
        real(rk),    intent( out )      :: clt             !> coefficient for low temperature stress
        real(rk),    intent( out )      :: tlt             !> threshold for low temperature stress (k)
        real(rk),    intent( out )      :: dtlt            !> delta threshold low temperature stress (k)
        real(rk),    intent( out )      :: chw             !> coefficient for high wind stress
        real(rk),    intent( out )      :: thw             !> threshold for high wind stress (m/s)
        real(rk),    intent( out )      :: dthw            !> delta threshold high wind stress (m/s)
!> \}

!> \defgroup bioparm_local_vars local variables
!! \brief local variables for parameter assignment
!! \{
        real(rk) :: ef1,ef2,ef3,ef4,ef5,ef6,ef7    !> plant emission factors (ef) (ug/m2 hr)
        real(rk) :: ef8,ef9,ef10,ef11,ef12,ef13    !> plant emission factors (ef) (ug/m2 hr)
        real(rk) :: ef14,ef15                      !> plant emission factors (ef) (ug/m2 hr)
!> \}

!> \defgroup bioparm_isop_params isoprene parameters
!! \brief plant-dependent emission capacity factors for isoprene from tables 2-3 of guenther et al. (2012)
!! \{

        !> \brief needleleaf evergreen temperate tree isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef1_isop    =  600.0_rk
        !> \brief needleleaf evergreen boreal tree isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef2_isop    =  3000.0_rk
        !> \brief needleleaf deciduous boreal tree isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef3_isop    =  1.0_rk
        !> \brief broadleaf evergreen tropical tree isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef4_isop    =  7000.0_rk
        !> \brief broadleaf evergreen temperate tree isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef5_isop    =  10000.0_rk
        !> \brief broadleaf deciduous tropical tree isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef6_isop    =  7000.0_rk
        !> \brief broadleaf deciduous temperate tree isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef7_isop    =  10000.0_rk
        !> \brief broadleaf deciduous boreal tree isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef8_isop    =  11000.0_rk
        !> \brief broadleaf evergreen temperate shrub isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef9_isop    =  2000.0_rk
        !> \brief broadleaf deciduous temperate shrub isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef10_isop   =  4000.0_rk
        !> \brief broadleaf deciduous boreal shrub isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef11_isop   =  4000.0_rk
        !> \brief arctic c3 grass isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef12_isop   =  1600.0_rk
        !> \brief cool c3 grass isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef13_isop   =  800.0_rk
        !> \brief warm c4 grass isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef14_isop   =  200.0_rk
        !> \brief crop1 isoprene ef (μg/m²/hr)
        real(rk),          parameter     :: ef15_isop   =  1.0_rk

        !> \brief isoprene leaf age factor for new foliage (table 4 of guenther et al., 2012)
        real(rk),          parameter     :: anew_isop  = 0.05_rk
        !> \brief isoprene leaf age factor for growing foliage (table 4 of guenther et al., 2012)
        real(rk),          parameter     :: agro_isop  = 0.6_rk
        !> \brief isoprene leaf age factor for mature foliage (table 4 of guenther et al., 2012)
        real(rk),          parameter     :: amat_isop  = 1.0_rk
        !> \brief isoprene leaf age factor for old/senescing foliage (table 4 of guenther et al., 2012)
        real(rk),          parameter     :: aold_isop  = 0.9_rk

!> \}

! plant-dependent emissions capacity/factors (efs) for myrcene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_myrc    =  70.0_rk      ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_myrc    =  70.0_rk      ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_myrc    =  60.0_rk      ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_myrc    =  80.0_rk      ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_myrc    =  30.0_rk      ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_myrc    =  80.0_rk      ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_myrc    =  30.0_rk      ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_myrc    =  30.0_rk      ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_myrc    =  30.0_rk      ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_myrc   =  50.0_rk      ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_myrc   =  30.0_rk      ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_myrc   =  0.3_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_myrc   =  0.3_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_myrc   =  0.3_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_myrc   =  0.3_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for myrcene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_myrc  = 2.0_rk
        real(rk),          parameter     :: agro_myrc  = 1.8_rk
        real(rk),          parameter     :: amat_myrc  = 1.0_rk
        real(rk),          parameter     :: aold_myrc  = 1.05_rk

! plant-dependent emissions capacity/factors (efs) for sabinene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_sabi    =  70.0_rk      ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_sabi    =  70.0_rk      ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_sabi    =  40.0_rk      ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_sabi    =  80.0_rk      ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_sabi    =  50.0_rk      ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_sabi    =  80.0_rk      ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_sabi    =  50.0_rk      ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_sabi    =  50.0_rk      ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_sabi    =  50.0_rk      ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_sabi   =  70.0_rk      ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_sabi   =  50.0_rk      ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_sabi   =  0.7_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_sabi   =  0.7_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_sabi   =  0.7_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_sabi   =  0.7_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for sabinene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_sabi  = 2.0_rk
        real(rk),          parameter     :: agro_sabi  = 1.8_rk
        real(rk),          parameter     :: amat_sabi  = 1.0_rk
        real(rk),          parameter     :: aold_sabi  = 1.05_rk

! plant-dependent emissions capacity/factors (efs) for limonene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_limo    =  100.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_limo    =  100.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_limo    =  130.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_limo    =  80.0_rk      ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_limo    =  80.0_rk      ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_limo    =  80.0_rk      ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_limo    =  80.0_rk      ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_limo    =  80.0_rk      ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_limo    =  60.0_rk      ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_limo   =  100.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_limo   =  60.0_rk      ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_limo   =  0.7_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_limo   =  0.7_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_limo   =  0.7_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_limo   =  0.7_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for limonene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_limo  = 2.0_rk
        real(rk),          parameter     :: agro_limo  = 1.8_rk
        real(rk),          parameter     :: amat_limo  = 1.0_rk
        real(rk),          parameter     :: aold_limo  = 1.05_rk

! plant-dependent emissions capacity/factors (efs) for 3-carene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_care    =  160.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_care    =  160.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_care    =  80.0_rk      ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_care    =  40.0_rk      ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_care    =  30.0_rk      ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_care    =  40.0_rk      ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_care    =  30.0_rk      ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_care    =  30.0_rk      ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_care    =  30.0_rk      ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_care   =  100.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_care   =  30.0_rk      ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_care   =  0.3_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_care   =  0.3_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_care   =  0.3_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_care   =  0.3_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for 3-carene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_care  = 2.0_rk
        real(rk),          parameter     :: agro_care  = 1.8_rk
        real(rk),          parameter     :: amat_care  = 1.0_rk
        real(rk),          parameter     :: aold_care  = 1.05_rk

! plant-dependent emissions capacity/factors (efs) for t-beta-ocimene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_ocim    =  70.0_rk      ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_ocim    =  70.0_rk      ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_ocim    =  60.0_rk      ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_ocim    =  150.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_ocim    =  120.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_ocim    =  150.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_ocim    =  120.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_ocim    =  120.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_ocim    =  90.0_rk      ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_ocim   =  150.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_ocim   =  90.0_rk      ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_ocim   =  2.0_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_ocim   =  2.0_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_ocim   =  2.0_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_ocim   =  2.0_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for t-beta-ocimene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_ocim  = 2.0_rk
        real(rk),          parameter     :: agro_ocim  = 1.8_rk
        real(rk),          parameter     :: amat_ocim  = 1.0_rk
        real(rk),          parameter     :: aold_ocim  = 1.05_rk

! plant-dependent emissions capacity/factors (efs) for beta-pinene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_bpin    =  300.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_bpin    =  300.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_bpin    =  200.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_bpin    =  120.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_bpin    =  130.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_bpin    =  120.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_bpin    =  130.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_bpin    =  130.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_bpin    =  100.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_bpin   =  150.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_bpin   =  100.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_bpin   =  1.5_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_bpin   =  1.5_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_bpin   =  1.5_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_bpin   =  1.5_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for beta-pinene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_bpin  = 2.0_rk
        real(rk),          parameter     :: agro_bpin  = 1.8_rk
        real(rk),          parameter     :: amat_bpin  = 1.0_rk
        real(rk),          parameter     :: aold_bpin  = 1.05_rk


! plant-dependent emissions capacity/factors (efs) for alpha-pinene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_apin    =  500.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_apin    =  500.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_apin    =  510.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_apin    =  600.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_apin    =  400.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_apin    =  600.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_apin    =  400.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_apin    =  400.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_apin    =  200.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_apin   =  300.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_apin   =  200.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_apin   =  2.0_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_apin   =  2.0_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_apin   =  2.0_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_apin   =  2.0_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for alpha-pinene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_apin  = 2.0_rk
        real(rk),          parameter     :: agro_apin  = 1.8_rk
        real(rk),          parameter     :: amat_apin  = 1.0_rk
        real(rk),          parameter     :: aold_apin  = 1.05_rk

! plant-dependent emissions capacity/factors (efs) for other monoterpenes (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
! ! other monoterpenes category (34 compounds):  see table 1 of guenther et al. (2012)
        real(rk),          parameter     :: ef1_mono    =  180.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_mono    =  180.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_mono    =  170.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_mono    =  150.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_mono    =  150.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_mono    =  150.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_mono    =  150.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_mono    =  150.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_mono    =  110.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_mono   =  200.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_mono   =  110.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_mono   =  5.0_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_mono   =  5.0_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_mono   =  5.0_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_mono   =  5.0_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for other monoterpenes as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_mono  = 2.0_rk
        real(rk),          parameter     :: agro_mono  = 1.8_rk
        real(rk),          parameter     :: amat_mono  = 1.0_rk
        real(rk),          parameter     :: aold_mono  = 1.05_rk

! plant-dependent emissions capacity/factors (efs) for alpha-farnesene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_farn    =  40.0_rk      ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_farn    =  40.0_rk      ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_farn    =  40.0_rk      ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_farn    =  60.0_rk      ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_farn    =  40.0_rk      ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_farn    =  60.0_rk      ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_farn    =  40.0_rk      ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_farn    =  40.0_rk      ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_farn    =  40.0_rk      ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_farn   =  40.0_rk      ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_farn   =  40.0_rk      ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_farn   =  3.0_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_farn   =  3.0_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_farn   =  3.0_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_farn   =  4.0_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for alpha-farnesene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_farn  = 0.4_rk
        real(rk),          parameter     :: agro_farn  = 0.6_rk
        real(rk),          parameter     :: amat_farn  = 1.0_rk
        real(rk),          parameter     :: aold_farn  = 0.95_rk

! plant-dependent emissions capacity/factors (efs) for beta-caryophyllene (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_cary    =  80.0_rk      ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_cary    =  80.0_rk      ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_cary    =  80.0_rk      ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_cary    =  60.0_rk      ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_cary    =  40.0_rk      ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_cary    =  60.0_rk      ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_cary    =  40.0_rk      ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_cary    =  40.0_rk      ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_cary    =  50.0_rk      ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_cary   =  50.0_rk      ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_cary   =  50.0_rk      ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_cary   =  1.0_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_cary   =  1.0_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_cary   =  1.0_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_cary   =  4.0_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for beta-caryophyllene as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_cary  = 0.4_rk
        real(rk),          parameter     :: agro_cary  = 0.6_rk
        real(rk),          parameter     :: amat_cary  = 1.0_rk
        real(rk),          parameter     :: aold_cary  = 0.95_rk

! plant-dependent emissions capacity/factors (efs) for other sesquieterpenes (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
! other sesquiterpenes category (30 compounds):  see table 1 of guenther et al. (2012)
        real(rk),          parameter     :: ef1_sesq    =  120.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_sesq    =  120.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_sesq    =  120.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_sesq    =  120.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_sesq    =  100.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_sesq    =  120.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_sesq    =  100.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_sesq    =  100.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_sesq    =  100.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_sesq   =  100.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_sesq   =  100.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_sesq   =  1.0_rk       ! arctic c3 grass
        real(rk),          parameter     :: ef13_sesq   =  1.0_rk       ! cool c3 grass
        real(rk),          parameter     :: ef14_sesq   =  1.0_rk       ! warm c4 grass
        real(rk),          parameter     :: ef15_sesq   =  1.0_rk       ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for other sesquieterpenes as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_sesq  = 0.4_rk
        real(rk),          parameter     :: agro_sesq  = 0.6_rk
        real(rk),          parameter     :: amat_sesq  = 1.0_rk
        real(rk),          parameter     :: aold_sesq  = 0.95_rk

! plant-dependent emissions capacity/factors (efs) for 232-mbo (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_mbol    =  700.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_mbol    =  60.0_rk      ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_mbol    =  0.01_rk      ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_mbol    =  0.01_rk      ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_mbol    =  0.01_rk      ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_mbol    =  0.01_rk      ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_mbol    =  0.01_rk      ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_mbol    =  2.0_rk       ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_mbol    =  0.01_rk      ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_mbol   =  0.01_rk      ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_mbol   =  0.01_rk      ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_mbol   =  0.01_rk      ! arctic c3 grass
        real(rk),          parameter     :: ef13_mbol   =  0.01_rk      ! cool c3 grass
        real(rk),          parameter     :: ef14_mbol   =  0.01_rk      ! warm c4 grass
        real(rk),          parameter     :: ef15_mbol   =  0.01_rk      ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for 232-mbo as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_mbol  = 0.05_rk
        real(rk),          parameter     :: agro_mbol  = 0.6_rk
        real(rk),          parameter     :: amat_mbol  = 1.0_rk
        real(rk),          parameter     :: aold_mbol  = 0.9_rk

! plant-dependent emissions capacity/factors (efs) for methanol (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_meth    =  900.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_meth    =  900.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_meth    =  900.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_meth    =  500.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_meth    =  900.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_meth    =  500.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_meth    =  900.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_meth    =  900.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_meth    =  900.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_meth   =  900.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_meth   =  900.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_meth   =  500.0_rk     ! arctic c3 grass
        real(rk),          parameter     :: ef13_meth   =  500.0_rk     ! cool c3 grass
        real(rk),          parameter     :: ef14_meth   =  500.0_rk     ! warm c4 grass
        real(rk),          parameter     :: ef15_meth   =  900.0_rk     ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for methanol as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_meth  = 3.5_rk
        real(rk),          parameter     :: agro_meth  = 3.0_rk
        real(rk),          parameter     :: amat_meth  = 1.0_rk
        real(rk),          parameter     :: aold_meth  = 1.2_rk

! plant-dependent emissions capacity/factors (efs) for acetone (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_acet    =  240.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_acet    =  240.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_acet    =  240.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_acet    =  240.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_acet    =  240.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_acet    =  240.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_acet    =  240.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_acet    =  240.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_acet    =  240.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_acet   =  240.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_acet   =  240.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_acet   =  80.0_rk      ! arctic c3 grass
        real(rk),          parameter     :: ef13_acet   =  80.0_rk      ! cool c3 grass
        real(rk),          parameter     :: ef14_acet   =  80.0_rk      ! warm c4 grass
        real(rk),          parameter     :: ef15_acet   =  80.0_rk      ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for acetone as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_acet  = 1.0_rk
        real(rk),          parameter     :: agro_acet  = 1.0_rk
        real(rk),          parameter     :: amat_acet  = 1.0_rk
        real(rk),          parameter     :: aold_acet  = 1.0_rk

! plant-dependent emissions capacity/factors (efs) for carbon monoxide (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
        real(rk),          parameter     :: ef1_co      =  600.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_co      =  600.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_co      =  600.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_co      =  600.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_co      =  600.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_co      =  600.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_co      =  600.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_co      =  600.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_co      =  600.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_co     =  600.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_co     =  600.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_co     =  600.0_rk     ! arctic c3 grass
        real(rk),          parameter     :: ef13_co     =  600.0_rk     ! cool c3 grass
        real(rk),          parameter     :: ef14_co     =  600.0_rk     ! warm c4 grass
        real(rk),          parameter     :: ef15_co     =  600.0_rk     ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for carbon monoxide as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_co  = 1.0_rk
        real(rk),          parameter     :: agro_co  = 1.0_rk
        real(rk),          parameter     :: amat_co  = 1.0_rk
        real(rk),          parameter     :: aold_co  = 1.0_rk

! plant-dependent emissions capacity/factors (efs) for bidi voc species (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
! bidirectional voc (5 compounds): see table 1 of guenther et al. (2012)
        real(rk),          parameter     :: ef1_bvoc    =  500.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_bvoc    =  500.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_bvoc    =  500.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_bvoc    =  500.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_bvoc    =  500.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_bvoc    =  500.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_bvoc    =  500.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_bvoc    =  500.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_bvoc    =  500.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_bvoc   =  500.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_bvoc   =  500.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_bvoc   =  80.0_rk      ! arctic c3 grass
        real(rk),          parameter     :: ef13_bvoc   =  80.0_rk      ! cool c3 grass
        real(rk),          parameter     :: ef14_bvoc   =  80.0_rk      ! warm c4 grass
        real(rk),          parameter     :: ef15_bvoc   =  80.0_rk      ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for bidi voc species as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_bvoc  = 1.0_rk
        real(rk),          parameter     :: agro_bvoc  = 1.0_rk
        real(rk),          parameter     :: amat_bvoc  = 1.0_rk
        real(rk),          parameter     :: aold_bvoc  = 1.0_rk

! plant-dependent emissions capacity/factors (efs) for stress vocs (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
! stress voc (15 compounds): see table 1 of guenther et al. (2012)
        real(rk),          parameter     :: ef1_svoc    =  300.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_svoc    =  300.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_svoc    =  300.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_svoc    =  300.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_svoc    =  300.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_svoc    =  300.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_svoc    =  300.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_svoc    =  300.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_svoc    =  300.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_svoc   =  300.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_svoc   =  300.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_svoc   =  300.0_rk     ! arctic c3 grass
        real(rk),          parameter     :: ef13_svoc   =  300.0_rk     ! cool c3 grass
        real(rk),          parameter     :: ef14_svoc   =  300.0_rk     ! warm c4 grass
        real(rk),          parameter     :: ef15_svoc   =  300.0_rk     ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for stress vocs as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_svoc  = 1.0_rk
        real(rk),          parameter     :: agro_svoc  = 1.0_rk
        real(rk),          parameter     :: amat_svoc  = 1.0_rk
        real(rk),          parameter     :: aold_svoc  = 1.0_rk

! plant-dependent emissions capacity/factors (efs) for other vocs (tables 2-3 of guenther et al., 2012) (ug/m2 hr)
! other voc (49 compounds): see table 1 of guenther et al. (2012)
        real(rk),          parameter     :: ef1_ovoc    =  140.0_rk     ! needleleaf evergreen temperate tree
        real(rk),          parameter     :: ef2_ovoc    =  140.0_rk     ! needleleaf evergreen boreal tree
        real(rk),          parameter     :: ef3_ovoc    =  140.0_rk     ! needleleaf deciduous boreal tree
        real(rk),          parameter     :: ef4_ovoc    =  140.0_rk     ! broadleaf evergreen tropical tree
        real(rk),          parameter     :: ef5_ovoc    =  140.0_rk     ! broadleaf evergreen temperate tree
        real(rk),          parameter     :: ef6_ovoc    =  140.0_rk     ! broadleaf deciduous tropical tree
        real(rk),          parameter     :: ef7_ovoc    =  140.0_rk     ! broadleaf deciduous temperate tree
        real(rk),          parameter     :: ef8_ovoc    =  140.0_rk     ! broadleaf deciduous boreal tree
        real(rk),          parameter     :: ef9_ovoc    =  140.0_rk     ! broadleaf evergreen temperate shrub
        real(rk),          parameter     :: ef10_ovoc   =  140.0_rk     ! broadleaf deciduous temperate shrub
        real(rk),          parameter     :: ef11_ovoc   =  140.0_rk     ! broadleaf deciduous boreal shrub
        real(rk),          parameter     :: ef12_ovoc   =  140.0_rk     ! arctic c3 grass
        real(rk),          parameter     :: ef13_ovoc   =  140.0_rk     ! cool c3 grass
        real(rk),          parameter     :: ef14_ovoc   =  140.0_rk     ! warm c4 grass
        real(rk),          parameter     :: ef15_ovoc   =  140.0_rk     ! crop1

!empirical factors or coefficients for: growing, mature, and old/senescing foliage, for other vocs as per table 4 of guenther et al., 2012
        real(rk),          parameter     :: anew_ovoc  = 1.0_rk
        real(rk),          parameter     :: agro_ovoc  = 1.0_rk
        real(rk),          parameter     :: amat_ovoc  = 1.0_rk
        real(rk),          parameter     :: aold_ovoc  = 1.0_rk

! species-dependent parameterized canopy model parameters (table 4 of guenther et al., 2012)
        real(rk),          parameter     :: ldf_isop         =  1.0_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_isop        =  0.13_rk    !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_isop         =  95.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_isop         =  2.0_rk     !empirical coefficient
        real(rk),          parameter     :: caq_isop         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_isop         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_isop        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_isop         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_isop         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_isop        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_isop         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_isop         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_isop        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_isop         =  1.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_isop         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_isop        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_myrc         =  0.6_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_myrc        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_myrc         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_myrc         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_myrc         =  5.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_myrc         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_myrc        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_myrc         =  5.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_myrc         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_myrc        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_myrc         =  5.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_myrc         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_myrc        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_myrc         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_myrc         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_myrc        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_sabi         =  0.6_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_sabi        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_sabi         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_sabi         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_sabi         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_sabi         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_sabi        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_sabi         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_sabi         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_sabi        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_sabi         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_sabi         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_sabi        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_sabi         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_sabi         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_sabi        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_limo         =  0.2_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_limo        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_limo         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_limo         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_limo         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_limo         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_limo        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_limo         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_limo         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_limo        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_limo         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_limo         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_limo        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_limo         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_limo         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_limo        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_care         =  0.2_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_care        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_care         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_care         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_care         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_care         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_care        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_care         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_care         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_care        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_care         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_care         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_care        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_care         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_care         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_care        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_ocim         =  0.8_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_ocim        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_ocim         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_ocim         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_ocim         =  5.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_ocim         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_ocim        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_ocim         =  5.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_ocim         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_ocim        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_ocim         =  5.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_ocim         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_ocim        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_ocim         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_ocim         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_ocim        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_bpin         =  0.2_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_bpin        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_bpin         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_bpin         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_bpin         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_bpin         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_bpin        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_bpin         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_bpin         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_bpin        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_bpin         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_bpin         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_bpin        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_bpin         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_bpin         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_bpin        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_apin         =  0.6_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_apin        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_apin         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_apin         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_apin         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_apin         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_apin        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_apin         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_apin         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_apin        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_apin         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_apin         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_apin        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_apin         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_apin         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_apin        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_mono         =  0.4_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_mono        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_mono         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_mono         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_mono         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_mono         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_mono        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_mono         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_mono         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_mono        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_mono         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_mono         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_mono        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_mono         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_mono         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_mono        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_farn         =  0.5_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_farn        =  0.17_rk    !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_farn         =  130.0_rk   !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_farn         =  2.37_rk    !empirical coefficient
        real(rk),          parameter     :: caq_farn         =  5.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_farn         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_farn        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_farn         =  5.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_farn         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_farn        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_farn         =  5.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_farn         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_farn        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_farn         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_farn         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_farn        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_cary         =  0.5_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_cary        =  0.17_rk    !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_cary         =  130.0_rk   !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_cary         =  2.37_rk    !empirical coefficient
        real(rk),          parameter     :: caq_cary         =  5.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_cary         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_cary        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_cary         =  5.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_cary         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_cary        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_cary         =  5.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_cary         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_cary        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_cary         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_cary         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_cary        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_sesq         =  0.5_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_sesq        =  0.17_rk    !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_sesq         =  130.0_rk   !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_sesq         =  2.37_rk    !empirical coefficient
        real(rk),          parameter     :: caq_sesq         =  5.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_sesq         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_sesq        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_sesq         =  5.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_sesq         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_sesq        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_sesq         =  5.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_sesq         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_sesq        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_sesq         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_sesq         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_sesq        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_mbol         =  1.0_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_mbol        =  0.13_rk    !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_mbol         =  95.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_mbol         =  2.0_rk     !empirical coefficient
        real(rk),          parameter     :: caq_mbol         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_mbol         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_mbol        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_mbol         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_mbol         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_mbol        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_mbol         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_mbol         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_mbol        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_mbol         =  1.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_mbol         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_mbol        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_meth         =  0.8_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_meth        =  0.08_rk    !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_meth         =  60.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_meth         =  1.6_rk     !empirical coefficient
        real(rk),          parameter     :: caq_meth         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_meth         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_meth        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_meth         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_meth         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_meth        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_meth         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_meth         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_meth        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_meth         =  1.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_meth         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_meth        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_acet         =  0.2_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_acet        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_acet         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_acet         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_acet         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_acet         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_acet        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_acet         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_acet         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_acet        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_acet         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_acet         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_acet        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_acet         =  1.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_acet         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_acet        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_co           =  1.0_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_co          =  0.08_rk    !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_co           =  60.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_co           =  1.6_rk     !empirical coefficient
        real(rk),          parameter     :: caq_co           =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_co           =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_co          =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_co           =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_co           =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_co          =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_co           =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_co           =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_co          =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_co           =  1.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_co           =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_co          =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_bvoc         =  0.8_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_bvoc        =  0.13_rk    !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_bvoc         =  95.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_bvoc         =  2.0_rk     !empirical coefficient
        real(rk),          parameter     :: caq_bvoc         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_bvoc         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_bvoc        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_bvoc         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_bvoc         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_bvoc        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_bvoc         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_bvoc         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_bvoc        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_bvoc         =  1.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_bvoc         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_bvoc        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_svoc         =  0.8_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_svoc        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_svoc         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_svoc         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_svoc         =  5.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_svoc         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_svoc        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_svoc         =  5.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_svoc         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_svoc        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_svoc         =  5.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_svoc         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_svoc        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_svoc         =  5.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_svoc         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_svoc        =  8.0_rk     !delta threshold for high wind stress (m/s)

        real(rk),          parameter     :: ldf_ovoc         =  0.2_rk     !light-dependent fraction
        real(rk),          parameter     :: beta_ovoc        =  0.1_rk     !empirical coefficient for temperature dependence of light-independent fraction
        real(rk),          parameter     :: ct1_ovoc         =  80.0_rk    !activation energy (kj/mol)
        real(rk),          parameter     :: ceo_ovoc         =  1.83_rk    !empirical coefficient
        real(rk),          parameter     :: caq_ovoc         =  1.0_rk     !coefficient for poor air quality stress
        real(rk),          parameter     :: taq_ovoc         =  20.0_rk    !threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: dtaq_ovoc        =  30.0_rk    !delta threshold for poor air quality stress (ppm-hours)
        real(rk),          parameter     :: cht_ovoc         =  1.0_rk     !coefficient for high temperature stress
        real(rk),          parameter     :: tht_ovoc         =  313.15_rk  !threshold for high temperature stress (k)
        real(rk),          parameter     :: dtht_ovoc        =  8.0_rk     !delta threshold for high temperature stress (k)
        real(rk),          parameter     :: clt_ovoc         =  1.0_rk     !coefficient for low temperature stress
        real(rk),          parameter     :: tlt_ovoc         =  283.15_rk  !threshold for low temperature stress (k)
        real(rk),          parameter     :: dtlt_ovoc        =  8.0_rk     !delta threshold for low temperature stress (k)
        real(rk),          parameter     :: chw_ovoc         =  1.0_rk     !coefficient for high wind stress
        real(rk),          parameter     :: thw_ovoc         =  12.0_rk    !threshold for high wind stress (m/s)
        real(rk),          parameter     :: dthw_ovoc        =  8.0_rk     !delta threshold for high wind stress (m/s)

! set tree and species dependent coefficients
        if (emi_ind .eq. 1 ) then
            ldf  = ldf_isop
            beta = beta_isop
            ct1  = ct1_isop
            ceo  = ceo_isop
            ef1  = ef1_isop
            ef2  = ef2_isop
            ef3  = ef3_isop
            ef4  = ef4_isop
            ef5  = ef5_isop
            ef6  = ef6_isop
            ef7  = ef7_isop
            ef8  = ef8_isop
            ef9  = ef9_isop
            ef10 = ef10_isop
            ef11 = ef11_isop
            ef12 = ef12_isop
            ef13 = ef13_isop
            ef14 = ef14_isop
            ef15 = ef15_isop
            anew  = anew_isop
            agro  = agro_isop
            amat  = amat_isop
            aold  = aold_isop
            caq   = caq_isop
            taq   = taq_isop
            dtaq  = dtaq_isop
            cht   = cht_isop
            tht   = tht_isop
            dtht  = dtht_isop
            clt   = clt_isop
            tlt   = tlt_isop
            dtlt  = dtlt_isop
            chw   = chw_isop
            thw   = thw_isop
            dthw  = dthw_isop
        else if (emi_ind .eq. 2 ) then
            ldf  = ldf_myrc
            beta = beta_myrc
            ct1 = ct1_myrc
            ceo = ceo_myrc
            ef1  = ef1_myrc
            ef2  = ef2_myrc
            ef3  = ef3_myrc
            ef4  = ef4_myrc
            ef5  = ef5_myrc
            ef6  = ef6_myrc
            ef7  = ef7_myrc
            ef8  = ef8_myrc
            ef9  = ef9_myrc
            ef10 = ef10_myrc
            ef11 = ef11_myrc
            ef12 = ef12_myrc
            ef13 = ef13_myrc
            ef14 = ef14_myrc
            ef15 = ef15_myrc
            anew  = anew_myrc
            agro  = agro_myrc
            amat  = amat_myrc
            aold  = aold_myrc
            caq   = caq_myrc
            taq   = taq_myrc
            dtaq  = dtaq_myrc
            cht   = cht_myrc
            tht   = tht_myrc
            dtht  = dtht_myrc
            clt   = clt_myrc
            tlt   = tlt_myrc
            dtlt  = dtlt_myrc
            chw   = chw_myrc
            thw   = thw_myrc
            dthw  = dthw_myrc
        else if (emi_ind .eq. 3 ) then
            ldf  = ldf_sabi
            beta = beta_sabi
            ct1 = ct1_sabi
            ceo = ceo_sabi
            ef1  = ef1_sabi
            ef2  = ef2_sabi
            ef3  = ef3_sabi
            ef4  = ef4_sabi
            ef5  = ef5_sabi
            ef6  = ef6_sabi
            ef7  = ef7_sabi
            ef8  = ef8_sabi
            ef9  = ef9_sabi
            ef10 = ef10_sabi
            ef11 = ef11_sabi
            ef12 = ef12_sabi
            ef13 = ef13_sabi
            ef14 = ef14_sabi
            ef15 = ef15_sabi
            anew  = anew_sabi
            agro  = agro_sabi
            amat  = amat_sabi
            aold  = aold_sabi
            caq   = caq_sabi
            taq   = taq_sabi
            dtaq  = dtaq_sabi
            cht   = cht_sabi
            tht   = tht_sabi
            dtht  = dtht_sabi
            clt   = clt_sabi
            tlt   = tlt_sabi
            dtlt  = dtlt_sabi
            chw   = chw_sabi
            thw   = thw_sabi
            dthw  = dthw_sabi
        else if (emi_ind .eq. 4 ) then
            ldf  = ldf_limo
            beta = beta_limo
            ct1 = ct1_limo
            ceo = ceo_limo
            ef1  = ef1_limo
            ef2  = ef2_limo
            ef3  = ef3_limo
            ef4  = ef4_limo
            ef5  = ef5_limo
            ef6  = ef6_limo
            ef7  = ef7_limo
            ef8  = ef8_limo
            ef9  = ef9_limo
            ef10 = ef10_limo
            ef11 = ef11_limo
            ef12 = ef12_limo
            ef13 = ef13_limo
            ef14 = ef14_limo
            ef15 = ef15_limo
            anew  = anew_limo
            agro  = agro_limo
            amat  = amat_limo
            aold  = aold_limo
            caq   = caq_limo
            taq   = taq_limo
            dtaq  = dtaq_limo
            cht   = cht_limo
            tht   = tht_limo
            dtht  = dtht_limo
            clt   = clt_limo
            tlt   = tlt_limo
            dtlt  = dtlt_limo
            chw   = chw_limo
            thw   = thw_limo
            dthw  = dthw_limo
        else if (emi_ind .eq. 5 ) then
            ldf  = ldf_care
            beta = beta_care
            ct1 = ct1_care
            ceo = ceo_care
            ef1  = ef1_care
            ef2  = ef2_care
            ef3  = ef3_care
            ef4  = ef4_care
            ef5  = ef5_care
            ef6  = ef6_care
            ef7  = ef7_care
            ef8  = ef8_care
            ef9  = ef9_care
            ef10 = ef10_care
            ef11 = ef11_care
            ef12 = ef12_care
            ef13 = ef13_care
            ef14 = ef14_care
            ef15 = ef15_care
            anew  = anew_care
            agro  = agro_care
            amat  = amat_care
            aold  = aold_care
            caq   = caq_care
            taq   = taq_care
            dtaq  = dtaq_care
            cht   = cht_care
            tht   = tht_care
            dtht  = dtht_care
            clt   = clt_care
            tlt   = tlt_care
            dtlt  = dtlt_care
            chw   = chw_care
            thw   = thw_care
            dthw  = dthw_care
        else if (emi_ind .eq. 6 ) then
            ldf  = ldf_ocim
            beta = beta_ocim
            ct1 = ct1_ocim
            ceo = ceo_ocim
            ef1  = ef1_ocim
            ef2  = ef2_ocim
            ef3  = ef3_ocim
            ef4  = ef4_ocim
            ef5  = ef5_ocim
            ef6  = ef6_ocim
            ef7  = ef7_ocim
            ef8  = ef8_ocim
            ef9  = ef9_ocim
            ef10 = ef10_ocim
            ef11 = ef11_ocim
            ef12 = ef12_ocim
            ef13 = ef13_ocim
            ef14 = ef14_ocim
            ef15 = ef15_ocim
            anew  = anew_ocim
            agro  = agro_ocim
            amat  = amat_ocim
            aold  = aold_ocim
            caq   = caq_ocim
            taq   = taq_ocim
            dtaq  = dtaq_ocim
            cht   = cht_ocim
            tht   = tht_ocim
            dtht  = dtht_ocim
            clt   = clt_ocim
            tlt   = tlt_ocim
            dtlt  = dtlt_ocim
            chw   = chw_ocim
            thw   = thw_ocim
            dthw  = dthw_ocim
        else if (emi_ind .eq. 7 ) then
            ldf  = ldf_bpin
            beta = beta_bpin
            ct1 = ct1_bpin
            ceo = ceo_bpin
            ef1  = ef1_bpin
            ef2  = ef2_bpin
            ef3  = ef3_bpin
            ef4  = ef4_bpin
            ef5  = ef5_bpin
            ef6  = ef6_bpin
            ef7  = ef7_bpin
            ef8  = ef8_bpin
            ef9  = ef9_bpin
            ef10 = ef10_bpin
            ef11 = ef11_bpin
            ef12 = ef12_bpin
            ef13 = ef13_bpin
            ef14 = ef14_bpin
            ef15 = ef15_bpin
            anew  = anew_bpin
            agro  = agro_bpin
            amat  = amat_bpin
            aold  = aold_bpin
            caq   = caq_bpin
            taq   = taq_bpin
            dtaq  = dtaq_bpin
            cht   = cht_bpin
            tht   = tht_bpin
            dtht  = dtht_bpin
            clt   = clt_bpin
            tlt   = tlt_bpin
            dtlt  = dtlt_bpin
            chw   = chw_bpin
            thw   = thw_bpin
            dthw  = dthw_bpin
        else if (emi_ind .eq. 8 ) then
            ldf  = ldf_apin
            beta = beta_apin
            ct1 = ct1_apin
            ceo = ceo_apin
            ef1  = ef1_apin
            ef2  = ef2_apin
            ef3  = ef3_apin
            ef4  = ef4_apin
            ef5  = ef5_apin
            ef6  = ef6_apin
            ef7  = ef7_apin
            ef8  = ef8_apin
            ef9  = ef9_apin
            ef10 = ef10_apin
            ef11 = ef11_apin
            ef12 = ef12_apin
            ef13 = ef13_apin
            ef14 = ef14_apin
            ef15 = ef15_apin
            anew  = anew_apin
            agro  = agro_apin
            amat  = amat_apin
            aold  = aold_apin
            caq   = caq_apin
            taq   = taq_apin
            dtaq  = dtaq_apin
            cht   = cht_apin
            tht   = tht_apin
            dtht  = dtht_apin
            clt   = clt_apin
            tlt   = tlt_apin
            dtlt  = dtlt_apin
            chw   = chw_apin
            thw   = thw_apin
            dthw  = dthw_apin
        else if (emi_ind .eq. 9 ) then
            ldf  = ldf_mono
            beta = beta_mono
            ct1 = ct1_mono
            ceo = ceo_mono
            ef1  = ef1_mono
            ef2  = ef2_mono
            ef3  = ef3_mono
            ef4  = ef4_mono
            ef5  = ef5_mono
            ef6  = ef6_mono
            ef7  = ef7_mono
            ef8  = ef8_mono
            ef9  = ef9_mono
            ef10 = ef10_mono
            ef11 = ef11_mono
            ef12 = ef12_mono
            ef13 = ef13_mono
            ef14 = ef14_mono
            ef15 = ef15_mono
            anew  = anew_mono
            agro  = agro_mono
            amat  = amat_mono
            aold  = aold_mono
            caq   = caq_mono
            taq   = taq_mono
            dtaq  = dtaq_mono
            cht   = cht_mono
            tht   = tht_mono
            dtht  = dtht_mono
            clt   = clt_mono
            tlt   = tlt_mono
            dtlt  = dtlt_mono
            chw   = chw_mono
            thw   = thw_mono
            dthw  = dthw_mono
        else if (emi_ind .eq. 10 ) then
            ldf  = ldf_farn
            beta = beta_farn
            ct1 = ct1_farn
            ceo = ceo_farn
            ef1  = ef1_farn
            ef2  = ef2_farn
            ef3  = ef3_farn
            ef4  = ef4_farn
            ef5  = ef5_farn
            ef6  = ef6_farn
            ef7  = ef7_farn
            ef8  = ef8_farn
            ef9  = ef9_farn
            ef10 = ef10_farn
            ef11 = ef11_farn
            ef12 = ef12_farn
            ef13 = ef13_farn
            ef14 = ef14_farn
            ef15 = ef15_farn
            anew  = anew_farn
            agro  = agro_farn
            amat  = amat_farn
            aold  = aold_farn
            caq   = caq_farn
            taq   = taq_farn
            dtaq  = dtaq_farn
            cht   = cht_farn
            tht   = tht_farn
            dtht  = dtht_farn
            clt   = clt_farn
            tlt   = tlt_farn
            dtlt  = dtlt_farn
            chw   = chw_farn
            thw   = thw_farn
            dthw  = dthw_farn
        else if (emi_ind .eq. 11 ) then
            ldf  = ldf_cary
            beta = beta_cary
            ct1 = ct1_cary
            ceo = ceo_cary
            ef1  = ef1_cary
            ef2  = ef2_cary
            ef3  = ef3_cary
            ef4  = ef4_cary
            ef5  = ef5_cary
            ef6  = ef6_cary
            ef7  = ef7_cary
            ef8  = ef8_cary
            ef9  = ef9_cary
            ef10 = ef10_cary
            ef11 = ef11_cary
            ef12 = ef12_cary
            ef13 = ef13_cary
            ef14 = ef14_cary
            ef15 = ef15_cary
            anew  = anew_cary
            agro  = agro_cary
            amat  = amat_cary
            aold  = aold_cary
            caq   = caq_cary
            taq   = taq_cary
            dtaq  = dtaq_cary
            cht   = cht_cary
            tht   = tht_cary
            dtht  = dtht_cary
            clt   = clt_cary
            tlt   = tlt_cary
            dtlt  = dtlt_cary
            chw   = chw_cary
            thw   = thw_cary
            dthw  = dthw_cary
        else if (emi_ind .eq. 12 ) then
            ldf  = ldf_sesq
            beta = beta_sesq
            ct1 = ct1_sesq
            ceo = ceo_sesq
            ef1  = ef1_sesq
            ef2  = ef2_sesq
            ef3  = ef3_sesq
            ef4  = ef4_sesq
            ef5  = ef5_sesq
            ef6  = ef6_sesq
            ef7  = ef7_sesq
            ef8  = ef8_sesq
            ef9  = ef9_sesq
            ef10 = ef10_sesq
            ef11 = ef11_sesq
            ef12 = ef12_sesq
            ef13 = ef13_sesq
            ef14 = ef14_sesq
            ef15 = ef15_sesq
            anew  = anew_sesq
            agro  = agro_sesq
            amat  = amat_sesq
            aold  = aold_sesq
            caq   = caq_sesq
            taq   = taq_sesq
            dtaq  = dtaq_sesq
            cht   = cht_sesq
            tht   = tht_sesq
            dtht  = dtht_sesq
            clt   = clt_sesq
            tlt   = tlt_sesq
            dtlt  = dtlt_sesq
            chw   = chw_sesq
            thw   = thw_sesq
            dthw  = dthw_sesq
        else if (emi_ind .eq. 13 ) then
            ldf  = ldf_mbol
            beta = beta_mbol
            ct1 = ct1_mbol
            ceo = ceo_mbol
            ef1  = ef1_mbol
            ef2  = ef2_mbol
            ef3  = ef3_mbol
            ef4  = ef4_mbol
            ef5  = ef5_mbol
            ef6  = ef6_mbol
            ef7  = ef7_mbol
            ef8  = ef8_mbol
            ef9  = ef9_mbol
            ef10 = ef10_mbol
            ef11 = ef11_mbol
            ef12 = ef12_mbol
            ef13 = ef13_mbol
            ef14 = ef14_mbol
            ef15 = ef15_mbol
            anew  = anew_mbol
            agro  = agro_mbol
            amat  = amat_mbol
            aold  = aold_mbol
            caq   = caq_mbol
            taq   = taq_mbol
            dtaq  = dtaq_mbol
            cht   = cht_mbol
            tht   = tht_mbol
            dtht  = dtht_mbol
            clt   = clt_mbol
            tlt   = tlt_mbol
            dtlt  = dtlt_mbol
            chw   = chw_mbol
            thw   = thw_mbol
            dthw  = dthw_mbol
        else if (emi_ind .eq. 14 ) then
            ldf  = ldf_meth
            beta = beta_meth
            ct1 = ct1_meth
            ceo = ceo_meth
            ef1  = ef1_meth
            ef2  = ef2_meth
            ef3  = ef3_meth
            ef4  = ef4_meth
            ef5  = ef5_meth
            ef6  = ef6_meth
            ef7  = ef7_meth
            ef8  = ef8_meth
            ef9  = ef9_meth
            ef10 = ef10_meth
            ef11 = ef11_meth
            ef12 = ef12_meth
            ef13 = ef13_meth
            ef14 = ef14_meth
            ef15 = ef15_meth
            anew  = anew_meth
            agro  = agro_meth
            amat  = amat_meth
            aold  = aold_meth
            caq   = caq_meth
            taq   = taq_meth
            dtaq  = dtaq_meth
            cht   = cht_meth
            tht   = tht_meth
            dtht  = dtht_meth
            clt   = clt_meth
            tlt   = tlt_meth
            dtlt  = dtlt_meth
            chw   = chw_meth
            thw   = thw_meth
            dthw  = dthw_meth
        else if (emi_ind .eq. 15 ) then
            ldf  = ldf_acet
            beta = beta_acet
            ct1 = ct1_acet
            ceo = ceo_acet
            ef1  = ef1_acet
            ef2  = ef2_acet
            ef3  = ef3_acet
            ef4  = ef4_acet
            ef5  = ef5_acet
            ef6  = ef6_acet
            ef7  = ef7_acet
            ef8  = ef8_acet
            ef9  = ef9_acet
            ef10 = ef10_acet
            ef11 = ef11_acet
            ef12 = ef12_acet
            ef13 = ef13_acet
            ef14 = ef14_acet
            ef15 = ef15_acet
            anew  = anew_acet
            agro  = agro_acet
            amat  = amat_acet
            aold  = aold_acet
            caq   = caq_acet
            taq   = taq_acet
            dtaq  = dtaq_acet
            cht   = cht_acet
            tht   = tht_acet
            dtht  = dtht_acet
            clt   = clt_acet
            tlt   = tlt_acet
            dtlt  = dtlt_acet
            chw   = chw_acet
            thw   = thw_acet
            dthw  = dthw_acet
        else if (emi_ind .eq. 16 ) then
            ldf  = ldf_co
            beta = beta_co
            ct1 = ct1_co
            ceo = ceo_co
            ef1  = ef1_co
            ef2  = ef2_co
            ef3  = ef3_co
            ef4  = ef4_co
            ef5  = ef5_co
            ef6  = ef6_co
            ef7  = ef7_co
            ef8  = ef8_co
            ef9  = ef9_co
            ef10 = ef10_co
            ef11 = ef11_co
            ef12 = ef12_co
            ef13 = ef13_co
            ef14 = ef14_co
            ef15 = ef15_co
            anew  = anew_co
            agro  = agro_co
            amat  = amat_co
            aold  = aold_co
            caq   = caq_co
            taq   = taq_co
            dtaq  = dtaq_co
            cht   = cht_co
            tht   = tht_co
            dtht  = dtht_co
            clt   = clt_co
            tlt   = tlt_co
            dtlt  = dtlt_co
            chw   = chw_co
            thw   = thw_co
            dthw  = dthw_co
        else if (emi_ind .eq. 17 ) then
            ldf  = ldf_bvoc
            beta = beta_bvoc
            ct1 = ct1_bvoc
            ceo = ceo_bvoc
            ef1  = ef1_bvoc
            ef2  = ef2_bvoc
            ef3  = ef3_bvoc
            ef4  = ef4_bvoc
            ef5  = ef5_bvoc
            ef6  = ef6_bvoc
            ef7  = ef7_bvoc
            ef8  = ef8_bvoc
            ef9  = ef9_bvoc
            ef10 = ef10_bvoc
            ef11 = ef11_bvoc
            ef12 = ef12_bvoc
            ef13 = ef13_bvoc
            ef14 = ef14_bvoc
            ef15 = ef15_bvoc
            anew  = anew_bvoc
            agro  = agro_bvoc
            amat  = amat_bvoc
            aold  = aold_bvoc
            caq   = caq_bvoc
            taq   = taq_bvoc
            dtaq  = dtaq_bvoc
            cht   = cht_bvoc
            tht   = tht_bvoc
            dtht  = dtht_bvoc
            clt   = clt_bvoc
            tlt   = tlt_bvoc
            dtlt  = dtlt_bvoc
            chw   = chw_bvoc
            thw   = thw_bvoc
            dthw  = dthw_bvoc
        else if (emi_ind .eq. 18 ) then
            ldf  = ldf_svoc
            beta = beta_svoc
            ct1 = ct1_svoc
            ceo = ceo_svoc
            ef1  = ef1_svoc
            ef2  = ef2_svoc
            ef3  = ef3_svoc
            ef4  = ef4_svoc
            ef5  = ef5_svoc
            ef6  = ef6_svoc
            ef7  = ef7_svoc
            ef8  = ef8_svoc
            ef9  = ef9_svoc
            ef10 = ef10_svoc
            ef11 = ef11_svoc
            ef12 = ef12_svoc
            ef13 = ef13_svoc
            ef14 = ef14_svoc
            ef15 = ef15_svoc
            anew  = anew_svoc
            agro  = agro_svoc
            amat  = amat_svoc
            aold  = aold_svoc
            caq   = caq_svoc
            taq   = taq_svoc
            dtaq  = dtaq_svoc
            cht   = cht_svoc
            tht   = tht_svoc
            dtht  = dtht_svoc
            clt   = clt_svoc
            tlt   = tlt_svoc
            dtlt  = dtlt_svoc
            chw   = chw_svoc
            thw   = thw_svoc
            dthw  = dthw_svoc
        else   ! emi_ind = 19
            ldf  = ldf_ovoc
            beta = beta_ovoc
            ct1 = ct1_ovoc
            ceo = ceo_ovoc
            ef1  = ef1_ovoc
            ef2  = ef2_ovoc
            ef3  = ef3_ovoc
            ef4  = ef4_ovoc
            ef5  = ef5_ovoc
            ef6  = ef6_ovoc
            ef7  = ef7_ovoc
            ef8  = ef8_ovoc
            ef9  = ef9_ovoc
            ef10 = ef10_ovoc
            ef11 = ef11_ovoc
            ef12 = ef12_ovoc
            ef13 = ef13_ovoc
            ef14 = ef14_ovoc
            ef15 = ef15_ovoc
            anew  = anew_ovoc
            agro  = agro_ovoc
            amat  = amat_ovoc
            aold  = aold_ovoc
            caq   = caq_ovoc
            taq   = taq_ovoc
            dtaq  = dtaq_ovoc
            cht   = cht_ovoc
            tht   = tht_ovoc
            dtht  = dtht_ovoc
            clt   = clt_ovoc
            tlt   = tlt_ovoc
            dtlt  = dtlt_ovoc
            chw   = chw_ovoc
            thw   = thw_ovoc
            dthw  = dthw_ovoc
        end if

        if (lu_opt .eq. 0 .or. lu_opt .eq. 1) then !viirs or modis  lu types

! simple megan (table 3 in guenther et al., 2012) pft to viirs/modis vtype mapping
            if (vtype .eq. 1) then !viirs cat 1 evergreen needleleaf
                !--> average needleleaf evergreen temperate tree and needleleaf evergreen boreal tree

                ef = (ef1+ef2)/2.0_rk
                !set pft dependent a and b coefficients for cumulative root depth fraction (zeng 2001)
                !see table 2 for igpb classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 6.706_rk
                rootb = 2.175_rk

            else if (vtype .eq. 2) then !viirs/modis cat 2 evergreen broadleaf
                !--> average broadleaf evergreen tropical tree and broadleaf evergreen temperate tree

                ef = (ef4+ef5)/2.0_rk
                !set pft dependent a and b coefficients for cumulative root depth fraction (zeng 2001)
                !see table 2 for igpb classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 7.344_rk
                rootb = 1.303_rk

            else if (vtype .eq. 3) then !viirs/modis cat 3 deciduous needleaf
                !--> average needleleaf deciduous boreal tree

                ef = ef3
                !set pft dependent a and b coefficients for cumulative root depth fraction (zeng 2001)
                !see table 2 for igpb classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 7.066_rk
                rootb = 1.953_rk

            else if (vtype .eq. 4) then !viirs/modis cat 4 deciduous broadleaf
                !--> average broadleaf deciduous tropical tree, broadleaf deciduous temperate tree,
                ! and broadleaf deciduous boreal tree

                ef = (ef6+ef7+ef8)/3.0_rk
                !set pft dependent a and b coefficients for cumulative root depth fraction (zeng 2001)
                !see table 2 for igpb classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 5.990_rk
                rootb = 1.955_rk

            else if (vtype .eq. 5) then !viirs/modis cat 5 mixed forests
                !--> avearge of all above ef1-ef8 pfts.

                ef = (ef1+ef2+ef3+ef4+ef5+ef6+ef7+ef8)/8.0_rk
                !set pft dependent a and b coefficients for cumulative root depth fraction (zeng 2001)
                !see table 2 for igpb classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 4.453_rk
                rootb = 1.631_rk

            else if (vtype .ge. 6 .and. vtype .le. 7) then !viirs/modis cat 6-7 closed/open shrublands
                !--> avearge broadleaf evergreen temperate shrub, broadleaf deciduous temperate shrub,
                ! and broadleaf deciduous boreal shrub

                ef = (ef9+ef10+ef11)/3.0_rk
                !set pft dependent a and b coefficients for cumulative root depth fraction (zeng 2001)
                !see table 2 for igpb classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = (6.326_rk + 7.718_rk)/2.0_rk
                rootb = (1.567_rk + 1.262_rk)/2.0_rk

            else if (vtype .ge. 8 .and. vtype .le. 11) then !viirs/modis cat 8-10 savannas and grasslands
                !--> avearge arctic c3 grass, cool c3 grass, warm c4 grass), intent(out) ef;
```




<hr>



### variable emi\_ind 

```Fortran
integer, intent(in) emi_ind;
```




<hr>



### variable emissions 

```Fortran
integer, intent(in) emissions;
```




<hr>



### variable et 

```Fortran
integer, intent(in) et;
```




<hr>



### variable for 

```Fortran
integer, intent(in) for;
```




<hr>



### variable from 

```Fortran
integer, intent(in) from;
```




<hr>



### variable grid 

```Fortran
integer, intent(in) grid;
```




<hr>



### variable index 

```Fortran
integer, intent(in) index;
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



### variable mapped 

```Fortran
real(rk), intent(out) mapped;
```




<hr>



### variable massman 

```Fortran
integer, intent(in) massman;
```




<hr>



### variable model 

```Fortran
integer, intent(in) model;
```




<hr>



### variable out 

```Fortran
real(rk), intent(out) out;
```




<hr>



### variable to 

```Fortran
integer, intent(in) to;
```




<hr>



### variable type 

```Fortran
integer, intent(in) type;
```




<hr>



### variable vegetation 

```Fortran
integer, intent(in) vegetation;
```




<hr>



### variable vtype 

```Fortran
integer, intent(in) vtype;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_bioparm_mod.F90`

