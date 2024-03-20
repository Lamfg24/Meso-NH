# Meso-NH

## Abstract of DELBEKE's thesis

Atmospheric aerosols are solid or liquid particles suspended in the air. Aerosols inter-
act directly with radiation in the atmosphere by absorbing or scattering radiations. On
the other hand, aerosols can be activated under proper thermodynamic conditions to form
cloud droplets. Therefore, aerosol properties can also determine cloud physical features
and thus indirectly affect atmospheric radiation profile. These direct and indirect effects
of aerosols on radiation have been demonstrated as a critical factor in regional and global
climate.

The focus of this thesis is the aerosol impact on the radiatively important low clouds
over southern West Africa. Africa is the largest source of mineral dust and biomass burn-
ing aerosols. Dust from the Saharan and Sahel regions and biomass burning (BB) aerosols
from areas in equatorial and central Africa often coexist in southern West Africa during
the monsoon season, coinciding with low-level stratiform clouds that cover a large area
frequently. Therefore, the interaction between aerosols in various types and the low clouds
could potentially impact on the regional weather and climate feature.

The first objective of this thesis is to understand the climatological features of aerosols
in Africa in particular of West Africa. For this purpose, the main characteristics of various
types of aerosols in the region including distribution and major transport paths have been
developed using different datasets and by applying machine learning methods.

Then the study is to focus on examining the impact of aerosols in different types on the
life cycle of low clouds using data obtained from the DACCIWA field campaign and the
meso-scale model Meso-NH. The 3 July 2016 case of DACCIWA campaign has been firstly
simulated in detail based on observed aerosol and cloud profiles with an ultra-fine resolu-
tion Large-Eddy simulation configuration. The model has successfully reproduced the life
cyle of observed cloud alongside evolution of all the major cloud macro and microphysical
features. The results of this reference run have advanced the understanding of the life cycle
as well as forcing particularly related to aerosol behind this cycle of the observed cloud
system. Simulations with different aerosol modules have also demonstrated the importance
of using interactive aerosol-chemistry model in studying the detailed impacts of aerosols
on the life cycle of studied low level cloud.

Lastly, several carefully designed sensitivity simulations have been conducted with dif-
ferent aerosol size distributions and chemical compositions to understand more precisely

the direct and indirect effect of different types of aerosols on the cloud life cycle. It is found
that aerosol size distribution has a clear impact on cloud droplet number and mean droplet
size and affects cloud reflectivity. However, the variation in cloud reflectivity induced by
different aerosol profiles is found to be not the only factor in determining the incoming
solar radiation at the ground and thus for the cloud life cycle after the sunrise. Instead,
the difference in cloud fraction brought by dry-air entrainment from above and thus the
speed of consequent evaporation – also influenced by aerosol concentration, plays a more
critical role. Thus, cloud with a higher number concentration and smaller size of cloud
droplets are found to evaporate more easily and hence impose a lower cloud fraction. In
addition, the sensitivity runs, including versus excluding aerosol direct radiative effects,
have also demonstrated the impacts specifically of solar absorption by black carbon on the
cloud life cycle. The semi-direct effect resulting from an excessive atmospheric heating by
black carbon in the modeled cases is found to lower the cloud top as well as the liquid
water path reducing surface incoming solar radiation and dry entrainment and increasing
the cloud fraction.

## Simulation of 3 July 2016

### Introduction

![](/Pictures/LLCS_3_july_meas.png)

### Meso-NH

Meso-NH is a French non hydrostatic atmospheric research model used in research studies
ranging from synoptic to turbulent scales developed by the Centre National de Recherches
Météorologiques (CNRM) and the Laboratoire d’Aérologie (LAERO) (Lafore et al., 1998;
Lac et al., 2018). The model can simulate real (realistic meteorological situations) or ide-
alized (simulations starting from highly idealized conditions, such as an atmosphere in
geostrophic equilibrium, with simple vertical and horizontal structure, and simple orog-
raphy, or no orography at all) cases which can be used from 1D to 3D. The model uses
three prognostic variables like three velocity terms (u,v,w), the potential temperature and
mixing ratios of some species like vapor, cloud droplets, raindroplets, ice crystals, snow,
graupel and hail. Turbulent kinetic energy (TKE) is also used and other reactive and
passive scalars like aerosols and hydrometeors concentration or chemical species. Meso-NH
conserves the grid nesting technique in order to simulate scale interactions. The software
is maintained and improved by computer and research team from CNRM and LAERO.
Code and scripts are written in Fortan 90 and shell respectively and use makefiles. Paral-
lelization has been applied on much of the soft using the Message Passing Interface (MPI).
Meso-NH is organized in three blocks corresponding to three steps. The first step is the
preparation of the simulation when the user has to compose initial fields for idealized or
real conditions and also the preparation of the nesting application. The second step is
the running of the simulation starting with an initialization step and then the simulation
is integrated. A post-processing step is able to compute supplementary diagnostic fields.
Meso-NH needs atmospheric initial and boundary conditions. These so called Large-Scale
(LS) fields are used to initialize prognostic variables and force them at lateral boundaries
with time-evolving fields. For real cases, AROME, ARPEGE, ECMWF and also GFS suits
data suits well. For ideal case, initial vertical profiles derived usually from radiosounding
data are then interpolated horizontally and vertically onto the Meso-NH grid. The forcing
methods implemented are from geostrophic winds to large-scale thermodynamical tenden-
cies.

### Using only LIMA

### Using ORILAM 

#### REF 

#### POL

#### CLEAN
