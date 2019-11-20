This is a growing collection (work in progress) of Julia script to process ENZO simulated sky models (in FITS files) and produce long
lightcones which can be further processed using TRECS (Bonaldi et al. https://github.com/abonaldi/TRECS )to generate realistic population
of radio sources starting from galaxy positions.


generate_TRECS_public.jl generates galaxy distributions in TRECS-like format
 
generate_web_TRECS_share.jl generates similar data structures, but for the diffuse radio emission from cosmic shocks
generate_web_TRECS_v1.2.jl  is a bug-fix and Julia v.1.2 version adaptation by T. Hodgson

Both routines  are serial and typically become memory bound when large z>0.5-1.0 lightcones are generated, depending on the needed resolution, Field of View or sensitivity. 
