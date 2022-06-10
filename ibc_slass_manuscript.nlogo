;; IBC-mammal: Dynamic, allometric, individual-based, spatial community model for cental-place forager mammals in heterogenous landscapes
;; copyright Leonna Szangolies, Florian Jeltsch, Marie-Sophie Rohwäder, University of Potsdam

extensions [profiler]                               ;;extension for tracking model run and debugging

globals [s1mass s1massdev  s2mass s2massdev s3mass s3massdev s4mass s4massdev s5mass s5massdev s6mass s6massdev
         s7mass s7massdev s8mass s8massdev s9mass s9massdev s0mass s0massdev s0shelt s1shelt s2shelt s3shelt s4shelt
         s5shelt s6shelt s7shelt s8shelt s9shelt s0food s1food s2food s3food s4food s5food s6food s7food
         s8food s9food s0fortype s1fortype s2fortype s3fortype s4fortype s5fortype s6fortype s7fortype s8fortype
         s9fortype feed_matrix feed_struct r_search feed shelter_need maxrad minfeed clump cover spec_color
         movecost foodshare feed2 day fail fail2 fail3 repro foodmort sp0 sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8 sp9 maxmass
         shannon actcov mass_order repetitions n_immi hr_try_juv max_nat_disp mortbold focal_spec spec_num
         stepstones density_dep small_habitat dynamic total_cover_constant cover_large cover_small disturb bold_die]
patches-own [habitat patchfeed patchquali spec_id spec_list eaten sss0 sss1 sss2 sss3 sss4 sss5 sss6 sss7 sss8 sss9 save patch_spec_num]
turtles-own [mass species shelter food_pref for_type maxhr hrsize lococost feedrate age average_age preg
             pregcount age_first_repro gest_period lact_period sex order stime hunger core patches_in_radius ss_use forage_large forage_small]

to setup
  clear-all
  reset-ticks
  set-parameters
  ifelse total_cover_constant
  [set cover_large cover * (1 - perc_small) set cover_small (cover * perc_small)]
  [set cover_large cover set cover_small (size_small * (perc_small * cover))]
  setup-patches-large
  setup-patches-small
  setup-turtles
  ask patches
  [
    set pcolor scale-color black habitat 1 0
  ]
  find-hr
end

;---------------------------------------------

to set-parameters                                      ;; characterize species und basic parameters
  set maxmass 0.1                                      ;; upper limit of mean body mass
  set s0mass 0.1 * maxmass                             ;; mean body mass of species 0 in kg - normal distributed
  set s0massdev 0.2 * s0mass                           ;; std dev of body mass species 0 - normal distribution
  set s1mass 0.2 * maxmass                             ;; mean body mass of species 0 - normal distributed
  set s1massdev 0.2 * s1mass                           ;; std dev of body mass species 0 - normal distribution
  set s2mass 0.3 * maxmass                             ;; mean body mass of species 0 - normal distributed
  set s2massdev 0.2 * s2mass                           ;; std dev of body mass species 0 - normal distribution
  set s3mass 0.4 * maxmass                             ;; mean body mass of species 0 - normal distributed
  set s3massdev 0.2 * s3mass                           ;; std dev of body mass species 0 - normal distribution
  set s4mass 0.5 * maxmass                             ;; mean body mass of species 0 - normal distributed
  set s4massdev 0.2 * s4mass                           ;; std dev of body mass species 0 - normal distribution
  set s5mass 0.6 * maxmass                             ;; mean body mass of species 0 - normal distributed
  set s5massdev 0.2 * s5mass                           ;; std dev of body mass species 0 - normal distribution
  set s6mass 0.7 * maxmass                             ;; mean body mass of species 0 - normal distributed
  set s6massdev 0.2 * s6mass                           ;; std dev of body mass species 0 - normal distribution
  set s7mass 0.8 * maxmass                             ;; mean body mass of species 0 - normal distributed
  set s7massdev 0.2 * s7mass                           ;; std dev of body mass species 0 - normal distribution
  set s8mass 0.9 * maxmass                             ;; mean body mass of species 0 - normal distributed
  set s8massdev 0.2 * s8mass                           ;; std dev of body mass species 0 - normal distribution
  set s9mass maxmass                                   ;; mean body mass of species 0 - normal distributed
  set s9massdev 0.2 * s9mass                           ;; std dev of body mass species 0 - normal distribution

  set s0shelt 1                                     ;; index 0..1 of shelter use/need
  set s1shelt 1                                     ;; index 0..1 of shelter use/need
  set s2shelt 1                                     ;; index 0..1 of shelter use/need
  set s3shelt 1                                     ;; index 0..1 of shelter use/need
  set s4shelt 1                                     ;; index 0..1 of shelter use/need
  set s5shelt 1                                     ;; index 0..1 of shelter use/need
  set s6shelt 1                                     ;; index 0..1 of shelter use/need
  set s7shelt 1                                     ;; index 0..1 of shelter use/need
  set s8shelt 1                                     ;; index 0..1 of shelter use/need
  set s9shelt 1                                     ;; index 0..1 of shelter use/need

  set s0food "omni"                                    ;; food preference; not distinguished yet
  set s1food "omni"                                    ;; food preference
  set s2food "omni"                                    ;; food preference
  set s3food "omni"                                    ;; food preference
  set s4food "omni"                                    ;; food preference
  set s5food "omni"                                    ;; food preference
  set s6food "omni"                                    ;; food preference
  set s7food "omni"                                    ;; food preference
  set s8food "omni"                                    ;; food preference
  set s9food "omni"                                    ;; food preference

  set s0fortype "central"                              ;; type of forage movement; not distinguished yet
  set s1fortype "central"                              ;; type of forage movement
  set s2fortype "central"                              ;; type of forage movement
  set s3fortype "central"                              ;; type of forage movement
  set s4fortype "central"                              ;; type of forage movement
  set s5fortype "central"                              ;; type of forage movement
  set s6fortype "central"                              ;; type of forage movement
  set s7fortype "central"                              ;; type of forage movement
  set s8fortype "central"                              ;; type of forage movement
  set s9fortype "central"                              ;; type of forage movement

  set actcov cover                                     ;; actual cover after changes
  set fail 0                                           ;; counting failing attempts by offspring to establish homerange
  set fail2 0                                          ;; counting cases with food shortage in old home range
  set fail3 0                                          ;; counting death events caused by food shortage in attempts to re-establish homerange next day
  set repro 0                                          ;; counting overall annual reproduction
  set foodmort 0                                       ;; mortality caused by starving

;;parameters to be changed
  set cover total_cover * 100 * 100                    ;; xx % of shrub cover * nr of patches; xx * 100 * 100
  set feed_matrix 0.0                                  ;; no food prod in matrix
  set feed_struct 6.85 ;1.7                            ;; 10% of 68.5g dry mass/ grid cell*day Buchmann et al. 2011
  set clump 0.99                                       ;; clumping index default 0.999=low fragmented, szenarien 0.99=med fragmented und 0.9=highly fragmented
  set mortbold 0.333 / 10000                           ;;daily mortality of bold indivdiuals: 74.1 => 20%, 3.35 => 10% , 1.67 => 5%, 0.333 => 1% per month ;;richtig??!!
  set mass_order 1;0.75                                ;; % of turtles that are ordered acording to body mass, remaining turtles are ordered randomly
  set repetitions 20                                   ;; nr of repetitions with new landscape but same landscape features and parameters
  set n_immi 0                                         ;;nr of random immigrants per timestep
  set hr_try_juv 10                                    ;;nr of attempts of juveniles to find own homerange
  set focal_spec 10                                    ;; fcoal species of which daily homeranges are shown (cumulative)
  set spec_num 10                                      ;; initially 10 species present
  set bold_die [0 0 0 0 0 0 0 0 0 0]                   ;; individuals dying because of boldness
  set stepstones False                                 ;; scneario definitions
  set small_habitat False ;True
  set density_dep False
  set dynamic False
  set total_cover_constant True ;False
  if scenario = "base" [set small_habitat True]
  if scenario = "stepping_stone" [set stepstones True]
  if scenario = "increasing_base" [set total_cover_constant False set small_habitat True]
  if scenario = "increasing_reduced" [set total_cover_constant False]
end
;---------------------------------------------

to setup-patches-large                                  ;; creates 4 large habitat patches
  let cov 0
  ask patches   [
      set patchquali 0
      set habitat 0
      set spec_id -1
      set spec_list [ ]
      set eaten 0
   ]
  ask n-of 4 patches with [habitat = 0] [set habitat 1 set patchquali 1 set cov cov + 1]

  while [cov < cover_large]
  [
   ask patches with [habitat = 1]
   [
    if ((any? neighbors with [habitat = 0]) and (cov < cover_large))
      [ask one-of neighbors with [habitat = 0] [set habitat 1 set patchquali 1 set cov cov + 1]]
   ]
  ]
end

to setup-patches-small                                   ;; creates the defined amount of small habitat patches
  if dynamic
  [
    ask patches with [patchquali = 2]
    [
      set habitat 0
      set patchquali 0
    ]
  ]

  ask n-of ceiling (cover_small / size_small) patches with [habitat = 0 and (count neighbors with [habitat = 0]) = 8]
  [
    ask (n-of size_small patches in-radius ceiling (size_small / 8 )) with [habitat = 0] [set habitat 1 set patchquali 2]
  ]
  ask patches[set save (sum [habitat] of neighbors > 6)]

  resource-update
end
;..............................................

to resource-update                                    ;;daily update of food resources, random variation normal distributed
  ask patches
   [
      ifelse habitat = 0 [set patchfeed feed_matrix] [set patchfeed random-normal feed_struct 0.7]
   ]
end
;---------------------------------------------

to setup-turtles                                       ;; create xx initial individuals, location at 0,0; identifies species, mass, ...
  crt 1000 [set age random 100 ]
                                                       ;; scenario 10 species, species identity 0 .. 9
                                                       ;; allometric frequency distribution with mass^-1.5, see Buchmann et al. 2012 Ecography
  let s0m s0mass ^ (-1.5)
  let s1m s1mass ^ (-1.5)
  let s2m s2mass ^ (-1.5)
  let s3m s3mass ^ (-1.5)
  let s4m s4mass ^ (-1.5)
  let s5m s5mass ^ (-1.5)
  let s6m s6mass ^ (-1.5)
  let s7m s7mass ^ (-1.5)
  let s8m s8mass ^ (-1.5)
  let s9m s9mass ^ (-1.5)
  let mass_sum s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m + s8m + s9m

  ask turtles
  [
   let zuf random-float 1.0

   ifelse zuf < s0m / mass_sum
   [ set species 0 ]
    [ ifelse zuf < (s0m + s1m) / mass_sum [set species 1 ]
      [ ifelse zuf < (s0m + s1m + s2m) / mass_sum [set species 2 ]
        [ ifelse zuf < (s0m + s1m + s2m + s3m) / mass_sum [set species 3 ]
          [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m) / mass_sum [set species 4 ]
            [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m) / mass_sum [set species 5 ]
              [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m) / mass_sum [set species 6 ]
                [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m) / mass_sum [set species 7 ]
                  [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m + s8m) / mass_sum [set species 8 ]
                    [set species 9] ;;last else
                  ];;ifelse8
                ];;ifelse7
              ];;ifelse6
            ];;ifelse5
          ];; ifelse4
        ] ;; ifelse3
      ] ;; ifelse2
    ] ;;ifelse1

    specify_turtle
  ] ;;ask turtles
end
;.............................................

to specify_turtle                                               ;;characterize and parameterize individuals
  set hrsize 0
  let shy 0
  let zuf1 random-float 1.0
  set color (species * 10 + 5)
  ifelse (zuf1 < bold_prob) [set shy 0][set shy 1]              ;; shy=0 means bold individuals

    if (species = 0) [set mass random-normal s0mass s0massdev   ;; body mass normal distributed
      set shelter s0shelt * shy
      set food_pref s0food
      set for_type s0fortype]
    if (species = 1) [set mass random-normal s1mass s1massdev
      set shelter s1shelt * shy
      set food_pref s1food
      set for_type s1fortype]
    if (species = 2) [set mass random-normal s2mass s2massdev
      set shelter s2shelt * shy
      set food_pref s2food
      set for_type s2fortype]
    if (species = 3) [set mass random-normal s3mass s3massdev
      set shelter s3shelt * shy
      set food_pref s3food
      set for_type s3fortype]
    if (species = 4) [set mass random-normal s4mass s4massdev
      set shelter s4shelt * shy
      set food_pref s4food
      set for_type s4fortype]
    if (species = 5) [set mass random-normal s5mass s5massdev
      set shelter s5shelt * shy
      set food_pref s5food
      set for_type s5fortype]
    if (species = 6) [set mass random-normal s6mass s6massdev
      set shelter s6shelt * shy
      set food_pref s6food
      set for_type s6fortype]
    if (species = 7) [set mass random-normal s7mass s7massdev
      set shelter s7shelt * shy
      set food_pref s7food
      set for_type s7fortype]
    if (species = 8) [set mass random-normal s8mass s8massdev
      set shelter s8shelt * shy
      set food_pref s8food
      set for_type s8fortype]
    if (species = 9) [set mass random-normal s9mass s9massdev
      set shelter s9shelt * shy
      set food_pref s9food
      set for_type s9fortype]

    if mass < 0 [set mass 0 die]

    set preg 0
    set pregcount 0
    set average_age 1766.53 * (mass ^ 0.21)                        ;; average lifespan: allometric formula of mammals after Hamilton et al 2011 [g, days]
    set age_first_repro 293.17 * (mass ^ 0.27)                     ;; age of first reproduction: allometric fromula after Hamilton et al 2011
    set gest_period 64.14 * (mass ^ 0.24)                          ;; gestation period allometric after Hamilton et al 2011
    set lact_period  57.16 * (mass ^ 0.22)                         ;; lactation period allometric after Hamilton et al 2011
    set stime 9.3 * mass ^ 0.44                                    ;; allometric survival time when fasting after Lindstedt and Boyce 1985
    set hunger 0
    set ss_use False
    set forage_small 0
    set forage_large 0

    let zuf random-float 1.0
    ifelse zuf < 0.5 [set sex "male"] [set sex "female"]

    calc-maxhr
    calc-lococost
    calc-feedrate
end
;---------------------------------------------

to calc-maxhr                                                      ;; calculate max home range after Kelt & Van Vuren 2001 in ha, see Buchmann et al 2011
  let maxhr1 56.23 * (mass ^ 0.91)                                 ;; max hr for herbivores and omnivores, larger one is used
  let maxhr2 47.863 * (mass ^ 1.18)
  ifelse maxhr2 > maxhr1 [set maxhr maxhr2] [set maxhr maxhr1]
  set maxhr maxhr * 10000                                          ;; in m2
  set maxhr sqrt (maxhr / pi)                                      ;; radius of maxhr in m
  set maxhr maxhr / 10                                             ;; radius in patch length (=10m)
end
;---------------------------------------------

to calc-lococost                                                   ;; calculate movement costs, specific for food type
  set lococost 10.7 * (mass ^ 0.68)                                ;; costs for mammals in J/m; mass in kg after Calder 1996
  set lococost lococost / 10000                                    ;; costs in g dry biomass/m; after Nagy'99 p.263 Buchmann 2011
  set lococost lococost * 10                                       ;; lococost in patchlength (= 10m)
end
;---------------------------------------------

to calc-feedrate                                                   ;; caluclate daily feeding rate, specific for food type
 set feedrate 55.1 * (mass ^ 0.74)                                 ;; daily feeding rate of mammals in dry g/day, mass in kg, after Nagy 2001
end
;----------------------------------------------

to change-order                                                    ;; order x% of turtles according to body mass, remaining turtles are ordered randomly
  ask turtles
  [
     let zuf random-float 1.0
     ifelse zuf < mass_order [set order mass] [set order random-float maxmass]
  ]
end
;-----------------------------------------------

to one-step                                                        ;; one time step

  if ( debug = 1 ) and ( ticks > 360 ) [
    profiler:start
  ]

  set day ticks mod 360
  set fail 0
  set fail2 0
  set fail3 0
  set repro 0
  ask patches
  [
    set pcolor scale-color green habitat 1 0
    set spec_list [ ]
    set eaten 0
  ]
  immigrate                                                        ;; function: immigration
  ask turtles
  [
   set age age + 1
  ]

  ifelse dynamic [setup-patches-small] [resource-update]
  ;resource-update                                                 ;; function: update resources

  check-hr                                                         ;; function: check if old homerange is still O.K., adapt if possible
  change-order                                                     ;; function: order mass_order% of turtles according to body mass
  foreach sort-on [ (- order)] turtles                             ;; offspring of heaviest turtles check hr first for mass_oder% of turtles
  [ ?1 -> ask ?1                                                   ;; ask turtles
   [
   if sex = "female"
    [
     ifelse preg = 1
      [set pregcount pregcount + 1]
      [if age > age_first_repro [ set preg 1 ] ]                   ;; getting pregnant deterministically => age allometric
     if pregcount >  gest_period + lact_period                     ;; days of gestation/pregnancy PLUS days of lactation==> allometric after Hamilton et al 2011
      [ offspring set  pregcount 0   set preg 0 ]
    ] ;; end if sex = female
   mort                                                            ;; mortality function
   ] ;;end ?1 -> ask ?1
   ] ;; end foreach sort
  patch-use
  output

  if ( debug = 1 ) and ( ticks > 360 ) [
    profiler:stop
    print profiler:report
    profiler:reset
  ]
end
;-----------------------------------------------

to go                                                              ;; to go function, multiple time steps
  one-step
  tick
end
;---------------------------------------------

to offspring                                                         ;; offspring, inherits everything from mother
  let young round (2.24 * (mass ^ (-0.13)))                          ;; allomteric litter size after Hamilton et al. 2011: 5.5 * (mass ^ (-0.13))
  set repro repro + young                                            ;; to count overall annual reproduction
  hatch young
    [
     set age 0
     let zuf random-float 1.0
     ifelse zuf < 0.5 [set sex "male"] [set sex "female"]
     set shape "bug"
     set color (species * 10 + 5)
     specify_turtle
     find-hr-offspring                                               ;; note: to offspring takes place after gestation AND lactation!
    ]
end
;.............................................

to mort                                                              ;; mortality function related to life span, normally distributed and boldness
  if shelter = 0 and random-float 1 < mortbold [set bold_die replace-item species bold_die ((item species bold_die) + 1) die]    ;; additional mortality for bold indivdiuals (bold=> shelter= 0)
  if age > random-normal average_age ( 0.1 * average_age ) [ die ]   ;; average life span is allometric, see above
  if density_dep [
    let this_spec species
    let spec_ind count turtles with [species = this_spec]
    let capac 3000 * ( mass ^ (-0.75))
    if random-float 1.0 < (spec_ind / capac) [die]
  ]
end
;.............................................

to immigrate                                                         ;; immigratin function
 crt n_immi                                                          ;; n_immi immigrants per day, random species, random age < 100 days, random location
  [
    set species random 10
    set age random 100
    specify_turtle
    move-to one-of patches with [habitat > 0]                        ;; random search for hr core cell ]
    set r_search 0
    set patches_in_radius []
    set core patchquali
    set maxrad maxhr
    while [r_search < maxrad]
     [
       set r_search r_search + 1
       set patches_in_radius lput (patches in-radius (r_search) with [patchfeed > 0 and eaten = 0]) patches_in_radius
     ]
  ]
end
;.............................................

to find-hr                                                         ;; find suitable homerange for initial distribution
  ask turtles
  [
  set maxrad maxhr
  set minfeed feedrate
  set shelter_need shelter
  set spec_color (species * 10 + 5)
  set movecost lococost
  set foodshare (mass / 0.001) ^ (-0.25)                           ;; allometric share of available food see Buchmann 2011
  let spec species
  let success 0
  let try 0

  while [success = 0 and try < 100]                                ;; 100 attempts for each mammal of inital community to find suitable hr
    [
      set try try + 1
      set r_search 0
      set feed 0
      set patches_in_radius []
      ifelse small_habitat [move-to one-of patches with [habitat > 0]][move-to one-of patches with [patchquali = 1]]                    ;; random search for potential hr core cell
      set core patchquali
      set feed patchfeed * foodshare                               ;; calculate food of core cell
      set eaten 1
      if feed >= minfeed                                           ;; if enough food in core cell - end search
       [
         set success 1
         set eaten 1                                               ;; 1 indicates patch as hr-patch forlater food reduction
         set spec_id spec                                          ;; identify patch as part of hr of species
         set spec_list lput spec_id spec_list                      ;; add species to patch-specific species list
       ]
      while [r_search < maxrad]
     [
       set r_search r_search + 1
       set patches_in_radius lput (patches in-radius (r_search) with [patchfeed > 0 and eaten = 0]) patches_in_radius

      if feed < minfeed;and feed < minfeed]                        ;; if not enough food in core cell - search in neighborhood
      [
       ask item (r_search - 1) patches_in_radius
        [
           set eaten 2                                             ;; 2 indicates potential use as hr-patch for later food reduction
           ifelse save                                             ;; edge effects: patches 'in the open' are less frequently visited or shorter time => reduced food intake
            [set feed feed + patchfeed * foodshare - 2 * movecost * distance myself]
            [set feed feed + (patchfeed * foodshare - 2 * movecost * distance myself) * (1 - shelter_need)]
        ]  ;; end ask-patches in r
       ]  ;; end while r-search
      ]
     ifelse feed >= minfeed
       [
          ;let patchbefore sum [patchfeed] of patches
          set success 1
          set patchfeed patchfeed * (1 - foodshare)
          set eaten 0
          ask patches in-radius (r_search) with [eaten = 2]
           [
            ;if eaten = 2 [set eaten 1]                             ;; only used patches
            ;if eaten = 1
             ;[
              ifelse save
                [set patchfeed patchfeed * (1 - foodshare)]        ;; reduce remaining food in patch
                [set patchfeed patchfeed * (1 - foodshare * (1 - shelter_need))]
              set spec_id spec                                     ;; identify patch as part of hr of species
              set spec_list lput spec_id spec_list                 ;; add species to patch-specific species list
              if spec = focal_spec [set pcolor spec_color]         ;; show only hr of focal species
              set eaten 0
             ;] ;; end if eaten = 1
           ] ;; end ask parches in-radius
          let patchafter sum [patchfeed] of patches
          ;print feed
          ;print patchbefore - patchafter
       ] ;; end ifelse feed >=minfeed cond1
       [
          ask patches in-radius (r_search) [set eaten 0]            ;; set back to not-eaten
       ] ;; end ifelse feed >=minfeed cond2

    ] ;; end while success = 0

    if success = 1 [set hrsize 2 * r_search]
    if success = 0 [die]

  ] ;; end ask turtles
end
;----------------------------------------------

to check-hr                                                              ;; check if existing homerange is still sufficiant => food is reduced!
 ask turtles
  [
     let zuf random-float 1.0
     ifelse zuf < 0.0 [set order mass] [set order random-float maxmass]  ;; optional: arrange order so that x% of heavier animals can check hr first, default x=0.0
  ]
  foreach sort-on [ (- order)] turtles
  [ ?1 -> ask ?1
   [
    if age > 0                                                           ;; only for non-offspring
    [
     let spec species
     set maxrad maxhr
     let act_feedrate feedrate
     set spec_color (species * 10 + 5)
     if (preg = 1)                                                       ;; gestating or lactating females
      [
       ifelse (pregcount < gest_period)                                  ;; gestating
        [
         set act_feedrate feedrate * 1.2                                 ;; increased feeding rate for gestating females
        ]
      [
        if (pregcount < gest_period + lact_period)
         [
          set act_feedrate feedrate * 1.8                                ;; increased feeing rate for lactating females
         ]
      ] ;;end ifelse pregcount
     ] ;; end preg=1

    set minfeed act_feedrate
    set shelter_need shelter
    set movecost lococost
    set foodshare (mass / 0.001) ^ (-0.25)                               ;; allometric share of available food see Buchmann 2011
    let success 0
    set r_search 0
    set feed 0
    set feed patchfeed * foodshare                                       ;; calculate food of core cell
    set patchfeed patchfeed * (1 - foodshare)                            ;; reduce remaining food in core patch
    set eaten 1
    if patchfeed < 0 [set patchfeed 0]
    if feed >= minfeed                                                   ;; if enough food in core cell - end
       [
         set success 1
         set hunger 0
         set pcolor spec_color
         set spec_id spec
         set spec_list lput spec_id spec_list
       ]
    if dynamic
    [
    set patches_in_radius []
    while [r_search < maxrad]
     [
       set r_search r_search + 1
       set patches_in_radius lput (patches in-radius (r_search) with [patchfeed > 0 and eaten = 0]) patches_in_radius
     ]
    set r_search 0
    ]
    while [r_search < maxrad and feed < minfeed]                         ;; if not enough food in core cell - search in neighborhood
      [
       set r_search r_search + 1
       ask item (r_search - 1) patches_in_radius
        [
         if spec = focal_spec [set pcolor spec_color]                    ;; show only hr of focal species
         set eaten 2
         set spec_id spec
         set spec_list lput spec_id spec_list
         ifelse save                                                     ;; edge effect: patches 'in the open' are less frequently or for shorter time visited
           [
             set feed feed + patchfeed * foodshare - 2 * movecost * distance myself
             set patchfeed patchfeed * (1 - foodshare)
           ]
           [
             set feed feed + (patchfeed * foodshare - 2 * movecost * distance myself) * (1 - shelter_need)
             set patchfeed patchfeed * (1 - foodshare * (1 - shelter_need))  ;;reduce food in patch
           ]
         if patchfeed < 0 [set patchfeed 0]
        ] ;; end ask patches in-radius
       ] ;; end while r-search
        ask patches in-radius r_search [set eaten 0]

        if feed >= minfeed
          [
           set success 1
           set hunger 0
          ]

        ifelse success = 1
          [
           set hrsize 2 * r_search                                           ;; diameter of hr
          ]
          [
           set hrsize 2 * maxhr                                              ;; max. diameter of hr
          ]
         if success = 0                                                      ;; increased mortality with food shortage
          [
           set hunger hunger + 1                                             ;; add one day of food shortage
           set foodmort (feed / minfeed)                                     ;; mortality depends on relative food shortage
           if random-float 1.0 < (1 - foodmort) * (hunger / stime)           ;; actual food shortage * (days of starving/allometric fasting endurance)
            [
             set fail3 fail3 + 1
             die                                                             ;; individual dies
            ]
          ]  ;; end if success=0
      ]   ;; end if age >0
     ]   ;;end ask turtles
  ]   ;; end foreach
end
;----------------------------------------------

to find-hr-offspring                                                         ;; offspring searches for new homerange
  set maxrad maxhr
  set minfeed feedrate;;
  set shelter_need shelter
  set spec_color (species * 10 + 5)
  set movecost lococost
  set foodshare (mass / 0.001) ^ (-0.25)                                     ;;allometric share of available food see Buchmann 2011
  set max_nat_disp 3.31 * (mass ^ 0.65)                                      ;; allometric maximum natal dispersal after Sutherland et al (2000) in km
  let spec species
  let success 0
  let try 0
  let xx xcor
  let yy ycor

  while [success = 0 and try < hr_try_juv]                                   ;; hr_try_juv attempts for each mammal to find suitable hr
    [
      setxy xx yy
      set try try + 1
      set patches_in_radius []
      ifelse small_habitat
          [move-to one-of patches in-radius (100 * max_nat_disp) with [habitat > 0]]     ;; random search for hr core cell in maximum natal dispersal distance in grid cell length = 10m
          [move-to one-of patches in-radius (100 * max_nat_disp) with [patchquali = 1]]
      set core patchquali
      set r_search 0
      set feed patchfeed * foodshare                                         ;; take food of core cell
      set eaten 1
      if feed >= minfeed                                                     ;; if enough food in core cell - end
       [
         set success 1
         set pcolor spec_color
         set eaten 1                                                         ;; 1 indicates patch as hr-patch for later food reduction
         set spec_id spec                                                    ;; identify patch as part of hr of species
         set spec_list lput spec_id spec_list                                ;; add species to patch-specific species list
       ]
       while [r_search < maxrad]
     [
       set r_search r_search + 1
       set patches_in_radius lput (patches in-radius (r_search) with [patchfeed > 0 and eaten = 0]) patches_in_radius

      if feed < minfeed;and feed < minfeed]                                  ;; if not enough food in core cell - search in neighborhood
      [
       ask item (r_search - 1) patches_in_radius
        [
          set eaten 2                                                        ;;  2 indicates potential use as hr-patch for later food reduction
          ifelse save                                                        ;; edge effect: patches 'in the open' are less frequently visited or fot shorter time
            [set feed feed + patchfeed * foodshare - 2 * movecost * distance myself]
            [set feed feed + (patchfeed * foodshare - 2 * movecost * distance myself) * (1 - shelter_need)]
        ] ;; end ask patches in radius
       ] ;; end while r-search
      ]

      ifelse feed >= minfeed
       [
         set success 1
         set patchfeed patchfeed * (1 - foodshare)
         set eaten 0
         ask patches in-radius (r_search) with [eaten = 2]
           [
              ifelse save
                [set patchfeed patchfeed * (1 - foodshare)]                  ;; reduce remaining food in patch
                [set patchfeed patchfeed * (1 - foodshare * (1 - shelter_need))]
              set spec_id spec                                               ;; identify patch as part of hr of species
              set spec_list lput spec_id spec_list                           ;; add species to patch-specific species list
              if spec = focal_spec [set pcolor spec_color]                   ;; show only hr of focal species
              set eaten 0
             ;; end if eaten = 1
           ] ;; end ask patches in-radius
         ] ;; end ifelse feed >=minfeed cond 1
         [
          ask patches in-radius (r_search) [set eaten 0]                     ;; set back to not-eaten
         ] ;;end ifelse feed >=minfeed cond 2
    ] ;; end while success = 0

  if stepstones[                                                             ;; if stepstone scenario: go to random small patch and search for home range again from there
  if success = 0 and any? patches in-radius (100 * max_nat_disp) with [patchquali = 2][
    move-to one-of patches in-radius (100 * max_nat_disp) with [patchquali = 2]
    set try 0
    set xx xcor
    set yy ycor

  while [success = 0 and try < hr_try_juv]                                   ;; hr_try_juv attempts for each mammal to find suitable hr
    [
      setxy xx yy
      set try try + 1
      set patches_in_radius []
      ifelse small_habitat
          [move-to one-of patches in-radius (100 * max_nat_disp) with [habitat > 0]]     ;; random search for hr core cell in maximum natal dispersal distance in grid cell length = 10m
          [move-to one-of patches in-radius (100 * max_nat_disp) with [patchquali = 1]]  ;; random search for hr core cell in maximum natal dispersal distance in grid cell length = 10m
      set core patchquali
      set r_search 0
      set feed patchfeed * foodshare                                         ;; take food of core cell
      set eaten 1
      if feed >= minfeed                                                     ;; if enough food in core cell - end
       [
         set success 1
         set pcolor spec_color
         set eaten 1                                                         ;; 1 indicates patch as hr-patch for later food reduction
         set spec_id spec                                                    ;; identify patch as part of hr of species
         set spec_list lput spec_id spec_list                                ;; add species to patch-specific species list
       ]
       while [r_search < maxrad]
     [
       set r_search r_search + 1
       set patches_in_radius lput (patches in-radius (r_search) with [patchfeed > 0 and eaten = 0]) patches_in_radius

      if feed < minfeed;and feed < minfeed]                                  ;; if not enough food in core cell - search in neighborhood
      [
       ask item (r_search - 1) patches_in_radius
        [
          set eaten 2                                                        ;;  2 indicates potential use as hr-patch for later food reduction
          ifelse save                                                        ;; edge effect: patches 'in the open' are less frequently visited or fot shorter time
            [set feed feed + patchfeed * foodshare - 2 * movecost * distance myself]
            [set feed feed + (patchfeed * foodshare - 2 * movecost * distance myself) * (1 - shelter_need)]
        ] ;; end ask patches in radius
       ] ;; end while r-search
      ]

      ifelse feed >= minfeed
       [
         set ss_use True
         set success 1
         set patchfeed patchfeed * (1 - foodshare)
         set eaten 0
         ask patches in-radius (r_search) with [eaten = 2]
           [
              ifelse save
                [set patchfeed patchfeed * (1 - foodshare)]                  ;; reduce remaining food in patch
                [set patchfeed patchfeed * (1 - foodshare * (1 - shelter_need))]
              set spec_id spec                                               ;; identify patch as part of hr of species
              set spec_list lput spec_id spec_list                           ;; add species to patch-specific species list
              if spec = focal_spec [set pcolor spec_color]                   ;; show only hr of focal species
              set eaten 0
             ;; end if eaten = 1
           ] ;; end ask patches in-radius
         ] ;; end ifelse feed >=minfeed cond 1
         [
          ask patches in-radius (r_search) [set eaten 0]                     ;; set back to not-eaten
         ] ;;end ifelse feed >=minfeed cond 2
    ]
   ]
  ]
   if success = 1 [set hrsize 2 * r_search]
   if success = 0 [set fail fail + 1]
   if success = 0 [ die ]                                                    ;; unsuccessful turtles die or emigrate from area
end
;----------------------------------------------

to patch-use                                                                 ;; document patch types in home range
  ask turtles
  [
    set forage_large count (patches in-radius (hrsize / 2) with [patchquali = 1])
    set forage_small count (patches in-radius (hrsize / 2) with [patchquali = 2])
  ]
end
;----------------------------------------------

to output                                                                    ;; document shannon diversity
  let specs [species] of turtles
  let specs_unique remove-duplicates specs
  set spec_num length specs_unique
  set shannon 0
  let spec 0
  while [spec < 9]
  [
    let p count turtles with [species = spec] / count turtles
    if p > 0 [set shannon shannon + (p * ln(p))]
    set spec spec + 1
  ]
  set shannon (- shannon)
  ask patches
  [
    let patch_unique remove-duplicates spec_list
    set patch_spec_num length patch_unique
  ]

end
;----------------------------------------------

to analyze                                                                     ;; histograms
  set-current-plot "homerange distribution" histogram [hrsize] of turtles with [hrsize > 0]
  set-current-plot "species richness" histogram [species] of turtles with [hrsize > 0]
  set-current-plot "body mass distribution" histogram [mass * 1000] of turtles with [hrsize > 0]
end
@#$#@#$#@
GRAPHICS-WINDOW
213
15
689
492
-1
-1
4.634
1
6
1
1
1
0
1
1
1
0
100
0
100
1
1
1
ticks
30.0

BUTTON
29
17
109
51
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
29
55
111
89
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
11
87
211
237
homerange distribution
hr-size/ diameter in grid cells
number
0.0
20.0
0.0
10.0
true
false
"" ""
PENS
"pen-0" 1.0 1 -7500403 true "histogram [hrsize] of turtles" "histogram [hrsize] of turtles"

PLOT
12
243
212
393
species richness
species ID
nr individuals
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -12895429 true "histogram [species] of turtles ;with [hrsize > 0]" "histogram [species] of turtles ;with [hrsize > 0]"

PLOT
12
399
212
549
body mass distribution
body mass / g
individuals
0.0
250.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "histogram [mass] of turtles with [hrsize > 0]" "histogram [mass * 1000] of turtles with [hrsize > 0]"

PLOT
698
14
1315
167
number of animals (total, gestating & lactating females)
time
number
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"gestating/ lactating" 1.0 0 -16777216 true "" "plot count turtles with [ preg = 1 ]"
"all" 1.0 0 -7500403 true "" "plot count turtles"

PLOT
698
173
1315
323
Number of offspring (total, failed, successful first year)
time
number
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"offspring" 1.0 0 -16777216 true "" "plot repro"
"failed" 1.0 0 -5298144 true "" "plot fail"
"age 1" 1.0 0 -8990512 true "" "plot count turtles with [ age = 0 ]"

PLOT
698
329
1316
492
Animals failing to re-establish homerange
number
time
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot fail2"
"pen-1" 1.0 0 -1184463 true "" "plot fail3"

TEXTBOX
424
506
591
526
1km
11
0.0
1

PLOT
410
541
1107
758
Species numbers
time
number of ind.
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"sp0" 1.0 0 -16777216 true "" "plot count turtles with [species = 0]"
"sp1" 1.0 0 -7500403 true "" "plot count turtles with [species = 1]"
"sp2" 1.0 0 -2674135 true "" "plot count turtles with [species = 2]"
"sp3" 1.0 0 -955883 true "" "plot count turtles with [species = 3]"
"sp4" 1.0 0 -6459832 true "" "plot count turtles with [species = 4]"
"sp5" 1.0 0 -1184463 true "" "plot count turtles with [species = 5]"
"sp6" 1.0 0 -10899396 true "" "plot count turtles with [species = 6]"
"sp7" 1.0 0 -13840069 true "" "plot count turtles with [species = 7]"
"sp8" 1.0 0 -14835848 true "" "plot count turtles with [species = 8]"
"sp9" 1.0 0 -11221820 true "" "plot count turtles with [species = 9]"

SLIDER
1390
19
1562
52
perc_small
perc_small
0
1
0.2
0.05
1
NIL
HORIZONTAL

SLIDER
1391
68
1563
101
size_small
size_small
1
100
1.0
1
1
NIL
HORIZONTAL

SLIDER
1391
120
1563
153
bold_prob
bold_prob
0
1
0.25
0.05
1
NIL
HORIZONTAL

SLIDER
1389
193
1561
226
total_cover
total_cover
0
1
0.05
0.05
1
NIL
HORIZONTAL

MONITOR
1116
541
1196
586
NIL
spec_num
17
1
11

CHOOSER
1391
243
1581
288
scenario
scenario
"base" "reduced" "stepping_stone" "increasing_base" "increasing_reduced"
0

CHOOSER
1393
307
1531
352
debug
debug
0 1
0

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="daily" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10950"/>
    <metric>ticks</metric>
    <metric>mean [hrsize] of turtles with [species = 0]</metric>
    <metric>mean [hrsize] of turtles with [species = 1]</metric>
    <metric>mean [hrsize] of turtles with [species = 2]</metric>
    <metric>mean [hrsize] of turtles with [species = 3]</metric>
    <metric>mean [hrsize] of turtles with [species = 4]</metric>
    <metric>mean [hrsize] of turtles with [species = 5]</metric>
    <metric>mean [hrsize] of turtles with [species = 6]</metric>
    <metric>mean [hrsize] of turtles with [species = 7]</metric>
    <metric>mean [hrsize] of turtles with [species = 8]</metric>
    <metric>mean [hrsize] of turtles with [species = 9]</metric>
    <metric>spec_num</metric>
    <metric>shannon</metric>
    <metric>mean [length spec_list] of patches with [patchquali = 1]</metric>
    <metric>mean [length spec_list] of patches with [patchquali = 2]</metric>
    <metric>standard-deviation [length spec_list] of patches with [patchquali = 1]</metric>
    <metric>standard-deviation [length spec_list] of patches with [patchquali = 2]</metric>
    <metric>min [length spec_list] of patches with [patchquali = 1]</metric>
    <metric>min [length spec_list] of patches with [patchquali = 2]</metric>
    <metric>max [length spec_list] of patches with [patchquali = 1]</metric>
    <metric>max [length spec_list] of patches with [patchquali = 2]</metric>
    <metric>(count turtles with [age = 0]) / repro</metric>
    <metric>count turtles with [species = 0 and shelter = 0]</metric>
    <metric>count turtles with [species = 1 and shelter = 0]</metric>
    <metric>count turtles with [species = 2 and shelter = 0]</metric>
    <metric>count turtles with [species = 3 and shelter = 0]</metric>
    <metric>count turtles with [species = 4 and shelter = 0]</metric>
    <metric>count turtles with [species = 5 and shelter = 0]</metric>
    <metric>count turtles with [species = 6 and shelter = 0]</metric>
    <metric>count turtles with [species = 7 and shelter = 0]</metric>
    <metric>count turtles with [species = 8 and shelter = 0]</metric>
    <metric>count turtles with [species = 9 and shelter = 0]</metric>
    <metric>count turtles with [species = 0]</metric>
    <metric>count turtles with [species = 1]</metric>
    <metric>count turtles with [species = 2]</metric>
    <metric>count turtles with [species = 3]</metric>
    <metric>count turtles with [species = 4]</metric>
    <metric>count turtles with [species = 5]</metric>
    <metric>count turtles with [species = 6]</metric>
    <metric>count turtles with [species = 7]</metric>
    <metric>count turtles with [species = 8]</metric>
    <metric>count turtles with [species = 9]</metric>
    <metric>mean [mass] of turtles</metric>
    <metric>mean [hrsize] of turtles</metric>
    <metric>mean [forage_small] of turtles with [species = 0]</metric>
    <metric>mean [forage_small] of turtles with [species = 1]</metric>
    <metric>mean [forage_small] of turtles with [species = 2]</metric>
    <metric>mean [forage_small] of turtles with [species = 3]</metric>
    <metric>mean [forage_small] of turtles with [species = 4]</metric>
    <metric>mean [forage_small] of turtles with [species = 5]</metric>
    <metric>mean [forage_small] of turtles with [species = 6]</metric>
    <metric>mean [forage_small] of turtles with [species = 7]</metric>
    <metric>mean [forage_small] of turtles with [species = 8]</metric>
    <metric>mean [forage_small] of turtles with [species = 9]</metric>
    <metric>count turtles with [core = 1 and species = 0]</metric>
    <metric>count turtles with [core = 1 and species = 1]</metric>
    <metric>count turtles with [core = 1 and species = 2]</metric>
    <metric>count turtles with [core = 1 and species = 3]</metric>
    <metric>count turtles with [core = 1 and species = 4]</metric>
    <metric>count turtles with [core = 1 and species = 5]</metric>
    <metric>count turtles with [core = 1 and species = 6]</metric>
    <metric>count turtles with [core = 1 and species = 7]</metric>
    <metric>count turtles with [core = 1 and species = 8]</metric>
    <metric>count turtles with [core = 1 and species = 9]</metric>
    <metric>count turtles with [core = 2 and species = 0]</metric>
    <metric>count turtles with [core = 2 and species = 1]</metric>
    <metric>count turtles with [core = 2 and species = 2]</metric>
    <metric>count turtles with [core = 2 and species = 3]</metric>
    <metric>count turtles with [core = 2 and species = 4]</metric>
    <metric>count turtles with [core = 2 and species = 5]</metric>
    <metric>count turtles with [core = 2 and species = 6]</metric>
    <metric>count turtles with [core = 2 and species = 7]</metric>
    <metric>count turtles with [core = 2 and species = 8]</metric>
    <metric>count turtles with [core = 2 and species = 9]</metric>
    <metric>count turtles with [ss_use and species = 0]</metric>
    <metric>count turtles with [ss_use and species = 1]</metric>
    <metric>count turtles with [ss_use and species = 2]</metric>
    <metric>count turtles with [ss_use and species = 3]</metric>
    <metric>count turtles with [ss_use and species = 4]</metric>
    <metric>count turtles with [ss_use and species = 5]</metric>
    <metric>count turtles with [ss_use and species = 6]</metric>
    <metric>count turtles with [ss_use and species = 7]</metric>
    <metric>count turtles with [ss_use and species = 8]</metric>
    <metric>count turtles with [ss_use and species = 9]</metric>
    <enumeratedValueSet variable="perc_small">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="size_small">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bold_prob">
      <value value="0.25"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_bold" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="18250"/>
    <metric>spec_num</metric>
    <enumeratedValueSet variable="total_cover">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perc_small">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bold_prob">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="size_small">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_base" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10950"/>
    <exitCondition>spec_num &lt; 5</exitCondition>
    <metric>ticks</metric>
    <metric>spec_num</metric>
    <metric>count turtles with [species = 0]</metric>
    <metric>count turtles with [species = 1]</metric>
    <metric>count turtles with [species = 2]</metric>
    <metric>count turtles with [species = 3]</metric>
    <metric>count turtles with [species = 4]</metric>
    <metric>count turtles with [species = 5]</metric>
    <metric>count turtles with [species = 6]</metric>
    <metric>count turtles with [species = 7]</metric>
    <metric>count turtles with [species = 8]</metric>
    <metric>count turtles with [species = 9]</metric>
    <metric>mean [mass] of turtles</metric>
    <metric>mean [hrsize] of turtles</metric>
    <metric>mean [forage_small] of turtles with [species = 0]</metric>
    <metric>mean [forage_small] of turtles with [species = 1]</metric>
    <metric>mean [forage_small] of turtles with [species = 2]</metric>
    <metric>mean [forage_small] of turtles with [species = 3]</metric>
    <metric>mean [forage_small] of turtles with [species = 4]</metric>
    <metric>mean [forage_small] of turtles with [species = 5]</metric>
    <metric>mean [forage_small] of turtles with [species = 6]</metric>
    <metric>mean [forage_small] of turtles with [species = 7]</metric>
    <metric>mean [forage_small] of turtles with [species = 8]</metric>
    <metric>mean [forage_small] of turtles with [species = 9]</metric>
    <metric>count turtles with [core = 1 and species = 0]</metric>
    <metric>count turtles with [core = 1 and species = 1]</metric>
    <metric>count turtles with [core = 1 and species = 2]</metric>
    <metric>count turtles with [core = 1 and species = 3]</metric>
    <metric>count turtles with [core = 1 and species = 4]</metric>
    <metric>count turtles with [core = 1 and species = 5]</metric>
    <metric>count turtles with [core = 1 and species = 6]</metric>
    <metric>count turtles with [core = 1 and species = 7]</metric>
    <metric>count turtles with [core = 1 and species = 8]</metric>
    <metric>count turtles with [core = 1 and species = 9]</metric>
    <metric>count turtles with [core = 2 and species = 0]</metric>
    <metric>count turtles with [core = 2 and species = 1]</metric>
    <metric>count turtles with [core = 2 and species = 2]</metric>
    <metric>count turtles with [core = 2 and species = 3]</metric>
    <metric>count turtles with [core = 2 and species = 4]</metric>
    <metric>count turtles with [core = 2 and species = 5]</metric>
    <metric>count turtles with [core = 2 and species = 6]</metric>
    <metric>count turtles with [core = 2 and species = 7]</metric>
    <metric>count turtles with [core = 2 and species = 8]</metric>
    <metric>count turtles with [core = 2 and species = 9]</metric>
    <metric>count turtles with [ss-use and species = 0]</metric>
    <metric>count turtles with [ss-use and species = 1]</metric>
    <metric>count turtles with [ss-use and species = 2]</metric>
    <metric>count turtles with [ss-use and species = 3]</metric>
    <metric>count turtles with [ss-use and species = 4]</metric>
    <metric>count turtles with [ss-use and species = 5]</metric>
    <metric>count turtles with [ss-use and species = 6]</metric>
    <metric>count turtles with [ss-use and species = 7]</metric>
    <metric>count turtles with [ss-use and species = 8]</metric>
    <metric>count turtles with [ss-use and species = 9]</metric>
    <enumeratedValueSet variable="perc_small">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="size_small">
      <value value="1"/>
      <value value="4"/>
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bold_prob">
      <value value="0"/>
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.15"/>
      <value value="0.25"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_cover">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario">
      <value value="&quot;base&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_reduced" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10950"/>
    <exitCondition>spec_num &lt; 5</exitCondition>
    <metric>ticks</metric>
    <metric>spec_num</metric>
    <metric>count turtles with [species = 0]</metric>
    <metric>count turtles with [species = 1]</metric>
    <metric>count turtles with [species = 2]</metric>
    <metric>count turtles with [species = 3]</metric>
    <metric>count turtles with [species = 4]</metric>
    <metric>count turtles with [species = 5]</metric>
    <metric>count turtles with [species = 6]</metric>
    <metric>count turtles with [species = 7]</metric>
    <metric>count turtles with [species = 8]</metric>
    <metric>count turtles with [species = 9]</metric>
    <metric>mean [mass] of turtles</metric>
    <metric>mean [hrsize] of turtles</metric>
    <metric>mean [forage_small] of turtles with [species = 0]</metric>
    <metric>mean [forage_small] of turtles with [species = 1]</metric>
    <metric>mean [forage_small] of turtles with [species = 2]</metric>
    <metric>mean [forage_small] of turtles with [species = 3]</metric>
    <metric>mean [forage_small] of turtles with [species = 4]</metric>
    <metric>mean [forage_small] of turtles with [species = 5]</metric>
    <metric>mean [forage_small] of turtles with [species = 6]</metric>
    <metric>mean [forage_small] of turtles with [species = 7]</metric>
    <metric>mean [forage_small] of turtles with [species = 8]</metric>
    <metric>mean [forage_small] of turtles with [species = 9]</metric>
    <metric>count turtles with [core = 1 and species = 0]</metric>
    <metric>count turtles with [core = 1 and species = 1]</metric>
    <metric>count turtles with [core = 1 and species = 2]</metric>
    <metric>count turtles with [core = 1 and species = 3]</metric>
    <metric>count turtles with [core = 1 and species = 4]</metric>
    <metric>count turtles with [core = 1 and species = 5]</metric>
    <metric>count turtles with [core = 1 and species = 6]</metric>
    <metric>count turtles with [core = 1 and species = 7]</metric>
    <metric>count turtles with [core = 1 and species = 8]</metric>
    <metric>count turtles with [core = 1 and species = 9]</metric>
    <metric>count turtles with [core = 2 and species = 0]</metric>
    <metric>count turtles with [core = 2 and species = 1]</metric>
    <metric>count turtles with [core = 2 and species = 2]</metric>
    <metric>count turtles with [core = 2 and species = 3]</metric>
    <metric>count turtles with [core = 2 and species = 4]</metric>
    <metric>count turtles with [core = 2 and species = 5]</metric>
    <metric>count turtles with [core = 2 and species = 6]</metric>
    <metric>count turtles with [core = 2 and species = 7]</metric>
    <metric>count turtles with [core = 2 and species = 8]</metric>
    <metric>count turtles with [core = 2 and species = 9]</metric>
    <metric>count turtles with [ss_use and species = 0]</metric>
    <metric>count turtles with [ss_use and species = 1]</metric>
    <metric>count turtles with [ss_use and species = 2]</metric>
    <metric>count turtles with [ss_use and species = 3]</metric>
    <metric>count turtles with [ss_use and species = 4]</metric>
    <metric>count turtles with [ss_use and species = 5]</metric>
    <metric>count turtles with [ss_use and species = 6]</metric>
    <metric>count turtles with [ss_use and species = 7]</metric>
    <metric>count turtles with [ss_use and species = 8]</metric>
    <metric>count turtles with [ss_use and species = 9]</metric>
    <enumeratedValueSet variable="perc_small">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="size_small">
      <value value="1"/>
      <value value="4"/>
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bold_prob">
      <value value="0"/>
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.15"/>
      <value value="0.25"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_cover">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario">
      <value value="&quot;reduced&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_stepping" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10950"/>
    <exitCondition>spec_num &lt; 5</exitCondition>
    <metric>ticks</metric>
    <metric>spec_num</metric>
    <metric>count turtles with [species = 0]</metric>
    <metric>count turtles with [species = 1]</metric>
    <metric>count turtles with [species = 2]</metric>
    <metric>count turtles with [species = 3]</metric>
    <metric>count turtles with [species = 4]</metric>
    <metric>count turtles with [species = 5]</metric>
    <metric>count turtles with [species = 6]</metric>
    <metric>count turtles with [species = 7]</metric>
    <metric>count turtles with [species = 8]</metric>
    <metric>count turtles with [species = 9]</metric>
    <metric>mean [mass] of turtles</metric>
    <metric>mean [hrsize] of turtles</metric>
    <metric>mean [forage_small] of turtles with [species = 0]</metric>
    <metric>mean [forage_small] of turtles with [species = 1]</metric>
    <metric>mean [forage_small] of turtles with [species = 2]</metric>
    <metric>mean [forage_small] of turtles with [species = 3]</metric>
    <metric>mean [forage_small] of turtles with [species = 4]</metric>
    <metric>mean [forage_small] of turtles with [species = 5]</metric>
    <metric>mean [forage_small] of turtles with [species = 6]</metric>
    <metric>mean [forage_small] of turtles with [species = 7]</metric>
    <metric>mean [forage_small] of turtles with [species = 8]</metric>
    <metric>mean [forage_small] of turtles with [species = 9]</metric>
    <metric>count turtles with [core = 1 and species = 0]</metric>
    <metric>count turtles with [core = 1 and species = 1]</metric>
    <metric>count turtles with [core = 1 and species = 2]</metric>
    <metric>count turtles with [core = 1 and species = 3]</metric>
    <metric>count turtles with [core = 1 and species = 4]</metric>
    <metric>count turtles with [core = 1 and species = 5]</metric>
    <metric>count turtles with [core = 1 and species = 6]</metric>
    <metric>count turtles with [core = 1 and species = 7]</metric>
    <metric>count turtles with [core = 1 and species = 8]</metric>
    <metric>count turtles with [core = 1 and species = 9]</metric>
    <metric>count turtles with [core = 2 and species = 0]</metric>
    <metric>count turtles with [core = 2 and species = 1]</metric>
    <metric>count turtles with [core = 2 and species = 2]</metric>
    <metric>count turtles with [core = 2 and species = 3]</metric>
    <metric>count turtles with [core = 2 and species = 4]</metric>
    <metric>count turtles with [core = 2 and species = 5]</metric>
    <metric>count turtles with [core = 2 and species = 6]</metric>
    <metric>count turtles with [core = 2 and species = 7]</metric>
    <metric>count turtles with [core = 2 and species = 8]</metric>
    <metric>count turtles with [core = 2 and species = 9]</metric>
    <metric>count turtles with [ss_use and species = 0]</metric>
    <metric>count turtles with [ss_use and species = 1]</metric>
    <metric>count turtles with [ss_use and species = 2]</metric>
    <metric>count turtles with [ss_use and species = 3]</metric>
    <metric>count turtles with [ss_use and species = 4]</metric>
    <metric>count turtles with [ss_use and species = 5]</metric>
    <metric>count turtles with [ss_use and species = 6]</metric>
    <metric>count turtles with [ss_use and species = 7]</metric>
    <metric>count turtles with [ss_use and species = 8]</metric>
    <metric>count turtles with [ss_use and species = 9]</metric>
    <enumeratedValueSet variable="perc_small">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="size_small">
      <value value="1"/>
      <value value="4"/>
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bold_prob">
      <value value="0"/>
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.15"/>
      <value value="0.25"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_cover">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario">
      <value value="&quot;stepping_stone&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_increasing_base" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10950"/>
    <exitCondition>spec_num &lt; 5</exitCondition>
    <metric>ticks</metric>
    <metric>spec_num</metric>
    <metric>count turtles with [species = 0]</metric>
    <metric>count turtles with [species = 1]</metric>
    <metric>count turtles with [species = 2]</metric>
    <metric>count turtles with [species = 3]</metric>
    <metric>count turtles with [species = 4]</metric>
    <metric>count turtles with [species = 5]</metric>
    <metric>count turtles with [species = 6]</metric>
    <metric>count turtles with [species = 7]</metric>
    <metric>count turtles with [species = 8]</metric>
    <metric>count turtles with [species = 9]</metric>
    <metric>mean [mass] of turtles</metric>
    <metric>mean [hrsize] of turtles</metric>
    <metric>mean [forage_small] of turtles with [species = 0]</metric>
    <metric>mean [forage_small] of turtles with [species = 1]</metric>
    <metric>mean [forage_small] of turtles with [species = 2]</metric>
    <metric>mean [forage_small] of turtles with [species = 3]</metric>
    <metric>mean [forage_small] of turtles with [species = 4]</metric>
    <metric>mean [forage_small] of turtles with [species = 5]</metric>
    <metric>mean [forage_small] of turtles with [species = 6]</metric>
    <metric>mean [forage_small] of turtles with [species = 7]</metric>
    <metric>mean [forage_small] of turtles with [species = 8]</metric>
    <metric>mean [forage_small] of turtles with [species = 9]</metric>
    <metric>count turtles with [core = 1 and species = 0]</metric>
    <metric>count turtles with [core = 1 and species = 1]</metric>
    <metric>count turtles with [core = 1 and species = 2]</metric>
    <metric>count turtles with [core = 1 and species = 3]</metric>
    <metric>count turtles with [core = 1 and species = 4]</metric>
    <metric>count turtles with [core = 1 and species = 5]</metric>
    <metric>count turtles with [core = 1 and species = 6]</metric>
    <metric>count turtles with [core = 1 and species = 7]</metric>
    <metric>count turtles with [core = 1 and species = 8]</metric>
    <metric>count turtles with [core = 1 and species = 9]</metric>
    <metric>count turtles with [core = 2 and species = 0]</metric>
    <metric>count turtles with [core = 2 and species = 1]</metric>
    <metric>count turtles with [core = 2 and species = 2]</metric>
    <metric>count turtles with [core = 2 and species = 3]</metric>
    <metric>count turtles with [core = 2 and species = 4]</metric>
    <metric>count turtles with [core = 2 and species = 5]</metric>
    <metric>count turtles with [core = 2 and species = 6]</metric>
    <metric>count turtles with [core = 2 and species = 7]</metric>
    <metric>count turtles with [core = 2 and species = 8]</metric>
    <metric>count turtles with [core = 2 and species = 9]</metric>
    <metric>count turtles with [ss_use and species = 0]</metric>
    <metric>count turtles with [ss_use and species = 1]</metric>
    <metric>count turtles with [ss_use and species = 2]</metric>
    <metric>count turtles with [ss_use and species = 3]</metric>
    <metric>count turtles with [ss_use and species = 4]</metric>
    <metric>count turtles with [ss_use and species = 5]</metric>
    <metric>count turtles with [ss_use and species = 6]</metric>
    <metric>count turtles with [ss_use and species = 7]</metric>
    <metric>count turtles with [ss_use and species = 8]</metric>
    <metric>count turtles with [ss_use and species = 9]</metric>
    <enumeratedValueSet variable="perc_small">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="size_small">
      <value value="1"/>
      <value value="4"/>
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bold_prob">
      <value value="0"/>
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.15"/>
      <value value="0.25"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_cover">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario">
      <value value="&quot;increasing_base&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_increasing_reduced" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10950"/>
    <exitCondition>spec_num &lt; 5</exitCondition>
    <metric>ticks</metric>
    <metric>spec_num</metric>
    <metric>count turtles with [species = 0]</metric>
    <metric>count turtles with [species = 1]</metric>
    <metric>count turtles with [species = 2]</metric>
    <metric>count turtles with [species = 3]</metric>
    <metric>count turtles with [species = 4]</metric>
    <metric>count turtles with [species = 5]</metric>
    <metric>count turtles with [species = 6]</metric>
    <metric>count turtles with [species = 7]</metric>
    <metric>count turtles with [species = 8]</metric>
    <metric>count turtles with [species = 9]</metric>
    <metric>mean [mass] of turtles</metric>
    <metric>mean [hrsize] of turtles</metric>
    <metric>mean [forage_small] of turtles with [species = 0]</metric>
    <metric>mean [forage_small] of turtles with [species = 1]</metric>
    <metric>mean [forage_small] of turtles with [species = 2]</metric>
    <metric>mean [forage_small] of turtles with [species = 3]</metric>
    <metric>mean [forage_small] of turtles with [species = 4]</metric>
    <metric>mean [forage_small] of turtles with [species = 5]</metric>
    <metric>mean [forage_small] of turtles with [species = 6]</metric>
    <metric>mean [forage_small] of turtles with [species = 7]</metric>
    <metric>mean [forage_small] of turtles with [species = 8]</metric>
    <metric>mean [forage_small] of turtles with [species = 9]</metric>
    <metric>count turtles with [core = 1 and species = 0]</metric>
    <metric>count turtles with [core = 1 and species = 1]</metric>
    <metric>count turtles with [core = 1 and species = 2]</metric>
    <metric>count turtles with [core = 1 and species = 3]</metric>
    <metric>count turtles with [core = 1 and species = 4]</metric>
    <metric>count turtles with [core = 1 and species = 5]</metric>
    <metric>count turtles with [core = 1 and species = 6]</metric>
    <metric>count turtles with [core = 1 and species = 7]</metric>
    <metric>count turtles with [core = 1 and species = 8]</metric>
    <metric>count turtles with [core = 1 and species = 9]</metric>
    <metric>count turtles with [core = 2 and species = 0]</metric>
    <metric>count turtles with [core = 2 and species = 1]</metric>
    <metric>count turtles with [core = 2 and species = 2]</metric>
    <metric>count turtles with [core = 2 and species = 3]</metric>
    <metric>count turtles with [core = 2 and species = 4]</metric>
    <metric>count turtles with [core = 2 and species = 5]</metric>
    <metric>count turtles with [core = 2 and species = 6]</metric>
    <metric>count turtles with [core = 2 and species = 7]</metric>
    <metric>count turtles with [core = 2 and species = 8]</metric>
    <metric>count turtles with [core = 2 and species = 9]</metric>
    <metric>count turtles with [ss_use and species = 0]</metric>
    <metric>count turtles with [ss_use and species = 1]</metric>
    <metric>count turtles with [ss_use and species = 2]</metric>
    <metric>count turtles with [ss_use and species = 3]</metric>
    <metric>count turtles with [ss_use and species = 4]</metric>
    <metric>count turtles with [ss_use and species = 5]</metric>
    <metric>count turtles with [ss_use and species = 6]</metric>
    <metric>count turtles with [ss_use and species = 7]</metric>
    <metric>count turtles with [ss_use and species = 8]</metric>
    <metric>count turtles with [ss_use and species = 9]</metric>
    <enumeratedValueSet variable="perc_small">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="size_small">
      <value value="1"/>
      <value value="4"/>
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bold_prob">
      <value value="0"/>
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.15"/>
      <value value="0.25"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="total_cover">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario">
      <value value="&quot;increasing_reduced&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
