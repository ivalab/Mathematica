(*--Defines the constants that allow one to specify the loading
  of certain objects.--*)

<<ivamatica/Basic/objects.m
<<ivamatica/Math/binary.m;

oPLEN = 11;

oEUCLIDEAN = oInit[Binary,oPLEN,2^0];
oBUNDLES = oInit[Binary,oPLEN,2^1];
oMANIFOLDS = oInit[Binary,oPLEN,2^2];
oVECTORS = oInit[Binary,oPLEN,2^3];
oTANGENTS = oInit[Binary,oPLEN,2^4];
oTBUNDLES = oInit[Binary, oPLEN,2^5];
oTANGENTMANIFOLDS = oInit[Binary, oPLEN, 2^6];

oLIEGROUPS = oInit[Binary,oPLEN,2^7];
oTLIEGROUPS = oInit[Binary,oPLEN,2^8];
oPRINCIPAL = oInit[Binary,oPLEN,2^9];
oTPRINCIPAL = oInit[Binary,oPLEN,2^10];

oALL = oInit[Binary, oPLEN, Sum[ 2^k , {k, 0, oPLEN-1}] ];
oNONE = oInit[Binary, oPLEN, 0 ];

EUCLIDEAN = oEUCLIDEAN;
MANIFOLDS = EUCLIDEAN+oMANIFOLDS;
BUNDLES = EUCLIDEAN+oBUNDLES;
VECTORS = EUCLIDEAN+oVECTORS;
TANGENTS = BUNDLES+VECTORS+oTANGENTS;
TBUNDLES = TANGENTS + oTBUNDLES;
TANGENTMANIFOLDS = TANGENTS+MANIFOLDS+oTANGENTMANIFOLDS;

LIEGROUPS = MANIFOLDS+VECTORS+oLIEGROUPS;
TLIEGROUPS = TANGENTMANIFOLDS + LIEGROUPS + oTLIEGROUPS;
PRINCIPAL = LIEGROUPS+BUNDLES+oPRINCIPAL;
TPRINCIPAL = PRINCIPAL+TLIEGROUPS+TBUNDLES+oTPRINCIPAL;

ALL = oALL;
NONE = oNONE;


LoadPackages[pnum_] := Module[{loaded},

  loaded = {};
  If[ pnum == oNONE, loaded = "No packages loaded."];

  If[ pnum && oEUCLIDEAN , 
    <<ivamatica/DiffGeometry/euclidean.m; 
    loaded = StringJoin[loaded, "Loaded Euclidean package \n" ]];

  If[ pnum && oBUNDLES , 
    <<ivamatica/DiffGeometry/bundles.m; 
    loaded = StringJoin[loaded, "Loaded Fiber Bundles package \n"]];

  If[ pnum && oMANIFOLDS , 
    <<ivamatica/DiffGeometry/manifolds.m; 
    loaded = StringJoin[loaded, "Loaded Manifolds package \n" ]];

  If[ pnum && oVECTORS , 
    <<ivamatica/DiffGeometry/vectors.m; 
    loaded = StringJoin[loaded, "Loaded Vectors package \n"] ];

  If[ pnum && oTANGENTS , 
    <<ivamatica/DiffGeometry/tangents.m; 
    loaded = StringJoin[loaded, "Loaded Tangents package \n"] ];

  If[ pnum && oTBUNDLES , 
    <<ivamatica/DiffGeometry/tbundles.m; 
    loaded = StringJoin[loaded, "Loaded Tangent Fiber Bundles package \n"] ];

  If[ pnum && oTANGENTMANIFOLDS , 
    <<ivamatica/DiffGeometry/tmanifolds.m; 
    loaded = StringJoin[loaded, "Loaded Tangent Manifolds package \n"] ];

  If[ pnum && oLIEGROUPS , 
    <<ivamatica/GeoMechanics/liegroups.m; 
    loaded = StringJoin[loaded, "Loaded Lie Groups package \n"] ];

  If[ pnum && oTLIEGROUPS , 
    <<ivamatica/GeoMechanics/tliegroups.m; 
    loaded = StringJoin[loaded, "Loaded Tangent Lie Groups package \n"] ];

  If[ pnum && oPRINCIPAL , 
    <<ivamatica/GeoMechanics/principal.m; 
    loaded = StringJoin[loaded, "Loaded Principal Bundles package \n"] ];

  If[ pnum && oTPRINCIPAL , 
    <<ivamatica/GeoMechanics/tprincipal.m; 
    loaded = StringJoin[loaded, "Loaded Tangent Principal Bundles package \n"] ];
  loaded
 ];


DoneLoading := Module[{},

  Remove[oPLEN, oEUCLIDEAN, oMANIFOLDS, oVECTORS, oBUNDLES,
    oTANGENTS, oTBUNDLES, oTANGENTMANIFOLDS, oLIEGROUPS, oTLIEGROUPS,
    oPRINCIPAL, oTPRINCIPAL, oALL, oNONE];

  Remove[EUCLIDEAN, MANIFOLDS, VECTORS, BUNDLES, TANGENTS,
    TBUNDLES, TANGENTMANIFOLDS, LIEGROUPS, TLIEGROUPS, PRINCIPAL, 
    TPRINCIPAL, ALL, NONE];

  "Removed package loading variables."
];

