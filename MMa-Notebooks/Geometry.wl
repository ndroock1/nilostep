(* Wolfram Language Package *)


(* Created by the Wolfram Workbench 03/7-'18 *)
(* Version X.1 05/09-'23 *)


BeginPackage["Geometry`"]
(* Exported symbols added here with SymbolName::usage *)


(* CONSTANTS *)
FACE:={Red,1};
EDGE:={Black,1,.0125,2\[Pi],2\[Pi]};

POINT:=0;
POLYGON:=1;
CIRCLE:=2;
LINE:=3;
DISK:=4;

FILLEDCURVE:=5;
HYPERBOLICLINE:=6

kanweg:=999


(* DATA STRUCTURES *)
pnt = {Repeated[Except[_Complex, _?NumericQ], {2}]};
pts = {RepeatedNull[pnt]};
nls = {Repeated[_?NumericQ]};

typ = {_?IntegerQ | _Real, {RepeatedNull[pnt]}};
face = {_RGBColor | _GrayLevel, _?IntegerQ | _Real};
edge = {_RGBColor | _GrayLevel, _?IntegerQ | _Real, _?IntegerQ | _Real, _?NumericQ, _?NumericQ};
clr = {typ, face};
edg = {{typ, face}, edge};



(* GEOMETRIES *)
(* ========== *)


(* POINT : Type = 0 *)
ClearAll[gPoint]
Options[gPoint] = {"Color" -> Red, "Opacity" -> 1, 
   "PointSize" -> 0.0125};
gPoint[p : pnt, 
  OptionsPattern[]] := {{{POINT, {p}}, {OptionValue[Color], 
    OptionValue[Opacity]}}, {GrayLevel[0], 1, OptionValue[PointSize], 
   2 Pi, 2 Pi}}
gPoint[p : pts, OptionsPattern[]] := 
 p /. {q : pnt :> 
    gPoint[q, Color -> OptionValue[Color], 
     Opacity -> OptionValue[Opacity], 
     PointSize -> OptionValue[PointSize]]}
gPoint[p_?NumericQ, 
  OptionsPattern[]] := {{{POINT, {ReIm[N[p]]}}, {OptionValue[Color], 
    OptionValue[Opacity]}}, {GrayLevel[0], 1, OptionValue[PointSize], 
   2 Pi, 2 Pi}}
gPoint[p : nls, OptionsPattern[]] := 
 p /. {q_?NumericQ :> 
    gPoint[q, Color -> OptionValue[Color], 
     Opacity -> OptionValue[Opacity], 
     PointSize -> OptionValue[PointSize]]}
gPoint[OptionsPattern[]] := 
 gPoint[0, Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  PointSize -> OptionValue[PointSize]]


(* POLYGON : Type = 1 *)
ClearAll[gPolygon];
Options[gPolygon] = {"FaceColor" -> Red, "FaceOpacity" -> 1, 
   "EdgeColor" -> Black, "EdgeOpacity" -> 1, "Thickness" -> 0.0125, 
   "Dashing" -> {2 Pi, 2 Pi}};
gPolygon[p : pts, 
  OptionsPattern[]] := {{{POLYGON, p}, {OptionValue[FaceColor], 
    OptionValue[FaceOpacity]}}, {OptionValue[EdgeColor], 
   OptionValue[EdgeOpacity], OptionValue[Thickness], 
   OptionValue[Dashing][[1]], OptionValue[Dashing][[2]]}}
gPolygon[p : nls, 
  OptionsPattern[]] := {{{POLYGON, ReIm[p]}, {OptionValue[FaceColor], 
    OptionValue[FaceOpacity]}}, {OptionValue[EdgeColor], 
   OptionValue[EdgeOpacity], OptionValue[Thickness], 
   OptionValue[Dashing][[1]], OptionValue[Dashing][[2]]}}
gPolygon[OptionsPattern[]] := gPolygon[{I, 0, 1},
  "FaceColor" -> OptionValue[FaceColor], 
  "FaceOpacity" -> OptionValue[FaceOpacity], 
  "EdgeColor" -> OptionValue[EdgeColor], 
  "EdgeOpacity" -> OptionValue[EdgeOpacity],
  "Thickness" -> OptionValue[Thickness], 
  "Dashing" -> {OptionValue[Dashing][[1]], OptionValue[Dashing][[2]]}]


(* CIRCLE : Type = 2 *)
ClearAll[gCircle]
Options[gCircle] = {"Color" -> Red, "Opacity" -> 1, 
   "Thickness" -> 0.0125, "Dashing" -> {2 Pi, 2 Pi}};
gCircle[cntr : pnt, r_?NumericQ, alfa1_?NumericQ, alfa2_?NumericQ, 
  OptionsPattern[]] := (ParametricPlot[
     {cntr[[1]] + r Cos[t], cntr[[2]] + r Sin[t]}, {t, alfa1, alfa2}, 
     PlotStyle -> {OptionValue[Color], 
       Opacity -> OptionValue[Opacity], 
       Thickness -> OptionValue[Thickness], 
       Dashing -> {OptionValue[Dashing][[1]], 
         OptionValue[Dashing][[2]]}}] // InputForm)[[1, 1, 1, 1, 1]]
gCircle[cntr : pnt, r_?NumericQ, OptionsPattern[]] := {
  {
   {
    CIRCLE, {cntr, {cntr[[1]], cntr[[2]] + r}}
    },
   {
    OptionValue[Color], OptionValue[Opacity]
    }
   },
  {
   GrayLevel[0], 1, OptionValue[Thickness], OptionValue[Dashing][[1]],
    OptionValue[Dashing][[2]]
   }
  }
gCircle[cntr : pnt, OptionsPattern[]] := 
 gCircle[cntr, 1, Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]
gCircle[cntr_?NumericQ, r_?NumericQ, OptionsPattern[]] := {
  {
   {
    CIRCLE, {ReIm[N[cntr]], ReIm[N[cntr + r]]}
    },
   {
    OptionValue[Color], OptionValue[Opacity]
    }
   },
  {
   GrayLevel[0], 1, OptionValue[Thickness], OptionValue[Dashing][[1]],
    OptionValue[Dashing][[2]]
   }
  }
gCircle[cntr_?NumericQ, OptionsPattern[]] := 
 gCircle[cntr, 1, Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]
gCircle[OptionsPattern[]] := 
 gCircle[0, 1, Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  Thickness -> OptionValue[Thickness], 
  Dashing -> OptionValue[Dashing]]


(* LINE : Type = 3 *)
ClearAll[gLine]
Options[gLine] = {"Color" -> Red, "Opacity" -> 1, 
   "Thickness" -> 0.0125, "Dashing" -> {2 Pi, 2 Pi}};
gLine[p : pnt, q : pnt, OptionsPattern[]] := {
  {
   {
    LINE, {p, q}
    },
   {
    OptionValue[Color], OptionValue[Opacity]
    }
   },
  {
   OptionValue[Color], OptionValue[Opacity], OptionValue[Thickness], 
   OptionValue[Dashing][[1]], OptionValue[Dashing][[2]]
   }
  }
gLine[p : pnt, OptionsPattern[]] := 
 gLine[p, {p[[1]] + 1, p[[2]] + 1}, Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]
gLine[p_?NumericQ, q_?NumericQ, OptionsPattern[]] := 
 gLine[ReIm[p], ReIm[q], Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]
gLine[p_?NumericQ, OptionsPattern[]] := 
 gLine[ReIm[p], ReIm[p] + {1, 1}, Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]
gLine[OptionsPattern[]] := 
 gLine[{0, 0}, {1, 1}, Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  Thickness -> OptionValue[Thickness], 
  Dashing -> OptionValue[Dashing]]


(* DISK : Type = 4 *)
ClearAll[gDisk];
Options[gDisk] = {"FaceColor" -> Red, "FaceOpacity" -> 1, 
   "EdgeColor" -> Black, "EdgeOpacity" -> 1, "Thickness" -> 0.0125, 
   "Dashing" -> {2 Pi, 2 Pi}};
gDisk[cntr : pnt, r_?NumericQ, 
  OptionsPattern[]] := {{{DISK, {cntr, {cntr[[1]], 
      cntr[[2]] + r}}}, {OptionValue[FaceColor], 
    OptionValue[FaceOpacity]}}, {OptionValue[EdgeColor], 
   OptionValue[EdgeOpacity], OptionValue[Thickness], 
   OptionValue[Dashing][[1]], OptionValue[Dashing][[2]]}}
gDisk[cntr : pnt, OptionsPattern[]] := 
 gDisk[cntr, 1, FaceColor -> OptionValue[FaceColor], 
  FaceOpacity -> OptionValue[FaceOpacity], 
  EdgeColor -> OptionValue[EdgeColor], 
  EdgeOpacity -> OptionValue[EdgeOpacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]
gDisk[cntr_?NumericQ, r_?NumericQ, OptionsPattern[]] := 
 gDisk[ReIm[cntr], r, FaceColor -> OptionValue[FaceColor], 
  FaceOpacity -> OptionValue[FaceOpacity], 
  EdgeColor -> OptionValue[EdgeColor], 
  EdgeOpacity -> OptionValue[EdgeOpacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]
gDisk[cntr_?NumericQ, OptionsPattern[]] := 
 gDisk[ReIm[cntr], 1, FaceColor -> OptionValue[FaceColor], 
  FaceOpacity -> OptionValue[FaceOpacity], 
  EdgeColor -> OptionValue[EdgeColor], 
  EdgeOpacity -> OptionValue[EdgeOpacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]
gDisk[OptionsPattern[]] := 
 gDisk[0, 1, FaceColor -> OptionValue[FaceColor], 
  FaceOpacity -> OptionValue[FaceOpacity], 
  EdgeColor -> OptionValue[EdgeColor], 
  EdgeOpacity -> OptionValue[EdgeOpacity], 
  Thickness -> OptionValue[Thickness], Dashing -> OptionValue[Dashing]]


(* HYPERBOLIC LINE UPPER HALF-PLANE : via gCircle *)
Options[gHLineHP] = {"Color" -> Red, "Opacity" -> 1, 
   "Thickness" -> 0.0125, "Dashing" -> {2 Pi, 2 Pi}};
gHLineHP[p : pnt, q : pnt, OptionsPattern[]] := Module[
   {A, B, c1, c2, R, a1, a2},
   p1 = p[[1]];
   p2 = p[[2]];
   q1 = q[[1]];
   q2 = q[[2]];
   A = -(q1 - p1)/(q2 - p2);
   B = (p2 + q2)/2 + ((q1 - p1)/(q2 - p2)) (p1 + q1)/2;
   c1 = -B/A;
   c2 = 0;
   R = EuclideanDistance[p, {c1, c2}];
   a1 = ArcTan[q1 - c1, q2];
   a2 = ArcTan[p1 - c1, p2];
   gCircle[{c1, c2}, R, a1, a2, Color -> OptionValue[Color], 
    Opacity -> OptionValue[Opacity], 
    Thickness -> OptionValue[Thickness], 
    Dashing -> {OptionValue[Dashing][[1]], OptionValue[Dashing][[2]]} ]
   ] /; p[[2]] != q[[2]]
   
   
(* HYPERBOLIC LINE POINCARE DISK : via gCircle *)
ClearAll[gHLinePD]
Options[gHLinePD] = {"Color" -> Red, "Opacity" -> 1, 
   "Thickness" -> 0.0125, "Dashing" -> {2 Pi, 2 Pi}};
gHLinePD[p : pnt, q : pnt, OptionsPattern[]] := Module[
  {mPP, mQQ, Ap, Bp, Aq, Bq, sol, x, y, cntr, rad, ang1, ang2},
  mPP = (p + cnj3[mob[p, {{0, 1}, {1, 0}}]])/2;
  mQQ = (q + cnj3[mob[q, {{0, 1}, {1, 0}}]])/2;
  Ap = (-mPP[[1]] + p[[1]])/(mPP[[2]] - p[[2]]);
  Bp = mPP[[2]] - Ap mPP[[1]];
  Aq = (-mQQ[[1]] + q[[1]])/(mQQ[[2]] - q[[2]]);
  Bq = mQQ[[2]] - Aq mQQ[[1]];
  sol = Solve[{y == Ap x + Bp, y == Aq x + Bq}, {x, y}];
  x = sol[[1, 1, 2]];
  y = sol[[1, 2, 2]];
  cntr = {x, y};
  rad = EuclideanDistance[cntr, p];
  ang1 = ArcTan[p[[1]] - x, p[[2]] - y];
  ang2 = ArcTan[q[[1]] - x, q[[2]] - y];
  (*Print[{cntr, rad, ang1, ang2}];
  Print[OptionValue[Color]];*)
  gCircle[cntr, rad, ang1, ang2, Color -> OptionValue[Color], 
  Opacity -> OptionValue[Opacity], 
  Thickness -> OptionValue[Thickness], 
  Dashing -> OptionValue[Dashing]]//toGL
  ]


(* GEOMETRIC TRANSFORMATIONS *)
tra3[figs_,t_]:=figs/.f:pnt:> TranslationTransform[t][f]
rot3[figs_,{a_,p_}]:=figs/.f:pnt:> RotationTransform[a,p][f]
ref3[figs_,{v_,p_}]:=figs/.f:pnt:> ReflectionTransform[v,p][f]
cnj3[figs_]:=ref3[figs,{{0,1},{0,0}}];
sca3[figs_,{s_,p_}]:=figs/.f:pnt :> ScalingTransform[s,p][f]

trf3[figs_,trf_]:=figs/.f:pnt:> trf[f]//N

mid[fig_]:={(Max[Cases[{fig},pnt,Infinity][[All,1]]]+Min[Cases[{fig},pnt,Infinity][[All,1]]])/2,(Max[Cases[{fig},pnt,Infinity][[All,2]]]+Min[Cases[{fig},pnt,Infinity][[All,2]]])/2}
rad[fig_]:=Max[{Max[Cases[{fig},pnt,Infinity][[All,1]]]-Min[Cases[{fig},pnt,Infinity][[All,1]]],Max[Cases[{fig},pnt,Infinity][[All,2]]]-Min[Cases[{fig},pnt,Infinity][[All,2]]]}]
sca2pnt[fig_,{ce_,ra_}]:=tra3[sca3[tra3[fig,-mid[fig]],{{(2 ra)/rad[fig],(2 ra)/rad[fig]},{0,0}}],ce]

complex2Euclidean[g_]:=ComplexExpand[Through[{Re,Im}[(g[#])&[x+I*y]]]]
mob[p_,mat_]:=Function[{p1,p2, mob},complex2Euclidean[(mob[[1,1]] #+mob[[1,2]])/(mob[[2,1]] #+mob[[2,2]])&]/.{x:>p1,y:>p2}][p[[1]],p[[2]], mat]
mob3[pts_,mat_]:=pts/.p:pnt:> mob[p,mat]

lattice[{{p_,q_},w_,h_}]:=Flatten[Table[i p + j q,{i,w},{j,h}],1]



sym4rot[figke_] := Module[{rotp},
  rotp = {1, 1};
  tra3[Scale[{figke, rot3[figke, {\[Pi]/2, rotp}], 
     rot3[figke, {\[Pi], rotp}], 
     rot3[figke, {3 \[Pi] /2, rotp}]}, .5], {-.5, -.5}]
  ]
  
sym4ref[figke_] := Module[{rotp},
  rotp = {1, 1};
  tra3[Scale[{figke, ref3[figke, {{1, 0}, rotp}], 
     rot3[figke, {\[Pi], rotp}], 
     ref3[figke, {{0, 1}, rotp}]}, .5], {-.5, -.5}]
  ]

sym4ovr[figke_] := Module[{rotp},
  rotp = {1, 1};
  tra3[Scale[
    {figke, ref3[figke, {{1, 0}, rotp}], tra3[figke, {1, 1}], 
     ref3[tra3[figke, {1, 1}], {{1, 0}, rotp}]},
    .5], {-.5, -.5}]]


(* *)

(* GEOMETRIES *)
pGon[n_Integer] := Module[{},Table[{Cos[(2 \[Pi] k)/n + \[Pi]/n], Sin[(2 \[Pi] k)/n + \[Pi]/n]},{k, 0, n - 1}]]

newGeo[data : edg] := data;
newGeo[data : clr] := {{data}, EDGE};
newGeo[data : typ] := {{data, FACE}, EDGE};
newGeo[data_Integer] := {{{1, pGon[data]}, FACE},EDGE}
newGeo[n_Integer, col:(_RGBColor | _GrayLevel)] := repCol[newGeo[n], {col}][[1]]
newGeo[n_Integer, col:(_RGBColor | _GrayLevel), alfa_Real] := setVal[newGeo[n, col], {{1, 2, 2}, {alfa}}][[1]]

cCircle[p:pnt, n_?NumericQ] := {{{2, {p,{p[[1]]+n,p[[2]]}}},FACE}, EDGE}


(* PROPERTIES *)
repCol[data_, cols_] := MapIndexed[ReplacePart[#1, {1, 2, 1} -> cols[[#2]][[1]]] &, Cases[{data}, edg, Infinity]]
repLineCol[data_, cols_] := MapIndexed[ReplacePart[#1, {2, 1} -> cols[[#2]][[1]]] &, Cases[{data}, edg, Infinity]]

setVal[data_, {pos_, val_}]:=MapIndexed[ReplacePart[#1, pos -> val[[#2]][[1]]] &, Cases[{data}, edg, Infinity]]
setEdges[data_, edges_] := MapIndexed[ReplacePart[#1, {2} -> edges[[#2[[1]]]]] &, Cases[{data}, edg, Infinity]]
repLineAlfa[data_, alfa_] := MapIndexed[ReplacePart[#1, {2, 2} -> alfa[[#2]][[1]]] &, Cases[{data}, edg, Infinity]]

(* Chain of connected nGons *)
toUnit[{p:pnt,q:pnt}]:=Composition[ScalingTransform[{1/Norm[q-p],1/Norm[q-p]},{0,0}],RotationTransform[-ArcTan[(q-p)[[1]],(q-p)[[2]]],{0,0}],TranslationTransform[-p]]
toVector[{p:pnt,q:pnt}]:=InverseFunction[toUnit[{q,p}]]
toEdges[data_]:=data/.geo:edg:>Flatten[Map[Function[k,{ReplacePart[geo,{1,1}->{3,k}]}],Partition[geo[[1,1,2]],2,1,{1,1}]],1]
nGonAtEdge3[fig_,{m_,n_}]:=trf3[trf3[newGeo[n],toUnit[toEdges[newGeo[n]][[1,1,1,2]]]//N],toVector[toEdges[fig][[m,1,1,2]]]]

(* Envelope of Lines *)
newEnvelope[points_, corner_, numlines_,col_,alfa_]:=repLineAlfa[repLineCol[newEnvelope[points, corner, numlines], Array[col&, (numlines+2)(Length[points]-2)]],Array[alfa&, (numlines+2)(Length[points]-2)]]

newEnvelope[points_, corner_, numlines_,col_]:=repLineCol[newEnvelope[points, corner, numlines], Array[col&, (numlines+2)(Length[points]-2)]]
newEnvelope[points_, corner_, numlines_] := Module[{poly},
  poly = RotateLeft[points, corner - 1];
  Join[{newGeo[{3, {poly[[1]], poly[[2]]}}],
    Flatten[
     Map[
      Function[line,
       Map[Function[k,
         newGeo[{3, {poly[[1]], (1 - k/(numlines + 1)) line[[1]] + 
             k/(numlines + 1) line[[2]]}}]
         ], Range[numlines + 1]]  
       ],
      Rest[Partition[poly, 2, 1]]
      ]
     , 1]}]]

popGeo3[p : pnt, q : pnt, n_] := popGeo3[p, q, {n, n}]
popGeo3[p : pnt, q : pnt, {n_, m_}] := Module[
  {
   d1 = (q[[1]] - p[[1]])/n,
   d2 = (q[[2]] - p[[2]])/m
   },
  Flatten[Table[
    {{x, y}, {x + d1, y}, {x + d1, y + d2}, {x, y + d2}},
    {x, p[[1]], q[[1]] - d1, d1},
    {y, p[[2]], q[[2]] - d2, d2}]
   , 1]
  ]
  
matrixToGeo3[mat_, at_, dim_] := 
 Module[{dimnew = Dimensions[mat]*dim},
  (*Print[mat];
  Print[at];
  Print[dimnew];
  Print[dimnew[[1]]];
  Print[dimnew[[2]]];
  Print["======"];*)
  Cases[
   Map[
    If[Head[#[[1]]] === Integer,
      {
       at + ((#[[2]] - {1, 1}))/dimnew,
       at + ((#[[2]] - {1, 1}))/dimnew + {1/dimnew[[1]], 0},
       at + ((#[[2]] - {1, 1}))/dimnew + {1/dimnew[[1]], 
         1/dimnew[[2]]},
       at + ((#[[2]] - {1, 1}))/dimnew + {0, 1/dimnew[[2]]}
       },
      matrixToGeo3[#[[1]], at + ((#[[2]] - {1, 1}))/dimnew, 
       Dimensions[mat]*dim]
      ] &,
    Flatten[MapIndexed[{#1, #2} &, mat, {2}], 1], 1],
   {Repeated[{_, _}, 4]}, Infinity](*/.p:pnt\[RuleDelayed]{p[[2]],p[[
  1]]}*)
  ]



(* GEOMETRY UTILITIES *)
pntc = {Repeated[pnt,{3}]};
pts2Circle[v : pntc] := Module[{p, q, r, centre, radius},
  p=v[[1]];
  q=v[[2]];
  r=v[[3]];
  sol = LinearSolve[
    {{2 p[[1]], 2 p[[2]], -1},
     {2 q[[1]], 2 q[[2]], -1},
     {2 r[[1]], 2 r[[2]], -1}
     },
    {p[[1]]^2 + p[[2]]^2,
     q[[1]]^2 + q[[2]]^2,
     r[[1]]^2 + r[[2]]^2
     }];
  centre = {sol[[1]], sol[[2]]};
  radius = Sqrt[sol[[1]]^2 + sol[[2]]^2 - sol[[3]]];
  cCircle[centre, radius]
]


ClearAll[lFun1]
lFun1[point1 : pnt, point2 : pnt] := Module[{a, b},
  {a, b} = 
   LinearSolve[{{point1[[1]], 1}, {point2[[1]], 1}}, {point1[[2]], 
     point2[[2]]}];
  Print[a];
  Print[b];
  Function[x, a x + b]
  ]
lFun1[point1 : pnt, r_Real] := Module[{a, b},
  {a, b} = {r, point1[[2]] - r point1[[1]]};
  Print[a];
  Print[b];
  Function[x, r x + b]
  ]


(* MAP DATA TO GL *)
toGL[data_] := 
 data /. fig : edg :> {fig[[1, 2, 1]], Opacity[fig[[1, 2, 2]]], Which[
     fig[[1, 1, 1]] == POINT, 
     	{PointSize[fig[[2, 3 ]]], 
     	Point[fig[[1, 1, 2]]]},
     fig[[1, 1, 1]] == POLYGON, 
     	{EdgeForm[{fig[[2, 1]], Opacity[fig[[2, 2]]], Thickness[fig[[2, 3]]], Dashing[{fig[[2, 4]], fig[[2, 5]]}]}], 
     	Polygon[fig[[1, 1, 2]]]},
     fig[[1, 1, 1]] == CIRCLE, 
     	{Thickness[fig[[2, 3]]], Dashing[{fig[[2, 4]], fig[[2, 5]]}], 
     	Circle[fig[[1, 1, 2, 1]], EuclideanDistance[fig[[1, 1, 2, 1]], fig[[1, 1, 2, 2]]]]},
     fig[[1, 1, 1]] == LINE, 
     	{fig[[2, 1]], Opacity[fig[[2, 2]]], Thickness[fig[[2, 3]]], Dashing[{fig[[2, 4]], fig[[2, 5]]}], 
     	Line[fig[[1, 1, 2]]]},
     fig[[1, 1, 1]] == DISK, 
     	{EdgeForm[{fig[[2, 1]], Opacity[fig[[2, 2]]], Thickness[fig[[2, 3]]], Dashing[{fig[[2, 4]], fig[[2, 5 ]]}]}], 
      	Disk[fig[[1, 1, 2, 1]], EuclideanDistance[fig[[1, 1, 2, 1]], fig[[1, 1, 2, 2]]]]},
     fig[[1, 1, 1]] == FILLEDCURVE, 
     	{EdgeForm[Directive[{fig[[2, 1]], Opacity[fig[[2, 2]]], Thickness[fig[[2, 3]]], Dashing[{fig[[2, 4]], fig[[2, 5]]}]}]], 
      	FilledCurve[BezierCurve[fig[[1, 1, 2]]]]},
     fig[[1, 1, 1]] == HYPERBOLICLINE, 
     	Map[With[{tmp = hLine32[#]}, ParametricPlot[{tmp[[2]] Cos[t] + tmp[[1, 1]], tmp[[2]] Sin[t] + tmp[[1, 2]]}, {t, tmp[[3]], tmp[[4]]}, 
          PlotStyle -> Directive[fig[[2, 1]], Opacity[fig[[2, 2]]], Thickness[fig[[2, 3]]], Dashing[{fig[[2, 4]], fig[[2, 5]]}]]][[1, 1, 3]]] &, 
      	Partition[fig[[1, 1, 2]], 2, 1]]
     ]}
toGL3D[x_] := (toGL[x] /. p : pnt :> Append[p, 0])



(* IMAGE PROCESSING *)
overlayColorsFromImage3[data_, {imgdata_, dw_}] := Module[{getColor3},
  getColor3[poly_ ] :=
   Module[{
     iw = Dimensions[imgdata][[2]],
     ih = Dimensions[imgdata][[1]],
     dh
	},
    dh = (ih/iw) dw;
    Mean[Flatten[
      Table[
       imgdata[[
         ih - Mod[y - 1, ih],
         Mod[x - 1, iw] + 1
         ]],
       {x,
        Ceiling[Min[poly[[All, 1]]]*iw/dw + 1],
        Ceiling[Max[poly[[All, 1]]]*iw/dw]
        },
       {y,
        Ceiling[Min[poly[[All, 2]]]*ih/dh + 1],
        Ceiling[Max[poly[[All, 2]]]*ih/dh]}
       ], 1]]];
  Map[
   Function[k, 
    Fold[setVal, k[[1]], {{{1, 2, 1}, {k[[2]]}}, {{2, 1}, {k[[2]]}},{{2, 2}, {0.000625}}}]],
    
   Map[
    {#, RGBColor[getColor3[#[[1, 1, 2]]]]} &,
    Cases[{data}, edg, Infinity]]]]




(* MINI-APPS *)
chairTiling[fig_] := Module[
  	{nPos, traR, traU, out},
  	traR = 2 Max[Cases[fig, pnt, Infinity][[All, 1]]];
  	traU = 2 Max[Cases[fig, pnt, Infinity][[All, 2]]];
  	nPos = {Max[Cases[fig, pnt, Infinity][[All, 1]]]/2, Max[Cases[fig, pnt, Infinity][[All, 1]]]/2};
	out = {
	fig, 
	tra3[fig, nPos], 
    tra3[ref3[fig, {{1, 0}, {0, 0}}], {traR, 0}], 
    tra3[ref3[fig, {{0, 1}, {0, 0}}], {0, traU}]
    }
]
chair := {{{1, {{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}}, {Cyan, 1}}, {Black, 1, .005, 2 \[Pi], 2 \[Pi]}}
showChair[n_] := toGL[Nest[chairTiling, chair, n]]




(* GUI-APPS *)

tilingdataset=Dataset[{
<|"tiling"->"3.3.4.3.4","geometry"->FoldList[nGonAtEdge3,newGeo[4],{{1,3},{2,3},{2,4},{2,3},{3,3}}], 
	"lattice"->{{Sqrt[2]+Sqrt[6]/2,-(1/2) Sqrt[2]},{1/2 Sqrt[2],Sqrt[2]+Sqrt[6]/2}}|>,
<|"tiling"->"3.6.3.6","geometry"->FoldList[nGonAtEdge3,newGeo[3],{{3,6},{4,3}}], 
	"lattice"->{{3.`,-1.732050807568877`},{3.`,1.7320508075688774`}}|>,
<|"tiling"->"4.8.8","geometry"->FoldList[nGonAtEdge3,newGeo[4],{{3,8}}], 
	"lattice"->{{2.414213562373095`,2.414213562373095`},{0,4.82842712474619`}}|>,
<|"tiling"->"4.6.12","geometry"->FoldList[nGonAtEdge3,newGeo[12],{{6,4},{2,6},{3,4},{3,6},{3,4}}],
	"lattice"->{{Sqrt[6],0},{Sqrt[3/2],3/Sqrt[2]}}|>,
<|"tiling"->"3.12.12","geometry"->{Fold[nGonAtEdge3,newGeo[12],{{3,3}}],FoldList[nGonAtEdge3,newGeo[12],{{9,3}}]},
	"lattice"->{{(2+2Sqrt[3])/(2Sqrt[2]),0},{(1+Sqrt[3])/(2 Sqrt[2]),(3+Sqrt[3])/(2 Sqrt[2])}}|>,
<|"tiling"->"3.4.6.4","geometry"->FoldList[nGonAtEdge3,newGeo[4],{{3,6},{3,4},{2,3},{2,4},{3,3}}],
	"lattice"->{{-((3+Sqrt[3])/Sqrt[2]),-((1+Sqrt[3])/Sqrt[2])},{0,Sqrt[2]+Sqrt[6]}}|>,
<|"tiling"->"3.3.3.4.4","geometry"->FoldList[nGonAtEdge3,newGeo[4],{{1,3},{2,3}}],
	"lattice"->{{Sqrt[2],0},{Sqrt[2]/2,Sqrt[3/2]+Sqrt[2]}}|>,
<|"tiling"->"3.3.3.3.6","geometry"->FoldList[nGonAtEdge3,newGeo[6],{{4,3},{3,3},{3,3},{3,3},{2,3},{3,3},{3,3},{2,3}}],
	"lattice"->{{(3Sqrt[3])/2,1/2},{Sqrt[3]/2,2.5}}|>,
<|"tiling"->"3","geometry"->FoldList[nGonAtEdge3,newGeo[3],{{1,3},{2,3},{2,3},{2,3},{2,3}}],
	"lattice"->{{3,0},{3/2,3/2 Sqrt[3]}}|>,
<|"tiling"->"4","geometry"->newGeo[4],
	"lattice"->{{Sqrt[2],0},{0,Sqrt[2]}}|>,
<|"tiling"->"6","geometry"->newGeo[6],
	"lattice"->{{Sqrt[3],0},{Sqrt[3]/2,3/2}}|>,
<|"tiling"->"L4.8.8","geometry"->{
	{{{1,{{0,0},{1.30656,0},{1.30656,1.30656}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
	{{{1,{{0,0},{1.30656,1.30656},{0,1.30656}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
	{{{1,{{0,0},{0,1.30656},{-1.30656,1.30656}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
	{{{1,{{0,0},{-1.30656,1.30656},{-1.30656,0}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
	{{{1,{{0,0},{-1.30656,0},{-1.30656,-1.30656}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
	{{{1,{{0,0},{-1.30656,-1.30656},{0,-1.30656}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
	{{{1,{{0,0},{0,-1.30656},{1.30656,-1.30656}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
	{{{1,{{0,0},{1.30656,-1.30656},{1.30656,0}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}}
	},
	"lattice"->{{2 1.3065629648763766`,0},{0,2 1.3065629648763766`}}|>,
<|"tiling"->"L3.6.3.6","geometry"->{{{{1,{{0,0},{1,0},{1/2,Sqrt[3]/2},{-(1/2),Sqrt[3]/2}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
		{{{1,{{0,0},{-(1/2),Sqrt[3]/2},{-1,0},{-1/2,-Sqrt[3]/2}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}},
		{{{1,{{0,0},{-1/2,-Sqrt[3]/2},{1/2,-Sqrt[3]/2},{1,0}}},{RGBColor[0, 1, 1],1}},{GrayLevel[0],1,0.005`,2 \[Pi],2 \[Pi]}}},
	"lattice"->{{3/2,Sqrt[3]/2},{3/2,-(Sqrt[3]/2)}}|>
}];
tilingData[a_,b_]:=Select[tilingdataset,#[[Key["tiling"]]]==a&][1,b]//Normal

showTiling[]:=DynamicModule[{
	col = Array[White &, 9], 
  	edg1 = Array[{GrayLevel[0], 1, 0.0005`, 2 \[Pi], 2 \[Pi]} &, 9]
	},
 Manipulate[
  rot3[
     sca2pnt[
      tra3[
         setEdges[
          repCol[
           tilingData[t, "geometry"],
           col],
          edg1],
         #] & /@ lattice[{tilingData[t, "lattice"], n, n}],
      {{.25, .25}, .25}],
     {rot, {0, 0}}] // toGL // gr,
  {t, tilingdataset[All, "tiling"] // Normal},
  {n, 1, 12, 1},
  {i, 1, 9, 1},
  {icol, White},
  {ecol, Black},
  {{ethi, 0.0005}, 0, .025, .0001},
  {{edas, 1}, 1, 1024, 1},
  {rot, 0, 2 \[Pi], .1309},
  Button["Update", col[[i]] = icol; 
   edg1[[i]] = {ecol, 1, ethi, 2 /edas \[Pi], 2 /(2 edas) \[Pi]}],
  ControlPlacement -> Right]
]


showLines1[] := Manipulate[
  Table[
     Fold[
     setVal,
      newGeo[{3,
        {
         {Sin[n r deg] Sin[m r deg],
          Cos[n r deg] Sin[m r deg]},
         {Sin[n (r + 1) deg ] Sin[m (r + 1) deg], 
          Cos[n (r + 1) deg ] Sin[m (r + 1) deg]}
         }
        }]
      , {{{2, 1}, {col}}, {{2, 3}, {0.0025}}}
      ]
     , {r, 0, k}] // toGL // gr,
  {{n, 3}, 1, 99, 1},
  {{m, 4}, 1, 99, 1},
  {{deg, 37 Degree}, 1 Degree, 360 Degree, 1 Degree},
  {{k, 360}, 1, 720, 1},
  {col, Red}, ControlPlacement -> Right]


showLines2[] := Manipulate[Table[
     Fold[setVal,
      newGeo[{3,
        {
         {
          Sin[n r deg] Cos[r deg],
          Sin[n r deg] Sin[r deg]},
         {
          Sin[n (r + 1) deg] Cos[(r + 1) deg],
          Sin[n (r + 1) deg] Sin[(r + 1) deg]}
         }
        }
       ], {{{2, 1}, {col}}, {{2, 3}, {0.0025}}}
      ],
     {r, 0, k}] // toGL // gr,
  {{n, 3/2}, 1, 99, .5},
  {{deg, 222 Degree}, 1 Degree, 360 Degree, 1 Degree},
  {{k, 180}, 1, 720, 1},
  {col, Gray}, ControlPlacement -> Right]


(* UTILITIES *)
tf:=TableForm
mf:=MatrixForm
grids[min_,max_]:=Join[Range[Ceiling[min],Floor[max]],Table[{j+.5,Dashed},{j,Round[min],Round[max-1],1}]]
gr:=Graphics
gr1[x_]:=Graphics[x, Axes->True,AxesOrigin->{0,0},GridLines->Automatic, AxesLabel -> {"x", "y"}]
gr2[x_]:=Graphics[x, Axes->True,AxesOrigin->{0,0},GridLines->grids, Background->Lighter[LightGray], AxesLabel -> {"x", "y"}]
gr3d[x_]:=Graphics3D[x, Boxed->False]
grD[x_] := 
 Graphics[{x, 
   gDisk[FaceColor -> Cyan, FaceOpacity -> 0.25, 
     Thickness -> 0.0025] // toGL}, Axes -> True, 
  AxesOrigin -> {0, 0}, GridLines -> Automatic]
si:=Simplify
im:=IdentityMatrix


sym2rot180::usage=
	"sym2rot180[fig, p]."
sym4rot90::usage=
	"sym4rot90[fig, p]."

hLine32::usage =
    "hLine32[{p1_,p2_}]."

colfig2Hypcoledges::usage =
	"colfig2Hypcoledges[data]."

g::usage =
    "g[a_]."

Begin["`Private`"]
(* Implementation of the package *)

(*Needs["Geometry`Datasets`"]*)
g[a_]:=2*fun[a]


(* Symmetries *)
sym2rot180[fig_,p_]:={fig,rot3[fig,{\[Pi],p}]}
sym4rot90[fig_,p_]:={fig,rot3[fig,{\[Pi]/2,p}],rot3[fig,{\[Pi],p}],rot3[fig,{3\[Pi]/2,p}]}

(* Parameters for hyperbolic line on Poincare disk *)
hLine32[{p1_,p2_}]:=Module[
{sol,center,radius,ext1, ext2, s1, s2,alfa},
sol=LinearSolve[{
{2p1[[1]],2p1[[2]],-1},
{2p2[[1]],2p2[[2]],-1},
{2 p1[[1]]/(p1[[1]]^2+p1[[2]]^2),2 p1[[2]]/(p1[[1]]^2+p1[[2]]^2),-1}},
{p1[[1]]^2+p1[[2]]^2,p2[[1]]^2+p2[[2]]^2,1/(p1[[1]]^2+p1[[2]]^2)}];
center={sol[[1]],sol[[2]]};
alfa=sol[[1]]+ I sol[[2]];
radius=Sqrt[sol[[1]]^2+sol[[2]]^2-sol[[3]]];
ext1=ArcTan[((p1-center)/radius)[[1]],((p1-center)/radius)[[2]]];
ext2=ArcTan[((p2-center)/radius)[[1]],((p2-center)/radius)[[2]]];
s1=If[ext2-ext1<\[Pi],ext1,ext1+2\[Pi]];
s2=If[ext1-ext2<\[Pi],ext2,ext2+2\[Pi]];
{center,radius,s1,s2,alfa}	
]

(* geometry to edges *)
colfig2Coledges[data_, type_]:=
   Flatten[Map[Function[k,{{type,k},#[[2]],#[[3]],#[[4]],#[[5]],#[[6]],#[[7]],#[[8]]}],Partition[#[[1,2]],2,1,{1,1}]]&/@Cases[{data},colfig2,\[Infinity]],1]
colfig2Coledges[data_]:=colfig2Coledges[data, 3]
colfig2Hypcoledges[data_]:=colfig2Coledges[data, 6]

End[]

EndPackage[]
