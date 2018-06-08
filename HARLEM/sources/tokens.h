/*! \file tokens.h
 
   Interger constants for  processing RASMOL-like text commands in HARLEM 
 
   derived from 
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6
 
   \author Igor Kurnikov
   \date 1997-2002

*/

/* Lexeme Tokens */
const int IdentTok    =  256;
const int NumberTok   =  257;
const int FloatTok    =  258;
const int StringTok   =  259;

/* Command Tokens */
const int AdviseTok   =  260;
const int BackboneTok =  261;
const int CartoonTok  =  262;
const int CentreTok    = 263;
const int ClipboardTok = 264;
const int ColourTok   =  265;
const int ConnectTok  =  266;
const int DashTok     =  267;
const int DefineTok   =  268;
const int DisplayTok  =  269;
const int EchoTok     =  270;
const int ExitTok     =  271;
const int HelpTok     =  272;
const int LabelTok    =  273;
const int LoadTok     =  274;
const int MonitorTok  =  275;
const int PrintTok    =  276;
const int QuitTok     =  277;
const int RefreshTok  =  278;
const int RenumTok    =  279;
const int ResetTok    =  280;
const int ResizeTok   =  281;
const int RestrictTok =  282;
const int RotateTok   =  283;
const int SaveTok     =  284;
const int ScriptTok   =  285;
const int SelectTok   =  286;
const int SetTok      =  287;
const int ShowTok     =  288;
const int SlabTok     =  289;
const int SourceTok   =  290;
const int SpacefillTok = 291;
const int StructureTok = 292;
const int SymmetryTok  = 293;
const int TraceTok     = 294;
const int TranslateTok = 295;
const int WaitTok      = 296;
const int WireframeTok = 297;
//const int WriteTok     = 298;
const int ZapTok       = 299;
const int ZoomTok      = 300;

/* Predicate Tokens */
#define IsPredTok(x)   (((x)>=310) && ((x)<=348))
#define PredTokOrd(x)  ((x)-310)
#define PredTokChr(x)  ((x)+310)

const int  AlphaTok     = 310;
const int  AminoTok     = 311;
const int  ATTok        = 312;
const int  BondedTok    = 313;
const int  CGTok        = 314;
const int  CystineTok   = 315;
const int  DNATok       = 316;
const int  HelixTok     = 317;
const int  HeteroTok    = 318;
const int  HydrogenTok  = 319;
const int  IonTok       = 320;
const int  LigandTok    = 321;
const int  MainChainTok = 322;
const int  NucleicTok   = 323;
const int  ProteinTok   = 324;
const int  PurineTok    = 325;
const int  PyrimidineTok =326;
const int  RNATok       = 327;
const int  SelectedTok  = 328;
const int  SheetTok     = 329;
const int  SidechainTok = 330;
const int  SolventTok   = 331;
const int  TurnTok      = 332;
const int  WaterTok     = 333;

const int  AcidicTok    = 334;
const int  AcyclicTok   = 335;
const int  AliphaticTok = 336;
const int  AromaticTok  = 337;
const int  BasicTok     = 338;
const int  BuriedTok    = 339;
const int  ChargedTok   = 340;
const int  CyclicTok    = 341;
const int  HydrophobicTok = 342;
const int  LargeTok     = 343;
const int  MediumTok    = 344;
const int  NeutralTok   = 345;
const int  PolarTok     = 346;
const int  SmallTok     = 347;
const int  SurfaceTok   = 348;


/* Property Tokens */
#define  IsPropTok(x)   (((x)>=350) && ((x)<=356))
const int  TemperatureTok = 350;
const int  RadiusTok    = 351;
const int  AtomNoTok    = 352;
const int  ElemNoTok    = 353;
const int  ModelTok     = 354;
const int  ResNoTok     = 355;
const int  ChemGrpTok   = 356;

/* File Format Tokens */
/* Warning! Tokens are related to Format values */
#define IsMoleculeFormat(x)  (((x)>=360) && ((x)<=377))

const int  PDBTok       = 360;
const int  MacroModelTok =361;
const int  GaussianTok  = 362;
const int  AlchemyTok   = 363;
const int  NMRPDBTok    = 364;
const int  CharmmTok    = 365;
const int  BiosymTok    = 366;
const int  MOPACTok     = 367;
const int  SHELXTok     = 368;
const int  Mol2Tok      = 369;
const int  FDATTok      = 370;
const int  MMDBTok      = 371;
const int  MDLTok       = 372;
const int  XYZTok       = 373;
const int  CIFTok       = 374;
const int  CEXTok       = 375;
const int  RWFTok       = 376;
const int  HarlemTok    = 377;

/* Raster Tokens */
#define IsImageFormat(x) (((x)>=1100) && ((x)<=1113))
const int  GIFTok       = 1100;
const int  PPMTok       = 1101;
const int  SUNTok       = 1102;
const int  SUNRLETok    = 1103;
const int  EPSFTok      = 1104;
const int  PICTTok      = 1105;
const int  IRISTok      = 1106;
const int  BMPTok       = 1107;
const int  MonoPSTok    = 1108;
const int  VectPSTok    = 1109;
const int  JPEGTok      = 1110;
const int  MolScriptTok = 1111;
const int  POVRayTok    = 1112;
const int  VRMLTok      = 1113;

/* Feature Tokens */
const int  AtomTok      = 390;
const int  BondTok      = 391;
const int  DotsTok      = 392;
const int  HBondTok     = 393;
const int  RibbonTok    = 394;
const int  SSBondTok    = 395;
const int  Ribbon1Tok   = 396;
const int  Ribbon2Tok   = 397;

/* Expression Tokens */
const int  TrueTok      = 400;
const int  FalseTok     = 401;
const int  AllTok       = 402;
const int  NoneTok      = 403;
const int  AndTok       = 404;
const int  OrTok        = 405;
const int  NotTok       = 406;
const int  WithinTok    = 407;
const int  XorTok       = 408;

/* Colour Tokens */
const int  BlueTok      = 410;
const int  BlueTintTok  = 411;
const int  BlackTok     = 412;
const int  BrownTok     = 413;
const int  CyanTok      = 414;
const int  GoldTok      = 415;
const int  GrayTok      = 416;
const int  GreenTok     = 417;
const int  GreenblueTok = 418;
const int  GreenTintTok = 419;
const int  HotPinkTok   = 420;
const int  MagentaTok   = 421;
const int  OrangeTok    = 422;
const int  PinkTok      = 423;
const int  PinkTintTok  = 424;
const int  PurpleTok    = 425;
const int  RedTok       = 426;
const int  RedorangeTok = 427;
const int  SeaTok       = 428;
const int  SkyTok       = 429;
const int  VioletTok    = 430;
const int  WhiteTok     = 431;
const int  YellowTok    = 432;
const int  YellowTintTok =433;

const int  CPKTok       = 434;
const int  ShapelyTok   = 435;
const int  ResidueTok   = 436;
const int  UserTok      = 437;
const int  GroupTok     = 438;
const int  ChainTok     = 439;
const int  TypeTok      = 440;
const int  PotentialTok = 441;
const int  ChargeTok    = 442;

/* Variable Tokens */
const int  AmbientTok   = 450;
const int  AxesTok      = 451;
const int  BackFadeTok  = 452;
const int  BackgroundTok = 453;
const int  BondModeTok  = 454;
const int  BoundBoxTok  = 455;
const int  DepthCueTok  = 456;
const int  FontSizeTok  = 457;
const int  HourGlassTok = 458;
const int  MenusTok     = 459;
const int  MouseTok     = 460;
const int  PickingTok   = 461;
const int  ShadowTok    = 462;
const int  SlabModeTok  = 463;
const int  SpecularTok  = 464;
const int  SpecPowerTok = 465;
const int  StrandsTok   = 466;
const int  TransparentTok = 467;
const int  UnitCellTok  = 468;

/* SlabMode Tokens */
const int  RejectTok    = 470;
const int  HalfTok      = 471;
const int  HollowTok    = 472;
const int  SolidTok     = 473;
const int  SectionTok   = 474;

/* MouseMode Tokens */
const int  RasMolTok    = 475;
const int  InsightTok   = 476;
const int  QuantaTok    = 477;
const int  SybylTok     = 478;

/* Information Tokens */
const int  InfoTok      = 480;
const int  SequenceTok  = 481;
const int  VersionTok   = 482;

/* Display Mode Tokens */
const int  NormalTok    = 485;
const int  StereoTok    = 486;
const int  MonoTok      = 487;
const int  HardwareTok  = 488;

/* Axis Tokens */
const int  XTok         = 490;
const int  YTok         = 491;
const int  ZTok         = 492;

/* Picking Tokens */
const int  IdentifyTok  = 495;
const int  DistanceTok  = 496;
const int  AngleTok     = 497;
const int  TorsionTok   = 498;

/* Misc Tokens */
const int  InLineTok    = 500;
const int  VDWTok       = 501;

/* Db Tokens */
const int  DbTok       = 1210;
const int  KeyTok      = 1220;
const int  ListTok     = 1230;
const int  OpenTok     = 1240;

/* IO Tokens */
const int  ConsoleTok  = 1320;
const int  OutputTok   = 1330;
const int  MolinfoTok  = 1340;
const int  DumpTok     = 1350;

/* Propagators Tokens      */
const int  RotAngleTok = 1400;
const int  ReadTok     = 1410;
const int  GrpContrTok = 1420;
const int  Idx2Tok     = 1430;
const int  Idx4Tok     = 1440;
const int  OperTok     = 1450;
const int  TestTok     = 1460;
const int  DaltonTok   = 1470;

/* ET Tokens   */
const int  ETTok       = 1600;
const int  AcceptorTok = 1610;
const int  DonorTok    = 1620;
const int  ChooseTok   = 1630;
const int  PathwaysTok = 1640;
const int  BestPathTok = 1650;
const int  CouplMapTok = 1660;
const int  ModeTok     = 1670;
const int  RunTok      = 1680;

/* Molecule structure tokens */
const int  GroupsTok   = 1700;

/* Quantum Chemical Tok */
const int  QChemTok    = 1810;
const int  MoleculeTok = 1820;
const int  BasisTok    = 1830;
const int  LocOrbTok   = 1840;
const int  InitTok     = 1850;
const int  MultTok     = 1860;  
const int  NameTok     = 1870;
const int  ForegroundTok = 1880;
const int  InpFileTok  = 1890;
const int  WaveFunTok  = 1900;
const int  MOTok       = 1910;

/* Molecular Manipulation Tokens */
const int  InterMolTok = 2110;
const int  MolMechTok  = 2120;
const int  CalculateTok = 2130;
const int  EnergyTok   = 2140;
const int  CoordTok    = 2150;
const int  UpdateTok   = 2160;
const int  ModuleTok   = 2170;
const int  ElStaticTok = 2180;
const int  SaveAmberParmTok = 2190;
const int  LoadAmberRestartTok = 2200;
const int  CombineBondingAtomsTok = 2210;
const int  AlignOverlapMolTok = 2220;

/* Edit Molecule Tokens */
const int  EditTok     = 2310;
const int  CreateTok   = 2320;
const int  DeleteTok   = 2330;

const int  OverlapTok  = 2410;
const int  OverlapMolTok  = 2420;



