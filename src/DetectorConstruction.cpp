#include "DetectorConstruction.h"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4RunManager.hh"
//#include "G4GeometryTolerance.hh"
//#include "G4GDMLParser.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RegionStore.hh"

#include "DetectorConstructionMessenger.h"
#include "DetectorConstants.h"

using namespace detector_constants;

DetectorConstruction* DetectorConstruction::fInstance = 0;

DetectorConstruction* DetectorConstruction::GetInstance()
{
    if (fInstance == 0)
    {
        fInstance = new DetectorConstruction();
    }

    return fInstance;
}


DetectorConstruction::DetectorConstruction()
:  G4VUserDetectorConstruction(), fRunNumber(0),  fLoadCADFrame(false), fLoadWrapping(true), fLoadModularLayer(true)
{

    InitializeMaterials();

    fMessenger = new DetectorConstructionMessenger(this);

}

DetectorConstruction::~DetectorConstruction()
{
     delete fMessenger;
}


void DetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}


G4VPhysicalVolume* DetectorConstruction::Construct() 
{

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();


    // world 
     worldSolid   = new G4Box("world", world_size[0], world_size[1], world_size[2]);
     worldLogical  = new G4LogicalVolume(worldSolid,air,"worldLogical");  
     worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);                
     
    // scintillators for standard setup; right now always loaded
     ConstructScintillators();

     if(fLoadModularLayer){
      ConstructScintillatorsModularLayer();
     }

     if(fLoadCADFrame)
     {
        ConstructFrameCAD();
     }

     if (fRunNumber == 3) {
         ConstructTargetRun3();
     }

     if (fRunNumber == 5) {
         ConstructTargetRun5();
     }

     if (fRunNumber == 6) {
         ConstructTargetRun6();
     }




    return worldPhysical;
}



void DetectorConstruction::LoadGeometryForRun(G4int nr)
{
    fRunNumber = nr;

     if (fRunNumber == 3 ||fRunNumber == 5 ||fRunNumber == 6 || fRunNumber == 0) {
        LoadFrame(true);  
     } else {
         G4Exception ("DetectorConstruction","DC02", FatalException, 
             " This run setup is not implemented ");    
     }

}


void DetectorConstruction::ConstructTargetRun6()
{

    ///  Parameters taken from http://koza.if.uj.edu.pl/petwiki/images/7/71/AnnihilationChamberRun6.pdf
    ///  below are stored main values. Please note that chamber is symmetric in z, therefore 
    ///  only half of the table is presented below.
    /// 
    ///    z       |  r Inner     | r Outer     | Comment
    ///   ---------|--------------|-------------|-------   
    ///    -31.0*cm|    0*cm      | 12.7*cm     |  chamber end
    ///    -29.0*cm|    0*cm      | 12.7*cm     |  start slope before endcup     
    ///    -28.9*cm|    0*cm      | 15*cm       |  stop slope before endcup 
    ///     -28*cm |   0*cm       | 15*cm       |  endcup
    ///   -27.9*cm |   12.2*cm    | 12.5*cm     |  end slope after endcup and start empty inside
    ///

   const double chamber_radius_inner = 12.2*cm;
   const double chamber_radius_outer = 12.5*cm;
   const double chamber_endcup_min = 12.7*cm; 
   const double chamber_endcup_max = 15.0*cm;

   const double xad_z_coverage = 27.5*cm; 
   const double xad_halfthickness = 0.2*cm;

   const double source_holder_connector_halfthickness = 1.125*mm;
   const double source_holder_connector_width = 4.5*mm;
   const double source_holder_connector_height = 44.5*mm;

   const double source_holder_radius_inner = 24*mm;
   const double source_holder_radius_outer = 31*mm;
   const double source_holder_radius_halfthickness = 1.78*mm; 

   const double kapton_foil_halfthickness = 0.1*cm; 


   G4RotationMatrix rot = G4RotationMatrix();

   G4double z[] =      {-31.0*cm, -29.0*cm, -28.9*cm, -28*cm, -27.9*cm, 27.9*cm, 28*cm, 28.9*cm, 29*cm,  31*cm }; 
   G4double rInner[] = { 0*cm,    0*cm,     0*cm,      0*cm,   chamber_radius_inner, chamber_radius_inner, 0*cm,  0*cm,    0*cm,   0*cm };
   G4double rOuter[] = { chamber_endcup_min, chamber_endcup_min, chamber_endcup_max, chamber_endcup_max ,  
       chamber_radius_outer, chamber_radius_outer, chamber_endcup_max, chamber_endcup_max,   chamber_endcup_min, chamber_endcup_min }; 

   G4Polycone* bigChamber = new G4Polycone("bigChamber",0*degree,360*degree, 10 , z, rInner, rOuter);
        

   G4LogicalVolume * bigChamber_logical = new G4LogicalVolume(bigChamber, plexiglassMaterial, "bigChamber_logical");

    G4VisAttributes* DetVisAtt =  new G4VisAttributes(G4Colour(0.9,0.9,.9));
    DetVisAtt->SetForceWireframe(true);
    DetVisAtt->SetForceSolid(true);
    bigChamber_logical->SetVisAttributes(DetVisAtt);


     G4ThreeVector loc = G4ThreeVector(0.0,0.0,0.0);
     G4Transform3D transform(rot,loc);
     new G4PVPlacement(transform,             //rotation,position
                       bigChamber_logical,            //its logical volume
                       "bigChamberGeom",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 


    G4Tubs* ringInner = new G4Tubs("ringInner",source_holder_radius_inner,source_holder_radius_outer,source_holder_radius_halfthickness,0*degree,360*degree);

    G4Box* conn = new G4Box("conn",source_holder_connector_height,source_holder_connector_width,source_holder_connector_halfthickness);
    G4LogicalVolume* conn_logical = new G4LogicalVolume(conn,plexiglassMaterial,"conn_logical");
    conn_logical->SetVisAttributes(DetVisAtt);

    G4ThreeVector loc2;
    G4Transform3D transform2;

    loc2 = G4ThreeVector(75.*mm,0.0,0.0);
    transform2 = G4Transform3D(rot,loc2);
    G4UnionSolid*  unionSolid =  new G4UnionSolid("c1", ringInner,conn,transform2);

    loc2 = G4ThreeVector(-53.*mm,53.*mm,0.0);
    transform2 = G4Transform3D(rot.rotateZ(-45*degree),loc2);
    unionSolid =  new G4UnionSolid("c2", unionSolid,conn,transform2);

    loc2 = G4ThreeVector(-53.*mm,-53.*mm,0.0);
    transform2 = G4Transform3D(rot.rotateZ(90*degree),loc2);
    unionSolid =  new G4UnionSolid("c3", unionSolid,conn,transform2);

    G4LogicalVolume* unionSolid_logical = new G4LogicalVolume(unionSolid,plexiglassMaterial,"union_logical");
    unionSolid_logical->SetVisAttributes(DetVisAtt);

    new G4PVPlacement(transform,             //rotation,position
                       unionSolid_logical,            //its logical volume
                       "bigChamberInnerStructure",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 




// ----------------------------------------------------------------------------------

    G4Tubs* xadFilling = new G4Tubs("xadFilling",chamber_radius_inner-0.1*cm-xad_halfthickness, chamber_radius_inner-0.1*cm, xad_z_coverage, 0*degree,360*degree);
    G4LogicalVolume* xadFilling_logical = new G4LogicalVolume(xadFilling,XADMaterial,"xadFilling_logical");
    G4VisAttributes* XADVisAtt =  new G4VisAttributes(G4Colour(0.2,0.3,.5));
    XADVisAtt->SetForceWireframe(true);
    XADVisAtt->SetForceSolid(true);

    xadFilling_logical->SetVisAttributes(XADVisAtt);

     new G4PVPlacement(transform,             //rotation,position
                       xadFilling_logical,            //its logical volume
                       "xadFillingGeom",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 



    // KAPTON foil   
// ----------------------------------------------------------------------------------
    G4Tubs* kaptonFilling = new G4Tubs("kaptonFilling",0.*cm, source_holder_radius_inner-0.01*mm , kapton_foil_halfthickness, 0*degree,360*degree);
    G4LogicalVolume* kaptonFilling_logical = new G4LogicalVolume(kaptonFilling,kapton,"kaptonFilling_logical");
    G4VisAttributes* kaptonVisAtt =  new G4VisAttributes(G4Colour(0.2,0.3,.5));
    kaptonVisAtt->SetForceWireframe(true);
    kaptonVisAtt->SetForceSolid(true);

    kaptonFilling_logical->SetVisAttributes(kaptonVisAtt);

     new G4PVPlacement(transform,             //rotation,position
                       kaptonFilling_logical,            //its logical volume
                       "kaptonFillingGeom",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 


}





void DetectorConstruction::ConstructTargetRun5()
{

    ///  Parameters taken from http://koza.if.uj.edu.pl/petwiki/index.php/Small_annihilation_chamber 
    ///  below are stored main values. Please note that chamber is symmetric in z, therefore 
    ///  only half of the table is presented below.
    ///  Chamber has four different outer radius 
    /// 
    ///    z       |  r Inner     | r Outer    | Comment
    ///   ---------|--------------|------------|-------   
    ///   -7.6*cm  |    0*cm      | 2.5*cm     | chamber end  
    ///   -7.0*cm  |    0*cm      | 2.5*cm     | start slope (R1->R2)
    ///   -6.9*cm  |    1.5*cm    | 2.0*cm     | end slope and start empty space
    ///   -4.3*cm  |    1.5*cm    | 2.0*cm     | start slope (R2->R3)
    ///    -4.2*cm |    1.5*cm    | 1.8*cm     | end slope
    ///    -2.7*cm |    1.5*cm    | 1.8*cm     | start slope (R3->R4)
    ///    -2.6*cm |    1.5*cm    | 1.57*cm    | end slope 
    ///

   const double chamber_radius_inner = 1.5*cm;
   const double chamber_radius_outer_1 = 2.5*cm;
   const double chamber_radius_outer_2 = 2.0*cm;
   const double chamber_radius_outer_3 = 1.8*cm;
   const double chamber_radius_outer_4 = 1.57*cm;

   const double xadFilling_halfthickness =  0.6*cm;


   G4RotationMatrix rot = G4RotationMatrix();

   G4double z[] = {-7.6*cm, -7.0*cm, -6.9*cm, -4.3*cm, -4.2*cm, -2.7*cm, -2.6*cm, 2.6*cm, 2.7*cm, 4.2*cm, 4.3*cm, 6.9*cm, 7.0*cm, 7.6*cm }; 
   G4double rInner[] = { 0*cm, 0*cm, chamber_radius_inner, chamber_radius_inner, chamber_radius_inner, chamber_radius_inner,
       chamber_radius_inner, chamber_radius_inner, chamber_radius_inner, chamber_radius_inner, chamber_radius_inner, chamber_radius_inner,  0*cm, 0*cm}; 
   G4double rOuter[] = {chamber_radius_outer_1, chamber_radius_outer_1, chamber_radius_outer_2, chamber_radius_outer_2, 
       chamber_radius_outer_3,chamber_radius_outer_3, chamber_radius_outer_4, chamber_radius_outer_4,
       chamber_radius_outer_3, chamber_radius_outer_3, chamber_radius_outer_2, chamber_radius_outer_2,
       chamber_radius_outer_1, chamber_radius_outer_1};

   G4Polycone* smallChamber = new G4Polycone("bigChamber",0*degree,360*degree, 14 , z, rInner, rOuter);
        

   G4LogicalVolume * smallChamber_logical = new G4LogicalVolume(smallChamber, smallChamberMaterial, "smallChamber_logical");

    G4VisAttributes* DetVisAtt =  new G4VisAttributes(G4Colour(0.9,0.9,.9));
    DetVisAtt->SetForceWireframe(true);
    DetVisAtt->SetForceSolid(true);
    smallChamber_logical->SetVisAttributes(DetVisAtt);


     G4ThreeVector loc = G4ThreeVector(0.0,0.0,0.0);
     G4Transform3D transform(rot,loc);
     new G4PVPlacement(transform,             //rotation,position
                       smallChamber_logical,            //its logical volume
                       "smallChamberGeom",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 




    G4Tubs* xadFilling = new G4Tubs("xadFilling",0*cm,  chamber_radius_inner-0.01*cm, xadFilling_halfthickness, 0*degree,360*degree);
    G4LogicalVolume* xadFilling_logical = new G4LogicalVolume(xadFilling,XADMaterial,"xadFilling_logical");
    G4VisAttributes* XADVisAtt =  new G4VisAttributes(G4Colour(0.2,0.3,.5));
    XADVisAtt->SetForceWireframe(true);
    XADVisAtt->SetForceSolid(true);

    xadFilling_logical->SetVisAttributes(XADVisAtt);

     new G4PVPlacement(transform,             //rotation,position
                       xadFilling_logical,            //its logical volume
                       "xadFillingGeom",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 


}

void DetectorConstruction::ConstructTargetRun3()
{
    ///  Please note that chamber is symmetric in z, therefore 
    ///  only half of the table is presented below.
    ///  Chamber has four different outer radius 
    /// 
    ///    z         |  r Inner     | r Outer    | Comment
    ///   -----------|--------------|------------|-------   
    ///   -37.00*cm  |    0*cm      | 3*cm       | chamber end  
    ///   -32.61*cm  |    0*cm      | 3*cm       | start endcup 
    ///   -32.60*cm  |    0*cm      | 10*cm      | end slope 
    ///   -31.10*cm  |    0*cm      | 10*cm      | start slope 
    ///   -31.00*cm  |    7.1*cm    | 7.5*cm     | end slope and start empty volume
    ///

   const double chamber_radius_inner = 7.1*cm;
   const double chamber_radius_outer = 7.5*cm;
   const double chamber_endcup_min = 3.0*cm; 
   const double chamber_endcup_max = 10.0*cm;


   G4RotationMatrix rot = G4RotationMatrix();

   G4double z[] = {-37*cm, -32.61*cm,-32.6*cm, -31.1*cm, -31*cm, 31*cm, 31.1*cm, 32.6*cm, 32.61*cm, 37*cm}; 
   G4double rInner[] = { 0*cm, 0*cm, 0*cm, 0*cm, chamber_radius_inner, chamber_radius_inner, 0*cm, 0*cm, 0*cm, 0*cm};
   G4double rOuter[] = { chamber_endcup_min, chamber_endcup_min, chamber_endcup_max, chamber_endcup_max,
       chamber_radius_outer, chamber_radius_outer, chamber_endcup_max, chamber_endcup_max, chamber_endcup_min,chamber_endcup_min};

   G4Polycone* bigChamber = new G4Polycone("bigChamber",0*degree,360*degree, 10 , z, rInner, rOuter);
        

   G4LogicalVolume * bigChamber_logical = new G4LogicalVolume(bigChamber, aluminiumMaterial, "bigChamber_logical");

    G4VisAttributes* DetVisAtt =  new G4VisAttributes(G4Colour(0.9,0.9,.9));
    DetVisAtt->SetForceWireframe(true);
    DetVisAtt->SetForceSolid(true);
    bigChamber_logical->SetVisAttributes(DetVisAtt);


     G4ThreeVector loc = G4ThreeVector(0.0,0.0,0.0);
     G4Transform3D transform(rot,loc);
     new G4PVPlacement(transform,             //rotation,position
                       bigChamber_logical,            //its logical volume
                       "bigChamberGeom",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 



    G4Tubs* ringInner = new G4Tubs("ringInner",15*mm,20.8*mm,0.8*mm,0*degree,360*degree);


    G4Box* conn = new G4Box("conn",25*mm,7.*mm,0.8*mm);
    G4LogicalVolume* conn_logical = new G4LogicalVolume(conn,aluminiumMaterial,"conn_logical");
    conn_logical->SetVisAttributes(DetVisAtt);


    G4ThreeVector loc2;
    G4Transform3D transform2;

    loc2 = G4ThreeVector(39.8*mm,0.0,0.0);
    transform2 = G4Transform3D(rot,loc2);
    G4UnionSolid*  unionSolid =  new G4UnionSolid("c1", ringInner,conn,transform2);

    loc2 = G4ThreeVector(-39.8*mm,0.0,0.0);
    transform2 = G4Transform3D(rot,loc2);
    unionSolid =  new G4UnionSolid("c2", unionSolid,conn,transform2);

    loc2 = G4ThreeVector(0.0,39.8*mm,0.0);
    transform2 = G4Transform3D(rot.rotateZ(90*degree),loc2);
    unionSolid =  new G4UnionSolid("c3", unionSolid,conn,transform2);

    loc2 = G4ThreeVector(0.0,-39.8*mm,0.0);
    transform2 = G4Transform3D(rot,loc2);
    unionSolid =  new G4UnionSolid("c4", unionSolid,conn,transform2);

    G4Tubs* ringOuter = new G4Tubs("ringOuter",60*mm,70*mm,0.8*mm,0*degree,360*degree);
    unionSolid =  new G4UnionSolid("c5", unionSolid,ringOuter);


    G4LogicalVolume* unionSolid_logical = new G4LogicalVolume(unionSolid,aluminiumMaterial,"union_logical");
    unionSolid_logical->SetVisAttributes(DetVisAtt);

    new G4PVPlacement(transform,             //rotation,position
                       unionSolid_logical,            //its logical volume
                       "bigChamberInnerStructure",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 


}

void DetectorConstruction::ConstructScintillators()
{
    // scintillator
    G4Box* scinBox = new G4Box("scinBox", scinDim[0]/2.0 ,scinDim[1]/2.0 , scinDim[2]/2.0 );
    scinLog = new G4LogicalVolume(scinBox, scinMaterial , "scinLogical");
    //G4VisAttributes* BoxVisAtt =  new G4VisAttributes(G4Colour(0.3,0.4,.9));
    G4VisAttributes* BoxVisAtt =  new G4VisAttributes(G4Colour(0.447059,0.623529,0.811765));
    BoxVisAtt->SetForceWireframe(true);
    BoxVisAtt->SetForceSolid(true);
    scinLog->SetVisAttributes(BoxVisAtt);

    G4Box* scinBoxFree = new G4Box("scinBoxFree", scinDim[0]/2.0+wrappingShift ,scinDim[1]/2.0+wrappingShift ,
            scinDim[2]/2.0 );
    G4Box* wrappingBox = new G4Box("wrappingBox", scinDim[0]/2.0+wrappingThickness,
            scinDim[1]/2.0+wrappingThickness , scinDim[2]/2.0-1*cm );
    G4LogicalVolume* wrappingLog; 

    G4VisAttributes* BoxVisAttWrapping =  new G4VisAttributes(G4Colour(0.447059,0.623529,0.811765));
    //G4VisAttributes* BoxVisAttWrapping =  new G4VisAttributes(G4Colour(0.4,0.4,.4));
    BoxVisAttWrapping->SetForceWireframe(true);
    BoxVisAttWrapping->SetForceSolid(true);



    G4int icopy = 1;

    for(int j=0;j<layers;j++){
        for(int i=0;i<nSegments[j];i++){
            G4double phi = i*2*M_PI/nSegments[j];
            G4double fi = M_PI/nSegments[j];


            if( j == 0 ){
                fi =0.;
            }

            G4RotationMatrix rot = G4RotationMatrix();
            rot.rotateZ(phi+fi);

            G4ThreeVector loc = G4ThreeVector(radius[j]*(cos(phi+fi)),radius[j]*(sin(phi+fi)),0.0);
            G4Transform3D transform(rot,loc);

            G4String name = "scin_"+G4UIcommand::ConvertToString(icopy);

            new G4PVPlacement(transform,             //rotation,position
                              scinLog,            //its logical volume
                              name,             //its name
                              worldLogical,      //its mother (logical) volume
                              true,                 //no boolean operation
                              icopy,                 //copy number
                              checkOverlaps);       // checking overlaps 


             if(fLoadWrapping) 
             {
             // wrapping 
             
             G4VSolid* unionSolid =  new G4SubtractionSolid("wrapping", wrappingBox, scinBoxFree);
             wrappingLog = new G4LogicalVolume(unionSolid, kapton , "wrappingLogical");
             wrappingLog->SetVisAttributes(BoxVisAttWrapping);

             G4String nameWrapping = "wrapping_"+G4UIcommand::ConvertToString(icopy);

             new G4PVPlacement(transform,             //rotation,position
                              wrappingLog,            //its logical volume
                              nameWrapping,             //its name
                              worldLogical,      //its mother (logical) volume
                              true,                 //no boolean operation
                              icopy,                 //copy number
                              checkOverlaps);       // checking overlaps 
              }

            icopy++;

        }
    }
}

void DetectorConstruction::ConstructScintillatorsModularLayer()
{

  //4th Layer : S. Sharma 20.06.2018
  G4Box* scinBoxInModule = new G4Box("scinBoxInModule", scinDim_inModule[0]/2.0, scinDim_inModule[1]/2.0, scinDim_inModule[2]/2.0);
  scinLogInModule = new G4LogicalVolume(scinBoxInModule, scinMaterial, "scinBoxInModule");
  G4VisAttributes* BoxVisAttI = new G4VisAttributes(G4Colour(0.105, 0.210, 0.210, 0.9));
  BoxVisAttI->SetForceWireframe(true);
  BoxVisAttI->SetForceSolid(true);
  scinLogInModule->SetVisAttributes(BoxVisAttI);

  // 13 straight scintillators in single module
  // radius of each scintillator in tabular form
  // taken from equation: radius0/cos(j*displacement) where
  // j -> number of scintillator 
  // radius0 = 38.186*cm
  // displacement -> Angular displacement (1.04 degree / 0.01815 radius - fixed; value provided by Sushil)
  const G4double radius_inner[13] = {38.416, 38.346, 38.289, 38.244, 38.212, 38.192,
    38.186, 38.192, 38.212, 38.244, 38.289, 38.346, 38.416};
  G4double phi1 = 0.0;
  G4int icopyI = 193; // sum of already constructed scintillators; 
  //for Framework newly inserted scintillators need to have a unique numbering 

  for(int i=0; i<modulesInModularLayer; i++) {
    G4double phi = (i*2*M_PI/modulesInModularLayer);
    // 13 centered modules
    for(int j=-6; j<7; j++) {
      // 0.01815 - Angular displacement of 1.04 degree
      phi1 = phi + j*0.01815;
      G4double radius1 = radius_inner[j+6]*cm;
      G4RotationMatrix rot = G4RotationMatrix();
      rot.rotateZ(phi);
      G4ThreeVector loc = G4ThreeVector(radius1*cos(phi1), radius1*sin(phi1), 0.0);
      G4Transform3D transform(rot, loc);
      G4String nameNewI = "scin_"+G4UIcommand::ConvertToString(icopyI+i*13+j+6);
      new G4PVPlacement(transform, scinLogInModule, nameNewI, worldLogical, true, icopyI+i*13+j+6, checkOverlaps);
      }
      }
 
}

void DetectorConstruction::InitializeMaterials()
{
    // define material
    G4NistManager* nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_AIR");   
    nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");   
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    nistManager->FindOrBuildMaterial("G4_Al");
    nistManager->FindOrBuildMaterial("G4_KAPTON");
    nistManager->FindOrBuildMaterial("G4_Galactic");
    nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");

    air =  new MaterialExtension("air",G4Material::GetMaterial("G4_AIR")); 
    vacuum =  new MaterialExtension("vacuum", G4Material::GetMaterial("G4_Galactic")); 
    scinMaterial      = new MaterialExtension("scinMaterial", G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE"));
    detectorMaterial  = new MaterialExtension("detectorMaterial",G4Material::GetMaterial("G4_Al"));
    kapton =  new MaterialExtension("kapton", G4Material::GetMaterial("G4_KAPTON")); 


    aluminiumMaterial = new MaterialExtension("bigChamberRun3", G4Material::GetMaterial("G4_Al"));
    aluminiumMaterial->AllowsAnnihilations(true); 
    aluminiumMaterial->Set3gProbability(foPsProbabilityAl); 
    aluminiumMaterial->SetoPsLifetime(fTauoPsAl); 

    plexiglassMaterial = new MaterialExtension("bigChamberRun6", G4Material::GetMaterial("G4_PLEXIGLASS"));
    plexiglassMaterial->AllowsAnnihilations(true); 
    plexiglassMaterial->Set3gProbability(foPsProbabilityAl); 
    plexiglassMaterial->SetoPsLifetime(fTauoPsAl); 


    smallChamberMaterial = new MaterialExtension("smallChamber", G4Material::GetMaterial("G4_Al"));
    smallChamberMaterial->AllowsAnnihilations(true); 
    smallChamberMaterial->Set3gProbability(foPsProbabilityAl); 
    smallChamberMaterial->SetoPsLifetime(fTauoPsAl); 

    //  /// https://www.sigmaaldrich.com/catalog/product/sigma/xad4
    XADMaterial = new MaterialExtension("XAD", G4Material::GetMaterial("G4_POLYSTYRENE"));
    XADMaterial->AllowsAnnihilations(true); 
    XADMaterial->Set3gProbability(foPsProbabilityAl); 
    XADMaterial->SetoPsLifetime(fTauoPsAl); 


}

void DetectorConstruction::ConstructFrameCAD()
{

     // in stl file the scintillator slots were made bigger !!!
     //             Length      width
     //    true      2.1 cm      0.9 cm
     //    used in stl       2.6 cm      1.7 cm

     CADMesh * mesh1 = new CADMesh((char*)"stl_geometry/Frame_JPET.stl" );
     mesh1->SetScale(mm);
     G4VSolid* cad_solid1 = mesh1->TessellatedMesh();

     G4LogicalVolume * cad_logical = new G4LogicalVolume(cad_solid1, detectorMaterial, "cad_logical");
 
     G4VisAttributes* DetVisAtt =  new G4VisAttributes(G4Colour(0.9,0.9,.9));
     DetVisAtt->SetForceWireframe(true);
     DetVisAtt->SetForceSolid(true);
     cad_logical->SetVisAttributes(DetVisAtt);

     G4RotationMatrix rot = G4RotationMatrix();
     rot.rotateY(90*deg);
     G4ThreeVector loc = G4ThreeVector(0*cm, 306.5*cm ,-23*cm);
     G4Transform3D transform(rot,loc);

     new G4PVPlacement(transform,             //rotation,position
                       cad_logical,            //its logical volume
                       "cadGeom",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 

}

G4int DetectorConstruction::ReturnNumberOfScintillators()
{
      if(fLoadModularLayer){
            return 504;
      } else {
            return 192;
      }
}

void DetectorConstruction::ConstructSDandField()
{
      if(!detectorSD.Get()){
        DetectorSD* det = new DetectorSD("/mydet/detector",ReturnNumberOfScintillators());
        detectorSD.Put(det);
      }
        G4SDManager::GetSDMpointer()->AddNewDetector(detectorSD.Get());
        SetSensitiveDetector(scinLog,detectorSD.Get());
        if(fLoadModularLayer) SetSensitiveDetector(scinLogInModule,detectorSD.Get());

}
