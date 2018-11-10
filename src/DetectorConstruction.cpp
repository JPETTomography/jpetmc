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
:  G4VUserDetectorConstruction(), fRunNumber(0),  fLoadCADFrame(false), fLoadWrapping(true)
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
     worldSolid   = new G4Box("world", world_hx, world_hy, world_hz);
     worldLogical  = new G4LogicalVolume(worldSolid,air,"worldLogical");  
     worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);                
     
    // scintillators for standard setup; right now always loaded
     ConstructScintillators();

     if(fLoadCADFrame)
     {
       // ConstructFrameCAD();
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
   G4RotationMatrix rot = G4RotationMatrix();

   G4double z[] =      {-31.0*cm, -29.0*cm, -28.9*cm, -28*cm, -27.9*cm, 27.9*cm, 28*cm, 28.9*cm, 29*cm,  31*cm }; 
   G4double rInner[] = { 0*cm,    0*cm,     0*cm,      0*cm,   12.2*cm, 12.2*cm, 0*cm,  0*cm,    0*cm,   0*cm };
   G4double rOuter[] = { 12.7*cm, 12.7*cm,  15*cm,    15*cm,   12.5*cm, 12.5*cm, 15*cm, 15*cm,   12.7*cm, 12.7*cm }; 

   G4Polycone* bigChamber = new G4Polycone("bigChamber",0*degree,360*degree, 10 , z, rInner, rOuter);
        

   G4LogicalVolume * bigChamber_logical = new G4LogicalVolume(bigChamber, bigChamberRun6Material, "bigChamber_logical");

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

    G4Tubs* ringInner = new G4Tubs("ringInner",24*mm,31*mm,1.78*mm,0*degree,360*degree);


    G4Box* conn = new G4Box("conn",44.5*mm,4.5*mm,1.125*mm);
    G4LogicalVolume* conn_logical = new G4LogicalVolume(conn,bigChamberRun6Material,"conn_logical");
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

    G4LogicalVolume* unionSolid_logical = new G4LogicalVolume(unionSolid,bigChamberRun6Material,"union_logical");
    unionSolid_logical->SetVisAttributes(DetVisAtt);

    new G4PVPlacement(transform,             //rotation,position
                       unionSolid_logical,            //its logical volume
                       "bigChamberInnerStructure",             //its name
                       worldLogical,      //its mother (logical) volume
                       true,                 //no boolean operation
                       0,                 //copy number
                       checkOverlaps);       // checking overlaps 




// ----------------------------------------------------------------------------------

    G4Tubs* xadFilling = new G4Tubs("xadFilling",12.*cm, 12.19*cm, 27.5*cm, 0*degree,360*degree);
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
    G4Tubs* kaptonFilling = new G4Tubs("kaptonFilling",0.*cm, 23.99*mm, 0.1*cm, 0*degree,360*degree);
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
   G4RotationMatrix rot = G4RotationMatrix();

   G4double z[] = {-7.6*cm, -7.0*cm, -6.9*cm, -4.3*cm, -4.2*cm, -2.7*cm, -2.6*cm, 2.6*cm, 2.7*cm, 4.2*cm, 4.3*cm, 6.9*cm, 7.0*cm, 7.6*cm }; 
   G4double rInner[] = { 0*cm, 0*cm, 1.5*cm, 1.5*cm, 1.5*cm, 1.5*cm , 1.5*cm,  1.5*cm, 1.5*cm, 1.5*cm, 1.5*cm , 1.5*cm,  0*cm, 0*cm};
   G4double rOuter[] = { 2.5*cm, 2.5*cm, 2.0*cm, 2.0*cm, 1.8*cm, 1.8*cm, 1.57*cm, 1.57*cm, 1.8*cm, 1.8*cm, 2.0*cm, 2.0*cm, 2.5*cm, 2.5*cm }; 

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


    G4Tubs* xadFilling = new G4Tubs("xadFilling",0*cm, 1.49*cm, 0.6*cm, 0*degree,360*degree);
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
   G4RotationMatrix rot = G4RotationMatrix();

   G4double z[] = {-37*cm, -32.61*cm,-32.6*cm, -31.1*cm, -31*cm, 31*cm, 31.1*cm, 32.6*cm, 32.61*cm, 37*cm}; 
   G4double rInner[] = { 0*cm, 0*cm, 0*cm, 0*cm, 7.1*cm, 7.1*cm, 0*cm, 0*cm, 0*cm, 0*cm};
   G4double rOuter[] = { 3*cm, 3*cm, 10*cm, 10*cm, 7.5*cm, 7.5*cm, 10*cm, 10*cm, 3*cm, 3*cm}; 

   // wartości bazujące na wykresach technicznych - nie odpowiadają rzeczywistości 
   //G4double rInner[] = { 0*cm, 0*cm, 0*cm, 0*cm, 9.1*cm, 9.1*cm, 0*cm, 0*cm, 0*cm, 0*cm};
   //G4double rOuter[] = { 3*cm, 3*cm, 12*cm, 12*cm, 9.5*cm, 9.5*cm, 12*cm, 12*cm, 3*cm, 3*cm}; 

   G4Polycone* bigChamber = new G4Polycone("bigChamber",0*degree,360*degree, 10 , z, rInner, rOuter);
        

   G4LogicalVolume * bigChamber_logical = new G4LogicalVolume(bigChamber, bigChamberRun3Material, "bigChamber_logical");

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
    G4LogicalVolume* conn_logical = new G4LogicalVolume(conn,bigChamberRun3Material,"conn_logical");
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

    //G4Tubs* ringOuter = new G4Tubs("ringOuter",80*mm,95*mm,0.8*mm,0*degree,360*degree);
    G4Tubs* ringOuter = new G4Tubs("ringOuter",60*mm,70*mm,0.8*mm,0*degree,360*degree);
    unionSolid =  new G4UnionSolid("c5", unionSolid,ringOuter);


    G4LogicalVolume* unionSolid_logical = new G4LogicalVolume(unionSolid,bigChamberRun3Material,"union_logical");
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
    G4Box* scinBox = new G4Box("scinBox", scinDim_x/2.0 ,scinDim_y/2.0 , scinDim_z/2.0 );
    scinLog = new G4LogicalVolume(scinBox, scinMaterial , "scinLogical");
    //G4VisAttributes* BoxVisAtt =  new G4VisAttributes(G4Colour(0.3,0.4,.9));
    G4VisAttributes* BoxVisAtt =  new G4VisAttributes(G4Colour(0.447059,0.623529,0.811765));
    BoxVisAtt->SetForceWireframe(true);
    BoxVisAtt->SetForceSolid(true);
    scinLog->SetVisAttributes(BoxVisAtt);

    G4Box* scinBoxFree = new G4Box("scinBoxFree", scinDim_x/2.0+wrappingShift ,scinDim_y/2.0+wrappingShift ,
            scinDim_z/2.0 );
    G4Box* wrappingBox = new G4Box("wrappingBox", scinDim_x/2.0+wrappingThickness,
            scinDim_y/2.0+wrappingThickness , scinDim_z/2.0-1*cm );
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

    //air = G4Material::GetMaterial("G4_AIR"); 
    //scinMaterial      =G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    //detectorMaterial  =G4Material::GetMaterial("G4_Al");
    //kapton  =G4Material::GetMaterial("G4_KAPTON");

    air =  new MaterialExtension("air",G4Material::GetMaterial("G4_AIR")); 
    vacuum =  new MaterialExtension("vacuum", G4Material::GetMaterial("G4_Galactic")); 
    scinMaterial      = new MaterialExtension("scinMaterial", G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE"));
    detectorMaterial  = new MaterialExtension("detectorMaterial",G4Material::GetMaterial("G4_Al"));
    kapton =  new MaterialExtension("kapton", G4Material::GetMaterial("G4_KAPTON")); 


    bigChamberRun3Material = new MaterialExtension("bigChamberRun3", G4Material::GetMaterial("G4_Al"));
    bigChamberRun3Material->AllowsAnnihilations(true); 
    bigChamberRun3Material->Set3gProbability(foPsProbabilityAl); 
    bigChamberRun3Material->SetoPsLifetime(fTauoPsAl); 

    bigChamberRun6Material = new MaterialExtension("bigChamberRun6", G4Material::GetMaterial("G4_PLEXIGLASS"));
    bigChamberRun6Material->AllowsAnnihilations(true); 
    bigChamberRun6Material->Set3gProbability(foPsProbabilityAl); 
    bigChamberRun6Material->SetoPsLifetime(fTauoPsAl); 



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

void DetectorConstruction::ConstructSDandField()
{
        if(!detectorSD.Get()){
        DetectorSD* det = new DetectorSD("/mydet/detector");
        detectorSD.Put(det);
        }
        G4SDManager::GetSDMpointer()->AddNewDetector(detectorSD.Get());
        SetSensitiveDetector(scinLog,detectorSD.Get());

}
