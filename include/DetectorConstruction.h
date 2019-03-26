#ifndef DetectorConstruction_h 
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorSD.h"

#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4NistManager.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SDManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4GeometryManager.hh"

#include "CADMesh.hh"
#include <string>

#include "MaterialExtension.h"

#include "G4Cache.hh"

class DetectorConstructionMessenger;


const G4double world_hx = 1.0*m; ///< max world size
const G4double world_hy = 1.0*m; ///< max world size
const G4double world_hz = 1.0*m; ///< max world size
const G4int layers = 3; ///< number of simulated laters in detector
const G4double scinDim_x = 1.9*cm; ///<  X dimension of simulated strip
const G4double scinDim_y = 0.7*cm; ///<  Y dimension of simulated strip
const G4double scinDim_z = 50.0*cm; ///<  Z dimension of simulated strip

const G4double wrappingThickness = 2*25.4*pow(10,-6)*m; ///<  total width of used wrapping 
const G4double wrappingShift = 1*pow(10,-5)*m; ///<  free space between wrapping and scinitlator 



const G4double radius[layers] = {42.5*cm,46.75*cm,57.5*cm}; ///< layer radius up to center of the strip
const int nSegments[layers] = {48,48,96}; ///< number of segments in each layer
//const int nSegments[layers] = {2,2,2}; ///< number of segments in each layer
const  G4bool checkOverlaps = false; ///< debugging purpose 

const int extraLayers = 2;
const int nSegmentsExtraLayers[2] = {96,96};
const G4double radiusExtraLayers[2] = {509*mm, 533*mm};

/**
* \class DetectorConstruction
* \brief creating detector; can read the CAD geometry 
*/
class DetectorConstruction : public G4VUserDetectorConstruction
{
    public:
    static DetectorConstruction* GetInstance(); ///< only single instance can exist        
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField(); 

    void LoadFrame(G4bool tf){fLoadCADFrame=tf;};  
    void LoadGeometryForRun(G4int nr);

    G4int GetRunNumber(){return fRunNumber;};

    void UpdateGeometry();


    private:
    static G4ThreadLocal G4bool fConstructedSDandField;
    static DetectorConstruction* fInstance;

    /// standard constructor
    DetectorConstruction();
    virtual ~DetectorConstruction();

    DetectorConstructionMessenger* fMessenger;


    /// load from NIST database
    void InitializeMaterials();
    ///  function is loading detector elements from CAD files
    void ConstructFrameCAD();
    //// create only scintillators; dimensions are right now fixed in code
    void ConstructScintillators();
    //// create target used in run3 - big chamber no XAD material inside
    void ConstructTargetRun3();
    //// create target used in run5 - small chamber + XAD material inside 
    void ConstructTargetRun5();
    //// create target used in run6 - big chamber + XAD material inside tube
    void ConstructTargetRun6();

    G4int fRunNumber; ///< corresponds to JPET measurements; run 0 = user setup 

    G4bool fLoadCADFrame; ///< load frame from cad file
    G4bool fLoadWrapping; ///< load scintillator wrapping


    G4Box* worldSolid;
    G4LogicalVolume* worldLogical;
    G4VPhysicalVolume* worldPhysical;

    MaterialExtension* vacuum; 
    MaterialExtension* air; 
    MaterialExtension* scinMaterial;    
    MaterialExtension* detectorMaterial;

    MaterialExtension* plexiglass;
    MaterialExtension* kapton;
    MaterialExtension* aluminiumMaterial;
    MaterialExtension* plexiglassMaterial;
    MaterialExtension* smallChamberMaterial;
    MaterialExtension* XADMaterial;

    G4LogicalVolume * scinLog;
    G4Cache<DetectorSD*>  detectorSD;



};


#endif
