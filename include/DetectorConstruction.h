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



// debugging parameters
const  G4bool checkOverlaps = false; ///< debugging purpose 


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
    /// Modular layer (known as 4th layer); 24 modules filled with scintillators
    void ConstructModularLayer(bool tf){fLoadModularLayer=tf;}
    G4int ReturnNumberOfScintillators();



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
    //// construct modular layer inserted into detector (refered as 4th layer) 
    void ConstructScintillatorsModularLayer();
    //// create target used in run3 - big chamber no XAD material inside
    void ConstructTargetRun3();
    //// create target used in run5 - small chamber + XAD material inside 
    void ConstructTargetRun5();
    //// create target used in run6 - big chamber + XAD material inside tube
    void ConstructTargetRun6();
    G4int fRunNumber; ///< corresponds to JPET measurements; run 0 = user setup 

    G4bool fLoadCADFrame; ///< load frame from cad file
    G4bool fLoadWrapping;
    G4bool fLoadModularLayer;



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
    G4LogicalVolume * scinLogInModule;
    G4Cache<DetectorSD*>  detectorSD;



};


#endif
