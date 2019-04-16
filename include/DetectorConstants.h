#ifndef DetectorConstants_h 
#define DetectorConstants_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"


namespace detector_constants
{
      const G4double world_size[3] = {1.0*m, 1.0*m, 1.0*m}; ///< max world size

      const G4int layers = 3; ///< number of simulated laters in detector
      const G4double scinDim[3] = {1.9*cm, 0.7*cm, 50.0*cm}; ///<  X dimension of simulated strip

      const G4double radius[layers] = {42.5*cm,46.75*cm,57.5*cm}; ///< layer radius up to center of the strip
      const int nSegments[layers] = {48,48,96}; ///< number of segments in each layer

      const int extraLayers = 2;
      const int nSegmentsExtraLayers[2] = {96,96};
      const G4double radiusExtraLayers[2] = {509*mm, 533*mm};

      const G4double wrappingThickness = 2*25.4*pow(10,-6)*m; ///<  total width of used wrapping 
      const G4double wrappingShift = 1*pow(10,-5)*m; ///<  free space between wrapping and scinitlator 


      // parameters for modular layer (known as 4th layer)
      const G4double scinDim_inModule[3] = {2.5*cm, 0.6*cm, 50.0*cm}; ///<  X dimension of simulated strip

      const int modulesInModularLayer = 24;


}


#endif
