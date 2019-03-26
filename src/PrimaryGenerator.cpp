#include "PrimaryGenerator.h"
#include "G4SystemOfUnits.hh"

#include "G4HadPhaseSpaceGenbod.hh"
#include "globals.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh" 


#include "G4TransportationManager.hh" 
#include "G4Navigator.hh"
#include "Randomize.hh"

#include "G4PhysicalConstants.hh"
#include "G4Circle.hh"

#include "VtxInformation.h"
#include "PrimaryParticleInformation.h"

#include "DetectorConstruction.h"


PrimaryGenerator::PrimaryGenerator()
    : G4VPrimaryGenerator()
{
}

PrimaryGenerator::~PrimaryGenerator()
{}

void PrimaryGenerator::GenerateEvtChamberWithSodiumAndPorousMaterial(G4Event* event, G4double maxXhalf, G4double maxYhalf, G4double maxZhalf)
{

    G4ThreeVector vtxPosition; // 2g/3g
    G4ThreeVector promptVtxPosition
         =  VertexUniformInCylinder(1.5*cm,0.2*cm);
    G4double ratio3g;
    G4double lifetime3g; 
    std::tie(vtxPosition, ratio3g, lifetime3g) = GetVerticesDistribution(maxXhalf,maxYhalf,maxZhalf);

    G4PrimaryVertex* vertex;
    vertex = new G4PrimaryVertex();
    VtxInformation* info = new VtxInformation();
    vertex->SetUserInformation(info);
    vertex->SetPosition(vtxPosition.x(),vtxPosition.y(),vtxPosition.z());


    G4PrimaryVertex* vertexPrompt = new G4PrimaryVertex() ;
    VtxInformation* infoPrompt = new VtxInformation();
    vertexPrompt->SetUserInformation(infoPrompt);
    vertexPrompt->SetPosition(promptVtxPosition.x(),promptVtxPosition.y(),promptVtxPosition.z());


    G4double lifetime;
    G4double promptLifetime = G4RandExponential::shoot(3.7*ps);


    if( ratio3g > G4UniformRand() )
    { 
        lifetime = G4RandExponential::shoot(lifetime3g);
        //generate 3g
        GenerateThreeGammaVertex(vertex);
        info->SetThreeGammaGen(true);

    } else {
        //generate 2g
        lifetime =G4RandExponential::shoot(fTauBulk);
        GenerateTwoGammaVertex(vertex);
        info->SetTwoGammaGen(true);
    }
    info->SetLifetime(lifetime/ps);
    info->SetVtxPosition(vtxPosition.x()/cm,vtxPosition.y()/cm,vtxPosition.z()/cm);
    event->AddPrimaryVertex(vertex);


    // add sodium prompt gamma
    GeneratePromptGammaSodium(vertexPrompt);
    infoPrompt->SetPromptGammaGen(true);
    infoPrompt->SetLifetime(promptLifetime);
    infoPrompt->SetVtxPosition(promptVtxPosition.x()/cm,promptVtxPosition.y()/cm,promptVtxPosition.z()/cm);

    event->AddPrimaryVertex(vertexPrompt);

}


void PrimaryGenerator::GenerateBeam(BeamParams* beamParams, G4Event* event)
{

      G4ThreeVector vtxCoor = beamParams->GetVtx();
    G4PrimaryVertex* vertex = new G4PrimaryVertex(vtxCoor,0);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefinition = particleTable->FindParticle("gamma");

    const G4double ene = beamParams->GetEne();
    G4ThreeVector momentum = beamParams->GetMomentum();


    G4double px = ene * momentum.x();
    G4double py = ene * momentum.y();
    G4double pz = ene * momentum.z();

    G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition,
            px, py, pz,ene);

    vertex->SetPrimary(particle1);
    VtxInformation* infoPrompt = new VtxInformation();
    vertex->SetUserInformation(infoPrompt);
    infoPrompt->SetPromptGammaGen(true);
    infoPrompt->SetVtxPosition(vtxCoor.x()/cm,vtxCoor.y()/cm,vtxCoor.z()/cm);

    event->AddPrimaryVertex(vertex);
}


void PrimaryGenerator::GenerateIsotope(SourceParams* sourceParams, G4Event* event)
{

    G4ThreeVector vtxPosition; 

    if (sourceParams->GetShape() == "cylinder")
    {
        vtxPosition =  VertexUniformInCylinder(sourceParams->GetShapeDim(0),sourceParams->GetShapeDim(1))+sourceParams->GetShapeCenterPosition();
    }

    G4PrimaryVertex* vertex = new G4PrimaryVertex(vtxPosition,0);
    VtxInformation* info = new VtxInformation();
    vertex->SetUserInformation(info);
    vertex->SetPosition(vtxPosition.x()*cm,vtxPosition.y()*cm,vtxPosition.z()*cm);
    //G4cout << vtxPosition.x() << " cm "  << vtxPosition.x()*cm << " *cm " << vtxPosition.x()/cm << " /cm \n";


    G4double lifetime = G4RandExponential::shoot(fTauBulk);


    if ( sourceParams->GetGammasNumber() == 1 )
    { 
        GeneratePromptGammaSodium(vertex);
        info->SetPromptGammaGen(true);
        info->SetLifetime(lifetime);

    } else if ( sourceParams->GetGammasNumber() == 2 )   {
        //generate 2g
        GenerateTwoGammaVertex(vertex);
        info->SetTwoGammaGen(true);
    } else if ( sourceParams->GetGammasNumber() == 3 )   {
        //generate 3g
        GenerateThreeGammaVertex(vertex);
        info->SetThreeGammaGen(true);
    } else {
        G4Exception("PrimaryGenerator","PG01",FatalException,
                " program does not know how many gamma quanta need to be simulated ");
    }

    info->SetLifetime(lifetime/ps);
    info->SetVtxPosition(vtxPosition.x(),vtxPosition.y(),vtxPosition.z());
    event->AddPrimaryVertex(vertex);

}

void PrimaryGenerator::GenerateNema(G4int nemaPoint, G4Event* event)
{

    /** NEMA positions - coordinate system centered in the middle of the detector
     * for( z={0,3/4L}){
     *  x = 1, 10, 20 cm
     *  y =0 
     *  }
     *
     *  20      3       6
     *  10      2       5
     *  1       1       4
     *  z ------0------3/4L ------
     */

	G4double x_creation = 0.0; 
	G4double z_creation = 0.0;

	    if(nemaPoint>3){
	        z_creation = z_creation + scinDim_z*3/8/cm;
	    }
	
	    if(nemaPoint==1 || nemaPoint == 4){
	        x_creation = x_creation + 1.0;
	    }
	
	    if(nemaPoint==2 || nemaPoint == 5){
	        x_creation = x_creation + 10.0;
	    }
	    if(nemaPoint==3 || nemaPoint == 6){
	        x_creation = x_creation + 20.0;
	    }



    G4ThreeVector vtxPosition
      =  VertexUniformInCylinder(0.1*mm,0.1*mm)+G4ThreeVector(x_creation,0.,z_creation);


    G4PrimaryVertex* vertex = new G4PrimaryVertex(vtxPosition,0);
    VtxInformation* info = new VtxInformation();
    vertex->SetUserInformation(info);
    vertex->SetPosition(vtxPosition.x()*cm,vtxPosition.y()*cm,vtxPosition.z()*cm);

    G4double lifetime = G4RandExponential::shoot(fTauBulk);
    GenerateTwoGammaVertex(vertex);
    info->SetTwoGammaGen(true);
    info->SetLifetime(lifetime/ps);
    info->SetVtxPosition(vtxPosition.x(),vtxPosition.y(),vtxPosition.z());
    event->AddPrimaryVertex(vertex);



    // add sodium prompt gamma
    G4double promptLifetime = G4RandExponential::shoot(3.7*ps);
    G4PrimaryVertex* vertexPrompt = new G4PrimaryVertex(vtxPosition,promptLifetime) ;
    VtxInformation* infoPrompt = new VtxInformation();
    vertexPrompt->SetUserInformation(infoPrompt);
    vertexPrompt->SetPosition(vtxPosition.x()*cm,vtxPosition.y()*cm,vtxPosition.z()*cm);


    GeneratePromptGammaSodium(vertexPrompt);
    infoPrompt->SetPromptGammaGen(true);
    infoPrompt->SetLifetime(promptLifetime/ps);
    infoPrompt->SetVtxPosition(vtxPosition.x(),vtxPosition.y(),vtxPosition.z());
    event->AddPrimaryVertex(vertexPrompt);


}

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{
    // create vertex of 2g/ 3g and if needed add de-excitation gamma quanta to this vertex
    G4double time = 0*s;


    //G4ThreeVector* boost = new G4ThreeVector(0,0,0);
    //HepBoost* boost = new HepBost(Hep3Vector direction(1,0,0),0.1);

    //G4ThreeVector vtxPosition;
    //int a;
    //std::tie(vtxPosition, a) = GetVerticesDistribution();

    //G4PrimaryVertex* vertex = new G4PrimaryVertex(VertexUniformInCylinder(25*cm,20*cm), time);
    G4PrimaryVertex* vertex = new G4PrimaryVertex(VertexUniformInCylinder(1*mm,1*mm), time);
    const G4ThreeVector bo(0.1,0,0);

    //    GeneratePromptGammaSodium(vertex);
    GenerateTwoGammaVertex(vertex);
    GeneratePromptGammaSodium(vertex);

    //GenerateThreeGammaVertex(vertex);

    event->AddPrimaryVertex(vertex);


}


G4ThreeVector PrimaryGenerator::VertexUniformInCylinder(G4double rIn, G4double zmax)
{
    //vertex A uniform on a cylinder
    //  
    //const G4double rSquare = 144*cm;
    //const G4double zmax = 34*cm;
    //
    //G4double r = std::sqrt(rSquare*G4UniformRand());
    G4double r = std::sqrt(pow(rIn,2)*G4UniformRand());

    G4double alpha = twopi*G4UniformRand();     //alpha uniform in (0, 2*pi)
    G4double ux = std::cos(alpha);
    G4double uy = std::sin(alpha);
    G4double z = zmax*(2*G4UniformRand() - 1);  //z uniform in (-zmax, +zmax)
    G4ThreeVector positionA(r*ux,r*uy,z);

    return positionA;
}



std::tuple<G4ThreeVector,G4double,G4double> PrimaryGenerator::GetVerticesDistribution(G4double maxXhalf, G4double maxYhalf, G4double maxZhalf)
{

    G4bool lookForVtx = false;
    G4ThreeVector myPoint(0,0,0);
    MaterialExtension* mat;

    // annihilation will occure only in materials where it was allowed @see MaterialExtension
    // annihilation rate 2g/3g also depends on the material type
    // now assumed equal distribution in the target - this may be modified in the future
    while (!lookForVtx)
    {
        G4double x_tmp = maxXhalf*(2*G4UniformRand() - 1)*cm;
        G4double y_tmp = maxYhalf*(2*G4UniformRand() - 1)*cm;
        G4double z_tmp = maxZhalf*(2*G4UniformRand() - 1)*cm;

        myPoint.setX(x_tmp);
        myPoint.setY(y_tmp);
        myPoint.setZ(z_tmp);

        G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
        mat = (MaterialExtension*)theNavigator->LocateGlobalPointAndSetup(myPoint)->GetLogicalVolume()->GetMaterial()  ; 
        lookForVtx = mat->IsTarget();

    };

    G4double ratio3g = mat->Get3gFraction();
    G4double  lifetime3g = mat->GetoPsLifetime();


    return std::make_tuple(myPoint,ratio3g,lifetime3g);
}



//void PrimaryGenerator::GenerateTwoGammaVertex(G4PrimaryVertex* vertex )
//{
//    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//    G4ParticleDefinition* particleDefinition = particleTable->FindParticle("gamma");
//
//
//    G4double mass = 1022*keV;
//    std::vector<G4double> mass_secondaries = {0., 0.};
//
//    std::vector<G4LorentzVector> out;
//
//    G4HadPhaseSpaceGenbod* phS = new G4HadPhaseSpaceGenbod();
//    phS->Generate(mass,mass_secondaries,out);
//
//    // boost gamma quanta
//    for (int i=0; i<2; i++)
//    {
//        // out[i].boost(BoostAxis,0.1);
//
//        G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition,
//                out[i].px(),out[i].py(),out[i].pz(),out[i].e());
//        PrimaryParticleInformation* info = new PrimaryParticleInformation();
//        info->SetGammaMultiplicity(2);
//        info->SetIndex(i+1);
//        particle1->SetUserInformation(info);
//
//        vertex->SetPrimary(particle1);
//    }
//
//}


void PrimaryGenerator::GenerateTwoGammaVertex(G4PrimaryVertex* vertex )
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefinition = particleTable->FindParticle("gamma");


    Double_t mass_secondaries[2] = {0., 0.};

    TGenPhaseSpace event;
    TLorentzVector vec_pozytonium(0.0,0.0,0.0,1022);
    Bool_t test =  event.SetDecay(vec_pozytonium, 2, mass_secondaries);
    if( !test){
        std::cout   << "error: generate_gamma : createTwoEvts:" << test << std::endl;  
    }


    event.Generate();   
    G4PrimaryParticle* particle[2];

    for (int i=0; i<2; i++)
    {
        TLorentzVector * out = event.GetDecay(i);

        particle[i] = new G4PrimaryParticle(particleDefinition,
                out->Px()*keV,out->Py()*keV,out->Pz()*keV,out->E()*keV);


        PrimaryParticleInformation* info = new PrimaryParticleInformation();
        info->SetGammaMultiplicity(2);
        info->SetIndex(i+1);
        particle[i]->SetUserInformation(info);
        vertex->SetPrimary(particle[i]);
    }

      G4ThreeVector gammaMom = particle[0]->GetMomentum();
      G4ThreeVector polarization1 = gammaMom.orthogonal().unit();
      polarization1 = polarization1.rotate( twopi * G4UniformRand() , gammaMom);
      particle[1]->SetPolarization( polarization1 );

      G4ThreeVector polarization2 = polarization1;

      polarization2 = polarization2.rotate( pi/2.0, gammaMom);
      particle[0]->SetPolarization( polarization2 );

}


void PrimaryGenerator::GenerateThreeGammaVertex(G4PrimaryVertex* vertex )
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefinition = particleTable->FindParticle("gamma");


    Double_t mass_secondaries[3] = {0., 0., 0.};

    TGenPhaseSpace event;
    TLorentzVector vec_pozytonium(0.0,0.0,0.0,1022);
    Bool_t test =  event.SetDecay(vec_pozytonium, 3, mass_secondaries);
    if( !test){
        std::cout   << "error: generate_gamma : createThreeEvts:" << test << std::endl;  
    }


    Double_t weight;
    Double_t weight_max= event.GetWtMax()*pow(10,5);
    Double_t rwt;  
    Double_t M_max = 7.65928*pow(10,-6);    
    do{
        weight = event.Generate();        
        weight = weight*calculate_mQED(511,event.GetDecay(0)->E(),event.GetDecay(1)->E(),event.GetDecay(2)->E());
        rwt = M_max*weight_max*(G4UniformRand());
    }while( rwt > weight );    

    event.Generate();   
    for (int i=0; i<3; i++)
    {
        TLorentzVector * out = event.GetDecay(i);

        G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition,
                out->Px()*keV,out->Py()*keV,out->Pz()*keV,out->E()*keV);

        PrimaryParticleInformation* info = new PrimaryParticleInformation();
        info->SetGammaMultiplicity(3);
        info->SetIndex(i+1);
        particle1->SetUserInformation(info);
        vertex->SetPrimary(particle1);
    }


}


void PrimaryGenerator::GeneratePromptGammaSodium(G4PrimaryVertex* vertex ) 
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefinition = particleTable->FindParticle("gamma");

    const G4double ene = 1.2770*MeV;

    double theta = 2 * M_PI * G4UniformRand();
    double phi = acos(1 - 2 * G4UniformRand());
    double px = ene * sin(phi) * cos(theta);
    double py = ene * sin(phi) * sin(theta);
    double pz = ene * cos(phi);

    G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition,
            px, py, pz,ene);
    PrimaryParticleInformation* info = new PrimaryParticleInformation();
    info->SetGammaMultiplicity(1);
    info->SetIndex(1);
    particle1->SetUserInformation(info);

    vertex->SetPrimary(particle1);
    //printf(" %f \n", sqrt(pow(px,2)+pow(py,2)+pow(pz,2)));

}


Double_t PrimaryGenerator::calculate_mQED(Double_t mass_e, Double_t w1, Double_t w2, Double_t w3){ 
    Double_t mQED = pow((mass_e-w1)/(w2*w3),2) + pow((mass_e-w2)/(w1*w3),2) + pow((mass_e-w3)/(w1*w2),2); 
    return mQED;
}
