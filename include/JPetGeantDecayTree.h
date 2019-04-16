#ifndef JPetGeantDecayTree_h
#define JPetGeantDecayTree_h

#include "TObject.h"
#include <vector>
#include "TVector3.h" 


/**
 * \class JPetGeantDecayTree
 * \brief store decay tree structures (in form of vertices and tracks) about preformed simulation
 * Class is not yet implemented
 *
 */
class JPetGeantDecayTree : public TObject 
{
 
    public:
        JPetGeantDecayTree();
        ~JPetGeantDecayTree();

    private:


     ClassDef(JPetGeantDecayTree,1)

};

#endif
