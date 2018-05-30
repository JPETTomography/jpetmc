#include "SourceParams.hh"

SourceParams::SourceParams()
    : howManyGammas(1), shape("cylinder"),shapeCenterPosition(0,0,0)
{}


SourceParams::~SourceParams()
{}

void SourceParams::SetShape(G4String& newShape)
{
    if (std::find(std::begin(allowedShapes), std::end(allowedShapes), newShape) != std::end(allowedShapes))
    {
        shape = newShape;
    } else {
        G4Exception("SourceParams","SP01",FatalException,   
                "shape is not allowed");
    }

}
