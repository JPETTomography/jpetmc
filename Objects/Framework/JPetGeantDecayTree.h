/**
 *  @copyright Copyright 2020 The J-PET Monte Carlo Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file JPetGeantDecayTree.h
 */

#ifndef JPETGEANTDECAYTREE_H
#define JPETGEANTDECAYTREE_H 1

#include <TObject.h>
#include <TVector3.h>
#include <iostream>
#include <vector>
#include <map>

/**
 * @class JPetGeantDecayTree
 * @brief Class stores decay tree structures (in form of nodes and tracks)
 */

enum InteractionType
{
  kPrimaryGamma,
  kScattActivePart,
  kScattNonActivePart,
  kSecondaryPart,
  kUnknownInteractionType
};

enum DecayChannel { 
  Para2G, Direct2G, Ortho2G, Para3G, Direct3G, Ortho3G, Unknown
};

struct Branch {
  Branch() {};
  Branch(int trackID, int primaryBranch);
  int fTrackID = -1;             //ID of the track corresponding to this branch
  std::vector<int> fNodeIDs;    //container for all of the nodes
  std::vector<InteractionType> fInteractionType;
  int fPrimaryBranchID = -1;       //-1 for branch coming from primary photon, primary branchId otherwise
  
  void AddNodeID(int nodeID, InteractionType interactionType);
  int GetTrackID() const { return fTrackID; };
  int GetPrimaryNodeID() const { return fNodeIDs[0]; };
  int GetLastNodeID() const { return fNodeIDs[fNodeIDs.size()-1]; };
  int GetPrimaryBranchID() const { return fPrimaryBranchID; };
  int GetPreviousNodeID(int nodeID) const;
  InteractionType GetInteractionType(int nodeID) const;
};

class JPetGeantDecayTree : public TObject
{

public:
  JPetGeantDecayTree();
  ~JPetGeantDecayTree();
  
  void Clean();
  void ClearVectors();
  
  void SetEventNumber(int eventID) { fEventID = eventID; };
  void SetDecayChannel(DecayChannel decayChannel) { fDecayChannel = decayChannel; };
  int FindPrimaryPhoton(int nodeID);
  void AddNodeToBranch(int nodeID, int trackID, InteractionType interactionType);
  Branch GetBranch(unsigned trackID) const;
  int GetEventNumber() { return fEventID; };
  DecayChannel GetDecayChannel() { return fDecayChannel; };

private:
  int fEventID;
  DecayChannel fDecayChannel;
  std::vector<Branch> fBranches;
  std::map<int, int> fTrackBranchConnection;
     
  ClassDef(JPetGeantDecayTree,3)

};

#endif /* !JPETGEANTDECAYTREE_H */
