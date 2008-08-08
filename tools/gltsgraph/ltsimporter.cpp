#include "ltsimporter.h"
#include <mcrl2/lts/liblts.h>
#include <map>

using namespace mcrl2::lts;

Graph* LTSImporter::importFile(std::string fn)
{
  Graph* result = new Graph();

  lts fileLTS;
  if(fileLTS.read_from(fn))
  {
    unsigned int initialState = fileLTS.initial_state(); 
   
    std::map<unsigned int, State*> numsToStates;
    
    for(state_iterator si = fileLTS.get_states(); si.more(); ++si)
    {
      unsigned int stNum = *si;
      State* s = new State(stNum,
                        stNum == initialState);
      result->addState(s);

      std::pair<unsigned int, State*> pNumToState(stNum, s);
      numsToStates.insert(pNumToState);

      if(s->isInitialState())
      {
        result->setInitialState(s);
      }
      
      // Generate a random position (x, y) for this state
      int x = static_cast<int> (
                (rand() / static_cast<float>(RAND_MAX) - .5) * 2000
              );
      int y = static_cast<int> (
                (rand() / static_cast<float>(RAND_MAX) - .5) * 2000
              );
      s->setX(x);
      s->setY(y);
    }

    for(transition_iterator ti = fileLTS.get_transitions(); ti.more(); ++ti)
    {
      unsigned int idFrom, idTo; 
      State *stFrom, *stTo;

      std::string label = fileLTS.label_value_str(
                            ti.label());
      idFrom = ti.from();
      idTo = ti.to();

      stFrom = numsToStates.find(idFrom)->second;
      stTo = numsToStates.find(idTo)->second;

      Transition* t = new Transition(stFrom, stTo, label);

      if (idFrom != idTo)
      {
        // t is a transition from stFrom to stTo.
        stFrom->addOutTransition(t);
        stTo->addInTransition(t);
      }
      else
      {
        // t is a self-loop, stFrom and stTo are the same
        stFrom->addSelfLoop(t);
      }
    }

    int numStates = fileLTS.num_states();
    int numLabels = fileLTS.num_labels();
    int numTrans = fileLTS.num_transitions();

    result->setInfo(initialState, numStates, numTrans, numLabels); 
  }
  else
  {
    //(CT) TODO: Error handling
  }

  return result;
}
