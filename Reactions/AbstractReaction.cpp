
#include "AbstractReaction.hpp"

double AbstractReaction::GetRateConstant()
{
    return mRateConstant;
}

void AbstractReaction::SetReactionName(string reaction_name)
{
    mReactionName = move(reaction_name);
}

string AbstractReaction::GetReactionName()
{
    return mReactionName;
}
