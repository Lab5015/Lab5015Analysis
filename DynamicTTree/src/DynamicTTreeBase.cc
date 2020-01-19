#include "interface/DynamicTTreeBase.h"

std::map<std::type_index, const char*> DynamicTTreeBase::type_map={
    {typeid(bool), "/O"},
    {typeid(char), "/C"},
    {typeid(short int), "/S"},
    {typeid(unsigned short int), "/s"},
    {typeid(int), "/I"},
    {typeid(unsigned int), "/i"},
    {typeid(float), "/F"},
    {typeid(double), "/D"},
    {typeid(long int), "/L"},
    {typeid(unsigned long int), "/l"}
};
            
bool DynamicTTreeBase::NextEntry(long int entry)
{
    if(entry > -1)
        currentEntry_ = entry;

    ++currentEntry_;
    if(currentEntry_ < tree_->GetEntriesFast())
    {
        tree_->GetEntry(currentEntry_);
        return true;
    }

    currentEntry_=-1;
    return false;
}
