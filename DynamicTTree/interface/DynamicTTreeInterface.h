#ifdef DYNAMIC_TREE_NAME

#ifndef DATA_TABLE
#define DATA_TABLE
#endif
#ifndef DATA_VECT_TABLE
#define DATA_VECT_TABLE
#endif
#ifndef DATA_CLASS_TABLE
#define DATA_CLASS_TABLE
#endif

class DYNAMIC_TREE_NAME: public DynamicTTreeBase
{
public:
    //---dynamical data memeber name
    //---basic types
#define DATA(t, name) argument_type<void(t)>::type name;
    DATA_TABLE                                                          
#undef DATA
    //---c arrays
#define DATA(t, name, size) argument_type<void(t)>::type* name;
    DATA_VECT_TABLE                                                          
#undef DATA
    //---c++ classes
#define DATA(t, name, ...) argument_type<void(t __VA_ARGS__)>::type* name;
    DATA_CLASS_TABLE                                                          
#undef DATA    
    
    //---ctors---
    //---costructor for a new TTree
    DYNAMIC_TREE_NAME(const char* name="", const char* title=""):
        DynamicTTreeBase()
        {
            //---create new TTree
            tree_ = new TTree(name, title);

            //---create branches
            //---basic types
            std::string leaf;
#define DATA(t, name)                                                   \
            name=0;                                                     \
            leaf = std::string(#name)+type_map[typeid(argument_type<void(t)>::type)]; \
            tree_->Branch(#name, &name, leaf.c_str()); 
DATA_TABLE                                                          
#undef DATA
    //---c arrays
#define DATA(t, name, size)                                             \
    name = new argument_type<void(t)>::type[size]();                    \
    leaf = std::string(#name)+"["+#size+"]"+type_map[typeid(argument_type<void(t)>::type)]; \
    tree_->Branch(#name, name, leaf.c_str());                                  
DATA_VECT_TABLE                                                          
#undef DATA
    //---c++ classes
#define DATA(t, name, ...)                       \
    name=new argument_type<void(t __VA_ARGS__)>::type(); \
    tree_->Branch(#name, &name);
DATA_CLASS_TABLE                                                          
#undef DATA
    
    }

    //---costructor for already existing TChain
    DYNAMIC_TREE_NAME(TChain* t):
        DynamicTTreeBase()
        {
            //---save TTree ptr
            tree_ = t;

            //---set branches
            //---basic types
#define DATA(t, name) name=0; tree_->SetBranchAddress(#name, &name);
DATA_TABLE                                                          
#undef DATA

    //---get first entry in case c-arrays range depends on one of the previous variables
    tree_->GetEntry(0);

    //---c array
#define DATA(t, name, size) name=new argument_type<void(t)>::type[size](); tree_->SetBranchAddress(#name, name);
DATA_VECT_TABLE                                                          
#undef DATA
    //---c++ classes
#define DATA(t, name, ...) name=new argument_type<void(t __VA_ARGS__)>::type(); tree_->SetBranchAddress(#name, &name);
DATA_CLASS_TABLE                                                          
#undef DATA
    
    }

    //---costructor for already existing TTree
    DYNAMIC_TREE_NAME(TTree* t):
        DynamicTTreeBase()
        {
            //---save TTree ptr
            tree_ = t;

            //---set branches
            //---basic types
#define DATA(t, name) name=0; tree_->SetBranchAddress(#name, &name);
DATA_TABLE                                                          
#undef DATA

    //---get first entry in case c-arrays range depends on one of the previous variables
    tree_->GetEntry(0);

    //---c array
#define DATA(t, name, size) name=new argument_type<void(t)>::type[size](); tree_->SetBranchAddress(#name, name);
DATA_VECT_TABLE                                                          
#undef DATA
    //---c++ classes
#define DATA(t, name, ...) name=new argument_type<void(t __VA_ARGS__)>::type(); tree_->SetBranchAddress(#name, &name);
DATA_CLASS_TABLE                                                          
#undef DATA
    
    }
    
    //---dtor---
    ~DYNAMIC_TREE_NAME() {};
};

#endif
