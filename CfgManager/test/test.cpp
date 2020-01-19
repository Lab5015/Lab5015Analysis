#include "interface/CfgManagerT.h"

int main(int argc, char* argv[])
{
    CfgManager cfg(argv[1]);

    cfg.Print();

    std::string test_string = cfg.GetOpt<std::string>("test.stringa");
    std::cout << test_string << std::endl;
    
    return 0;
}
