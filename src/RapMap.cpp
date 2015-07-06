#include <iostream>

int rapMapIndex(int argc, char* argv[]);
int rapMapMap(int argc, char* argv[]);

int main(int argc, char* argv[]) {
    if (std::string(argv[1]) == "index") {
        return rapMapIndex(argc - 2, &argv[2]);
    } else if (std::string(argv[1]) == "map") {
        return rapMapMap(argc - 2, &argv[2]);      
    } else {
        std::cerr << "the command " << argv[1]
                  << " is not yet implemented\n";
        return 1;
    }
    return 0;
}
