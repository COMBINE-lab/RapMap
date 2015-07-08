#include <iostream>
#include <vector>

int rapMapIndex(int argc, char* argv[]);
int rapMapMap(int argc, char* argv[]);

int main(int argc, char* argv[]) {

    std::vector<char*> args;
    args.push_back(argv[0]);
    for (int i = 2; i < argc; ++i) {
        args.push_back(argv[i]);
    }

    if (std::string(argv[1]) == "index") {
        return rapMapIndex(argc - 1, args.data()); 
    } else if (std::string(argv[1]) == "map") {
        return rapMapMap(argc - 1, args.data()); 
    } else {
        std::cerr << "the command " << argv[1]
                  << " is not yet implemented\n";
        return 1;
    }
    return 0;
}
