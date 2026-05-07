#include "fastcv/cli.h"

int main(int argc, char* argv[]) {
    std::vector<std::string> args(argv + 1, argv + argc);
    return fastcv::run_cli(args);
}
