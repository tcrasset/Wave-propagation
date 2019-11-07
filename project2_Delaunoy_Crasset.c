#include <stdio.h>
#include <stdlib.h>

int main(int argc, char const* argv[]) {
    char* parameter_file;
    char* map_file;
    unsigned int scheme;

    assert(argc == 4);  // Check number of arguments
    //Scan arguments
    sscanf(argv[1], "%u", &parameter_file);
    sscanf(argv[2], "%u", &map_file);
    sscanf(argv[3], "%f", &scheme);

    //Check argument validity
    assert((scheme == 0) || (scheme == 1));

    // Explicit
    if (scheme == 0) {
        return 0;
    }
    // Implicit
    else {
        return 0;
    }
    /* code */
    return 0;
}
