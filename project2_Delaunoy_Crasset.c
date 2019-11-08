#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char const* argv[]) {

    // Check number of arguments
    assert(argc == 4);  
    const char* parameter_file = argv[1];
    const char* map_file = argv[2];
    const unsigned int scheme = atoi(argv[3]);
    
    //Check argument validity
    assert((scheme == 0) || (scheme == 1));

    // Explicit
    if (scheme == 0) {
        printf("Explicit");
        printf("%s %s %u", parameter_file, map_file, scheme);
        return 0;
    }
    // Implicit
    else {
        printf("Implicit");
        printf("%s %s %u", parameter_file, map_file, scheme);
        return 0;
    }
    /* code */
    return 0;
}
