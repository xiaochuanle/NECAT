#include <stdlib.h>

#include <stdio.h>

int main(int argc, char* argv[])
{
    const char* cmd = "rm ok.txt";
    int r = system(cmd);
    if (r != 0) {
        printf("command execute fail: %d\n", r);
    } else {
        printf("command execute successufally: %d\n", r);
    }
    return 0;
}
