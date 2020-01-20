#ifndef OC_SET_H
#define OC_SET_H

#include <pthread.h>

typedef struct{
    pthread_mutex_t* lock;
    void* set;
} OcSetInt;

#endif // OC_SET_H
