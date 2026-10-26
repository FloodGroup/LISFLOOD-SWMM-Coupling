#ifndef HYDRAULIC_RESTART_IO_H
#define HYDRAULIC_RESTART_IO_H
#include <stdio.h>
static int hydraulic_restart_bytes(FILE* f, void* data, size_t bytes, int reading)
{
    if (!bytes) return 0;
    return (reading ? fread(data,1,bytes,f) : fwrite(data,1,bytes,f)) != bytes;
}
#define HR_VALUE(value) do { if(hydraulic_restart_bytes(f,&(value),sizeof(value),reading)) return 1; } while(0)
#define HR_ARRAY(value,count) do { if(hydraulic_restart_bytes(f,(value),sizeof(*(value))*(count),reading)) return 1; } while(0)
#endif
