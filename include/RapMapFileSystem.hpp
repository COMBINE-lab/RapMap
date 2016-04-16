#ifndef __RAPMAP_FILESYSTEM_HPP__
#define __RAPMAP_FILESYSTEM_HPP__

namespace rapmap {
    namespace fs {
        // Taken from http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
        bool FileExists(const char *path);
        // Taken from http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
        bool DirExists(const char *path);
        void MakeDir(const char* path);
    }
}


#endif //__RAPMAP_FILESYSTEM_HPP__
