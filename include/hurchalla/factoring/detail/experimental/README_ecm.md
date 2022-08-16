
All the ECM headers and ECM C++ and C files in this folder are experimental.
It's important to note that these particular files are currently (temporarily) unavailable to you under any license.  You can contact me for express permission, but I expect I'll very soon offer them with a normal permissive license.

To use the experimental ECM algorithms, predefine the macro HURCHALLA_FACTORING_ALLOW_ECM_EXPERIMENTAL.  For example, if you want ECM and you are compiling with clang, you could invoke clang as follows:  
clang++ -DHURCHALLA_FACTORING_ALLOW_ECM_EXPERIMENTAL  ...more options and files...  
