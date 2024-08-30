
#ifndef VTK_META_EXPORT_H
#define VTK_META_EXPORT_H

#ifdef VTK_META_STATIC_DEFINE
#  define VTK_META_EXPORT
#  define VTK_META_NO_EXPORT
#else
#  ifndef VTK_META_EXPORT
#    ifdef vtkLocalExample_EXPORTS
        /* We are building this library */
#      define VTK_META_EXPORT 
#    else
        /* We are using this library */
#      define VTK_META_EXPORT 
#    endif
#  endif

#  ifndef VTK_META_NO_EXPORT
#    define VTK_META_NO_EXPORT 
#  endif
#endif

#ifndef VTK_META_DEPRECATED
#  define VTK_META_DEPRECATED __attribute__ ((__deprecated__))
#  define VTK_META_DEPRECATED_EXPORT VTK_META_EXPORT __attribute__ ((__deprecated__))
#  define VTK_META_DEPRECATED_NO_EXPORT VTK_META_NO_EXPORT __attribute__ ((__deprecated__))
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define VTK_META_NO_DEPRECATED
#endif



#endif
