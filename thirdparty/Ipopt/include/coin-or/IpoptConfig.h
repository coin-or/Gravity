/* src/Common/config_ipopt.h.  Generated from config_ipopt.h.in by configure.  */
/* src/Common/config_ipopt.h.in. */

#ifndef __CONFIG_IPOPT_H__
#define __CONFIG_IPOPT_H__

/* Version number of project */
#define IPOPT_VERSION "3.13.2"

/* Major Version number of project */
#define IPOPT_VERSION_MAJOR 3

/* Minor Version number of project */
#define IPOPT_VERSION_MINOR 13

/* Release Version number of project */
#define IPOPT_VERSION_RELEASE 2

/* Define to the debug sanity check level (0 is no test) */
#define IPOPT_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
#define IPOPT_VERBOSITY 0

/* Define to the C type corresponding to Fortran INTEGER */
#define IPOPT_FORTRAN_INTEGER_TYPE int

/* Library Visibility Attribute */
#define IPOPTAMPLINTERFACELIB_EXPORT __declspec(dllimport)

/* Library Visibility Attribute */
#define IPOPTLIB_EXPORT __declspec(dllimport)

/* for backward compatibility: will be removed in Ipopt 3.14 */
#define FORTRAN_INTEGER_TYPE  IPOPT_FORTRAN_INTEGER_TYPE
#define COIN_IPOPT_CHECKLEVEL IPOPT_CHECKLEVEL
#define COIN_IPOPT_VERBOSITY  IPOPT_VERBOSITY

#endif
