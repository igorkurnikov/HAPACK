#ifndef FC_FUNC_HEADER_INCLUDED
#define FC_FUNC_HEADER_INCLUDED

#if defined(_MSC_VER) || defined(_IFORT_NAMES)

#define FC_FUNC_GLOBAL(name) name##_
#define FC_FUNC_GLOBAL_(name) name##_
#define FC_FUNC_MODULE(mod_name,name) mod_name##_mp_##name##_
#define FC_FUNC_MODULE_(mod_name,name) mod_name##_mp_##name##_

#else

#define FC_FUNC_GLOBAL(name) name##_
#define FC_FUNC_GLOBAL_(name) name##_
#define FC_FUNC_MODULE(mod_name,name) __##mod_name##_MOD_##name
#define FC_FUNC_MODULE_(mod_name,name) __##mod_name##_MOD_##name

#endif

 
#endif
