#pragma once

#define VIRTUAL_VOID_MERGE_(a, b) a##b
#define VIRTUAL_VOID_LABEL_(a) VIRTUAL_VOID_MERGE_(unique_name_, a)
#define VIRTUAL_VOID_UNIQUE_NAME_ VIRTUAL_VOID_LABEL_(__LINE__)
#define VIRTUAL_VOID_UNIQUE_NAME VIRTUAL_VOID_UNIQUE_NAME_
#define __ VIRTUAL_VOID_UNIQUE_NAME_
